library(tidyverse)

# ============================================================
# 📥 Load data
# ============================================================

not_contextualized <- read_csv(
  "analysis/models/exp_constrained/reaction_fluxes_FBA_pFBA.csv"
)

# contextualized <- read_csv(
#   "analysis/models/GIMME/reaction_fluxes_FBA_pFBA.csv"
# )

contextualized <- read_csv(
  "analysis/models/GIMME_25perc/reaction_fluxes_FBA_pFBA.csv"
)

# ============================================================
# 🔍  Pre-join diagnostics: reaction counts per model/context
# ============================================================

# Count reactions in each dataset
not_ctx_counts <- not_contextualized %>%
  count(common_name, Condition, Window, name = "n_not_ctx")

ctx_counts <- contextualized %>%
  count(common_name, Condition, Window, name = "n_ctx")

# Combine counts and calculate differences
reaction_counts <- not_ctx_counts %>%
  full_join(ctx_counts, by = c("common_name", "Condition", "Window")) %>%
  mutate(
    diff = n_not_ctx - n_ctx,
    percent_removed = (diff / n_not_ctx) * 100
  )

# Count overlapping reactions (actual intersection)
overlap_counts <- not_contextualized %>%
  inner_join(
    contextualized,
    by = c("common_name", "Condition", "Window", "Reaction")
  ) %>%
  count(common_name, Condition, Window, name = "n_overlap")

# ============================================================
# 🔗 1. Join datasets for flux comparison
# ============================================================

comparison_df <- not_contextualized %>%
  inner_join(
    contextualized,
    by = c("common_name", "Condition", "Window", "Reaction"),
    suffix = c("_not_ctx", "_ctx")
  )

# ============================================================
# ⚖️  Compare fluxes (with tolerance)
# ============================================================

tolerance <- 1e-6

comparison_df <- comparison_df %>%
  mutate(
    FBA_diff  = FBA_Flux_not_ctx  - FBA_Flux_ctx,
    pFBA_diff = pFBA_Flux_not_ctx - pFBA_Flux_ctx,
    
    FBA_equal_tol  = abs(FBA_diff)  < tolerance,
    pFBA_equal_tol = abs(pFBA_diff) < tolerance
  )

# ============================================================
# 📊  Summary statistics
# ============================================================

match_summary <- comparison_df %>%
  summarise(
    total_reactions = n(),
    FBA_match_rate_tol  = mean(FBA_equal_tol, na.rm = TRUE),
    pFBA_match_rate_tol = mean(pFBA_equal_tol, na.rm = TRUE)
  )

match_summary

# ============================================================
# ⚠️  Extract differences (beyond numerical noise)
# ============================================================

differences <- comparison_df %>%
  filter(!FBA_equal_tol | !pFBA_equal_tol)

# ============================================================
# 📈  Distribution of flux differences
# ============================================================

summary(comparison_df$FBA_diff)
summary(comparison_df$pFBA_diff)

plot_df_clean <- comparison_df %>%
  mutate(
    FBA_diff_clean  = ifelse(abs(FBA_diff)  < 1e-6, 0, FBA_diff),
    pFBA_diff_clean = ifelse(abs(pFBA_diff) < 1e-6, 0, pFBA_diff)
  ) %>%
  select(FBA_diff_clean, pFBA_diff_clean) %>%
  pivot_longer(
    cols = everything(),
    names_to = "method",
    values_to = "flux_diff"
  )

ggplot(plot_df_clean, aes(x = method, y = flux_diff, fill = method)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.1, size = 0.7) +
  theme_minimal() +
  labs(
    title = "Flux Differences",
    y = "Flux difference (non-contextualized − contextualized)"
  ) +
  theme(legend.position = "none")


plot_df_signal <- comparison_df %>%
  filter(abs(FBA_diff) > 1e-6 | abs(pFBA_diff) > 1e-6) %>%
  select(FBA_diff, pFBA_diff) %>%
  pivot_longer(
    cols = everything(),
    names_to = "method",
    values_to = "flux_diff"
  )

ggsave(
  filename = "figures/flux_diff_boxplot.png",
  width = 6,
  height = 4,
  dpi = 300
)

ggplot(plot_df_signal, aes(x = method, y = flux_diff, fill = method)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.1, size = 0.7) +
  theme_minimal() +
  labs(
    title = "Flux Differences (cleaned)",
    y = "Flux difference (non-contextualized − contextualized)"
  ) +
  theme(legend.position = "none")

comparison_df %>%
  count(FBA_diff, sort = TRUE) %>%
  head(10)
# ============================================================
# 🔝 2.  Top changing reactions (global)
# ============================================================

top_changes <- comparison_df %>%
  mutate(max_diff = pmax(abs(FBA_diff), abs(pFBA_diff))) %>%
  arrange(desc(max_diff)) %>%
  select(common_name, Condition, Window, Reaction,
         FBA_Flux_not_ctx, FBA_Flux_ctx, FBA_diff,
         pFBA_Flux_not_ctx, pFBA_Flux_ctx, pFBA_diff) %>%
  slice_head(n = 20)

top_changes


top_changes_filtered <- comparison_df %>%
  filter(abs(FBA_diff) > 1e-6 | abs(pFBA_diff) > 1e-6) %>%
  mutate(max_diff = pmax(abs(FBA_diff), abs(pFBA_diff))) %>%
  arrange(desc(max_diff)) %>%
  slice_head(n = 20)

top_changes_filtered

# ============================================================
# 🔁 Consistently changing reactions
# ============================================================

consistent_changes <- comparison_df %>%
  filter(!FBA_equal_tol | !pFBA_equal_tol) %>%
  count(Reaction, sort = TRUE)

consistent_changes

consistent_changes_weighted <- comparison_df %>%
  filter(!FBA_equal_tol | !pFBA_equal_tol) %>%
  group_by(Reaction) %>%
  summarise(
    n_changes = n(),
    mean_abs_diff = mean(pmax(abs(FBA_diff), abs(pFBA_diff))),
    .groups = "drop"
  ) %>%
  arrange(desc(n_changes), desc(mean_abs_diff))

consistent_changes_weighted

# ============================================================
# 🔁 3. Global scatter plots with metrics
# ============================================================

ggplot(comparison_df, aes(x = FBA_Flux_not_ctx, y = FBA_Flux_ctx)) +
  geom_point(alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  theme_minimal() +
  labs(
    title = "FBA Flux Comparison",
    x = "Not contextualized",
    y = "Contextualized"
  )

# ggplot(comparison_df, aes(x = FBA_Flux_not_ctx, y = FBA_Flux_ctx)) +
#   geom_point(alpha = 0.2) +
#   geom_abline(slope = 1, intercept = 0, color = "red") +
#   scale_x_continuous(trans = "pseudo_log") +
#   scale_y_continuous(trans = "pseudo_log") +
#   theme_minimal()

ggplot(comparison_df, aes(x = pFBA_Flux_not_ctx, y = pFBA_Flux_ctx)) +
  geom_point(alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, color = "blue") +
  theme_minimal() +
  labs(
    title = "pFBA Flux Comparison",
    x = "Not contextualized",
    y = "Contextualized"
  )

# ggplot(comparison_df, aes(x = pFBA_Flux_not_ctx, y = pFBA_Flux_ctx)) +
#   geom_point(alpha = 0.2) +
#   geom_abline(slope = 1, intercept = 0, color = "blue") +
#   scale_x_continuous(trans = "pseudo_log") +
#   scale_y_continuous(trans = "pseudo_log") +
#   theme_minimal()


# ============================================================
# 📊 4. Match rate per condition / window
# ============================================================

condition_summary <- comparison_df %>%
  group_by(common_name, Condition, Window) %>%
  summarise(
    n_reactions = n(),
    FBA_match_rate  = mean(FBA_equal_tol),
    pFBA_match_rate = mean(pFBA_equal_tol),
    mean_abs_FBA_diff  = mean(abs(FBA_diff)),
    mean_abs_pFBA_diff = mean(abs(pFBA_diff)),
    .groups = "drop"
  )

condition_summary

ggplot(condition_summary,
       aes(x = Window, y = FBA_match_rate, color = Condition, group = Condition)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ common_name) +
  theme_minimal() +
  labs(title = "FBA Match Rate Across Conditions")

ggplot(condition_summary,
       aes(x = Window, y = mean_abs_FBA_diff, color = Condition)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ common_name) +
  theme_minimal() +
  labs(title = "Mean Absolute FBA Difference")

