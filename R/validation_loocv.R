# libraries
source("R/plotting_utils.R", local = FALSE)

library(tidyverse)
library(janitor)
library(tidytext) 

#load csv for matcing metabolite with reaction_id

identifiers <- read_csv("models/original/exchange_igg_rxn_id.csv")
head(identifiers)



#load experimental data
measured <-read_csv("data/rates/rates_mM_gDCW_h_new.csv") %>%
  select(Condition, Window, Metabolite, rate_mM_gDCW_h, SD_mM_gDCW_h) %>%
  clean_names() %>%
  filter(!metabolite %in% c("LDH"))
head(measured)

# prepare measured data
measured_clean <- measured %>%
  filter(metabolite != "Titer") %>%   # just in case
  rename(label = metabolite)          # matches identifiers

#match fluxes with the experimntal rates
measured_mapped <- measured_clean %>%
  left_join(identifiers, by = "label")
#check
measured_mapped %>%
  filter(is.na(reaction_id)) %>%
  distinct(label)



#load predicted fluxes
# predicted <- read_csv("analysis/models/GIMME_25perc/reaction_fluxes_CV_ONLY_GIMME25perc.csv")
predicted <- read_csv("analysis/models/exp_constrained/reaction_fluxes_CV_ONLY_exp_constrained.csv")

head(predicted)

#prepare predicted data
predicted_clean <- predicted %>%
  filter(!Reaction %in% c("DM_igg[g]", "DM_igg_g")) %>%
  rename(reaction_id = Reaction)



#join measured with predicted
df <- predicted_clean %>%
  left_join(
    measured_mapped %>%
      select(condition, window, reaction_id, rate_m_m_g_dcw_h, sd_m_m_g_dcw_h),
    by = c("Condition" = "condition",
           "Window" = "window",
           "reaction_id")
  ) %>%
  select(model_name, common_name, Condition, Window, reaction_id, FBA_Flux, pFBA_Flux, rate_m_m_g_dcw_h, sd_m_m_g_dcw_h) %>%
  rename(
    measured_flux = rate_m_m_g_dcw_h,
    measured_sd = sd_m_m_g_dcw_h,
    predicted_pFBA_flux = pFBA_Flux,
    predicted_FBA_flux = FBA_Flux
  ) %>%
  mutate(Window = as.character(Window))

df %>%
  summarize(
    missing_measured = sum(is.na(measured_flux)),
    total = n()
  )
#R2 for a given FBA/pFBA flux

r2_by_rxn <- df %>%
  group_by(reaction_id, common_name, Condition) %>%
  summarize(
    r2_pFBA = cor(measured_flux, predicted_pFBA_flux, use = "complete.obs")^2,
    r2_FBA = cor(measured_flux, predicted_FBA_flux, use = "complete.obs")^2,
    n = n(),
    .groups = "drop"
  ) %>%
mutate(common_name = recode(
  common_name,
  "iCHO3K_cho_prod_generic_unblocked" = "iCHO3K"
  ),
  common_name = factor(
         common_name,
         levels = c(
           "iCHOv1",
           "iCHO2441",
           "iCHO3K")
       ))


#plot experimental vs predicted, facet by reaction id, add R2 text
color_mapping_condition <- c(
  "Constant" = "#E6641E",
  "Temp. shifted" = "#4B288C"
)

shape_mapping_window <- c(
  "1" = 16,
  "2" = 17,
  "3" = 15,
  "4" = 18,
  "5" = 8
)

# # only ichov1
# r2_plot <- r2_by_rxn %>%
#   filter(common_name == "iCHO2441") %>%
#   select(reaction_id, common_name, r2_pFBA)
# 
# df_plot <- df %>%
#   filter(common_name == "iCHO2441") %>%
#   left_join(r2_plot, by = c("reaction_id", "common_name"))


# ggplot(df_plot, aes(x = measured_flux, y = predicted_pFBA_flux)) +
#   geom_point(aes(color = Condition, shape = Window), 
#              size = 2.5, alpha = 0.8) +
#   geom_abline(slope = 1, intercept = 0, 
#               linetype = "dashed", color = "grey40") +
#   facet_wrap(~reaction_id, scales = "free") +
#   scale_color_manual(values = color_mapping_condition) +
#   
#   geom_text(
#     data = df_plot %>% distinct(reaction_id, r2_pFBA),
#     aes(label = paste0("R² = ", sprintf("%.2f", r2_pFBA))),
#     x = -Inf, y = Inf,
#     hjust = -0.1, vjust = 1.2,
#     size = 3.2,
#     inherit.aes = FALSE
#   ) +
#   
#   labs(
#     x = "Measured flux",
#     y = "Predicted flux (pFBA)",
#     color = "Condition",
#     shape = "Window"
#   ) +
#   
#   theme_minimal(base_size = 12) +
#   theme(
#     legend.position = "bottom",
#     strip.text = element_text(face = "bold"),
#     panel.grid.minor = element_blank(),
#     panel.spacing = unit(1, "lines")
#   )

# # barplot for the R2 only
# ggplot(r2_plot, aes(x = reorder(reaction_id, -r2_pFBA), y = r2_pFBA)) +
#   geom_col(fill = "steelblue") +
#   
#   labs(
#     x = "Reaction",
#     y = expression(R^2)
#   ) +
#   
#   theme_minimal(base_size = 12) +
#   theme(
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#     panel.grid.minor = element_blank()
#   )

r2_long <- r2_by_rxn %>%
  pivot_longer(
    cols = starts_with("r2_"),
    names_to = "method",
    values_to = "R2"
  ) %>%
  mutate(
    method = recode(method,
                    "r2_pFBA" = "pFBA",
                    "r2_FBA"  = "FBA"),
    reaction_id_ordered = reorder_within(reaction_id, R2, common_name)
  )


r2_long <- r2_long %>%
  mutate(
    reaction_clean = reaction_id %>%
      
      # remove prefix
      str_remove("^EX_") %>%
      
      # remove model suffix
      str_remove("___.*$") %>%
      
      # unify parentheses and underscores
      str_replace_all("\\(e\\)", "_e") %>%
      str_replace_all("_+", "_") %>%
      
      # REMOVE stereochemistry / compartment junk AFTER splitting
      str_remove("_.*$") %>%
      
      # final cleanup
      str_to_lower()
  )

r2_long %>%
  distinct(reaction_id, reaction_clean) %>%
  arrange(reaction_clean)

order_levels <- r2_long %>%
  filter(common_name == "iCHOv1", method == "pFBA") %>%
  arrange(desc(R2)) %>%
  pull(reaction_clean) %>%
  unique()

r2_long <- r2_long %>%
  mutate(reaction_clean = factor(reaction_clean, levels = order_levels),
         common_name = factor(
           common_name,
           levels = c(
             "iCHOv1",
             "iCHO2441",
             "iCHO3K"
           )
         ))

r2_long %>%
  filter(grepl("glc", reaction_id)) %>%
  distinct(reaction_id, reaction_clean)

ggplot(r2_long, aes(x = reaction_clean, y = R2, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~common_name, scales = "fixed") +
  
  labs(
    x = "Reaction",
    y = expression(R^2),
    fill = "Method"
  ) +
  scale_fill_manual(
    values = c(
      "pFBA" = "#4C78A8",
      "FBA" = "#59A14F"
    )
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 11)
  )

ggsave("figures/models/exp_constrained/R2_by_reaction_barplot.png",
       width = 300,
       height = 150,
       units = "mm",
       dpi = 300)

ggplot(r2_long, aes(x = method, y = R2)) +
  # Boxplot for all data (ignoring Condition)
  geom_boxplot(alpha = 0.2, outlier.shape = NA, 
               fill = "gray90", color = "gray50") +
  # Points colored by Condition, jittered
  geom_jitter(aes(color = Condition), 
              width = 0.15, alpha = 0.7, size = pt_to_geom(2)) +
  
  facet_wrap(~common_name) +
  
  labs(x = NULL, y = expression(R^2)) +
  scale_color_manual(values = color_mapping_condition) +
  
  theme_bw() +
  theme_tidy() +
  theme(axis.text.x = element_text(angle = 90))

ggsave("figures/models/exp_constrained/R2_distribution_boxplot_condition.png",
       width = 50,
       height = 35,
       units = "mm",
       dpi = 300)

ggplot(r2_by_rxn, aes(x = r2_FBA, y = r2_pFBA)) +
  geom_point(size = 2.5, alpha = 0.8, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  
  facet_wrap(~common_name) +
  
  labs(
    x = expression("FBA " * R^2),
    y = expression("pFBA " * R^2),
    title = "Method comparison per reaction"
  ) +
  
  
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 11),
    axis.text.x = element_text(angle = 90),
    
  )

ggsave("figures/models/exp_constrained/R2_method_comparison.png",
       width = 125,
       height = 75,
       units = "mm",
       dpi = 300)
