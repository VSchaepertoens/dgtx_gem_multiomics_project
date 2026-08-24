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

model_type <- "GIMME_tpm_25perc" #or GIMME_tpm_25perc or GIMME_25perc or exp_constrained

#load predicted fluxes
predicted <- read_csv(paste0("analysis/models/",model_type, "/reaction_fluxes_CV_ONLY_exp_constrained.csv"))

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
#R2 and MAE for a given FBA/pFBA flux

r2_by_rxn <- df %>%
  group_by(reaction_id, common_name, Condition) %>%
  summarize(
    r2_pFBA = cor(measured_flux, predicted_pFBA_flux, use = "complete.obs")^2,
    r2_FBA = cor(measured_flux, predicted_FBA_flux, use = "complete.obs")^2,
    n = n(),
    # MAE
    MAE_FBA = mean(
      abs(measured_flux - predicted_FBA_flux),
      na.rm = TRUE
    ),
    MAE_pFBA = mean(
      abs(measured_flux - predicted_pFBA_flux),
      na.rm = TRUE
    ),
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


# checking Asp NA ---------------------------------------------------------
predicted_ASP <- predicted_clean %>% filter(reaction_id == "EX_asp_L_e") 
# View(predicted_ASP)

ASP_r2 <- predicted_ASP %>%
  group_by(common_name, Condition) %>%
  summarize(
    r2_pFBA = cor(measured_flux, pFBA_Flux, use = "complete.obs")^2,
    r2_FBA = cor(measured_flux, FBA_Flux, use = "complete.obs")^2,
    n = n(),
    .groups = "drop"
  )

#For this model  ASP pFBA is always 0 → the correlation-based R² is undefined.

# plotting-------------------------------------------------------
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

#R2 as barplots ---------------------------------------------
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

ggsave("figures/models/GIMME_tpm_25perc/R2_distribution_boxplot_condition.png",
       width = 50,
       height = 35,
       units = "mm",
       dpi = 300)

#R2 as barplots Larissa theme ---------------------------------------------

# Bar plot settings
PLOT_FONT_SIZE  <- 12 #6
PANEL_LABEL_SIZE <- 12 #8
FONT_FAMILY     <- "Liberation Sans"

BAR_ALPHA       <- 0.25
ERRORBAR_ALPHA  <- 0.8
POINT_ALPHA     <- 0.6

POINT_SIZE      <- 1
LINE_WIDTH      <- 0.6

# Mean ± standard deviation
mean_sd <- function(x) {
  data.frame(
    y = mean(x, na.rm = TRUE),
    ymin = mean(x, na.rm = TRUE) - sd(x, na.rm = TRUE),
    ymax = mean(x, na.rm = TRUE) + sd(x, na.rm = TRUE)
  )
}

# Dodge settings
DODGE_WIDTH <- 0.9

ggplot(
  r2_long,
  aes(
    x = common_name,
    y = R2,
    fill = Condition,
    color = Condition
  )
) +
  
  # Individual observations, jittered within each Condition
  geom_point(
    aes(fill = Condition, color = Condition),
    alpha = POINT_ALPHA,
    size = POINT_SIZE,
    shape = 21,
    stroke = 0.5,
    position = position_jitterdodge(
      jitter.width = 0.1,
      dodge.width = DODGE_WIDTH
    ),
    show.legend = FALSE,
    na.rm = TRUE
  ) +
  
  # Boxplots, side-by-side by Condition 
  geom_boxplot( alpha = BAR_ALPHA, 
                outlier.shape = NA, 
                width = 0.7, 
                position = position_dodge(width = DODGE_WIDTH), 
                na.rm = TRUE ) +
  
  labs(
    x = "",
    y = bquote("Exchange rate" ~ R^2)
  ) +
  facet_wrap( ~method) +
  
  scale_fill_manual(values = color_condition_long) +
  scale_color_manual(values = color_condition_long) +
  
  theme_bw(
    base_size = PLOT_FONT_SIZE,
    base_family = FONT_FAMILY
  ) +
  
  theme_tidy() +
  
  theme(
    strip.text = element_text(size = PANEL_LABEL_SIZE),
    axis.title = element_text(size = PLOT_FONT_SIZE),
    axis.text = element_text(size = PLOT_FONT_SIZE),
    legend.title = element_text(size = PLOT_FONT_SIZE),
    legend.text = element_text(size = PLOT_FONT_SIZE),
    # axis.text.x = element_text(angle = 90)
  )

ggsave(paste0("figures/models/",model_type, "/R2_distribution_boxplot_condition_LH.png"),
       width = 50,
       height = 45,
       units = "mm",
       dpi = 300)

ggsave(paste0("figures/models/",model_type, "/poster/R2_distribution_boxplot_condition_LH.png"),
       width = 200,
       height = 150,
       units = "mm",
       dpi = 300)

#Summary table-----------------------------

# Create summary table
r2_summary <- r2_long %>%
  group_by(common_name, Condition, method) %>%
  summarise(
    n = sum(!is.na(R2)),
    mean_R2 = mean(R2, na.rm = TRUE),
    sd_R2 = sd(R2, na.rm = TRUE),
    median_R2 = median(R2, na.rm = TRUE),
    min_R2 = min(R2, na.rm = TRUE),
    max_R2 = max(R2, na.rm = TRUE),
    IQR_R2 = IQR(R2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    across(
      c(mean_R2, sd_R2, median_R2, min_R2, max_R2, IQR_R2),
      ~ round(.x, 2)
    )
  ) %>%
  arrange(common_name, Condition, method)


# Print summary table
r2_summary

# Save r2_summary as CSV 
write_csv( r2_summary, paste0("analysis/models/",model_type,"/r2_summary.csv") )

#MAE plot ---------------------------------------
mae_long <- r2_by_rxn %>%
  pivot_longer(
    cols = starts_with("MAE_"),
    names_to = "method",
    values_to = "MAE"
  ) %>%
  mutate(
    method = recode(method,
                    "MAE_pFBA" = "pFBA",
                    "MAE_FBA"  = "FBA"),
    reaction_id_ordered = reorder_within(reaction_id, MAE, common_name)
  )


mae_long <- mae_long %>%
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

mae_long %>%
  distinct(reaction_id, reaction_clean) %>%
  arrange(reaction_clean)

order_levels <- mae_long %>%
  filter(common_name == "iCHOv1", method == "pFBA") %>%
  arrange(desc(MAE)) %>%
  pull(reaction_clean) %>%
  unique()

mae_long <- mae_long %>%
  mutate(reaction_clean = factor(reaction_clean, levels = order_levels),
         common_name = factor(
           common_name,
           levels = c(
             "iCHOv1",
             "iCHO2441",
             "iCHO3K"
           )
         ))

mae_long %>%
  filter(grepl("glc", reaction_id)) %>%
  distinct(reaction_id, reaction_clean)

ggplot(
  mae_long,
  aes(
    x = common_name,
    y = MAE,
    fill = Condition,
    color = Condition
  )
) +
  
  # Individual observations, jittered within each Condition
  geom_point(
    aes(fill = Condition, color = Condition),
    alpha = POINT_ALPHA,
    size = POINT_SIZE,
    shape = 21,
    stroke = 0.5,
    position = position_jitterdodge(
      jitter.width = 0.1,
      dodge.width = DODGE_WIDTH
    ),
    show.legend = FALSE,
    na.rm = TRUE
  ) +
  
  # Boxplots, side-by-side by Condition 
  geom_boxplot( alpha = BAR_ALPHA, 
                outlier.shape = NA, 
                width = 0.7, 
                position = position_dodge(width = DODGE_WIDTH), 
                na.rm = TRUE ) +
  
  labs(
    x = "",
    y = "Mean absolute error"
  ) +
  facet_wrap( ~method) +
  
  scale_fill_manual(values = color_condition_long) +
  scale_color_manual(values = color_condition_long) +
  
  theme_bw(
    base_size = PLOT_FONT_SIZE,
    base_family = FONT_FAMILY
  ) +
  
  theme_tidy() +
  
  theme(
    strip.text = element_text(size = PANEL_LABEL_SIZE),
    axis.title = element_text(size = PLOT_FONT_SIZE),
    axis.text = element_text(size = PLOT_FONT_SIZE),
    legend.title = element_text(size = PLOT_FONT_SIZE),
    legend.text = element_text(size = PLOT_FONT_SIZE)
  )

ggsave(filename = paste0("figures/models/", model_type,"/exchange_rate_exp_pred_mae.png"),
       # device = cairo_pdf,
       height = 50,
       width = 140,
       units = "mm",
       dpi = 600,
       bg = "white")

ggsave(filename = paste0("figures/models/", model_type,"/poster/exchange_rate_exp_pred_mae.png"),
       # device = cairo_pdf,
       height = 150,
       width = 230,
       units = "mm",
       dpi = 600,
       bg = "white")

# Save r2_summary as CSV 
write_csv( mae_long,paste0("analysis/models/",model_type,"/mae_exchange_summary.csv" ))
# other plots -------------------------------------------------------------


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

ggsave("figures/models/exp_constrained/R2_distribution_boxplot_condition.png",
       width = 50,
       height = 35,
       units = "mm",
       dpi = 300)


