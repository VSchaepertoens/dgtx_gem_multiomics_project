source("R/plotting_utils.R", local = FALSE)
library(tidyverse)

measured <- read_csv("data/rates/growth_rate.csv") 

# GROWTH rates ------------------------------------------
## MEASURED ##
# m_growth_rates <- measured %>%
#   # filter(AA_meta == "Growth_rate") %>%
#   mutate(Type = "measured")
#   # select(!Condition)

predicted <- read_csv("analysis/models/exp_constrained/mus_FBA_pFBA.csv")
#predicted <- read_csv("analysis/models/GIMME/mus_FBA_pFBA_GIMME.csv")
# predicted <- read_csv("analysis/models/GIMME_tpm_25perc/mus_FBA_pFBA.csv")

color_mapping_condition <- c(
  "Constant" = "#E6641E",
  "Temp. shifted" = "#4B288C"
)

merged <- predicted %>%
  left_join(measured, by = c("Condition", "Window"))

df_long <- merged %>%
  pivot_longer(
    cols = c(FBA_mu, pFBA_mu),
    names_to = "prediction_type",
    values_to = "predicted_mu"
  ) %>%
  mutate(Window = as.character(Window),
         common_name = recode(
           common_name,
           "iCHO3K_cho_prod_generic_unblocked" = "iCHO3K"
         ),
         common_name = factor(
           common_name,
           levels = c("iCHOv1", "iCHO2441", "iCHO3K")
         ))

#plot
# Define abline parameters
abline_params <- list(
  geom_abline(intercept = 0, slope = 1, linewidth = 0.5),
  geom_abline(intercept = 0, slope = 0.75, linetype = 2, color = 'grey40'),
  geom_abline(intercept = 0, slope = 1.25, linetype = 2, color = 'grey40')
)

r2_df <- df_long %>%
  group_by(common_name, Condition) %>%
  summarize( r2 = round( cor( avg_growth_µ, predicted_mu, use = "complete.obs" )^2, 2 ), n = n(), .groups = "drop" ) %>%
  group_by(common_name) %>%
  mutate(
    label_y = c(0.20, 0.10)
  ) %>%
  ungroup()

r2_df

# Main plotting dataframe
df_plot <- df_long %>%
  filter(prediction_type == "pFBA_mu")

ggplot(df_plot, aes(x = avg_growth_µ, y = predicted_mu)) +
  abline_params +
  
  geom_point(aes(color = Condition, shape = Window), size = 1) +
  
  geom_errorbar(
    aes(
      xmin = avg_growth_µ - sd_growth_µ,
      xmax = avg_growth_µ + sd_growth_µ
    ),
    width = .00025,
    linewidth = .25,
    alpha = 0.5
  ) +
  
  # Add R² labels
  geom_text(
    data = r2_df,
    aes(
      label = paste0("R² = ", sprintf("%.2f", r2)),
      color = Condition
      # vjust = vjust_pos
    ),
    x = -Inf,
    y = Inf,
    hjust = -0.1,
    size = 2,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = color_mapping_condition) +
  
  scale_x_continuous(limits = c(0, 0.043)) +
  scale_y_continuous(limits = c(0, 0.043)) +
  coord_equal() +
  
  facet_wrap(~ common_name) +
  
  labs(
    x = expression(paste("Measured (Experimental) Rates [", h^-1, "]")),
    y = expression(paste("Predicted Rates [", h^-1, "]"))
  ) +
  theme_tidy() +
  # theme_bw() +
  theme(legend.position = "none")
  #       panel.grid.minor = element_blank(),
  #       axis.text = element_text(color = "black", size = 11),
  #       axis.title = element_text(color = "black", size = 11))

ggsave(filename = "figures/models/exp_constrained/growth_rate_exp_pred_legend.pdf",
       device = cairo_pdf,
       height = 3,
       width = 4,
       # units = "mm",
       dpi = 600,
       bg = "white")


# only plot R2 distribution -----------------------------------------------

ggplot(r2_df, aes(x = "Growth", y = r2, color = Condition)) +
  geom_point(alpha = 0.7, size = 2) +
  
  geom_text(
    data = r2_df,
    aes(
      y = label_y,
      label = sprintf("R² = %.2f", r2),
      color = Condition
    ),
    x = -Inf,
    hjust = -0.1,
    vjust = 0,
    inherit.aes = FALSE,
    size = pt_to_geom(6),
    show.legend = FALSE
  ) +
  
  facet_wrap(~common_name) +

  
  scale_y_continuous(limits = c(0, 1)) +
  
  scale_color_manual(values = color_mapping_condition) +
  
  theme_bw() +
  labs(x = NULL, 
       y = expression(R^2)) +
  theme_tidy() 


ggsave("figures/models/exp_constrained/R2_distribution_boxplot_growth rate.png",
       width = 50,
       height = 35,
       units = "mm",
       dpi = 300)


# R2 as barplots ----------------------------------------------------------

ggplot(r2_df, aes(x = common_name, y = r2, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), 
            width = 0.7) +
  scale_fill_manual(values = color_mapping_condition) +
  coord_cartesian(ylim = c(0.75, 1)) +
  labs(x = "", y = bquote("Growth rate" ~ R^2)) +
  theme_bw() +
  theme_tidy()


ggsave("figures/models/GIMME_tpm_25perc/R2_distribution_barplot_growth rate.png",
       width = 50,
       height = 35,
       units = "mm",
       dpi = 300)

#R2 as barplots Larissa theme ---------------------------------------------

# Bar plot settings
PLOT_FONT_SIZE  <- 6
PANEL_LABEL_SIZE <- 8
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
  r2_df,
  aes(
    x = common_name,
    y = r2,
    fill = Condition,
    color = Condition
  )
) +
  
  # Mean bars, side-by-side by Condition
  stat_summary(
    fun = mean,
    geom = "bar",
    alpha = BAR_ALPHA,
    width = 0.85,
    position = position_dodge(width = DODGE_WIDTH),
    na.rm = TRUE
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
  
  # Mean ± SD error bars, aligned with each Condition
  stat_summary(
    fun.data = mean_sd,
    geom = "errorbar",
    width = 0.2,
    linewidth = LINE_WIDTH,
    alpha = ERRORBAR_ALPHA,
    position = position_dodge(width = DODGE_WIDTH),
    na.rm = TRUE
  ) +
  coord_cartesian(ylim = c(0.75, 1)) +
  
  labs(
    x = "",
    y = bquote("Growth rate" ~ R^2)
  ) +
  
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

ggsave("figures/models/exp_constrained/R2_distribution_barplot_growth_rate_LH.png",
       width = 50,
       height = 35,
       units = "mm",
       dpi = 300)



# Save r2_summary as CSV 
write_csv( r2_df, "analysis/models/GIMME_tpm_25perc/r2_summary_growth.csv" )
# mean absolute error -----------------------------------------------------

error_df <- df_long %>%
  mutate(
    error = predicted_mu - avg_growth_µ,
    abs_error = abs(error),
    sq_error = error^2
  )

mae_summary <- error_df %>%
  group_by(common_name, prediction_type, Condition) %>%
  summarise(
    MAE = mean(abs_error, na.rm = TRUE),
    SD = sd(abs_error, na.rm = TRUE),
    n = n(),
    SE = SD / sqrt(n),
    .groups = "drop"
  )

ggplot(mae_summary,
       aes(x = common_name,
           y = MAE,
           fill = Condition)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_errorbar(
    aes(ymin = MAE - SE,
        ymax = MAE + SE),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  scale_y_continuous(limits = c(0, 0.01)) +
  facet_grid(~prediction_type) +
  scale_fill_manual(values = color_mapping_condition) +
  labs(
    x = "Model",
    y = "Mean absolute error",
    fill = "Condition"
  ) +

theme_tidy()

ggsave(filename = "figures/models/exp_constrained/growth_rate_exp_pred_mae.pdf",
       device = cairo_pdf,
       height = 3,
       width = 4,
       # units = "mm",
       dpi = 600,
       bg = "white")

# MAE barplot with LH theme



# Bar plot of MAE
ggplot(
  mae_summary,
  aes(
    x = common_name,
    y = MAE,
    fill = Condition,
    color = Condition
  )
) +
  
  # Mean MAE bars, side-by-side by Condition
  geom_col(
    alpha = BAR_ALPHA,
    width = 0.85,
    position = position_dodge(width = DODGE_WIDTH)
  ) +
  
  # MAE ± SE error bars
  geom_errorbar(
    aes(
      ymin = MAE - SE,
      ymax = MAE + SE
    ),
    width = 0.2,
    linewidth = LINE_WIDTH,
    alpha = ERRORBAR_ALPHA,
    position = position_dodge(width = DODGE_WIDTH)
  ) +
  
  # Set y-axis limits
  scale_y_continuous(
    limits = c(0, 0.01)
  ) +
  
  # Separate panels for prediction type
  facet_grid(~prediction_type) +
  
  # Condition colours
  scale_fill_manual(
    values = color_condition_long
  ) +
  scale_color_manual(
    values = color_condition_long
  ) +
  
  # Labels
  labs(
    x = "",
    y = "Mean absolute error",
    fill = "Condition",
    color = "Condition"
  ) +
  
  # Theme
  theme_bw(
    base_size = PLOT_FONT_SIZE,
    base_family = FONT_FAMILY
  ) +
  theme_tidy() +
  
  theme(
    strip.text = element_text(
      size = PANEL_LABEL_SIZE
    ),
    axis.title = element_text(
      size = PLOT_FONT_SIZE
    ),
    axis.text = element_text(
      size = PLOT_FONT_SIZE
    ),
    legend.title = element_text(
      size = PLOT_FONT_SIZE
    ),
    legend.text = element_text(
      size = PLOT_FONT_SIZE
    )
  )

ggsave(filename = "figures/models/exp_constrained/growth_rate_exp_pred_mae.pdf",
       device = cairo_pdf,
       height = 3,
       width = 4,
       # units = "mm",
       dpi = 600,
       bg = "white")

ggsave(filename = "figures/models/exp_constrained/growth_rate_exp_pred_mae.png",
       # device = cairo_pdf,
       height = 2,
       width = 3.5,
       # units = "mm",
       dpi = 600,
       bg = "white")

# Save r2_summary as CSV 
write_csv( mae_summary, "analysis/models/exp_constrained/mae_growth_summary.csv" )
