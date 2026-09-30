library(tidyverse)
library(ComplexHeatmap)
library(circlize)
source("R/utils/plotting_utils.R", local = FALSE)
library(tidytext)

flux_df <- read_csv("analysis/models/exp_constrained/reaction_fluxes_FBA_pFBA.csv")
flux_threshold <- 1e-9

# Classify reactions based on activity AND direction consistency
rxn_classes <- flux_df %>%
  mutate(
    active = abs(pFBA_Flux) > flux_threshold,
    direction = case_when(
      pFBA_Flux > flux_threshold ~ "positive",
      pFBA_Flux < -flux_threshold ~ "negative",
      TRUE ~ "inactive"
    )
  ) %>%
  group_by(common_name, Reaction) %>%
  summarise(
    n_contexts = n(),
    n_active = sum(active),
    n_positive = sum(direction == "positive"),
    n_negative = sum(direction == "negative"),
    direction_consistency = case_when(
      n_active == 0 ~ "Always inactive",
      n_positive == n_active ~ "Always active (positive only)",
      n_negative == n_active ~ "Always active (negative only)",
      n_positive > 0 & n_negative > 0 ~ "Direction switching",
      TRUE ~ "Context-specific (activity switching)"
    ),
    .groups = "drop"
  )

# Recode model name
rxn_classes <- rxn_classes %>%
  mutate(
    common_name = recode(
      common_name,
      "iCHO3K_cho_prod_generic_unblocked" = "iCHO3K"
    )
  )

# Barplot of reaction classes
class_summary <- rxn_classes %>%
  count(common_name, direction_consistency) %>%
  mutate(
    common_name = factor(
      common_name,
      levels = c(
        "iCHOv1",
        "iCHO2441",
        "iCHO3K"
      )
    )
  )

# Total reactions per model
totals <- class_summary %>%
  group_by(common_name) %>%
  summarise(total_rxns = sum(n), .groups = "drop")

# Calculate percentages relative to ALL reactions
class_summary <- class_summary %>%
  left_join(totals, by = "common_name") %>%
  mutate(percent = 100 * n / total_rxns)

# Extract inactive percentages for labels
inactive_labels <- class_summary %>%
  filter(direction_consistency == "Always inactive") %>%
  mutate(y_pos = 20,label = paste0("Always inactive\n", round(percent, 1), "%"))


# Keep only non-inactive classes for plotting
plot_df <- class_summary %>%
  filter(direction_consistency != "Always inactive")

pt_to_geom <- function(pt) pt / 2.845

#plot pie chart of class summary ----------------------
plot_class_pie <- function(model_name) {
  
  plot_data <- class_summary %>%
    filter(common_name == model_name) %>%
    mutate(
      direction_consistency = recode(
        direction_consistency,
        "Always active (negative only)" = "Active (negative)",
        "Always active (positive only)" = "Active (positive)"
      ),
      direction_consistency = factor(
        direction_consistency,
        levels = c(
          
          "Always inactive",
          "Active (positive)",
          "Active (negative)",
          "Direction switching"
        )
      ),
    ) %>%
    arrange(desc(direction_consistency)) %>%
    mutate(
      label = sprintf("%.1f%%", percent),
      ymax = cumsum(percent),
      ymin = lag(ymax, default = 0),
      midpoint = (ymax + ymin) / 2,
      angle = 90 - 360 * midpoint / 100,
      
      # Position labels outside the pie
      x_label = 1.65,
      hjust = ifelse(angle < -90, 1, 0)
    )
  
  ggplot(
    plot_data,
    aes(x = 1, y = percent, fill = direction_consistency)
  ) +
    
    # Pie
    geom_col(
      width = 1,
      color = "white",
      linewidth = 0
    ) +
    
    # Leader lines
    geom_segment(
      aes(
        x = 1.5,
        xend = x_label,
        y = midpoint,
        yend = midpoint
      ),
      linewidth = 0.1,
      inherit.aes = FALSE
    ) +
    
    # Percentage labels outside
    geom_text(
      aes(
        x = x_label + ifelse(hjust == 0, 0.04, -0.04),
        y = midpoint,
        label = label,
        hjust = hjust
      ),
      inherit.aes = FALSE,
      size = pt_to_geom(6)
    ) +
    
    coord_polar(theta = "y") +
    
    scale_fill_manual(
      values = c(
        "Always inactive" = "#BDBDBD",
        "Active (positive)" = "#78c679",
        "Active (negative)" = "#41ab5d",
        "Direction switching" = "#006837"
      )
    ) +
    
    # Make space around the pie
    xlim(0.5, 1.85) +
    
    labs(
      fill = "Reaction class"
    ) +
    
    theme_void() +
    theme(
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 6),
      legend.position = "right"
    ) 
    # theme_tidy()
}
# p_iCHOv1 <- plot_class_pie("iCHOv1")
# p_iCHO2441 <- plot_class_pie("iCHO2441")
p_iCHO3K <- plot_class_pie("iCHO3K")

# p_iCHOv1
# p_iCHO2441
p_iCHO3K


ggsave(
  filename = "figures/models/exp_constrained/reaction_activity_classes_iCHO3K_pie_legend.png",
  width = 45,
  height = 60,
  units = "mm",
  bg = "white",
  dpi = 600
)

# # plot general landscape plots-----------------------


# only filter iCHO3k
# Keep only iCHO3K and exclude inactive class
plot_df <- class_summary %>%
  filter(
    common_name == "iCHO3K",
    direction_consistency != "Always inactive"
  ) %>%
  mutate(
    direction_consistency = recode(
      direction_consistency,
      "Always active (negative only)" = "Active (negative)",
      "Always active (positive only)" = "Active (positive)"
    )
  )
# Plot
ggplot(
  plot_df,
  aes(
    x = direction_consistency,
    y = percent,
    group = direction_consistency
  )
) +
  geom_col(
    fill = "black",
    width = 0.9
  ) +
  geom_text(
    aes(label = paste0(round(percent, 1), "%")),
    vjust = -0.3,
    size = pt_to_geom(6)
  ) +
  scale_y_continuous(
    limits = c(0, 10)
  ) +
  labs(
    x = NULL,
    y = "% of all reactions"
  ) +
  scale_x_discrete(
    expand = c(0, 0)
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    )
  ) +
  theme_tidy()


ggsave(
  filename = "figures/models/exp_constrained/reaction_activity_classes_iCHO3K.png",
  width = 40,
  height = 60,
  units = "mm"
)


# # Stacked barplot (percent)
# ggplot(class_summary, aes(x = common_name, y = n, fill = direction_consistency)) +
#   geom_col(position = "fill") +
#   scale_y_continuous(labels = scales::percent) +
#   labs(x = NULL, y = "% reactions", fill = "Activity pattern") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# # Absolute count barplot
# ggplot(class_summary, aes(x = common_name, y = n, fill = direction_consistency)) +
#   geom_col() +
#   labs(x = NULL, y = "Number of reactions", fill = "Activity pattern") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Example reaction flux profiles for supplementary data -----------------------

# Select one example reaction from each class specifically for iCHO3K
example_rxns <- rxn_classes %>%
  filter(
    common_name == "iCHO3K",
    direction_consistency %in% c(
      "Always active (positive only)",
      "Always active (negative only)",
      "Direction switching"
    )
  ) %>%
  mutate(
    direction_consistency = recode(
      direction_consistency,
      "Always active (negative only)" = "Active (negative)",
      "Always active (positive only)" = "Active (positive)"
    )
  ) %>%
  group_by(direction_consistency) %>%
  slice_head(n = 20) %>%
  ungroup()


# Retrieve the flux profiles for the selected reactions
example_flux <- flux_df %>%
  filter(
    common_name == "iCHO3K_cho_prod_generic_unblocked",
    Reaction %in% example_rxns$Reaction
  ) %>%
  left_join(
    example_rxns %>%
      select(
        Reaction,
        direction_consistency
      ),
    by = "Reaction"
  ) %>%
  mutate(
  #   direction_consistency = recode(
  #     direction_consistency,
  #     "Always active (positive only)" = "Always positive",
  #     "Always active (negative only)" = "Always negative",
  #     "Direction switching" = "Direction switching"
  #   ),
    # direction_consistency = factor(
    #   direction_consistency,
    #   levels = c(
    #     "Always positive",
    #     "Always negative",
    #     "Direction switching"
    #   )
    # ),
    Condition = ifelse(str_detect(Condition, "Constant"), "CT", "TS")
  )
  


plot_data <- example_flux %>%
  filter(Reaction == "NH4t3r") %>% # "ACONTam" for positive, "NH4t3r" for direction switching, "ARTPLM2" for negative
  mutate(
    Window = factor(Window, levels = sort(unique(Window))),
    Condition = factor(
      Condition,
      levels = c("CT", "TS")
    )
  )

ggplot(
  plot_data,
  aes(
    x = Window,
    y = pFBA_Flux,
    group = 1
  )
) +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~ Condition) +
  labs(
    title = unique(plot_data$Reaction),
    subtitle = paste(
      "Direction:",
      unique(plot_data$direction_consistency)
    ),
    x = "Window",
    y = expression("pFBA Flux [mM gDCW"^-1~"h"^-1~"]")
  ) +
  theme_tidy()

ggsave(
  filename = paste0("figures/models/exp_constrained/reaction_activity_classes_", plot_data$Reaction[1],".png"),
  width = 60,
  height = 60,
  units = "mm",
  bg = "white"
)


# filter only active ------------------------------------------------------
# Get always-active reactions (positive OR negative, but consistent)
always_active_rxns <- rxn_classes %>%
  filter(str_detect(direction_consistency, "^Always active")) %>%
  pull(Reaction)

# Filter flux data to only these reactions
active_fluxes <- flux_df %>%
  filter(Reaction %in% always_active_rxns)

# Parse context information (assuming naming pattern)
active_fluxes <- active_fluxes %>%
  mutate(common_name = recode(
    common_name,
    "iCHO3K_cho_prod_generic_unblocked" = "iCHO3K"
  ),
    context_label = paste(Condition, Window, sep = "_")
  )

# Split by model and create individual heatmaps
models <- unique(active_fluxes$common_name)


# heatmap_flux_distr ------------------------------------------------------

# Create a list to store heatmaps
heatmap_list <- list()

# Prepare matrix for heatmap
# Option A: Absolute flux values
for(model in models) {

  # Filter and prepare matrix
  model_fluxes <- active_fluxes %>%
    filter(common_name == model)
  
  flux_matrix <- model_fluxes %>%
    select(Reaction, context_label, pFBA_Flux) %>%
    pivot_wider(id_cols = Reaction, names_from = context_label, values_from = pFBA_Flux) %>%
    column_to_rownames("Reaction") %>%
    as.matrix()
  
  # After creating flux_matrix, add this:
  # Remove rows with zero variance (identical values across all contexts)
  row_var <- apply(flux_matrix, 1, var, na.rm = TRUE)
  flux_matrix <- flux_matrix[row_var > 0 & !is.na(row_var), ]
  
  # Also remove any rows with NA/NaN/Inf
  flux_matrix <- flux_matrix[complete.cases(flux_matrix), ]
  
  # Reorder columns: Constant 1-5, then TS 1-5
  col_order_df <- data.frame(
    col_name = colnames(flux_matrix),
    condition = ifelse(str_detect(colnames(flux_matrix), "Constant"), "Constant", "TS"),
    window = as.numeric(str_extract(colnames(flux_matrix), "\\d+"))
  ) %>%
    arrange(condition, window)
  
  flux_matrix <- flux_matrix[, col_order_df$col_name]
  
  # Column annotation
  col_ha <- HeatmapAnnotation(
    Condition = ifelse(str_detect(colnames(flux_matrix), "Constant"), 
                       "CT", "TS"),
    Window = as.numeric(str_extract(colnames(flux_matrix), "\\d+")),
    col = list(
      Condition = c("CT" = "#E6641E", "TS" = "#4B288C"),
      Window = c("1" = "#fef0d9", "2" = "#fdcc8a", "3" = "#fc8d59", 
               "4" = "#d7301f", "5" = "#7f0000")
    )
  )
  
  # Get number of reactions
  n_reactions <- nrow(flux_matrix)
  
  # Create heatmap with adjusted color scale (use 95th percentile)
  flux_95 <- quantile(abs(flux_matrix[flux_matrix != 0]), 0.95, na.rm = TRUE)
  
  # flux_min <- min(flux_matrix, na.rm = TRUE)
  # flux_max <- max(flux_matrix, na.rm = TRUE)
  
  hm <- Heatmap(flux_matrix,
                name = "Flux",
                col = colorRamp2(
                  c(-flux_95, 0, flux_95),
                  # c(flux_min, 0, flux_max),
                  c("#2166ac", "white", "#b2182b")),
                cluster_rows = TRUE,
                cluster_columns = FALSE,
                show_row_names = FALSE,
                show_column_names = TRUE,
                row_title = paste(model, "-", n_reactions, "active reactions"),
                column_title = "Contexts",
                top_annotation = col_ha,
                heatmap_legend_param = list(title = "Flux", direction = "vertical")
  )
  
  # Save
  pdf(paste0("figures/models/exp_constrained/flux_distr_heatmap_", model, ".pdf"), width = 5, height = 4)
  draw(hm)
  dev.off()
}



# PCA plots ---------------------------------------------------------------
# pca_rxns <- rxn_classes %>%
#   filter(direction_consistency != "Always inactive") %>%
#   select(common_name, Reaction)
# 
# pca_fluxes <- flux_df %>%
#   mutate(
#     common_name = recode(
#       common_name,
#       "iCHO3K_cho_prod_generic_unblocked" = "iCHO3K"
#     )
#   ) %>%
#   inner_join(
#     pca_rxns,
#     by = c("common_name", "Reaction")
#   )
# 
# 
# flux_df_clean <- flux_df
#   # mutate(
#   #   common_name = recode(
#   #     common_name,
#   #     "iCHO3K_cho_prod_generic_unblocked" = "iCHO3K"
#   #   )
#   # )
# 
# # Create a list to store PCA results
# pca_results <- list()

# for(model in c("iCHO3K")) {
#   
#   message("Running PCA for: ", model)
#   
#   model_flux <- pca_fluxes %>%
#     filter(common_name == model)
#   
#   pca_matrix <- model_flux %>%
#     mutate(
#       sample_id = paste(Condition, Window, sep = "_")
#     ) %>%
#     select(
#       sample_id,
#       Condition,
#       Window,
#       Reaction,
#       pFBA_Flux
#     ) %>%
#     pivot_wider(
#       id_cols = c(sample_id, Condition, Window),
#       names_from = Reaction,
#       values_from = pFBA_Flux,
#       values_fill = 0
#     )
#   
#   metadata <- pca_matrix %>%
#     select(sample_id, Condition, Window)
#   
#   flux_mat <- pca_matrix %>%
#     select(-sample_id, -Condition, -Window) %>%
#     as.matrix()
#   
#   reaction_sd <- apply(flux_mat, 2, sd)
#   
#   flux_mat <- flux_mat[
#     ,
#     reaction_sd > 0 & !is.na(reaction_sd),
#     drop = FALSE
#   ]
#   
#   message("Reactions before variance filtering: ",
#           length(reaction_sd))
#   
#   message("Reactions retained for PCA: ",
#           ncol(flux_mat))
#   
#   pca_res <- prcomp(
#     flux_mat,
#     center = TRUE,
#     scale. = TRUE
#   )
#   
#   pca_df <- as.data.frame(pca_res$x) %>%
#     bind_cols(metadata) %>%
#     mutate(
#       Model = model,
#       Window = as.character(Window)
#     )
#   
#   var_exp <- summary(pca_res)$importance[2, ] * 100
#   
#   pca_results[[model]] <- list(
#     pca_df = pca_df,
#     var_exp = var_exp,
#     pca_res = pca_res
#   )
# }
# 
# # Combine all PCA dataframes
# all_pca_df <- bind_rows(lapply(pca_results, function(x) x$pca_df))
# 
# # # Set the order of models
# # all_pca_df <- all_pca_df %>%
# #   mutate(Model = factor(Model, levels = c("iCHOv1", "iCHO2441", "iCHO3K")))
# 
# # Color mapping for conditions
# color_mapping_condition <- c(
#   "Constant" = "#E6641E",
#   "Temp. shifted" = "#4B288C"
# )
# 
# # # Create PCA plot
# p_combined <- ggplot(
#   all_pca_df,
#   aes(x = PC1, y = PC2, color = Condition, shape = Window)
# ) +
#   geom_point(size = 3) +
#   scale_color_manual(values = color_mapping_condition) +
#   facet_wrap(~Model, scales = "free") +
#   labs(
#     x = "PC1",
#     y = "PC2"
#   ) +
#   theme_bw() +
#   theme(
#     legend.position = "bottom",
#     panel.grid.minor = element_blank(),
#     axis.text = element_text(color = "black", size = 11),
#     axis.title = element_text(color = "black", size = 11),
#     strip.background = element_rect(fill = "lightgray"),
#     strip.text = element_text(face = "bold", size = 11)
#   )
# 
# plot(p_combined)
# 
# # Save the plot
# ggsave("figures/models/exp_constrained/pca_combined_models.pdf", p_combined, width = 6.5, height = 3)

# # Also create individual plots for each model with variance in titles
# for(model in c("iCHO3K")) {
#   if(model %in% names(pca_results)) {
#     p_individual <- ggplot(
#       pca_results[[model]]$pca_df,
#       aes(x = PC1, y = PC2, color = Condition, shape = Window)
#     ) +
#       geom_point(size = 3) +
#       scale_color_manual(values = color_mapping_condition) +
#       labs(
#         title = paste0(model, "\n(PC1: ", round(pca_results[[model]]$var_exp[1], 1), "%, PC2: ", round(pca_results[[model]]$var_exp[2], 1), "%)"),
#         x = paste0("PC1 (", round(pca_results[[model]]$var_exp[1], 1), "%)"),
#         y = paste0("PC2 (", round(pca_results[[model]]$var_exp[2], 1), "%)")
#       ) +
#       # theme_bw() +
#       theme_tidy() +
#       theme(
#         legend.position = "bottom",
#         # panel.grid.minor = element_blank(),
#         # axis.text = element_text(color = "black", size = 11),
#         # axis.title = element_text(color = "black", size = 11),
#         # plot.title = element_text(hjust = 0.5, face = "bold", size = 11)
#       ) +
#       guides(
#         color = "none",
#         shape = "none"
#       )
#     
#     # Save individual plot
# 
#     
#     plot(p_individual)
#   }
# }
# 
# ggsave("figures/models/exp_constrained/pca_iCHO3k.png", 
#        p_individual, width = 2, height = 2)


##
# PCA of pFBA fluxes
# Model: iCHO3K
# 1. Prepare reactions for PCA


pca_rxns <- rxn_classes %>%
  filter(direction_consistency != "Always inactive") %>%
  select(common_name, Reaction)

# 2. Prepare flux data


pca_fluxes <- flux_df %>%
  mutate(
    common_name = recode(
      common_name,
      "iCHO3K_cho_prod_generic_unblocked" = "iCHO3K"
    )
  ) %>%
  inner_join(
    pca_rxns,
    by = c("common_name", "Reaction")
  )

# 3. Run PCA

pca_results <- list()

for (model in c("iCHO3K")) {
  
  message("Running PCA for: ", model)
  
  model_flux <- pca_fluxes %>%
    filter(common_name == model)
  
  # Create sample × reaction matrix
  pca_matrix <- model_flux %>%
    mutate(
      sample_id = paste(Condition, Window, sep = "_")
    ) %>%
    select(
      sample_id,
      Condition,
      Window,
      Reaction,
      pFBA_Flux
    ) %>%
    pivot_wider(
      id_cols = c(sample_id, Condition, Window),
      names_from = Reaction,
      values_from = pFBA_Flux,
      values_fill = 0
    )
  
  # Store sample metadata
  metadata <- pca_matrix %>%
    select(
      sample_id,
      Condition,
      Window
    )
  
  # Extract numerical flux matrix
  flux_mat <- pca_matrix %>%
    select(
      -sample_id,
      -Condition,
      -Window
    ) %>%
    as.matrix()
  
  # Remove reactions with zero variance
  reaction_sd <- apply(flux_mat, 2, sd)
  
  keep <- reaction_sd > 0 & !is.na(reaction_sd)
  
  flux_mat <- flux_mat[
    ,
    keep,
    drop = FALSE
  ]
  
  message(
    "Reactions before variance filtering: ",
    length(reaction_sd)
  )
  
  message(
    "Reactions retained for PCA: ",
    ncol(flux_mat)
  )
  
  # Standardize fluxes
  flux_scaled <- scale(flux_mat)
  
  # PCA
  pca_res <- prcomp(
    flux_mat,
    center = TRUE,
    scale. = TRUE
  )
  
  # PCA scores
  pca_df <- as.data.frame(pca_res$x) %>%
    bind_cols(metadata) %>%
    mutate(
      Model = model,
      Window = as.character(Window)
    )
  
  # Variance explained
  var_exp <- summary(pca_res)$importance[2, ] * 100
  
  # Store everything needed for downstream analysis
  pca_results[[model]] <- list(
    pca_df = pca_df,
    var_exp = var_exp,
    pca_res = pca_res,
    flux_mat = flux_mat,
    flux_scaled = flux_scaled,
    metadata = metadata
  )
}


# 4. Combine PCA results


all_pca_df <- bind_rows(
  lapply(
    pca_results,
    function(x) x$pca_df
  )
)


# 5. PCA plot

color_mapping_condition <- c(
  "Constant" = "#E6641E",
  "Temp. shifted" = "#4B288C"
)


p_combined <- ggplot(
  all_pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = Condition,
    shape = Window
  )
) +
  geom_point(size = 3) +
  scale_color_manual(
    values = color_mapping_condition
  ) +
  facet_wrap(
    ~Model,
    scales = "free"
  ) +
  labs(
    x = "PC1",
    y = "PC2"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text = element_text(
      color = "black",
      size = 11
    ),
    axis.title = element_text(
      color = "black",
      size = 11
    ),
    strip.background = element_rect(
      fill = "lightgray"
    ),
    strip.text = element_text(
      face = "bold",
      size = 11
    )
  )

plot(p_combined)


# Save combined PCA plot
ggsave(
  "figures/models/exp_constrained/pca_combined_models.pdf",
  p_combined,
  width = 6.5,
  height = 3
)


# 6. Individual PCA plot

for (model in c("iCHO3K")) {
  
  if (model %in% names(pca_results)) {
    
    p_individual <- ggplot(
      pca_results[[model]]$pca_df,
      aes(
        x = PC1,
        y = PC2,
        color = Condition,
        shape = Window
      )
    ) +
      geom_point(size = 3) +
      scale_color_manual(
        values = color_mapping_condition
      ) +
      labs(
        title = paste0(
          model,
          "\n(PC1: ",
          round(pca_results[[model]]$var_exp[1], 1),
          "%, PC2: ",
          round(pca_results[[model]]$var_exp[2], 1),
          "%)"
        ),
        x = paste0(
          "PC1 (",
          round(pca_results[[model]]$var_exp[1], 1),
          "%)"
        ),
        y = paste0(
          "PC2 (",
          round(pca_results[[model]]$var_exp[2], 1),
          "%)"
        )
      ) +
      theme_tidy() +
      theme(
        legend.position = "bottom"
      ) +
      guides(
        color = "none",
        shape = "none"
      )
    
    plot(p_individual)
  }
}


ggsave(
  "figures/models/exp_constrained/pca_iCHO3k.png",
  p_individual,
  width = 2,
  height = 2
)


# OUTLIER ANALYSIS


# 7. Identify Constant, Window 2

outlier_sample <- all_pca_df %>%
  filter(
    Condition == "Constant",
    Window == "2"
  ) %>%
  select(
    sample_id,
    Model,
    Condition,
    Window,
    PC1,
    PC2
  )

outlier_sample


# 8. Extract exact PCA objects

pca_res <- pca_results[["iCHO3K"]]$pca_res

flux_mat <- pca_results[["iCHO3K"]]$flux_mat

flux_scaled <- pca_results[["iCHO3K"]]$flux_scaled

metadata <- pca_results[["iCHO3K"]]$metadata


# 9. Identify Constant_2

constant2_index <- which(
  metadata$sample_id == "Constant_2"
)



# 10. Calculate reaction-level contributions to Constant_2 PC2


constant2_PC2 <- tibble(
  Reaction = colnames(flux_scaled),
  z_flux = as.numeric(
    flux_scaled[constant2_index, ]
  ),
  PC2_loading = pca_res$rotation[, "PC2"]
) %>%
  mutate(
    PC2_contribution = z_flux * PC2_loading
  )


# 11. Verify that reaction contributions reconstruct the PCA score


constant2_PC2_check <- constant2_PC2 %>%
  summarise(
    PCA_score = pca_res$x[
      constant2_index,
      "PC2"
    ],
    reconstructed_score = sum(
      PC2_contribution
    ),
    negative = sum(
      PC2_contribution[
        PC2_contribution < 0
      ]
    ),
    positive = sum(
      PC2_contribution[
        PC2_contribution > 0
      ]
    ),
    n_negative = sum(
      PC2_contribution < 0
    ),
    n_positive = sum(
      PC2_contribution > 0
    )
  )

constant2_PC2_check



# 12. Inspect the strongest PC2 contributions


top_PC2_contributions <- constant2_PC2 %>%
  mutate(
    abs_contribution = abs(
      PC2_contribution
    )
  ) %>%
  arrange(
    desc(abs_contribution)
  )

top_PC2_contributions %>%
  select(
    Reaction,
    z_flux,
    PC2_loading,
    PC2_contribution
  ) %>%
  slice_head(n = 30)



# 13. Examine how concentrated the PC2 contribution is


PC2_concentration <- constant2_PC2 %>%
  filter(
    PC2_contribution < 0
  ) %>%
  mutate(
    abs_contribution = abs(
      PC2_contribution
    )
  ) %>%
  arrange(
    desc(abs_contribution)
  ) %>%
  mutate(
    cumulative_abs = cumsum(
      abs_contribution
    ),
    cumulative_fraction =
      cumulative_abs /
      sum(abs_contribution)
  )


PC2_concentration %>%
  slice_head(n = 50) %>%
  select(
    Reaction,
    z_flux,
    PC2_loading,
    PC2_contribution,
    cumulative_fraction
  )


# Contribution captured by the top reactions
PC2_concentration %>%
  summarise(
    total_negative = sum(abs_contribution),
    top_10 = sum(
      head(
        abs_contribution,
        10
      )
    ),
    top_25 = sum(
      head(
        abs_contribution,
        25
      )
    ),
    top_50 = sum(
      head(
        abs_contribution,
        50
      )
    ),
    top_100 = sum(
      head(
        abs_contribution,
        100
      )
    )
  )



# CHARACTERIZING THE FLUX DIFFERENCES

# 14. Convert PCA matrix back to long format


pca_flux <- as.data.frame(
  flux_mat
) %>%
  mutate(
    sample_id = metadata$sample_id
  ) %>%
  relocate(
    sample_id
  )



# 15. Compare Constant_2 against all other samples


constant2_vs_rest <- pca_flux %>%
  pivot_longer(
    cols = -sample_id,
    names_to = "Reaction",
    values_to = "Flux"
  ) %>%
  group_by(
    Reaction
  ) %>%
  summarise(
    Constant2_flux =
      Flux[
        sample_id == "Constant_2"
      ][1],
    
    mean_other_flux =
      mean(
        Flux[
          sample_id != "Constant_2"
        ],
        na.rm = TRUE
      ),
    
    max_other_abs_flux =
      max(
        abs(
          Flux[
            sample_id != "Constant_2"
          ]
        ),
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  mutate(
    difference_from_rest =
      Constant2_flux -
      mean_other_flux
  )



# 16. Classify the type of flux difference


constant2_patterns <- constant2_vs_rest %>%
  mutate(
    pattern = case_when(
      
      Constant2_flux != 0 &
        max_other_abs_flux == 0 ~
        "Unique to Constant_2",
      
      Constant2_flux == 0 &
        max_other_abs_flux > 0 ~
        "Absent in Constant_2",
      
      TRUE ~
        "Shared / quantitative difference"
    )
  )


# Number of reactions in each category
constant2_patterns %>%
  count(pattern)



# 17. Largest actual flux differences


constant2_vs_rest %>%
  arrange(
    desc(abs(difference_from_rest))
  ) %>%
  slice_head(n = 30)



# 18. Reactions uniquely active in Constant_2


constant2_patterns %>%
  filter(
    pattern == "Unique to Constant_2"
  ) %>%
  arrange(
    desc(abs(Constant2_flux))
  )



# 19. Reactions absent specifically in Constant_2


constant2_patterns %>%
  filter(
    pattern == "Absent in Constant_2"
  ) %>%
  arrange(
    desc(max_other_abs_flux)
  ) %>%
  slice_head(n = 50)

# 

# # OPTIONAL: SPECIFIC REACTIONS IDENTIFIED DURING EXPLORATION

# 
# 
# # These were previously identified as uniquely active in Constant_2.
# # Keep this section only if you want to explicitly inspect them.
# 
# key_reactions <- c(
#   "HMR_1526",
#   "RE2410C",
#   "r1562"
# )
# 
# key_fluxes <- pca_flux %>%
#   select(
#     sample_id,
#     all_of(key_reactions)
#   )
# 
# key_fluxes
# 



# Re-run PCA excluding Constant_2


# Original PCA objects
pca_res_original <- pca_results[["iCHO3K"]]$pca_res
flux_mat_original <- pca_results[["iCHO3K"]]$flux_mat
metadata_original <- pca_results[["iCHO3K"]]$metadata



# 1. Remove Constant_2


keep_samples <- (
  metadata_original$sample_id != "Constant_2" &
  metadata_original$sample_id != "Temp. shifted_4" # next one which looks like 'outlier'
)

flux_mat_no_C2 <- flux_mat_original[
  keep_samples,
  ,
  drop = FALSE
]

metadata_no_C2 <- metadata_original[
  keep_samples,
  ,
  drop = FALSE
]



# 2. Recalculate reaction variance

# Some reactions may become constant after removing
# Constant_2, so we need to remove them before PCA.

reaction_sd_no_C2 <- apply(
  flux_mat_no_C2,
  2,
  sd
)

keep_reactions_no_C2 <- (
  reaction_sd_no_C2 > 0 &
    !is.na(reaction_sd_no_C2)
)

# 2a. Count reactions removed due to zero variance # 
n_removed_no_C2 <- sum( !keep_reactions_no_C2 ) 
n_retained_no_C2 <- sum( keep_reactions_no_C2 ) 
cat( "Reactions before filtering:", ncol(flux_mat_no_C2), "\n" ) 
cat( "Reactions removed due to zero variance:", n_removed_no_C2, "\n" ) 
cat( "Reactions retained:", n_retained_no_C2, "\n" )

# 3. Filter zero-variance reactions


flux_mat_no_C2_filtered <- flux_mat_no_C2[
  ,
  keep_reactions_no_C2,
  drop = FALSE
]



# 4. Run PCA


pca_res_no_C2 <- prcomp(
  flux_mat_no_C2_filtered,
  center = TRUE,
  scale. = TRUE
)



# 5. Extract PCA scores


pca_df_no_C2 <- as.data.frame(
  pca_res_no_C2$x
) %>%
  bind_cols(
    metadata_no_C2
  ) %>% 
  mutate( Window = as.character(Window) )



# 6. Calculate variance explained


var_exp_no_C2 <- (
  summary(pca_res_no_C2)$importance[2, ]
) * 100



# 7. Check the results


pca_df_no_C2

var_exp_no_C2[1:5]



# 8. Plot PCA without Constant_2


ggplot(
  pca_df_no_C2,
  aes(
    x = PC1,
    y = PC2,
    color = Condition,
    shape = Window,
    label = sample_id
  )
) +
  geom_point(
    size = 4
  ) +
  labs(
    x = paste0(
      "PC1 (",
      round(var_exp_no_C2[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(var_exp_no_C2[2], 1),
      "%)"
    )
  ) +
  # geom_text(
  #   vjust = -1,
  #   size = 3
  # ) +
  scale_color_manual(
    values = color_mapping_condition
  ) +
  labs(
    title = "PCA of pFBA fluxes excluding Constant_2",
    x = paste0(
      "PC1 (",
      round(var_exp_no_C2[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(var_exp_no_C2[2], 1),
      "%)"
    )
  ) +
  theme_tidy() +
  theme(
    legend.position = "bottom"
  ) 
    
ggsave(
  "figures/models/exp_constrained/pca_iCHO3k_wo_CT_w2.png",
  width = 6,
  height = 5,
  bg= "white"
)


# 9. Compare variance explained with original PCA


variance_comparison <- tibble(
  PC = paste0("PC", 1:5),
  Original = pca_results[["iCHO3K"]]$var_exp[1:5],
  Without_Constant_2 = var_exp_no_C2[1:5]
)

variance_comparison


#----Calculate pair-wise Spearmann correlations and plot as heatmap ------------------------

# --- Reuse the same reaction set as your PCA, filtered to iCHO3K
pca_rxns <- rxn_classes %>%
  filter(
    direction_consistency != "Always inactive",
    common_name == "iCHO3K"
  ) %>%
  select(Reaction) %>%
  pull()


# --- Flux matrix (samples x reactions) for iCHO3K only
sample_mat <- flux_df %>%
  filter(common_name == "iCHO3K_cho_prod_generic_unblocked") %>%
  filter(Reaction %in% pca_rxns) %>%
  mutate(
    sample_id = paste(Condition, Window, sep = "_")
  ) %>%
  select(sample_id, Reaction, pFBA_Flux) %>%
  pivot_wider(
    names_from  = Reaction,
    values_from = pFBA_Flux,
    values_fill = 0
  ) %>%
  column_to_rownames("sample_id") %>%
  as.matrix()

cat(
  "samples:", nrow(sample_mat),
  " reactions:", ncol(sample_mat), "\n"
)


# --- Drop all-zero reactions
#     (they only add shared zeros to Spearman)
keep <- colSums(sample_mat != 0) > 0

sample_mat <- sample_mat[
  ,
  keep,
  drop = FALSE
]

cat(
  "kept reactions:", ncol(sample_mat), "\n"
)


# Paired CT/TS Spearman correlation

# 1. Spearman correlation between samples


# --- Spearman correlation between samples ---
cor_mat <- cor(
  t(sample_mat),
  method = "spearman",
  use = "pairwise.complete.obs"
)

# Rename samples in cor_mat: Constant -> CT, Temp. shifted -> TS
new_names <- rownames(cor_mat) %>%
  gsub("^Constant_", "CT_", .) %>%
  gsub("^Temp\\. shifted_", "TS_", .)

rownames(cor_mat) <- new_names
colnames(cor_mat) <- new_names

cat(
  "samples:", nrow(sample_mat),
  " reactions:", ncol(sample_mat), "\n"
)

print(round(cor_mat, 3))



# 2. Metadata aligned to cor_mat


meta <- tibble(sample_id = rownames(cor_mat)) %>%
  separate(
    sample_id,
    into = c("Condition", "Window"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    Condition = factor(Condition, levels = c("CT", "TS")),
    Window = factor(
      Window,
      levels = sort(unique(as.numeric(Window)))
    )
  )



# 3. Order samples:
#
#    CT_1
#    TS_1
#    CT_2
#    TS_2
#    ...


meta <- meta %>%
  arrange(Window, Condition)

sample_order <- meta$sample_id



# 4. Reorder correlation matrix


cor_mat <- cor_mat[
  sample_order,
  sample_order,
  drop = FALSE
]


# Check alignment
stopifnot(
  all(meta$sample_id == rownames(cor_mat)),
  all(meta$sample_id == colnames(cor_mat))
)

print(
  meta %>%
    select(sample_id, Condition, Window)
)



# 5. Fixed colours


condition_col <- c(
  "CT" = "#E6641E",
  "TS" = "#4B288C"
)


# Fixed window colours.
# The colours remain consistent across figures.
window_col <- c(
  "1" = "#D9EAF7",
  "2" = "#C7E9C0",
  "3" = "#FDD0A2",
  "4" = "#DADAEB",
  "5" = "#F2F0F7"
)

# Keep only windows actually present
window_col <- window_col[
  names(window_col) %in% as.character(unique(meta$Window))
]



# 6. Heatmap annotations

ht_opt$simple_anno_size <- unit(1.5, "mm")

top_ha <- HeatmapAnnotation(
  
  Condition = meta$Condition,
  
  Window = factor(
    meta$Window,
    levels = sort(unique(meta$Window))
  ),
  
  col = list(
    Condition = condition_col,
    Window = window_col
  ),
  
  annotation_name_gp = gpar(
    fontsize = 6
  ),
  
  annotation_legend_param = list(
    
    Condition = list(
      title_gp = gpar(fontsize = 6),
      labels_gp = gpar(fontsize = 6)
    ),
    
    Window = list(
      title_gp = gpar(fontsize = 6),
      labels_gp = gpar(fontsize = 6)
    )
  )
)



# 7. Correlation heatmap

ht <- Heatmap(
  
  cor_mat,
  
  name = "Spearman\n\u03c1",
  
  col = colorRamp2(
    c(0.5, 0.75, 1),
    c(
      "white",
      "#F4A582",
      "#B2182B"
    )
  ),
  
  top_annotation = top_ha,
  
  # ----------------------------------------------------------
  # Text
  # ----------------------------------------------------------
  
  show_row_names = TRUE,
  show_column_names = TRUE,
  
  row_names_gp = gpar(
    fontsize = 6
  ),
  
  column_names_gp = gpar(
    fontsize = 6
  ),
  
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 6),
    labels_gp = gpar(fontsize = 6)
  ),
  
  # ----------------------------------------------------------
  # Split heatmap by WINDOW
  # ----------------------------------------------------------
  
  row_split = factor(
    meta$Window,
    levels = sort(unique(meta$Window))
  ),
  
  column_split = factor(
    meta$Window,
    levels = sort(unique(meta$Window))
  ),
  
  # Remove the "1", "2", "3", ... split labels
  row_title = NULL,
  column_title = NULL,
  
  # Preserve CT/TS order
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  border = TRUE
)



# 8. Draw
# --- Save heatmap as PNG ---
png(
  filename = "figures/models/exp_constrained/Spearman_heatmap_iCHO3K.png",
  width = 90,
  height = 60,
  units = "mm",
  res = 600
)

# draw(ht)



draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()

# # reaction-levels dynamics ------------------------------------------------
# For each model, get top reactions by trajectory class

traj_features <- active_fluxes %>%
  
  arrange(common_name, Condition, Window) %>%
  
  group_by(common_name, Reaction, Condition) %>%
  
  summarise(
    
    mean_flux =
      mean(pFBA_Flux),
    
    flux_range =
      max(pFBA_Flux) -
      min(pFBA_Flux),
    
    cumulative_travel =
      sum(abs(diff(pFBA_Flux))),
    
    net_change =
      last(pFBA_Flux) -
      first(pFBA_Flux),
    
    slope =
      coef(
        lm(
          pFBA_Flux ~ Window
        )
      )[2],
    
    .groups = "drop"
    
  )

traj_features <- traj_features %>% mutate(
  travel_ratio = cumulative_travel / (abs(net_change) + 1e-9)
)

# Define thresholds and create trajectory_class
slope_thresh <- 0.01
travel_thresh <- quantile(traj_features$cumulative_travel, 0.25, na.rm = TRUE)

traj_features <- traj_features %>%
  mutate(
    trajectory_class = case_when(
      cumulative_travel < travel_thresh ~ "Stable",
      slope > slope_thresh & travel_ratio < 2 ~ "Increasing",
      slope < -slope_thresh & travel_ratio < 2 ~ "Decreasing",
      travel_ratio >= 2 ~ "Transient",
      TRUE ~ "Other"
    )
  )

model_traj_summary <- list()

# for(model in unique(active_fluxes$common_name)) {
#   
#   # Get trajectory features for this model
#   model_traj <- traj_features %>%
#     filter(common_name == model)
#   
#   # For each trajectory class, get top reactions based on relevant metrics
#   top_reactions <- model_traj %>%
#     filter(!is.na(trajectory_class)) %>%  # Remove NA classes
#     group_by(trajectory_class) %>%
#     mutate(
#       # Define ranking criteria based on class
#       rank_metric = case_when(
#         trajectory_class == "Stable" ~ -cumulative_travel,  # Most stable (lowest travel)
#         trajectory_class == "Increasing" ~ slope,           # Steepest increase
#         trajectory_class == "Decreasing" ~ -slope,          # Steepest decrease  
#         trajectory_class == "Transient" ~ travel_ratio,     # Most transient
#         TRUE ~ cumulative_travel
#       )
#     ) %>%
#     arrange(trajectory_class, desc(rank_metric)) %>%
#     # Take up to 5 per class (or whatever is available)
#     slice_head(n = 5) %>%  # Changed from 10 to 5
#     ungroup() %>%
#     distinct(Reaction, trajectory_class, .keep_all = TRUE)  # Ensure unique reactions
#   
#   model_traj_summary[[model]] <- top_reactions
#   
#   # Get flux data for these top reactions
#   top_flux_data <- active_fluxes %>%
#     filter(common_name == model,
#            Reaction %in% top_reactions$Reaction) %>%
#     left_join(top_reactions %>% select(Reaction, trajectory_class, rank_metric), 
#               by = "Reaction",
#               relationship = "many-to-many")
#   
#   # Print summary for this model
#   cat("\n", paste(rep("=", 60), collapse = ""), "\n")
#   cat("MODEL:", model, "\n")
#   cat(paste(rep("=", 60), collapse = ""), "\n")
#   
#   class_counts <- top_reactions %>%
#     count(trajectory_class)
#   
#   for(i in 1:nrow(class_counts)) {
#     cat("\n", class_counts$trajectory_class[i], ":", 
#         class_counts$n[i], "reactions\n")
#   }
#   
#   # Create plots for this model
#   if(nrow(top_flux_data) > 0) {
#     
#     # Calculate reasonable height (max 20 inches, min 8)
#     n_unique_reactions <- length(unique(top_flux_data$Reaction))
#     plot_height <- min(20, max(8, n_unique_reactions / 2.5))  # Adjusted for fewer reactions
#     
#     # Plot 1: All top reactions faceted
#     p1 <- ggplot(top_flux_data, 
#                  aes(x = as.numeric(Window), 
#                      y = pFBA_Flux, 
#                      color = Condition,
#                      group = Condition)) +
#       geom_line(linewidth = 1) +
#       geom_point(size = 2) +
#       scale_color_manual(values = color_mapping_condition) +
#       facet_grid(Reaction ~ trajectory_class, 
#                  scales = "free_y",
#                  space = "free_y") +
#       labs(
#         title = paste(model, "- Top 5 Reactions per Trajectory Class"),
#         x = "Time Window",
#         y = "Flux Value",
#         color = "Condition"
#       ) +
#       theme_bw() +
#       theme(
#         strip.background.y = element_rect(fill = "lightblue"),
#         strip.text.y = element_text(face = "bold", size = 7),
#         strip.background.x = element_rect(fill = "lightgray"),
#         strip.text.x = element_text(face = "bold", size = 9),
#         axis.text = element_text(color = "black", size = 7),
#         axis.title = element_text(color = "black", size = 10),
#         plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
#         legend.position = "bottom",
#         panel.grid.minor = element_blank()
#       ) +
#       scale_x_continuous(breaks = 1:5)
#     
#     # Save
#     ggsave(paste0("figures/models/exp_constrained/trajectories/top5_trajectories_", model, ".pdf"),
#            p1, width = 14, height = plot_height, limitsize = FALSE)
#     
#     # Plot 2: Separate plot for each trajectory class
#     for(class_type in unique(top_flux_data$trajectory_class)) {
#       
#       class_data <- top_flux_data %>%
#         filter(trajectory_class == class_type)
#       
#       if(nrow(class_data) > 0) {
#         n_reactions_class <- length(unique(class_data$Reaction))
#         class_height <- min(15, max(5, ceiling(n_reactions_class/2) * 2))
#         
#         p2 <- ggplot(class_data, 
#                      aes(x = as.numeric(Window), 
#                          y = pFBA_Flux, 
#                          color = Condition,
#                          group = Condition)) +
#           geom_line(linewidth = 1) +
#           geom_point(size = 2.5) +
#           scale_color_manual(values = color_mapping_condition) +
#           facet_wrap(~Reaction, 
#                      scales = "free_y",
#                      ncol = 2) +
#           labs(
#             title = paste(model, "-", class_type, "Reactions"),
#             subtitle = paste("Top", n_reactions_class, "reactions"),
#             x = "Time Window",
#             y = "Flux Value",
#             color = "Condition"
#           ) +
#           theme_bw() +
#           theme(
#             strip.background = element_rect(fill = "lightgray"),
#             strip.text = element_text(face = "bold", size = 8),
#             axis.text = element_text(color = "black", size = 8),
#             axis.title = element_text(color = "black", size = 10),
#             plot.title = element_text(hjust = 0.5, face = "bold"),
#             plot.subtitle = element_text(hjust = 0.5),
#             legend.position = "bottom"
#           ) +
#           scale_x_continuous(breaks = 1:5)
#         
#         # Save
#         ggsave(paste0("figures/models/exp_constrained/trajectories/top5_", class_type, "_", model, ".pdf"),
#                p2, width = 10, height = class_height, limitsize = FALSE)
#         
#         print(p2)
#       }
#     }
#     
#     print(p1)
#   } else {
#     cat("\nNo top reactions found for model:", model, "\n")
#   }
# }

# # Create a combined summary table (using top 5 from your loop)
# combined_summary <- bind_rows(model_traj_summary) %>%
#   select(common_name, Reaction, trajectory_class, 
#          slope, cumulative_travel, travel_ratio, net_change, flux_range, rank_metric) %>%
#   arrange(common_name, trajectory_class, desc(rank_metric))
# 
# # Print the combined summary
# cat("COMBINED SUMMARY: Top 5 Reactions Across All Models\n")
# print(combined_summary)

# Save the summary table
#write_csv(combined_summary, "figures/top5_trajectories_summary.csv")

# # Optional: Create slope bar plots for the top reactions
# for(model in names(model_traj_summary)) {
#   
#   top_reactions <- model_traj_summary[[model]]
#   
#   if(nrow(top_reactions) > 0) {
#     # Create unique labels for reactions by adding trajectory class
#     top_reactions_for_plot <- top_reactions %>%
#       mutate(Reaction_label = paste0(Reaction, " (", trajectory_class, ")"))
#     
#     # Simple bar plot of slopes by class
#     p_slopes <- ggplot(top_reactions_for_plot, 
#                        aes(x = reorder(Reaction_label, slope), 
#                            y = slope, 
#                            fill = trajectory_class)) +
#       geom_col() +
#       coord_flip() +
#       labs(
#         title = paste(model, "- Slope Values for Top 5 Reactions"),
#         x = "Reaction",
#         y = "Slope (flux units per window)", #Temporal trend (Δflux/window) OR Linear trend (flux/window)
#         fill = "Trajectory Class",
#         subtitle = "Positive = increasing over time; Negative = decreasing over time; Calculated separately for Constant (orange) and TS (purple)"
#       ) +
#       theme_bw() +
#       theme(
#         axis.text.y = element_text(size = 6),
#         legend.position = "bottom"
#       ) +
#       geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")
#     
#     # Create directory if it doesn't exist
#     if(!dir.exists("figures/models/exp_constrained")) {
#       dir.create("figures/models/exp_constrained", recursive = TRUE)
#     }
#     
#     ggsave(paste0("figures/models/exp_constrained/top5_slopes_", model, ".pdf"),
#            p_slopes, width = 10, height = 8, limitsize = FALSE)
#   }
# }

# ADD NEW CLASSIFICATIONS FOR STABLE REACTIONS --------------------
# Calculate stability metrics from the original flux data
stability_classification <- active_fluxes %>%
  group_by(common_name, Reaction, Condition) %>%
  summarise(
    # Check if flux is constant across windows within each condition
    flux_range_within_condition = max(pFBA_Flux) - min(pFBA_Flux),
    is_constant_within_condition = flux_range_within_condition < 1e-6,  # Threshold for constant
    .groups = "drop"
  ) %>%
  group_by(common_name, Reaction) %>%
  summarise(
    # Check if constant across conditions (comparing Constant vs TS at same window)
    # AND constant across windows
    constant_across_windows_and_conditions = all(is_constant_within_condition),
    constant_across_windows_only = any(is_constant_within_condition) & !all(is_constant_within_condition),
    .groups = "drop"
  )

# Join back to traj_features
traj_features <- traj_features %>%
  left_join(stability_classification, by = c("common_name", "Reaction"))

# Update trajectory_class with new categories
traj_features <- traj_features %>%
  mutate(
    trajectory_class = case_when(
      # New: Constant across everything
      constant_across_windows_and_conditions ~ "Stable (no change)",
      # New: Constant within each condition but differs between conditions
      constant_across_windows_only ~ "Stable within conditions",
      # Original classifications
      cumulative_travel < travel_thresh ~ "Stable",
      slope > slope_thresh & travel_ratio < 2 ~ "Increasing",
      slope < -slope_thresh & travel_ratio < 2 ~ "Decreasing",
      travel_ratio >= 2 ~ "Transient",
      TRUE ~ "Other"
    )
  )

# Add a separate stability type column
traj_features <- traj_features %>%
  mutate(
    stability_type = case_when(
      constant_across_windows_and_conditions ~ "Stable (no change across windows or conditions)",
      constant_across_windows_only ~ "Stable within conditions (differs between CT/TS)",
      TRUE ~ "Variable"
    )
  )


traj_features %>%
  count(common_name, trajectory_class) %>%
  arrange(common_name, trajectory_class) %>%
  print()
# 
# # Create summary plots for the new classifications
# p_stability <- ggplot(traj_features, 
#                       aes(x = common_name, fill = stability_type)) +
#   geom_bar(position = "fill") +
#   scale_fill_manual(
#     values = c(
#       "Stable (no change across windows or conditions)" = "#A6D854",  # pastel green
#       "Stable within conditions (differs between CT/TS)" = "#8DA0CB", # pastel blue
#       "Variable" = "#FC8D62"                                          # pastel orange
#     )
#   ) +
#   scale_y_continuous(labels = scales::percent) +
#   labs(
#     title = "Reaction Stability Patterns",
#     x = NULL,
#     y = "% of active reactions",
#     fill = "Stability Pattern"
#   ) +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.position = "bottom",
#     legend.direction = "vertical"
#   )
# 
# plot(p_stability)
# 
# ggsave("figures/models/exp_constrained/reaction_stability_patterns.pdf", p_stability, width = 8, height = 6)

traj_class_summary <- traj_features %>%
  mutate(
    # common_name = recode(
    #   common_name,
    #   "iCHO3K_cho_prod_generic_unblocked" = "iCHO3K"
    # ),
    common_name = factor(
      common_name,
      levels = c(
        "iCHOv1",
        "iCHO2441",
        "iCHO3K"
      )
    ),
    trajectory_class = factor(
      trajectory_class,
      levels = c(
        "Stable (no change)",
        "Stable within conditions",
        "Increasing",
        "Decreasing",
        "Transient",
        "Other"
      )
    )
  ) %>%
  count(common_name, trajectory_class)

# p_traj_classes <- ggplot(
#   traj_class_summary,
#   aes(
#     x = common_name,
#     y = n,
#     fill = trajectory_class
#   )
# ) +
#   geom_col(position = "fill") +
#   scale_y_continuous(
#     labels = scales::percent
#   ) +
#   scale_fill_manual(
#     values = c(
#       "Stable (no change)" = "#009E73",
#       "Stable within conditions" = "#56B4E9",
#       "Increasing" = "#E69F00",
#       "Decreasing" = "#D55E00",
#       "Transient" = "#CC79A7",
#       "Other" = "grey80"
#     )
#   ) +
#   labs(
#     x = NULL,
#     y = "% reactions",
#     fill = "Trajectory class"
#   ) +
#   theme_bw() +
#   theme(
#     axis.text = element_text(
#       color = "black"
#     ),
#     axis.title = element_text(
#       color = "black"
#     ),
#     legend.position = "right"
#   )
#  plot(p_traj_classes)
# ggsave(
#   "figures/models/exp_constrained/trajectory_classes.pdf",
#   p_traj_classes,
#   width = 6,
#   height = 4
# )

#only iCHO3k

traj_class_summary_iCHO3K <- traj_class_summary %>%
  filter(common_name == "iCHO3K")

p_traj_classes <- ggplot(
  traj_class_summary_iCHO3K,
  aes(
    x = trajectory_class,
    y = n
  )
) +
  geom_col(fill = "black", width = 0.9) +
  geom_text(
    aes(label = n),
    vjust = -0.3,
    size = pt_to_geom(6)
  ) +
  labs(
    x = NULL,
    y = "Number of reactions"
  ) +
  scale_y_continuous(limits = c(0, 1200)) +
  
  theme_bw() +
  theme_tidy() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    # axis.text = element_text(color = "black"),
    # axis.title = element_text(color = "black"),
    # legend.position = "none"
  )

plot(p_traj_classes)

ggsave(
  filename = "figures/models/exp_constrained/trajectory_classes_iCHO3K.png",
  width = 80,
  height = 60,
  units = "mm"
)

# 27.08.2026-------------


# REACTION-LEVEL TEMPORAL DYNAMICS

#
# Input:
#   active_fluxes
#
# active_fluxes should contain the same reaction universe used
# for the PCA:
# reactions active in at least one experimental context.
#
# Trajectories are calculated separately for each:
#   model × reaction × condition
#
# Classes:
#   Stable
#   Increasing
#   Decreasing
#   Transient
#   Other
#


# 1. Calculate trajectory features


traj_features <- pca_fluxes %>%
  
  arrange(
    common_name,
    Reaction,
    Condition,
    Window
  ) %>%
  
  group_by(
    common_name,
    Reaction,
    Condition
  ) %>%
  
  summarise(
    
    # Mean flux across time windows
    mean_flux =
      mean(pFBA_Flux, na.rm = TRUE),
    
    # Total range across time windows
    flux_range =
      max(pFBA_Flux, na.rm = TRUE) -
      min(pFBA_Flux, na.rm = TRUE),
    
    # Total amount of flux movement across consecutive windows
    cumulative_travel =
      sum(abs(diff(pFBA_Flux))),
    
    # Net change from first to last window
    net_change =
      last(pFBA_Flux) -
      first(pFBA_Flux),
    
    # Linear temporal trend
    slope =
      coef(
        lm(
          pFBA_Flux ~ Window
        )
      )[2],
    
    .groups = "drop"
  )



# 2. Calculate travel ratio

#
# Approximately 1:
#   mostly monotonic trajectory
#
# > 1:
#   some reversals in direction
#
# >= 2:
#   strongly non-monotonic / transient behavior
#


traj_features <- traj_features %>%
  mutate(
    
    travel_ratio =
      cumulative_travel /
      (abs(net_change) + 1e-9)
    
  )



# 3. Define trajectory thresholds


# Flux range below this value is considered effectively constant.
# Adjust if necessary according to the numerical precision/scale
# of your flux data.

flux_stability_thresh <- 1e-6


# Minimum absolute temporal slope required to classify a reaction
# as increasing or decreasing.

slope_thresh <- 0.01



# 4. Classify trajectories


traj_features <- traj_features %>%
  mutate(
    
    trajectory_class = case_when(
      
      # No meaningful change across windows
      flux_range <= flux_stability_thresh ~
        "Stable",
      
      # Predominantly monotonic increase
      slope > slope_thresh &
        travel_ratio < 2 ~
        "Increasing",
      
      # Predominantly monotonic decrease
      slope < -slope_thresh &
        travel_ratio < 2 ~
        "Decreasing",
      
      # Strongly non-monotonic trajectory
      travel_ratio >= 2 ~
        "Transient",
      
      # Remaining cases with small/ambiguous trends
      TRUE ~
        "Other"
    )
  )



# 5. Check trajectory counts


traj_features %>%
  count(
    common_name,
    Condition,
    trajectory_class
  ) %>%
  arrange(
    common_name,
    Condition,
    trajectory_class
  ) %>%
  print()



# 6. Summarise trajectory classes by model and condition


traj_class_summary <- traj_features %>%
  
  mutate(
    
    common_name = factor(
      common_name,
      levels = c(
        "iCHOv1",
        "iCHO2441",
        "iCHO3K"
      )
    ),
    
    trajectory_class = factor(
      trajectory_class,
      levels = c(
        "Stable",
        "Increasing",
        "Decreasing",
        "Transient",
        "Other"
      )
    )
    
  ) %>%
  
  count(
    common_name,
    Condition,
    trajectory_class
  )



# 7. iCHO3K trajectory summary


traj_class_summary_iCHO3K <- traj_class_summary %>%
  filter(
    common_name == "iCHO3K"
  )



# 8. Plot trajectory classes separately for CT and TS

# Bar plot settings
PLOT_FONT_SIZE  <- 6 # was 6 for paper figures
PANEL_LABEL_SIZE <- 6 # was 8 for paper figures
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
p_traj_classes <- ggplot(
  traj_class_summary_iCHO3K,
  aes(
    x = trajectory_class,
    y = n,
    fill = Condition,
    color = Condition
  )
) +
  
  geom_col(
    #position = "dodge",
    #width = 0.8,
    linewidth = LINE_WIDTH,
    alpha = BAR_ALPHA,
    width = 0.85,
    position = position_dodge(width = DODGE_WIDTH),
    na.rm = TRUE
    
  ) +
  
  geom_text(
    aes(
      label = n
    ),
    position = position_dodge(width = DODGE_WIDTH),
    vjust = -0.3,
    size = pt_to_geom(6)
  ) +
  scale_fill_manual(values = color_condition_long) +
  scale_color_manual(values = color_condition_long) +  
  labs(
    x = NULL,
    y = "Number of reactions",
    fill = "Condition"
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(0, 0.08)
    )
  ) +
  
  theme_bw(
    base_size = PLOT_FONT_SIZE,
    base_family = FONT_FAMILY
  ) +  
  theme_tidy() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    )
  )


plot(p_traj_classes)



# Save trajectory-class plot


ggsave(
  "figures/models/exp_constrained/trajectory_classes_iCHO3K_new.png",
  p_traj_classes,
  width = 80,
  height = 60,
  units = "mm",
  limitsize = FALSE
)



# 9. Compare trajectory behavior between CT and TS

#
# For each reaction, determine whether its trajectory class
# is the same or different between conditions.
#
# Example:
#
# Reaction A:
#   CT = Increasing
#   TS = Increasing
#   -> Same
#
# Reaction B:
#   CT = Stable
#   TS = Increasing
#   -> Different
#


trajectory_comparison <- traj_features %>%
  
  select(
    common_name,
    Reaction,
    Condition,
    trajectory_class
  ) %>%
  
  pivot_wider(
    names_from = Condition,
    values_from = trajectory_class
  )



# 10. Identify whether CT and TS have the same trajectory


# This assumes the condition names are "CT" and "TS".

trajectory_comparison <- trajectory_comparison %>%
  mutate(
    
    trajectory_response = case_when(
      
      Constant == `Temp. shifted` ~
        "Same trajectory",
      
      TRUE ~
        "Different trajectory"
    )
  )



# 11. Count same vs different trajectories


trajectory_comparison %>%
  count(
    common_name,
    trajectory_response
  ) %>%
  arrange(
    common_name,
    trajectory_response
  ) %>%
  print()



# 12. Detailed CT vs TS trajectory combinations

#
# This tells us exactly how reactions change their trajectory
# behavior between conditions.
#
# Example:
#
# Stable      -> Increasing
# Increasing  -> Transient
# Decreasing  -> Stable
#


trajectory_comparison %>%
  count(
    common_name,
    Constant,
    `Temp. shifted`
  ) %>%
  arrange(
    common_name,
    desc(n)
  ) %>%
  print()



# 13. iCHO3K: same vs different trajectory


trajectory_comparison_iCHO3K <- trajectory_comparison %>%
  filter(
    common_name == "iCHO3K"
  )


trajectory_comparison_iCHO3K %>%
  count(
    trajectory_response
  ) %>%
  print()



# 14. Plot same vs different trajectories

trajectory_response_summary <- trajectory_comparison_iCHO3K %>%
  count(
    trajectory_response
  )


p_trajectory_response <- ggplot(
  trajectory_response_summary,
  aes(
    x = trajectory_response,
    y = n
  )
) +
  
  geom_col(
    fill = "black",
    width = 0.8
  ) +
  
  geom_text(
    aes(
      label = n
    ),
    vjust = -0.3,
    size = pt_to_geom(6)
  ) +
  
  labs(
    x = NULL,
    y = "Number of reactions"
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(0, 0.08)
    )
  ) +
  
  theme_bw() +
  
  theme_tidy() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    )
  )


plot(p_trajectory_response)



# Save


ggsave(
  "figures/models/exp_constrained/trajectory_response_CT_vs_TS_iCHO3K.pdf",
  p_trajectory_response,
  width = 7,
  height = 5,
  limitsize = FALSE
)



# 15. Save detailed trajectory tables


write_csv(
  traj_features,
  "figures/models/exp_constrained/trajectory_features.csv"
)


write_csv(
  trajectory_comparison,
  "figures/models/exp_constrained/trajectory_comparison_CT_vs_TS.csv"
)


# # Create a detailed summary of constant reactions
# constant_reactions_summary <- traj_features %>%
#   filter(stability_type != "Variable") %>%
#   select(common_name, Reaction, stability_type, trajectory_class,
#          slope, cumulative_travel, travel_ratio, 
#          mean_flux, flux_range) %>%
#   arrange(common_name, stability_type, Reaction)
# 
# cat("STABLE REACTION SUMMARY\n")
# print(constant_reactions_summary)
# 
# # Save to CSV
# #write_csv(constant_reactions_summary, "figures/constant_reactions_summary.csv")
# 
# # Visualize constant reactions across models
# p_constant_heatmap <- traj_features %>%
#   filter(stability_type != "Variable") %>%
#   group_by(common_name, stability_type) %>%
#   summarise(count = n(), .groups = "drop") %>%
#   ggplot(aes(x = common_name, y = stability_type, fill = count)) +
#   geom_tile() +
#   geom_text(aes(label = count), color = "white", fontface = "bold") +
#   scale_fill_gradient(low = "lightblue", high = "darkblue") +
#   labs(
#     title = "Number of Stable Reactions per Model",
#     x = NULL,
#     y = NULL,
#     fill = "Count"
#   ) +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     panel.grid = element_blank()
#   )
# plot(p_constant_heatmap)
# ggsave("figures/models/exp_constrained/constant_reactions_heatmap.pdf", p_constant_heatmap, width = 6, height = 4)
# 
# # condition differences delta flux ----------------------------------------
# delta_flux <- active_fluxes %>%
#   select(
#     common_name,
#     Reaction,
#     Window,
#     Condition,
#     pFBA_Flux
#   ) %>%
#   pivot_wider(
#     names_from = Condition,
#     values_from = pFBA_Flux
#   ) %>%
#   mutate(
#     delta_flux =
#       `Temp. shifted` - Constant,
#     
#     relative_delta =
#       delta_flux /
#       pmax(
#         abs(`Temp. shifted`),
#         abs(Constant),
#         1e-9
#       )
#   )
# 
# delta_summary <- delta_flux %>%
#   group_by(
#     common_name,
#     Reaction
#   ) %>%
#   summarise(
#     
#     mean_delta =
#       mean(delta_flux),
#     
#     max_delta =
#       max(abs(delta_flux)),
#     
#     cumulative_delta =
#       sum(abs(delta_flux)),
#     
#     mean_relative_delta =
#       mean(abs(relative_delta)),
#     
#     .groups = "drop"
#     
#   )
# 
# # Combine cumulative delta and mean delta
# delta_analysis <- delta_summary %>%
#   mutate(
#     # Classify by condition bias
#     condition_bias = case_when(
#       mean_delta > 0 ~ "Higher in TS",
#       mean_delta < 0 ~ "Higher in CT",
#       TRUE ~ "Similar"
#     ),
#     # Classify by effect magnitude (using cumulative delta)
#     effect_magnitude = case_when(
#       cumulative_delta > quantile(cumulative_delta, 0.75) ~ "Large",
#       cumulative_delta > quantile(cumulative_delta, 0.5) ~ "Medium",
#       TRUE ~ "Small"
#     ),
#     # Combined classification
#     combined_class = paste(effect_magnitude, condition_bias)
#   ) %>%
#   mutate(
#     common_name = factor(
#       common_name,
#       levels = c(
#         "iCHOv1",
#         "iCHO2441",
#         "iCHO3K"
#       ))
#   )
# 
# # Plot 1: Scatter plot - cumulative delta vs mean delta
# p_scatter <- ggplot(delta_analysis, 
#                     aes(x = mean_delta, 
#                         y = cumulative_delta, 
#                         color = condition_bias
#                         # alpha = effect_magnitude
#                         )) +
#   geom_point(size = 2) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
#   scale_color_manual(values = c("Higher in TS" = "#4B288C", 
#                                 "Higher in CT" = "#E6641E", 
#                                 "Similar" = "gray")) +
#   # scale_alpha_manual(values = c("Large" = 1, "Medium" = 0.7, "Small" = 0.3)) +
#   facet_wrap(~common_name, scales = "free") +
#   labs(
#     title = "Condition Effects: Magnitude vs Direction",
#     x = "Mean Delta (TS - CT) → Positive = Higher in TS",
#     y = "Cumulative |Delta| (Total Condition Effect)",
#     color = "Direction"
#     # alpha = "Magnitude"
#   ) +
#   theme_bw() +
#   theme(
#     legend.position = "bottom",
#     panel.grid.minor = element_blank(),
#     strip.background = element_rect(fill = "lightgray"),
#     strip.text = element_text(face = "bold")
#   )
# 
# print(p_scatter)
# ggsave("figures/models/exp_constrained/delta_scatter_cumulative_vs_mean.pdf", 
#        p_scatter, 
#        width = 6, 
#        height = 4)
# 
# 
# top20_by_cumulative <- delta_analysis %>%
#   group_by(common_name) %>%
#   slice_max(cumulative_delta, n = 30) %>%
#   ungroup() %>%
#   # reorder_within creates factor that orders within each group
#   mutate(Reaction = reorder_within(Reaction, cumulative_delta, common_name))
# 
# p_top20 <- ggplot(top20_by_cumulative, 
#                   aes(x = Reaction, 
#                       y = cumulative_delta, 
#                       fill = condition_bias)) +
#   geom_col() +
#   coord_flip() +
#   facet_wrap(~common_name, scales = "free_y") +
#   scale_x_reordered() +  # This cleans up the axis labels
#   scale_fill_manual(values = c("Higher in TS" = "#4B288C", 
#                                "Higher in CT" = "#E6641E", 
#                                "Similar" = "gray")) +
#   labs(
#     title = "Top 30 Reactions by Cumulative Condition Effect",
#     x = "Reaction",
#     y = "Cumulative |Delta| (sum of |TS - CT| across windows)",
#     fill = "Direction"
#   ) +
#   theme_bw() +
#   theme(
#     legend.position = "bottom",
#     axis.text.y = element_text(size = 7),
#     strip.background = element_rect(fill = "lightgray"),
#     strip.text = element_text(face = "bold")
#   )
# 
# print(p_top20)
# ggsave("figures/models/exp_constrained/delta_top30_cumulative.pdf", p_top20, width = 6, height = 5)
# 
# # # Plot 3: Binned summary - how many reactions fall into each quadrant?
# # delta_analysis %>%
# #   mutate(
# #     direction_binned = case_when(
# #       mean_delta > 1 ~ "Strongly higher in TS",
# #       mean_delta > 0 ~ "Slightly higher in TS",
# #       mean_delta > -1 ~ "Similar",
# #       mean_delta > -5 ~ "Slightly higher in CT",
# #       TRUE ~ "Strongly higher in CT"
# #     ),
# #     magnitude_binned = case_when(
# #       cumulative_delta > quantile(cumulative_delta, 0.9) ~ "Top 10%",
# #       cumulative_delta > quantile(cumulative_delta, 0.75) ~ "Top 25%",
# #       cumulative_delta > quantile(cumulative_delta, 0.5) ~ "Top 50%",
# #       TRUE ~ "Bottom 50%"
# #     )
# #   ) %>%
# #   count(common_name, direction_binned, magnitude_binned) %>%
# #   ggplot(aes(x = direction_binned, y = n, fill = magnitude_binned)) +
# #   geom_col(position = "dodge") +
# #   facet_wrap(~common_name) +
# #   coord_flip() +
# #   labs(
# #     title = "Distribution of Condition Effects",
# #     x = "Direction (TS vs CT)",
# #     y = "Number of Reactions",
# #     fill = "Effect Magnitude"
# #   ) +
# #   theme_bw() +
# #   theme(legend.position = "bottom") +
# #   scale_fill_brewer(palette = "Blues")
# # 
# # ggsave("figures/delta_distribution_summary.pdf", width = 10, height = 6)
# 
# # Plot 4: Heatmap of top reactions by BOTH cumulative and mean delta
# # Select reactions that are high in cumulative delta AND have strong direction
# top_strong_effects <- delta_analysis %>%
#   group_by(common_name) %>%
#   filter(
#     cumulative_delta > quantile(cumulative_delta, 0.75),
#     abs(mean_delta) > quantile(abs(mean_delta), 0.75)
#   ) %>%
#   slice_max(abs(mean_delta), n = 15)
# 
# # Get the delta_flux data for these reactions
# heatmap_strong <- delta_flux %>%
#   semi_join(top_strong_effects, by = c("common_name", "Reaction")) %>%
#   mutate(Reaction_label = paste0(Reaction, " (", 
#                                  ifelse(mean(delta_flux) > 0, "↑TS", "↑CT"), 
#                                  ")"))
# 
# ggplot(heatmap_strong, 
#        aes(x = Window, y = reorder(Reaction_label, delta_flux), fill = delta_flux)) +
#   geom_tile() +
#   facet_grid(common_name ~ ., scales = "free_y", space = "free_y") +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
#                        midpoint = 0,
#                        name = "ΔFlux\n(TS - CT)") +
#   labs(
#     title = "Reactions with Strong Condition Effects (High Magnitude & Clear Direction)",
#     x = "Time Window",
#     y = "Reaction"
#   ) +
#   theme_bw() +
#   theme(
#     axis.text.y = element_text(size = 7),
#     strip.background = element_rect(fill = "lightgray"),
#     strip.text = element_text(face = "bold")
#   )
# 
# ggsave("figures/delta_strong_effects_heatmap.pdf", width = 10, height = 12)
# 
# # Print summary statistics
# cat("CONDITION EFFECT SUMMARY\n")
# 
# delta_analysis %>%
#   group_by(common_name) %>%
#   summarise(
#     total_reactions = n(),
#     higher_in_TS = sum(mean_delta > 0),
#     higher_in_CT = sum(mean_delta < 0),
#     pct_higher_TS = round(100 * higher_in_TS / total_reactions, 1),
#     pct_higher_CT = round(100 * higher_in_CT / total_reactions, 1),
#     mean_cumulative_delta = round(mean(cumulative_delta), 2),
#     max_cumulative_delta = round(max(cumulative_delta), 2),
#     .groups = "drop"
#   ) %>%
#   print()
