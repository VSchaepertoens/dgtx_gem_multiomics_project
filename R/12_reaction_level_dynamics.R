library(tidyverse)
source("R/utils/plotting_utils.R", local = FALSE)
library(ComplexHeatmap)
library(circlize)  # for color mapping functions

type <- "exp_constrained" # or "GIMME_tpm_25perc" or "exp_constrained"

flux_df <- read_csv(paste0("analysis/models/",type,"/reaction_fluxes_FBA_pFBA.csv")) %>%
  mutate(common_name = recode(
    common_name,
    "iCHO3K_cho_prod_generic_unblocked" = "iCHO3K"
  )) %>%
  filter(common_name == "iCHO3K")


reactions <- c("HEX1", "PGI", "TKT2", "TALA", "TKT1", "RPE",
               "RPI", "GND", "PGL", "G6PDH2r",
               "PFK", "FBA", "TPI", "GAPD",
               "PGK",  "PGM", "ENO", "r0122", "PYK", "LDH_D","LDH_L",
               "ALATA_L", "PYRt2m", "PYRtm", "PCm", "CSm", "ACONTam", "ACONTbm", "ICDHyrm",
               "AKGDm","SUCOASm", "SUCOAS1m", "SUCD1m", "FUMm", "MDHm",
               "OCOAT1m", "OAADC", "LDH_Lm", "L_LACtcm", "L_LACt4r")


# reactions <- c("3SPYRSPm", "CITMCOALm", "DDPGAm", "LDH_Lm", "AGTim", "D3AIBTm", "ME1m", "ME2m", "MERCSULFTm",
#                "MESCOALm", "r0081", "r0157", "r1464", "SPAm")
# 
# reactions <- c("L_LACtcm", "L_LACtm", "LCADim")
# 
# reactions <- c("15KPROSTGF2t", "L_LACt2r", "L_LACt4r", "L_LACt4r_2", "L_LACtm", "LACLt", "LCADi", "OROTi",
#                "r2083", "r2090", "r2097", "r2105", "r2112", "r2117", "r2121", "r2127", "r2128", "r2129", "r2516")
# reactions <- c("igg_formation", "igg_hc", "igg_lc")
# 
# 
# reactions <- c("F6Tg", "F6T3g","F6T2g", "M1316M2g", "M16NT1g", "G14Tg",
#                "G14Tg_1", "G14Tg_new2", "G14Tg_new3", "M1316Mg",
#                "FUCtly", "FK", "F1PGT", "GDPFUCtg", "GFUCS", "HMR_0873",
#                "GMAND", "MAN1PT2", "r0208", "GALK", "UGLT", "UDPGALtg",
#                "UDPG4E", "GALU", "PGMT","PGI", "MAN6PI", "PMANM", "GALT")

# reactions <- c("F6Tg", "F6T3g","F6T2g", "M1316M2g", "M16NT1g", "G14Tg",
#                "G14Tg_1", "G14Tg_new2", "G14Tg_new3", "M1316Mg",
#                "FUCtly", "FK", "F1PGT", "GDPFUCtg", "GFUCS", "HMR_0873",
#                "GMAND", "MAN1PT2", "r0208", "GALK", "UGLT", "UDPGALtg",
#                "UDPG4E", "GALU", "PGMT","PGI", "MAN6PI", "PMANM", "GALT")
# reactions <- c("3SPYRSP", "ACNMLr", "ALR", "GapFill-R04941", "HMR_3859", "L_LACDcm", "MCPST", "ME2", "PGAM1_PhosHydro", "PYRSMCT1", "PYRt2p", "PYRt4_2",
#                "r0060", "r0122", "r0153", "r0160", "r0165", "r0280", "r0413", "r1464", "r2079", "r2086", "r2093", "r2094", "r2095", "r2096", "r2097",
#                "r2098", "r2099", "r2100", "RE2954C")

color_mapping_condition <- c(
  "Constant" = "#E6641E",
  "Temp. shifted" = "#4B288C"
)

flux_plot_df <- flux_df %>%
  filter(Reaction %in% reactions) %>%
  # mutate(
  #   Window = factor(Window, levels = 1:5),
  #   Reaction = factor(Reaction, levels = reactions)
  # ) %>%
  mutate(
    # Reverse pFBA flux sign for selected reactions
    pFBA_Flux = if_else(
      Reaction %in% c("PGK", "PGM", "ACONTbm", "RPI"),
      -pFBA_Flux,
      pFBA_Flux
    ),
    Window = factor(Window, levels = 1:5),
    Reaction = factor(Reaction, levels = reactions)
  )

# Create succoa -> succ for each Condition + Window
new_reaction <- flux_plot_df %>%
  filter(Reaction %in% c("SUCOAS1m", "OCOAT1m")) %>%
  group_by(model_name, common_name, Condition, Window) %>%
  summarise(
    FBA_LB    = -(sum(FBA_LB)),
    FBA_UB    = -(sum(FBA_UB)),
    FBA_Flux  = -(sum(FBA_Flux)),
    pFBA_LB   = -(sum(pFBA_LB)),
    pFBA_UB   = -(sum(pFBA_UB)),
    pFBA_Flux = -(sum(pFBA_Flux)),
    Reaction  = "succoa -> succ",
    .groups = "drop"
  )

# Remove the original reactions and add the new reaction
flux_plot_df <- flux_plot_df %>%
  filter(!Reaction %in% c("SUCOAS1m", "OCOAT1m")) %>%
  bind_rows(new_reaction)

reactions <- c(reactions, "succoa -> succ")

reaction_plots <- lapply(reactions, function(rxn) {
  
  ggplot(
    filter(flux_plot_df, Reaction == rxn),
    aes(
      x = Window,
      y = pFBA_Flux,
      colour = Condition,
      group = Condition
    )
  ) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 0.5) +
    labs(
      title = rxn,
      x = "",
      y = ""
      # x = "Window",
      # y = expression(paste("pFBA Flux [mM ", gDCW^{-1}, " ", h^{-1}, "]"))
    ) +
    scale_color_manual(values = color_mapping_condition) +
    theme_bw() +
    theme_tidy() +
    theme(
      plot.margin = margin(t = 1, 
                           r = 1, 
                           b = 0.1, 
                           l = 0.1, unit = "pt")
    )
  
})

names(reaction_plots) <- reactions

reaction_plots[["HEX1"]]

reaction_plots[["succoa -> succ"]]
reaction_plots[["PFK"]]
reaction_plots[["PGM"]]



reaction_plots[["PYRtm"]]
reaction_plots[["PYRt2m"]]
reaction_plots[["OAADC"]]




dir.create(paste0("figures/models/",type,"/reaction_fluxes"), recursive = TRUE, showWarnings = FALSE)

for (rxn in names(reaction_plots)) {
  
  ggsave(
    file.path(paste0("figures/models/",type, "/reaction_fluxes/",rxn, ".png")),
    reaction_plots[[rxn]],
    width = 20,
    height = 20,
    units = "mm",
    dpi = 300,
    bg = "white"
  )
  
  # ggsave(
  #   file.path("figures/models/exp_constrained/reaction_fluxes",
  #             paste0(rxn, ".pdf")),
  #   reaction_plots[[rxn]],
  #   width = 3,
  #   height = 2,
  #   device = cairo_pdf
  # )
}

out_dir <- file.path(
  "figures/models",
  type,
  "reaction_fluxes"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

for (rxn in names(reaction_plots)) {
  
  # Make a filesystem-safe filename
  safe_rxn <- rxn |>
    gsub(" -> ", "_to_", x = _) |>
    gsub("[^[:alnum:]_.-]", "_", x = _)
  
  ggsave(
    filename = file.path(out_dir, paste0(safe_rxn, ".png")),
    plot = reaction_plots[[rxn]],
    width = 20,
    height = 20,
    units = "mm",
    dpi = 300,
    bg = "white"
  )
}




# Heatmap -----------------------------------------------------------------


# Define reactions with their pathway groupings
# Define reactions with their pathway groupings
pathway_groups <- data.frame(
  Reaction = c(
    # Glycolysis
    "HEX1", "PGI", "PFK", "FBA", "TPI", "GAPD",
    "PGK", "PGM", "ENO", "PYK","r0122",
    
    # Pentose Phosphate Pathway
    "TKT2", "TALA", "TKT1", "RPE", "RPI", "GND", "PGL", "G6PDH2r",
    
    # Amino Acid Metabolism
    "ALATA_L",
    
    # Fermentation
    "LDH_D", "LDH_L",
    
    # Pyruvate Transport
    "PYRt2m", "PYRtm",
    
    # Anaplerotic Reactions
    "PCm", "OAADC",
    
    # TCA Cycle
    "CSm", "ACONTam", "ACONTbm", "ICDHyrm",
    "AKGDm", "SUCD1m", "SUCOAS1m", "succoa -> succ", "FUMm", "MDHm",
    
    # Lactate Metabolism and Transport
    "LDH_Lm", "L_LACtcm", "L_LACt4r"
  ),
  
  Pathway = c(
    # Glycolysis (10 reactions)
    rep("Glycolysis", 11),
    
    # Pentose Phosphate Pathway (8 reactions)
    rep("Pentose Phosphate", 8),
    
    # Amino Acid Metabolism (1 reaction)
    rep("Amino Acid", 1),
    
    # Fermentation (2 reactions)
    rep("Fermentation", 2),
    
    # Pyruvate Transport (2 reactions)
    rep("Pyruvate Transport", 2),
    
    # Anaplerotic Reactions (3 reactions)
    rep("Anaplerosis", 2),
    
    # TCA Cycle (10 reactions)
    rep("TCA", 10),
    
    # Lactate Metabolism and Transport (3 reactions)
    rep("Lactate Metabolism", 3)
  ),
  
  stringsAsFactors = FALSE
)

flux_cols <- c("CT_1", "CT_2", "CT_3", "CT_4", "CT_5",
               "TS_1", "TS_2", "TS_3", "TS_4", "TS_5")

# Prepare data for heatmap
flux_heatmap_df <- flux_plot_df %>%
  filter(Reaction %in% pathway_groups$Reaction) %>%
  mutate(
    Condition = recode(
      Condition,
      "Constant" = "CT",
      "Temp. shifted" = "TS"
    ),
    Column = paste(Condition, Window, sep = "_")
  ) %>%  
  select(Reaction, Column, pFBA_Flux) %>%
  pivot_wider(
    names_from = Column,
    values_from = pFBA_Flux
  ) %>%
  left_join(pathway_groups, by = "Reaction") %>%
  mutate(
    Reaction = recode(Reaction, "r0122" = "r0122/PYK3"),
    Reaction = factor(
      Reaction,
      levels = recode(pathway_groups$Reaction, "r0122" = "r0122/PYK3")
    )
  ) %>%
  arrange(Pathway, Reaction)

# # Create the new reaction
# new_reaction <- flux_heatmap_df %>%
#   filter(Reaction %in% c("SUCOAS1m", "OCOAT1m")) %>%
#   summarise(
#     across(all_of(flux_cols), ~ -sum(.x)),
#     Reaction = "succoa -> succ",
#     Pathway = "TCA"
#   )
# 
# # Remove the two original reactions and add the new one
# flux_heatmap_df <- flux_heatmap_df %>%
#   filter(!Reaction %in% c("SUCOAS1m", "OCOAT1m")) %>%
#   bind_rows(new_reaction)


# Extract matrix for heatmap
matrix_data <- flux_heatmap_df %>%
  select(-Reaction, -Pathway) %>%
  as.matrix()

rownames(matrix_data) <- flux_heatmap_df$Reaction

# Define column order: CT1-CT5, TS1-TS5
# column_order <- c(
#   paste("CT", 1:5, sep = "_"),
#   paste("TS", 1:5, sep = "_")
# )

column_order <- as.vector(rbind(
  paste0("CT_", 1:5),
  paste0("TS_", 1:5)
))

# Make sure all columns exist
column_order <- column_order[column_order %in% colnames(matrix_data)]

# Reorder matrix columns
matrix_data <- matrix_data[, column_order]

# For horizontal orientation, transpose the matrix
matrix_data_t <- t(matrix_data)


# # Color settings - quantile-based for both negative and positive
# neg_vals <- matrix_data[matrix_data < 0]
# pos_vals <- matrix_data[matrix_data > 0]
# # # Color settings - symmetric around zero
# max_abs <- max(abs(matrix_data), na.rm = TRUE)
# 
# # Define the range that should receive the full color gradient
# color_limit <- 0.1

# Define color palette for pathways
pathway_colors <- c(
  "Glycolysis"          = "#D73027",  # red
  "Pentose Phosphate"   = "#4575B4",  # blue
  "Amino Acid"          = "#1A9850",  # green
  "Fermentation"        = "#FDAE61",  # orange
  "TCA"                 = "#7B3294",  # purple
  "Gluconeogenesis"     = "#008C95",  # teal
  "Pyruvate Transport"  = "#666666",  # dark gray
  "Anaplerosis"         = "#A6611A",  # brown
  "Lactate Metabolism"  = "#C51B7D"   # magenta
)

# Column annotation for pathways (at bottom)
column_annotation_pathway <- HeatmapAnnotation(
  Pathway = flux_heatmap_df$Pathway,
  col = list(Pathway = pathway_colors),
  show_legend = TRUE,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 6, fontface = "bold"),
  show_annotation_name = TRUE,
  annotation_height = unit(0.3, "cm"),
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontsize = 6)
  )
)

# # Row annotation for conditions (at left)
# row_annotation_condition <- rowAnnotation(
#   Condition = ifelse(grepl("Constant", rownames(matrix_data_t)), "Constant", "Temp. shifted"),
#   Window = rep(1:5, 2),
#   col = list(
#     Condition = c("Constant" = "#E6641E", "Temp. shifted" = "#4B288C"),
#     Window = setNames(colorRampPalette(c("#FFF5F5", "#8B0000"))(5), 1:5)
#   ),
#   show_legend = TRUE,
#   annotation_name_side = "top",
#   annotation_name_gp = gpar(fontsize = 6, fontface = "bold"),
#   show_annotation_name = TRUE,
#   annotation_width = unit(0.1, "cm"),
#   annotation_legend_param = list(
#     title_gp = gpar(fontsize = 6, fontface = "bold"),
#     labels_gp = gpar(fontsize = 6)
#   )
# )



# Create horizontal heatmap
heatmap_plot <- Heatmap(
  matrix_data_t,
  name = "pFBA Flux\n[mM gDCW-1 h-1]",
  
  # Column settings (pathways are columns)
  column_split = flux_heatmap_df$Pathway,
  column_title = NULL,
  column_names_side = "bottom",
  column_names_gp = gpar(fontsize = 6),
  column_names_max_height = unit(8, "cm"),
  cluster_columns = TRUE,
  cluster_column_slices = FALSE,
  show_column_dend = TRUE,
  column_dend_side = "top",
  column_dend_height = unit(1, "cm"),
  
  # Row settings (conditions are rows)
  row_order = column_order,
  row_split = sub(".*_", "", rownames(matrix_data_t)),
  # row_split = ifelse(grepl("CT", rownames(matrix_data_t)), "CT", "TS"),
  row_title = NULL,
  cluster_rows = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 6),
  row_names_centered = FALSE,
  
  
  # # Color settings - shift midpoint to make low positive values more visible
  # col = colorRamp2(
  #   c(min(matrix_data, na.rm = TRUE), 
  #     0.1,  # Changed from 0 to 0.1 to shift the white point
  #     max(matrix_data, na.rm = TRUE)),
  #   c("blue", "white", "red")
  # ),
  
  # # Color settings with discrete breaks for better separation
  # col = colorRamp2(
  #   c(min(matrix_data, na.rm = TRUE),
  #     quantile(matrix_data, 0.2, na.rm = TRUE),
  #     quantile(matrix_data, 0.4, na.rm = TRUE),
  #     quantile(matrix_data, 0.6, na.rm = TRUE),
  #     quantile(matrix_data, 0.8, na.rm = TRUE),
  #     max(matrix_data, na.rm = TRUE)),
  #   c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#FDAE61", "#D73027")
  # ),
  
  
  # col = colorRamp2(
  #   c(min(matrix_data, na.rm = TRUE),
  #     quantile(neg_vals, 0.25, na.rm = TRUE),
  #     quantile(neg_vals, 0.5, na.rm = TRUE),
  #     quantile(neg_vals, 0.75, na.rm = TRUE),
  #     0,
  #     quantile(pos_vals, 0.25, na.rm = TRUE),
  #     quantile(pos_vals, 0.5, na.rm = TRUE),
  #     quantile(pos_vals, 0.75, na.rm = TRUE),
  #     max(matrix_data, na.rm = TRUE)),
  #   c("#00008B", "#4040FF", "#8080FF", "#E0E0FF", "white",
  #     "#FEE08B", "#FDAE61", "#D73027", "#8B0000")
  # ),
  
  # col = colorRamp2(
  #   c(-max_abs, 0, max_abs),
  #   c("#313695", "white", "#D73027")
  # ),
  
  # #option 1: symmetric continuous scale, ±0.2
  # col <- colorRamp2(
  #   c(-0.20, -0.10, -0.05, 0, 0.05, 0.10, 0.20),
  #   c(
  #     "#313695",
  #     "#4575B4",
  #     "#74ADD1",
  #     "white",
  #     "#FDAE61",
  #     "#D73027",
  #     "#A50026"
  #   )
  # ),
  
  #option 2 nonlinear set of breakpoints:
  col_fun <- colorRamp2(
    c(
      -0.20,
      -0.05,
      -0.01,
      0,
      0.01,
      0.05,
      0.10,
      0.20
    ),
    c(
      "#313695",
      "#4575B4",
      "#74ADD1",
      "white",
      "#FDAE61",
      "#F46D43",
      "#D73027",
      "#A50026"
    )
  ),
  
  # Other settings
  na_col = "lightgrey",
  border = TRUE,
  border_gp = gpar(col = "black", lwd = 0.5),
  
  # Heatmap body settings - max height 60mm, max width 180mm
  heatmap_width = unit(120, "mm"),  # Slightly less than max to accommodate annotations
  heatmap_height = unit(60, "mm"),  # Slightly less than max to accommodate annotations
  
  # Legend settings option 1.
  heatmap_legend_param = list(
    title = expression("pFBA Flux\n[mM gDCW"^-1~"h"^-1~"]"),
    at = c(-0.20, -0.10, -0.05, 0, 0.05, 0.10, 0.20),
    labels = c(
      "-0.20",
      "-0.10",
      "-0.05",
      "0",
      "0.05",
      "0.10",
      "0.20"
    ),

    # direction = "horizontal",

    title_gp = gpar(
      fontsize = 6,
      fontface = "bold"
    ),
    labels_gp = gpar(
      fontsize = 6
    ),
    grid_height = unit(2, "mm"),
    grid_width = unit(2, "mm")
  ),

  # #Legend settings option 2
  # heatmap_legend_param = list(
  #   title = expression("pFBA Flux\n[mM gDCW"^-1~"h"^-1~"]"),
  #   at = c(-0.20, -0.05, -0.01, 0, 0.01, 0.05, 0.10, 0.20),
  #   labels = c(
  #     "-0.20",
  #     "-0.05",
  #     "-0.01",
  #     "0",
  #     "0.01",
  #     "0.05",
  #     "0.10",
  #     "0.20"
  #   ),
  #   title_gp = gpar(fontsize = 6, fontface = "bold"),
  #   labels_gp = gpar(fontsize = 6),
  #   grid_height = unit(2, "mm"),
  #   grid_width = unit(2, "mm")
  # ),
  
  # Annotations
  bottom_annotation = column_annotation_pathway,
  # left_annotation = row_annotation_condition
)

# Draw the heatmap
draw(heatmap_plot, 
     heatmap_legend_side = "right",
     annotation_legend_side = "right")

# Save the heatmap as PDF with exact dimensions
pdf(paste0("figures/models/",type,"/reaction_fluxes/flux_heatmap_horizontal_opt1.pdf"), 
    width = 180/25.4,  # Convert mm to inches (180 mm)
    height = 60/25.4,  # Convert mm to inches (60 mm)
    pointsize = 6)
draw(heatmap_plot, 
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()

# Also save as PNG for quick viewing
png(paste0("figures/models/",type, "/reaction_fluxes/flux_heatmap_horizontal_opt1_orderBYwindows.png"),
    width = 180, height = 60, units = "mm", res = 300)
draw(heatmap_plot, 
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()
