library(tidyverse)
source("R/plotting_utils.R", local = FALSE)


flux_df <- read_csv("analysis/models/exp_constrained/reaction_fluxes_FBA_pFBA.csv") %>%
  mutate(common_name = recode(
    common_name,
    "iCHO3K_cho_prod_generic_unblocked" = "iCHO3K"
  )) %>%
  filter(common_name == "iCHO3K")


reactions <- c("HEX1", "PGI", "TALA", "TKT1", "RPE",
               "RPI", "GAUGE-R10221", "PGL", "G6PDH2r",
               "PFK", "FBA", "TPI", "GAPD", "ACYP", "PGCD", "PSERT",
               "PGK", "DPGM", "DPGASE", "BPGP", "ENO", "PYK", "LDH_D","LDH_L",
               "ALATA_L", "PYRt2m", "PCm", "CSm", "r2373", "HMR_4777",
               "SUCD1m", "FUMm", "MDHm")


color_mapping_condition <- c(
  "Constant" = "#E6641E",
  "Temp. shifted" = "#4B288C"
)

flux_plot_df <- flux_df %>%
  filter(Reaction %in% reactions) %>%
  mutate(
    Window = factor(Window, levels = 1:5),
    Reaction = factor(Reaction, levels = reactions)
  )

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
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    labs(
      title = rxn,
      x = "Window",
      y = "pFBA Flux"
    ) +
    scale_color_manual(values = color_mapping_condition) +
    theme_bw() +
    theme_tidy()
  
})

names(reaction_plots) <- reactions

reaction_plots[["HEX1"]]
reaction_plots[["PGI"]]
reaction_plots[["PFK"]]



dir.create("figures/models/exp_constrained/reaction_fluxes", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/models/exp_constrained/reaction_fluxes", recursive = TRUE, showWarnings = FALSE)

for (rxn in names(reaction_plots)) {
  
  ggsave(
    file.path("figures/models/exp_constrained/reaction_fluxes", paste0(rxn, ".png")),
    reaction_plots[[rxn]],
    width = 6,
    height = 4,
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    file.path("figures/models/exp_constrained/reaction_fluxes",
              paste0(rxn, ".pdf")),
    reaction_plots[[rxn]],
    width = 6,
    height = 4,
    device = cairo_pdf
  )
}
