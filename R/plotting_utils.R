
# =============================================================================
# Global plot settings
# =============================================================================

PLOT_FONT_SIZE <- 6
PANEL_LABEL_SIZE <- 8
FONT_FAMILY <- "Liberation Sans"

LINE_WIDTH <- 0.6
RIBBON_ALPHA <- 0.25
ERRORBAR_ALPHA <- 0.8
POINT_ALPHA <- 0.6
POINT_SIZE <- 1
PVAL_BRACKET_OFFSET <- 0.06
PVAL_TEXT_OFFSET <- 0.02
BRACKET_LINE_WIDTH <- 0.3

color_condition <- c(
  "CT" = "#D94801",
  "TS" = "#6A51A3"
)

color_condition_long <- c(
  "Constant" = "#D94801",
  "Temp. shifted" = "#6A51A3",
  "Temp. Shifted" = "#6A51A3"
)

event_color <- "#58A787"

theme_tidy <- function(
    base_size = PLOT_FONT_SIZE,
    base_family = FONT_FAMILY,
    text_colour = "black") {
  theme(
    plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"),
    plot.background = element_rect(fill = NA, colour = NA),
    legend.background = element_rect(fill = NA, colour = NA),
    legend.key = element_rect(fill = NA, colour = NA),
    strip.background = element_rect(fill = NA, colour = NA),
    panel.background = element_rect(fill = NA, colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks = element_line(linewidth = 0.3, colour = "black"),
    legend.position = "none",
    legend.direction = "horizontal",
    text = element_text(size = base_size, family = base_family, colour = text_colour),
    plot.title = element_text(size = base_size, family = base_family, colour = text_colour, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = base_size, family = base_family, colour = text_colour),
    axis.text = element_text(size = base_size, family = base_family, colour = text_colour),
    axis.title = element_text(size = base_size, family = base_family, colour = text_colour),
    legend.title = element_text(size = base_size, family = base_family, face = "bold", colour = text_colour),
    legend.text = element_text(size = base_size, family = base_family, colour = text_colour),
    strip.text = element_text(size = base_size, family = base_family, colour = text_colour)
  )
}

pt_to_geom <- function(pt) pt / 2.845