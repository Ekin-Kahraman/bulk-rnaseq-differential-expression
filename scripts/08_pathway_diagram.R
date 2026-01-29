#!/usr/bin/env Rscript
# Conceptual ISG signalling pathway diagram

library(ggplot2)
library(ggrepel)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

message("Generating pathway diagram...")

# Pathway components
nodes <- data.frame(
  label = c("SARS-CoV-2\nRNA", "DDX58\n(RIG-I)", "MAVS", "TBK1/IKKε", 
            "IRF3/7", "Type I IFN", "IFNAR", "JAK1/TYK2", "STAT1/2",
            "IFIT1", "IFIT2", "IFIT3", "OAS3", "CXCL10", "ISG15"),
  x = c(1, 1, 1, 1, 1, 3, 5, 5, 5, 7, 7, 7, 7, 7, 7),
  y = c(10, 8, 6, 4, 2, 2, 2, 4, 6, 9, 7.5, 6, 4.5, 3, 1.5),
  type = c("virus", "sensor", "adaptor", "kinase", "tf", "cytokine", 
           "receptor", "kinase", "tf", rep("effector", 6))
)

# Colour scheme
type_colors <- c(
  virus = "#c0392b", sensor = "#2980b9", adaptor = "#27ae60",
  kinase = "#8e44ad", tf = "#d35400", cytokine = "#16a085",
  receptor = "#2c3e50", effector = "#e74c3c"
)

# Connections
edges <- data.frame(
  x = c(1, 1, 1, 1, 1, 3, 5, 5, 5, 5, 5, 5, 5, 5),
  xend = c(1, 1, 1, 1, 3, 5, 5, 5, 7, 7, 7, 7, 7, 7),
  y = c(10, 8, 6, 4, 2, 2, 2, 4, 6, 6, 6, 6, 6, 6),
  yend = c(8, 6, 4, 2, 2, 2, 4, 6, 9, 7.5, 6, 4.5, 3, 1.5)
)

# Plot
ggplot() +
  geom_segment(data = edges, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
               color = "grey40", linewidth = 0.6) +
  geom_point(data = nodes, aes(x = x, y = y, fill = type), 
             shape = 21, size = 12, color = "white", stroke = 1.5) +
  geom_text(data = nodes, aes(x = x, y = y, label = label), 
            size = 2.5, fontface = "bold", lineheight = 0.85) +
  scale_fill_manual(values = type_colors, guide = "none") +
  annotate("text", x = 1, y = 11, label = "Viral Detection", fontface = "bold", size = 4) +
  annotate("text", x = 3, y = 3.2, label = "IFN Production", fontface = "bold", size = 4) +
  annotate("text", x = 5, y = 7.2, label = "IFN Signalling", fontface = "bold", size = 4) +
  annotate("text", x = 7, y = 10.2, label = "Antiviral Effectors", fontface = "bold", size = 4) +
  labs(title = "Interferon-Stimulated Gene (ISG) Signalling Cascade",
       subtitle = "Host antiviral response to SARS-CoV-2 infection") +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10),
    plot.margin = margin(20, 20, 20, 20),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  coord_cartesian(xlim = c(0, 8), ylim = c(0, 11.5))

ggsave("results/figures/pathway_diagram.png", width = 10, height = 8, dpi = 300, bg = "white")

message("Pathway diagram saved")