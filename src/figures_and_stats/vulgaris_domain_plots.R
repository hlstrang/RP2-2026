library(tidyverse)

setwd("Documents/year_5/rp2")
domains <- read_csv("written_report/vulgaris_all_domains.csv")
domains <- domains %>%
  mutate(
    domain_name = case_when(
      domain_name %in% c("Col1", "Col2") ~ "Triple Helix",
      TRUE ~ domain_name
    ))
hcol_types <- unique(domains$gene_id)
domains <- domains %>%
  mutate(y_pos = match(gene_id, hcol_types))

backbones <- domains %>%
  distinct(gene_id, seq_length, y_pos)

domain_colors <- domains %>%
  distinct(domain_name, colour) %>%
  deframe()

ggplot() +
  geom_segment(
    data = backbones,
    aes(x = 1, xend = seq_length, y = y_pos, yend = y_pos),
    color = "grey60", linewidth = 1.5
  ) +
  geom_rect(
    data = domains,
    aes(
      xmin = start, xmax = end,
      ymin = y_pos - 0.25, ymax = y_pos + 0.25,
      fill = domain_name
    ),
    color = "black", linewidth = 0.3
  ) +
  scale_fill_manual(values = domain_colors, name = "Domain Name") +
  scale_y_continuous(breaks = seq_along(hcol_types), labels = hcol_types) +
  labs(
    x = "Amino Acid Position (aa)",
    y = "Hcol Type",
    title = "Hcol Domain Architecture"
  ) +
  theme_basic() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right"
  )
