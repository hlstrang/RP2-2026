library(tidyverse)
library(emmeans)
library(glmmTMB)
library(car)
setwd("Documents/year_5/rp2")

df <- read_csv("collagenome_stats/collagenome_stats.csv")
df <- df %>%
  filter(hcol_type != "others")

domain_colors <- read_csv("written_report/vulgaris_all_domains.csv") %>%
  distinct(domain_name, colour) %>%
  pull(colour, name = domain_name)
domain_names_cleaned <- tolower(names(domain_colors))
names(domain_colors) <- domain_names_cleaned
domain_colors$col1 <- "forestgreen"
domain_colors$col2 <- "#8FD6A5"

## Median Lengths
structure_median <- df %>%
  select(hcol_type, species_class, median_sp_length, median_tspn_length,
        median_col1_length_agg,median_col2_length_agg,hmm_avg_len_COLFI_median,
         hmm_avg_len_WAP_median, hmm_avg_len_VWA_median, hmm_avg_len_C4_median) %>%
  pivot_longer(
    cols = -c(hcol_type, species_class),
    names_to = "region",
    values_to = "length"
  ) %>%
  mutate(length = ifelse(is.na(length) | is.nan(length), 0, length)) %>%
  mutate(region = gsub("median_|mean_|_length|_agg|hmm_avg_len_|_median|_mean", "", region))%>%
  mutate(region = tolower(region)) %>%
  mutate(
    region = factor(
      region,
      levels = c(
        "sp",
        "tspn",
        "vwa",
        "wap",
        "col1",
        "col2",
        "colf1",
        "c4"
      )
    )
  )

ggplot(structure_median, aes(y = hcol_type, x = length, fill = region)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.7) +
  facet_grid(species_class ~ ., scales = "free_y", space = "free_y") +
  theme_basic() +
  scale_fill_manual(
    values = domain_colors,
    guide = guide_legend(reverse = TRUE)
  ) +
  labs(
    title = "Domain Architecture",
    x = "Median Domain Instance Length (Amino Acids)",
    y = "Collagen Type",
    fill = "Domain Type"
  ) +
  theme(
    strip.text.y = element_text(angle = 0, face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

## Mean Lengths
structure_mean <- df %>%
  select(hcol_type, species_class, mean_sp_length, mean_tspn_length,
         mean_col1_length_agg, mean_col2_length_agg, hmm_avg_len_COLFI_mean,
         hmm_avg_len_WAP_mean, hmm_avg_len_VWA_mean, hmm_avg_len_C4_mean) %>%
  pivot_longer(
    cols = -c(hcol_type, species_class),
    names_to = "region",
    values_to = "length"
  ) %>%
  mutate(length = ifelse(is.na(length) | is.nan(length), 0, length)) %>%
  mutate(region = gsub("median_|mean_|_length|_agg|hmm_avg_len_|_median|_mean", "", region))%>%
  mutate(region = tolower(region)) %>%
  mutate(
    region = factor(
      region,
      levels = c(
        "sp",
        "tspn",
        "vwa",
        "wap",
        "col1",
        "col2",
        "colf1",
        "c4"
      )
    )
  )

ggplot(structure_mean, aes(y = hcol_type, x = length, fill = region)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.7) +
  facet_grid(species_class ~ ., scales = "free_y", space = "free_y") +
  theme_basic() +
  scale_fill_manual(
    values = domain_colors,
    guide = guide_legend(reverse = TRUE)
  ) +
  labs(
    title = "Domain Architecture Comparison",
    x = "Average Region Length (Amino Acids)",
    y = "Collagen Type",
    fill = "Domain / Region"
  ) +
  theme(
    strip.text.y = element_text(angle = 0, face = "bold"),
    panel.grid.minor = element_blank()
  )

