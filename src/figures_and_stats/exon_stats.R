library(tidyverse)
library(car)

setwd("Documents/year_5/rp2")

df <- read_csv("total_exon_stats.csv")

df <- df %>%
  mutate(
    intron_fraction = (total_intron_bp / locus_span_bp) * 100,
    exon_fraction = (total_exonic_bp / locus_span_bp) * 100
  )

df %>%
  group_by(species) %>%
  summarize(
    avg_exon = mean(exon_fraction, na.rm = TRUE),
    min        = min(exon_fraction, na.rm = TRUE),
    max        = max(exon_fraction, na.rm = TRUE),
    .groups = "drop"
  )
df %>%
  group_by(species) %>%
  summarize(
    avg_intron = mean(intron_fraction, na.rm = TRUE),
    min        = min(intron_fraction, na.rm = TRUE),
    max        = max(intron_fraction, na.rm = TRUE),
    .groups = "drop"
  )
df %>%
  group_by(species) %>%
  summarize(
    avg_ex_no = mean(n_exons, na.rm = TRUE),
    min        = min(n_exons, na.rm = TRUE),
    max        = max(n_exons, na.rm = TRUE),
    .groups = "drop"
  )
df %>%
  group_by(species) %>%
  summarize(
    avg_exon_bp = mean(mean_exon_len, na.rm = TRUE),
    min        = min(mean_exon_len, na.rm = TRUE),
    max        = max(mean_exon_len, na.rm = TRUE),
    .groups = "drop"
  )
df %>%
  group_by(species) %>%
  summarize(
    avg_locus_len = mean(locus_span_bp, na.rm = TRUE),
    min        = min(locus_span_bp, na.rm = TRUE),
    max        = max(locus_span_bp, na.rm = TRUE),
    .groups = "drop"
  )

df %>%
  group_by(hcol_type) %>%
  summarize(
    avg_exon = mean(exon_fraction, na.rm = TRUE),
    min        = min(exon_fraction, na.rm = TRUE),
    max        = max(exon_fraction, na.rm = TRUE),
    .groups = "drop"
  )
df %>%
  group_by(hcol_type) %>%
  summarize(
    avg_intron = mean(intron_fraction, na.rm = TRUE),
    min        = min(intron_fraction, na.rm = TRUE),
    max        = max(intron_fraction, na.rm = TRUE),
    .groups = "drop"
  )
df %>%
  group_by(hcol_type) %>%
  summarize(
    avg_ex_no = mean(n_exons, na.rm = TRUE),
    min        = min(n_exons, na.rm = TRUE),
    max        = max(n_exons, na.rm = TRUE),
    .groups = "drop"
  )
df %>%
  group_by(hcol_type) %>%
  summarize(
    avg_exon_bp = mean(mean_exon_len, na.rm = TRUE),
    min        = min(mean_exon_len, na.rm = TRUE),
    max        = max(mean_exon_len, na.rm = TRUE),
    .groups = "drop"
  )
df %>%
  group_by(hcol_type) %>%
  summarize(
    avg_locus_len = mean(locus_span_bp, na.rm = TRUE),
    min        = min(locus_span_bp, na.rm = TRUE),
    max        = max(locus_span_bp, na.rm = TRUE),
    .groups = "drop"
  )

df_long <- df %>%
  filter(hcol_type != "hcol2b") %>%
  select(species, hcol_type, exon_fraction, intron_fraction, locus_span_bp, n_exons) %>%
  pivot_longer(
    cols = c(exon_fraction, intron_fraction, n_exons, locus_span_bp),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(Metric = case_when(
    Metric == "exon_fraction" ~ "Exonic Fraction (%)",
    Metric == "intron_fraction" ~ "Intronic Fraction (%)",
    Metric == "n_exons" ~ "Number of Exons",
    Metric == "locus_span_bp" ~ "Length of Whole Locus"
  ),
  hcol_type = str_trim(as.character(hcol_type)),
  species = str_trim(as.character(species)))

species_map <- c(
  "h_octo" = "H. octoradiatus",
  "tripedalia" = "Tripedalia sp.",
  "h_viridissima" = "H. viridissima",
  "r_luteum" = "R. luteum"
)

plot_df <- df_long %>%
  mutate(species_clean = species_map[species])

fig <- ggplot(plot_df, aes(x = hcol_type, y = Value, fill = species_clean)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75, alpha = 0.6) +
  facet_wrap(~Metric, scales = "free_y") +
  scale_fill_manual(values = c("H. octoradiatus" = "purple", "R. luteum" = "deepskyblue",
                               "H. viridissima" = "olivedrab", "Tripedalia sp." = "red")) +
  theme_basic() +
  labs(x = "Collagen Type", y = NULL, fill = "Species") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6))
ggsave("exon_stats.png", fig, height = 9, width = 14)

model1a <- lm(exon_fraction ~ species, data = df)
summary(model1a)
anova(model1a)
model1b <- lm(exon_fraction ~ hcol_type, data = df)
summary(model1b)
anova(model1b)

model2a <- lm(intron_fraction ~ species, data = df)
summary(model2a)
model2b <- lm(intron_fraction ~ hcol_type, data = df)
summary(model2b)

model3a <- lm(locus_span_bp ~ species, data = df)
summary(model3a)
anova(model3a)
model3b <- lm(locus_span_bp ~ hcol_type, data = df)
summary(model3b)

model4a <- lm(n_exons ~ species, data = df)
summary(model4a)
anova(model4a)
model4b <- lm(n_exons ~ hcol_type, data = df)
summary(model4b)
anova(model4b)

model5a <- lm(mean_exon_len ~ species, data = df)
summary(model5a)
anova(model5a)
model5b <- lm(mean_exon_len ~ hcol_type, data = df)
summary(model5b)
anova(model5b)
