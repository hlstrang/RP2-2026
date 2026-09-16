library(tidyverse)
library(Biostrings)
library(FSA)
library(ggstatsplot)
library(statsExpressions)
library(ggpubr)
library(pheatmap)
library(flextable)
library(rstatix)

## Statistics
full_df <- read_csv("collagenome_stats/combined_collagenome.csv")
full_df <- full_df %>%
  mutate(hcol_type = tolower(hcol_type) %>% str_remove(" .*"))

hcol1_full_df <- full_df %>%
  filter(hcol_type == "hcol1") %>%
  select(species, class, accession, sequence)
names(hcol1_full_df$sequence) <- hcol1_full_df$species
aa_set <- AAStringSet(hcol1_full_df$sequence)
writeXStringSet(aa_set, filepath = "hcol1_seqs.fasta")

all_domains <- read_csv("collagenome_stats/collagenome_domains_long.csv")
all_domains <- all_domains %>%
  filter(hcol_type != "others")

hcol1 <- full_df %>%
  filter(hcol_type == "hcol1")
nrow(hcol1)
hcol2a <- full_df %>%
  filter(hcol_type == "hcol2a")
nrow(hcol2a)
hcol2b <- full_df %>%
  filter(hcol_type == "hcol2b")
nrow(hcol2b)
hcol3 <- full_df %>%
  filter(hcol_type == "hcol3")
nrow(hcol3)
hcol4 <- full_df %>%
  filter(hcol_type == "hcol4")
nrow(hcol4)
hcol5 <- full_df %>%
  filter(hcol_type == "hcol5")
nrow(hcol5)
hcol6 <- full_df %>%
  filter(hcol_type == "hcol6")
nrow(hcol6)
hcol7 <- full_df %>%
  filter(hcol_type == "hcol7")
nrow(hcol7)
hcol8 <- full_df %>%
  filter(hcol_type == "hcol8")
nrow(hcol8)

## Hcol1 Stats
hcol1 <- hcol1 %>%
  mutate(seq_length = nchar(sequence))
hcol1_accession_domains <- all_domains %>%
  filter(hcol_type == "hcol1") %>%
  group_by(accession, domain_name) %>%
  summarise(
    domain_count = n(),
    total_domain_length = sum(length, na.rm = TRUE),
    .groups = "drop"
  )
hcol1_domain_stats <- hcol1 %>%
  select(accession, species, class, sequence) %>%
  mutate(total_seq_length = nchar(sequence)) %>%
  left_join(hcol1_accession_domains, by = "accession") %>%
  mutate(
    domain_prop = total_domain_length / total_seq_length
  )
ggplot(data = hcol1, aes(x = seq_length)) +
  geom_bar()

hcol1_summary <- hcol1_domain_stats %>%
  group_by(domain_name) %>%
  summarise(
    n_sequences       = n_distinct(accession),
    med_domain_count = median(domain_count, na.rm = TRUE),
    med_total_length = median(total_domain_length, na.rm = TRUE),
    sd_total_length   = sd(total_domain_length, na.rm = TRUE),
    med_prop         = median(domain_prop, na.rm = TRUE),
    sd_prop           = sd(domain_prop, na.rm = TRUE),
    med_prop_pct     = sprintf("%.2f%%", med_prop * 100),
    .groups           = "drop"
  )

print(hcol1_summary)
median(hcol1_domain_stats$total_seq_length)

## Hcol2 Stats
hcol2 <- hcol2a %>%
  bind_rows(hcol2b) %>%
  mutate(seq_length = nchar(sequence))
hcol2_accession_domains <- all_domains %>%
  filter(hcol_type == "hcol2") %>%
  group_by(accession, domain_name) %>%
  summarise(
    domain_count = n(),
    total_domain_length = sum(length, na.rm = TRUE),
    .groups = "drop"
  )
hcol2_domain_stats <- hcol2 %>%
  select(accession, species, class, sequence) %>%
  mutate(total_seq_length = nchar(sequence)) %>%
  left_join(hcol2_accession_domains, by = "accession") %>%
  mutate(
    domain_prop = total_domain_length / total_seq_length
  )
ggplot(data = hcol2, aes(x = seq_length)) +
  geom_bar()

hcol2_summary <- hcol2_domain_stats %>%
  group_by(domain_name) %>%
  summarise(
    n_sequences       = n_distinct(accession),
    med_domain_count = median(domain_count, na.rm = TRUE),
    med_total_length = median(total_domain_length, na.rm = TRUE),
    sd_total_length   = sd(total_domain_length, na.rm = TRUE),
    med_prop         = median(domain_prop, na.rm = TRUE),
    sd_prop           = sd(domain_prop, na.rm = TRUE),
    med_prop_pct     = sprintf("%.2f%%", med_prop * 100),
    .groups           = "drop"
  )

print(hcol2_summary)
median(hcol2_domain_stats$total_seq_length)

## Hcol3 Stats
hcol3 <- hcol3 %>%
  filter(accession != "ICXV01221851.1") %>%
  filter(accession != "CDWRCV010000266.1") %>%
  mutate(seq_length = nchar(sequence))
hcol3_accession_domains <- all_domains %>%
  filter(hcol_type == "hcol3") %>%
  group_by(accession, domain_name) %>%
  summarise(
    domain_count = n(),
    total_domain_length = sum(length, na.rm = TRUE),
    .groups = "drop"
  )
hcol3_domain_stats <- hcol3 %>%
  select(accession, species, class, sequence) %>%
  mutate(total_seq_length = nchar(sequence)) %>%
  left_join(hcol3_accession_domains, by = "accession") %>%
  mutate(
    domain_prop = total_domain_length / total_seq_length
  )
ggplot(data = hcol3, aes(x = seq_length)) +
  geom_bar()

hcol3_summary <- hcol3_domain_stats %>%
  group_by(domain_name) %>%
  summarise(
    n_sequences       = n_distinct(accession),
    med_domain_count = median(domain_count, na.rm = TRUE),
    med_total_length = median(total_domain_length, na.rm = TRUE),
    sd_total_length   = sd(total_domain_length, na.rm = TRUE),
    med_prop         = median(domain_prop, na.rm = TRUE),
    sd_prop           = sd(domain_prop, na.rm = TRUE),
    med_prop_pct     = sprintf("%.2f%%", med_prop * 100),
    .groups           = "drop"
  )

print(hcol3_summary)
mean(hcol3_domain_stats$total_seq_length)

## Hcol4 Stats
hcol4 <- hcol4 %>%
  filter(accession != "GANB01000001.1") %>%
  filter(accession != "HAGZ01055176.1") %>%
  mutate(seq_length = nchar(sequence))
hcol4_accession_domains <- all_domains %>%
  filter(hcol_type == "hcol4") %>%
  group_by(accession, domain_name) %>%
  summarise(
    domain_count = n(),
    total_domain_length = sum(length, na.rm = TRUE),
    .groups = "drop"
  )
hcol4_domain_stats <- hcol4 %>%
  select(accession, species, class, sequence) %>%
  mutate(total_seq_length = nchar(sequence)) %>%
  left_join(hcol4_accession_domains, by = "accession") %>%
  mutate(
    domain_prop = total_domain_length / total_seq_length
  )
ggplot(data = hcol4, aes(x = seq_length)) +
  geom_bar()

hcol4_summary <- hcol4_domain_stats %>%
  group_by(domain_name) %>%
  summarise(
    n_sequences       = n_distinct(accession),
    med_domain_count = median(domain_count, na.rm = TRUE),
    med_total_length = median(total_domain_length, na.rm = TRUE),
    sd_total_length   = sd(total_domain_length, na.rm = TRUE),
    med_prop         = median(domain_prop, na.rm = TRUE),
    sd_prop           = sd(domain_prop, na.rm = TRUE),
    med_prop_pct     = sprintf("%.2f%%", med_prop * 100),
    .groups           = "drop"
  )

print(hcol4_summary)
mean(hcol4_domain_stats$total_seq_length)

## Hcol5 Stats
hcol5 <- hcol5 %>%
  mutate(seq_length = nchar(sequence))
hcol5_accession_domains <- all_domains %>%
  filter(hcol_type == "hcol5") %>%
  group_by(accession, domain_name) %>%
  summarise(
    domain_count = n(),
    total_domain_length = sum(length, na.rm = TRUE),
    .groups = "drop"
  )
hcol5_domain_stats <- hcol5 %>%
  select(accession, species, class, sequence) %>%
  mutate(total_seq_length = nchar(sequence)) %>%
  left_join(hcol5_accession_domains, by = "accession") %>%
  mutate(
    domain_prop = total_domain_length / total_seq_length
  )
ggplot(data = hcol5, aes(x = seq_length)) +
  geom_bar()

hcol5_summary <- hcol5_domain_stats %>%
  group_by(domain_name) %>%
  summarise(
    n_sequences       = n_distinct(accession),
    med_domain_count = median(domain_count, na.rm = TRUE),
    med_total_length = median(total_domain_length, na.rm = TRUE),
    sd_total_length   = sd(total_domain_length, na.rm = TRUE),
    med_prop         = median(domain_prop, na.rm = TRUE),
    sd_prop           = sd(domain_prop, na.rm = TRUE),
    med_prop_pct     = sprintf("%.2f%%", med_prop * 100),
    .groups           = "drop"
  )
print(hcol5_summary)
mean(hcol5_domain_stats$total_seq_length)

## Hcol6 Stats
hcol6 <- hcol6 %>%
  filter(accession != "CAXPCZ010003408.1") %>%
  filter(accession != "GFGU01247931.1") %>%
  mutate(seq_length = nchar(sequence))
hcol6_accession_domains <- all_domains %>%
  filter(hcol_type == "hcol6") %>%
  group_by(accession, domain_name) %>%
  summarise(
    domain_count = n(),
    total_domain_length = sum(length, na.rm = TRUE),
    .groups = "drop"
  )
hcol6_domain_stats <- hcol6 %>%
  select(accession, species, class, sequence) %>%
  mutate(total_seq_length = nchar(sequence)) %>%
  left_join(hcol6_accession_domains, by = "accession") %>%
  mutate(
    domain_prop = total_domain_length / total_seq_length
  )
ggplot(data = hcol6, aes(x = seq_length)) +
  geom_bar()

hcol6_summary <- hcol6_domain_stats %>%
  group_by(domain_name) %>%
  summarise(
    n_sequences       = n_distinct(accession),
    med_domain_count = median(domain_count, na.rm = TRUE),
    med_total_length = median(total_domain_length, na.rm = TRUE),
    sd_total_length   = sd(total_domain_length, na.rm = TRUE),
    med_prop         = median(domain_prop, na.rm = TRUE),
    sd_prop           = sd(domain_prop, na.rm = TRUE),
    med_prop_pct     = sprintf("%.2f%%", med_prop * 100),
    .groups           = "drop"
  )
print(hcol6_summary)
mean(hcol6_domain_stats$total_seq_length)

## Hcol7 Stats
hcol7 <- hcol7 %>%
  filter(!accession %in% c("GHUC01000492.1", "GAOL01025660.1", "GFAS01287089.1", "HAHB01067895.1")) %>%
  filter(species != "Podocoryna carnea") %>%
  mutate(seq_length = nchar(sequence))
hcol7_accession_domains <- all_domains %>%
  filter(hcol_type == "hcol7") %>%
  group_by(accession, domain_name) %>%
  summarise(
    domain_count = n(),
    total_domain_length = sum(length, na.rm = TRUE),
    .groups = "drop"
  )
hcol7_domain_stats <- hcol7 %>%
  select(accession, species, class, sequence) %>%
  mutate(total_seq_length = nchar(sequence)) %>%
  left_join(hcol7_accession_domains, by = "accession") %>%
  mutate(
    domain_prop = total_domain_length / total_seq_length
  )
ggplot(data = hcol7, aes(x = seq_length)) +
  geom_bar()

hcol7_summary <- hcol7_domain_stats %>%
  group_by(domain_name) %>%
  summarise(
    n_sequences       = n_distinct(accession),
    med_domain_count = median(domain_count, na.rm = TRUE),
    med_total_length = median(total_domain_length, na.rm = TRUE),
    sd_total_length   = sd(total_domain_length, na.rm = TRUE),
    med_prop         = median(domain_prop, na.rm = TRUE),
    sd_prop           = sd(domain_prop, na.rm = TRUE),
    med_prop_pct     = sprintf("%.2f%%", med_prop * 100),
    .groups           = "drop"
  )
print(hcol7_summary)
mean(hcol7_domain_stats$total_seq_length)

## Hcol8 Stats
hcol8 <- hcol8 %>%
  mutate(seq_length = nchar(sequence))
hcol8_accession_domains <- all_domains %>%
  filter(hcol_type == "hcol8") %>%
  group_by(accession, domain_name) %>%
  summarise(
    domain_count = n(),
    total_domain_length = sum(length, na.rm = TRUE),
    .groups = "drop"
  )
hcol8_domain_stats <- hcol8 %>%
  select(accession, species, class, sequence) %>%
  mutate(total_seq_length = nchar(sequence)) %>%
  left_join(hcol8_accession_domains, by = "accession") %>%
  mutate(
    domain_prop = total_domain_length / total_seq_length
  )
ggplot(data = hcol8, aes(x = seq_length)) +
  geom_bar()

hcol8_summary <- hcol8_domain_stats %>%
  group_by(domain_name) %>%
  summarise(
    n_sequences       = n_distinct(accession),
    med_domain_count = median(domain_count, na.rm = TRUE),
    med_total_length = median(total_domain_length, na.rm = TRUE),
    sd_total_length   = sd(total_domain_length, na.rm = TRUE),
    med_prop         = median(domain_prop, na.rm = TRUE),
    sd_prop           = sd(domain_prop, na.rm = TRUE),
    med_prop_pct     = sprintf("%.2f%%", med_prop * 100),
    .groups           = "drop"
  )
print(hcol8_summary)
mean(hcol8_domain_stats$total_seq_length)

## All Types
complete_summary <- bind_rows(
  hcol1 = hcol1_summary,
  hcol2 = hcol2_summary,
  hcol3 = hcol3_summary,
  hcol4 = hcol4_summary,
  hcol5 = hcol5_summary,
  hcol6 = hcol6_summary,
  hcol7 = hcol7_summary,
  hcol8 = hcol8_summary,
  .id = "collagen_type"
) %>%
  select(collagen_type, domain_name, med_domain_count, med_prop_pct)

seq_lengths <- bind_rows(
  hcol1_domain_stats %>% mutate(collagen_type = "hcol1"),
  hcol2_domain_stats %>% mutate(collagen_type = "hcol2"),
  hcol3_domain_stats %>% mutate(collagen_type = "hcol3"),
  hcol4_domain_stats %>% mutate(collagen_type = "hcol4"),
  hcol5_domain_stats %>% mutate(collagen_type = "hcol5"),
  hcol6_domain_stats %>% mutate(collagen_type = "hcol6"),
  hcol7_domain_stats %>% mutate(collagen_type = "hcol7"),
  hcol8_domain_stats %>% mutate(collagen_type = "hcol8")
) %>%
  group_by(collagen_type) %>%
  summarise(
    med_seq_length = round(median(total_seq_length, na.rm = TRUE),0),
    sd_seq_length = sd(total_seq_length, na.rm = TRUE),
    .groups = "drop"
  )

complete_summary <- complete_summary %>%
  left_join(seq_lengths, by = "collagen_type")

wide_summary <- complete_summary %>%
  pivot_wider(
    id_cols = c(collagen_type, med_seq_length),
    names_from = domain_name,
    values_from = med_prop_pct
  ) %>%
  arrange(collagen_type)
write_csv(wide_summary, "collagenome_stats/median_summary.csv")

ft_summary <- flextable(wide_summary) %>%
  theme_zebra() %>%
  autofit()
ft_summary <- align(ft_summary, align = "center", part = "all")
ft_summary <- set_header_labels(ft_summary,
                                collagen_type = "Collagen Type",
                                med_seq_length = "Median Length (aa)",
                                Signal_Peptide = "Signal Peptide")
ft_summary

all_domain_comparisons <- hcol1_domain_stats %>%
  mutate(collagen_type = "Hcol1") %>%
  bind_rows(
    hcol2_domain_stats %>% mutate(collagen_type = "Hcol2"),
    hcol3_domain_stats %>% mutate(collagen_type = "Hcol3"),
    hcol4_domain_stats %>% mutate(collagen_type = "Hcol4"),
    hcol5_domain_stats %>% mutate(collagen_type = "Hcol5"),
    hcol6_domain_stats %>% mutate(collagen_type = "Hcol6"),
    hcol7_domain_stats %>% mutate(collagen_type = "Hcol7"),
    hcol8_domain_stats %>% mutate(collagen_type = "Hcol8")
  )

get_pairwise_vs_hydrozoa <- function(data_df, collagen_type, domain_list) {

  results <- lapply(domain_list, function(dom) {
    domain_data <- data_df %>% filter(domain_name == dom)

    if(nrow(domain_data) >= 6 && length(unique(domain_data$class)) >= 2) {
      classes <- unique(domain_data$class)

      hydrozoa_comparisons <- lapply(classes, function(other_class) {
        if(other_class != "Hydrozoans") {
          d1 <- (domain_data %>% filter(class == "Hydrozoans"))$total_domain_length
          d2 <- (domain_data %>% filter(class == other_class))$total_domain_length

          wtest <- wilcox.test(d2, d1, exact = FALSE, alternative = "two.sided")

          n1 <- length(d1)
          n2 <- length(d2)
          N <- n1 + n2

          z_val <- qnorm(wtest$p.value / 2) * -1
          r_effect <- z_val / sqrt(N)

          med1 <- median(d1, na.rm = TRUE)
          med2 <- median(d2, na.rm = TRUE)
          pct_diff <- ((med2 - med1) / med1) * 100

          tibble(
            domain = dom,
            compared_class = other_class,
            baseline_class = "Hydrozoans",
            median_baseline = med1,
            median_compared = med2,
            raw_diff = med2 - med1,
            pct_diff = pct_diff,
            abs_pct = abs(pct_diff),
            p.raw = wtest$p.value,
            p.adj = wtest$p.value,
            z_value = z_val,
            r_effect = r_effect,
            effect_size_label = case_when(
              r_effect < 0.1 ~ "negligible",
              r_effect < 0.3 ~ "small",
              r_effect < 0.5 ~ "medium",
              TRUE ~ "large"
            )
          )
        }
      }) %>% purrr::discard(is.null)

      if(length(hydrozoa_comparisons) > 0) {
        effects_df <- bind_rows(hydrozoa_comparisons)

        sorted_idx <- order(effects_df$p.raw)
        adjusted_p <- p.adjust(effects_df$p.raw[sorted_idx], method = "holm")
        effects_df$p.adj[sorted_idx] <- adjusted_p

        effects_df %>%
          mutate(
            p.signif = case_when(
              p.adj < 0.001 ~ "***",
              p.adj < 0.01 ~ "**",
              p.adj < 0.05 ~ "*",
              TRUE ~ ""
            ),
            log_p = -log10(p.adj)
          )
      } else {
        tibble(
          domain = dom,
          compared_class = character(),
          baseline_class = character(),
          median_baseline = numeric(),
          median_compared = numeric(),
          raw_diff = numeric(),
          pct_diff = numeric(),
          abs_pct = numeric(),
          p.raw = numeric(),
          p.adj = numeric(),
          p.signif = character(),
          r_effect = numeric(),
          effect_size_label = character()
        )
      }
    } else {
      tibble(
        domain = dom,
        compared_class = character(),
        baseline_class = character(),
        median_baseline = numeric(),
        median_compared = numeric(),
        raw_diff = numeric(),
        pct_diff = numeric(),
        abs_pct = numeric(),
        p.raw = numeric(),
        p.adj = numeric(),
        p.signif = character(),
        r_effect = numeric(),
        effect_size_label = character()
      )
    }
  })

  bind_rows(results) %>%
    mutate(collagen_type = collagen_type)
}

all_vs_hydrozoa <- bind_rows(
  get_pairwise_vs_hydrozoa(hcol1_domain_stats, "Hcol1", c("Col1", "Col2", "COLFI", "Signal_Peptide")),
  get_pairwise_vs_hydrozoa(hcol2_domain_stats, "Hcol2", c("Col1", "Col2", "COLFI", "WAP", "Signal_Peptide")),
  get_pairwise_vs_hydrozoa(hcol3_domain_stats, "Hcol3", c("Col1", "Col2", "COLFI", "Signal_Peptide", "WAP")),
  get_pairwise_vs_hydrozoa(hcol4_domain_stats, "Hcol4", c("Col1", "Col2", "C4", "Signal_Peptide")),
  get_pairwise_vs_hydrozoa(hcol5_domain_stats, "Hcol5", c("Col1", "Col2", "COLFI", "Signal_Peptide")),
  get_pairwise_vs_hydrozoa(hcol6_domain_stats, "Hcol6", c("Col1", "Col2", "C4", "Signal_Peptide")),
  get_pairwise_vs_hydrozoa(hcol7_domain_stats, "Hcol7", c("Col1", "Col2", "COLFI", "Signal_Peptide", "TSPN")),
  get_pairwise_vs_hydrozoa(hcol8_domain_stats, "Hcol8", c("Col1", "Col2", "COLFI", "Signal_Peptide"))
)

sig_vs_hydrozoa <- all_vs_hydrozoa %>%
  filter(p.adj < 0.05 & !is.na(p.adj))

ft_summary <- flextable(sig_vs_hydrozoa) %>%
  theme_zebra() %>%
  autofit()
ft_summary <- align(ft_summary, align = "center", part = "all")
ft_summary

abbr_map <- c(
  "Cubozoans" = "Cubozoan",
  "Hydrozoans" = "Hydrozoan",
  "Scyphozoans" = "Scyphozoan",
  "Staurozoans" = "Staurozoan"
)

heatmap_data <- sig_vs_hydrozoa %>%
  mutate(
    short_class = abbr_map[compared_class],
    direction = ifelse(pct_diff > 0, "↑", "↓"),
    pct_label = paste0(direction, abs(round(pct_diff)), "%"),
    raw_label = paste0(raw_diff, " aa"),
    combined_label = paste0(pct_label)
  )

final_fig <- ggplot(heatmap_data, aes(x = domain, y = collagen_type, fill = pct_diff)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = combined_label), size = 2.8, lineheight = 0.75) +
  scale_fill_gradient2(
    low = "#e74c3c", mid = "#ecf0f1", high = "#2ecc71",
    midpoint = 0, limits = c(-35, 35),
    name = "% Change\nvs. Hydrozoans"
  ) +
  labs(
    x = "Domain",
    y = "Collagen Type"
  ) +
  facet_wrap(~short_class, ncol = 2, strip.position = "top") +
  theme_basic() +
  theme(
    legend.position = c(0.78, 0.2),
    legend.justification = c("center", "center"),
    legend.direction = "vertical",
    legend.box.margin = margin(-10, -10, -10, -10),
    legend.background = element_rect(fill = alpha("white", 0.8)),
    legend.key.width = unit(0.8, "cm"),
    legend.key.height = unit(3, "cm"),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.margin = margin(10, 40, 10, 10)
  ) +
  guides(fill = guide_colorbar(
    barwidth = 0.5,
    barheight = 8,
    title.position = "top",
    label.position = "right"
  ))
final_fig