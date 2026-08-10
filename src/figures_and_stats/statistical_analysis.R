library(tidyverse)
library(Biostrings)
library(FSA)
library(ggstatsplot)
library(statsExpressions)
library(ggpubr)
library(pheatmap)
library(flextable)
setwd("Documents/year_5/rp2")

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

col1_data <- hcol1_domain_stats %>% filter(domain_name == "Col1")
kruskal.test(domain_prop ~ class, data = col1_data)
col1 <- ggbetweenstats(
  data = col1_data,
  x = class,
  y = domain_prop,
  type = "nonparametric",
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  stats.test = FALSE,
  ggtheme = theme_basic()
) +
  labs(title = "Col1") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14),
    axis.text.y = element_text(size = 10)
  )
col1 <- col1 +
  labs(subtitle = NULL) +
  theme(plot.subtitle = element_blank())
col1

col2_data <- hcol1_domain_stats %>% filter(domain_name == "Col2")
kruskal.test(domain_prop ~ class, data = col2_data)

colfi_data <- hcol1_domain_stats %>% filter(domain_name == "COLFI")
kruskal.test(domain_prop ~ class, data = colfi_data)
colfi <- ggbetweenstats(
  data = colfi_data,
  x = class,
  y = domain_prop,
  type = "nonparametric",
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  stats.test = FALSE,
  ggtheme = theme_basic()
) +
  labs(title = "COLFI") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14),
    axis.text.y = element_text(size = 10)
  )
colfi <- colfi +
  labs(subtitle = NULL) +
  theme(plot.subtitle = element_blank())
colfi

sp_data <- hcol1_domain_stats %>% filter(domain_name == "Signal_Peptide")
kruskal.test(domain_prop ~ class, data = sp_data)

combined_plots <- combine_plots(
  list(col1, colfi),
  annotation.args = list(
    title = "Domain Property Differences Across Classes",
    caption = "Kruskal-Wallis tests: Col1 χ²=18.81, p=0.0003; COLFI χ²=30.08, p=1.33e-6"
  )
) +
  theme(plot.margin = margin(20, 20, 20, 20))
combined_plots

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

col1_data <- hcol2_domain_stats %>% filter(domain_name == "Col1")
kruskal.test(domain_prop ~ class, data = col1_data)
col1 <- ggbetweenstats(
  data = col1_data,
  x = class,
  y = domain_prop,
  type = "nonparametric",
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  stats.test = FALSE,
  ggtheme = theme_basic()
) +
  labs(title = "Col1") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14),
    axis.text.y = element_text(size = 10)
  )
col1 <- col1 +
  labs(subtitle = NULL) +
  theme(plot.subtitle = element_blank())
col1

col2_data <- hcol2_domain_stats %>% filter(domain_name == "Col2")
kruskal.test(domain_prop ~ class, data = col2_data)

colfi_data <- hcol2_domain_stats %>% filter(domain_name == "COLFI")
kruskal.test(domain_prop ~ class, data = colfi_data)
colfi <- ggbetweenstats(
  data = colfi_data,
  x = class,
  y = domain_prop,
  type = "nonparametric",
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  stats.test = FALSE,
  ggtheme = theme_basic()
) +
  labs(title = "COLFI") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14),
    axis.text.y = element_text(size = 10)
  )
colfi <- colfi +
  labs(subtitle = NULL) +
  theme(plot.subtitle = element_blank())
colfi

sp_data <- hcol2_domain_stats %>% filter(domain_name == "Signal_Peptide")
kruskal.test(domain_prop ~ class, data = sp_data)

wap_data <- hcol2_domain_stats %>% filter(domain_name == "WAP")
kruskal.test(domain_prop ~ class, data = wap_data)
wap <- ggbetweenstats(
  data = wap_data,
  x = class,
  y = domain_prop,
  type = "nonparametric",
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  stats.test = FALSE,
  ggtheme = theme_basic()
) +
  labs(title = "WAP") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14),
    axis.text.y = element_text(size = 10)
  )
wap <- wap +
  labs(subtitle = NULL) +
  theme(plot.subtitle = element_blank())
wap

combined_plots <- combine_plots(
  list(col1, colfi, wap),
  annotation.args = list(
    title = "Domain Property Differences Across Classes",
    caption = "Kruskal-Wallis tests: Col1 χ²=19.38, p=0.0003; COLFI χ²=19.26, p=2.41e-04; WAP χ²=15.21, p=0.00165"
  )
) +
  theme(plot.margin = margin(20, 20, 20, 20))
combined_plots

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

col1_data <- hcol3_domain_stats %>% filter(domain_name == "Col1")
kruskal.test(domain_prop ~ class, data = col1_data)
col1 <- ggbetweenstats(
  data = col1_data,
  x = class,
  y = domain_prop,
  type = "nonparametric",
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  stats.test = FALSE,
  ggtheme = theme_basic()
) +
  labs(title = "Col1") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14),
    axis.text.y = element_text(size = 10)
  )
col1 <- col1 +
  labs(subtitle = NULL) +
  theme(plot.subtitle = element_blank())
col1

col2_data <- hcol3_domain_stats %>% filter(domain_name == "Col2")
kruskal.test(domain_prop ~ class, data = col2_data)

colfi_data <- hcol3_domain_stats %>% filter(domain_name == "COLFI")
kruskal.test(domain_prop ~ class, data = colfi_data)
colfi <- ggbetweenstats(
  data = colfi_data,
  x = class,
  y = domain_prop,
  type = "nonparametric",
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  stats.test = FALSE,
  ggtheme = theme_basic()
) +
  labs(title = "COLFI") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14),
    axis.text.y = element_text(size = 10)
  )
colfi <- colfi +
  labs(subtitle = NULL) +
  theme(plot.subtitle = element_blank())
colfi

sp_data <- hcol3_domain_stats %>% filter(domain_name == "Signal_Peptide")
kruskal.test(domain_prop ~ class, data = sp_data)
sp <- ggbetweenstats(
  data = sp_data,
  x = class,
  y = domain_prop,
  type = "nonparametric",
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  stats.test = FALSE,
  ggtheme = theme_basic()
) +
  labs(title = "Signal Peptide") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14),
    axis.text.y = element_text(size = 10)
  )
sp <- sp +
  labs(subtitle = NULL) +
  theme(plot.subtitle = element_blank())
sp

wap_data <- hcol3_domain_stats %>% filter(domain_name == "WAP")
kruskal.test(domain_prop ~ class, data = wap_data)

vwa_data <- hcol3_domain_stats %>% filter(domain_name == "VWA")
kruskal.test(domain_prop ~ class, data = vwa_data)

combined_plots <- combine_plots(
  list(col1, colfi, sp),
  annotation.args = list(
    title = "Domain Property Differences Across Classes",
    caption = "Kruskal-Wallis tests: Col1 χ²=11.82, p=0.00801; COLFI χ²=10.32, p=0.0161; SP χ²=8.69, p=0.0337"
  )
) +
  theme(plot.margin = margin(20, 20, 20, 20))
combined_plots

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

col1_data <- hcol4_domain_stats %>% filter(domain_name == "Col1")
kruskal.test(domain_prop ~ class, data = col1_data)
col1 <- ggbetweenstats(
  data = col1_data,
  x = class,
  y = domain_prop,
  type = "nonparametric",
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  stats.test = FALSE,
  ggtheme = theme_basic()
) +
  labs(title = "Col1") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14),
    axis.text.y = element_text(size = 10)
  )
col1 <- col1 +
  labs(subtitle = NULL) +
  theme(plot.subtitle = element_blank())
col1

col2_data <- hcol4_domain_stats %>% filter(domain_name == "Col2")

c4_data <- hcol4_domain_stats %>% filter(domain_name == "C4")
kruskal.test(domain_prop ~ class, data = c4_data)

sp_data <- hcol4_domain_stats %>% filter(domain_name == "Signal_Peptide")
kruskal.test(domain_prop ~ class, data = sp_data)
sp <- ggbetweenstats(
  data = sp_data,
  x = class,
  y = domain_prop,
  type = "nonparametric",
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  stats.test = FALSE,
  ggtheme = theme_basic()
) +
  labs(title = "Signal Peptide") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14),
    axis.text.y = element_text(size = 10)
  )
sp <- sp +
  labs(subtitle = NULL) +
  theme(plot.subtitle = element_blank())
sp

combined_plots <- combine_plots(
  list(col1, sp),
  annotation.args = list(
    title = "Domain Property Differences Across Classes",
    caption = "Kruskal-Wallis tests: Col1 χ²=8.25, p=0.0393; SP χ²=15.57, p=0.0014"
  )
) +
  theme(plot.margin = margin(20, 20, 20, 20))
combined_plots

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

col1_data <- hcol5_domain_stats %>% filter(domain_name == "Col1")
kruskal.test(domain_prop ~ class, data = col1_data)
col1 <- ggbetweenstats(
  data = col1_data,
  x = class,
  y = domain_prop,
  type = "nonparametric",
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  stats.test = FALSE,
  ggtheme = theme_basic()
) +
  labs(title = "Col1") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14),
    axis.text.y = element_text(size = 10)
  )
col1 <- col1 +
  labs(subtitle = NULL) +
  theme(plot.subtitle = element_blank())
col1

col2_data <- hcol5_domain_stats %>% filter(domain_name == "Col2")
kruskal.test(domain_prop ~ class, data = col2_data)
col2 <- ggbetweenstats(
  data = col2_data,
  x = class,
  y = domain_prop,
  type = "nonparametric",
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  stats.test = FALSE,
  ggtheme = theme_basic()
) +
  labs(title = "Col2") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14),
    axis.text.y = element_text(size = 10)
  )
col2 <- col2 +
  labs(subtitle = NULL) +
  theme(plot.subtitle = element_blank())
col2

colfi_data <- hcol5_domain_stats %>% filter(domain_name == "COLFI")
kruskal.test(domain_prop ~ class, data = colfi_data)
colfi <- ggbetweenstats(
  data = colfi_data,
  x = class,
  y = domain_prop,
  type = "nonparametric",
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  stats.test = FALSE,
  ggtheme = theme_basic()
) +
  labs(title = "COLFI") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14),
    axis.text.y = element_text(size = 10)
  )
colfi <- colfi +
  labs(subtitle = NULL) +
  theme(plot.subtitle = element_blank())
colfi

sp_data <- hcol5_domain_stats %>% filter(domain_name == "Signal_Peptide")
kruskal.test(domain_prop ~ class, data = sp_data)
sp <- ggbetweenstats(
  data = sp_data,
  x = class,
  y = domain_prop,
  type = "nonparametric",
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  stats.test = FALSE,
  ggtheme = theme_basic()
) +
  labs(title = "Signal Peptide") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14),
    axis.text.y = element_text(size = 10)
  )
sp <- sp +
  labs(subtitle = NULL) +
  theme(plot.subtitle = element_blank())
sp

combined_plots <- combine_plots(
  list(col1, col2, colfi, sp),
  annotation.args = list(
    title = "Domain Property Differences Across Classes",
    caption = "Kruskal-Wallis tests: Col1 χ²=12.70, p=0.00037; Col2 χ²=5.981, p=0.0145; COLFI χ²=6.244, p=0.0125; SP χ²=10.01, p=0.00152"
  )
) +
  theme(plot.margin = margin(20, 20, 20, 20))
combined_plots

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

col1_data <- hcol6_domain_stats %>% filter(domain_name == "Col1")
kruskal.test(domain_prop ~ class, data = col1_data)
col1 <- ggbetweenstats(
  data = col1_data,
  x = class,
  y = domain_prop,
  type = "nonparametric",
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  stats.test = FALSE,
  ggtheme = theme_basic()
) +
  labs(title = "Col1") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14),
    axis.text.y = element_text(size = 10)
  )
col1 <- col1 +
  labs(subtitle = NULL) +
  theme(plot.subtitle = element_blank())
col1

col2_data <- hcol6_domain_stats %>% filter(domain_name == "Col2")
kruskal.test(domain_prop ~ class, data = col2_data)

c4_data <- hcol6_domain_stats %>% filter(domain_name == "C4")
kruskal.test(domain_prop ~ class, data = c4_data)
c4 <- ggbetweenstats(
  data = c4_data,
  x = class,
  y = domain_prop,
  type = "nonparametric",
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  stats.test = FALSE,
  ggtheme = theme_basic()
) +
  labs(title = "COLFI") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14),
    axis.text.y = element_text(size = 10)
  )
c4 <- c4 +
  labs(subtitle = NULL) +
  theme(plot.subtitle = element_blank())
c4

sp_data <- hcol6_domain_stats %>% filter(domain_name == "Signal_Peptide")
kruskal.test(domain_prop ~ class, data = sp_data)
sp <- ggbetweenstats(
  data = sp_data,
  x = class,
  y = domain_prop,
  type = "nonparametric",
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  stats.test = FALSE,
  ggtheme = theme_basic()
) +
  labs(title = "Signal Peptide") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14),
    axis.text.y = element_text(size = 10)
  )
sp <- sp +
  labs(subtitle = NULL) +
  theme(plot.subtitle = element_blank())
sp

combined_plots <- combine_plots(
  list(col1, c4, sp),
  annotation.args = list(
    title = "Domain Property Differences Across Classes",
    caption = "Kruskal-Wallis tests: Col1 χ²=20.8, p=0.00012; C4 χ²=18.83, p=0.0003; SP χ²=12.51, p=0.0058"
  )
) +
  theme(plot.margin = margin(20, 20, 20, 20))
combined_plots

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

col1_data <- hcol7_domain_stats %>% filter(domain_name == "Col1")
kruskal.test(domain_prop ~ class, data = col1_data)
col1 <- ggbetweenstats(
  data = col1_data,
  x = class,
  y = domain_prop,
  type = "nonparametric",
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  stats.test = FALSE,
  ggtheme = theme_basic()
) +
  labs(title = "Col1") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14),
    axis.text.y = element_text(size = 10)
  )
col1 <- col1 +
  labs(subtitle = NULL) +
  theme(plot.subtitle = element_blank())
col1

col2_data <- hcol7_domain_stats %>% filter(domain_name == "Col2")

colfi_data <- hcol7_domain_stats %>% filter(domain_name == "COLFI")
kruskal.test(domain_prop ~ class, data = colfi_data)

sp_data <- hcol7_domain_stats %>% filter(domain_name == "Signal_Peptide")
kruskal.test(domain_prop ~ class, data = sp_data)

tspn_data <- hcol7_domain_stats %>% filter(domain_name == "TSPN")
kruskal.test(domain_prop ~ class, data = tspn_data)
tspn <- ggbetweenstats(
  data = tspn_data,
  x = class,
  y = domain_prop,
  type = "nonparametric",
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  stats.test = FALSE,
  ggtheme = theme_basic()
) +
  labs(title = "TSPN") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14),
    axis.text.y = element_text(size = 10)
  )
tspn <- tspn +
  labs(subtitle = NULL) +
  theme(plot.subtitle = element_blank())
tspn

combined_plots <- combine_plots(
  list(col1, tspn),
  annotation.args = list(
    title = "Domain Property Differences Across Classes",
    caption = "Kruskal-Wallis tests: Col1 χ²=17.39, p=0.00059; TSPN χ²=10.65, p=0.0138"
  )
) +
  theme(plot.margin = margin(20, 20, 20, 20))
combined_plots

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

col1_data <- hcol8_domain_stats %>% filter(domain_name == "Col1")
kruskal.test(domain_prop ~ class, data = col1_data)
col1 <- ggbetweenstats(
  data = col1_data,
  x = class,
  y = domain_prop,
  type = "nonparametric",
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  stats.test = FALSE,
  ggtheme = theme_basic()
) +
  labs(title = "Col1") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14),
    axis.text.y = element_text(size = 10)
  )
col1 <- col1 +
  labs(subtitle = NULL) +
  theme(plot.subtitle = element_blank())
col1

col2_data <- hcol8_domain_stats %>% filter(domain_name == "Col2")
kruskal.test(domain_prop ~ class, data = col2_data)

colfi_data <- hcol8_domain_stats %>% filter(domain_name == "COLFI")
kruskal.test(domain_prop ~ class, data = colfi_data)
colfi <- ggbetweenstats(
  data = colfi_data,
  x = class,
  y = domain_prop,
  type = "nonparametric",
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  stats.test = FALSE,
  ggtheme = theme_basic()
) +
  labs(title = "COLFI") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14),
    axis.text.y = element_text(size = 10)
  )
colfi <- colfi +
  labs(subtitle = NULL) +
  theme(plot.subtitle = element_blank())
colfi

sp_data <- hcol8_domain_stats %>% filter(domain_name == "Signal_Peptide")
kruskal.test(domain_prop ~ class, data = sp_data)
sp <- ggbetweenstats(
  data = sp_data,
  x = class,
  y = domain_prop,
  type = "nonparametric",
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  stats.test = FALSE,
  ggtheme = theme_basic()
) +
  labs(title = "Signal Peptide") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14),
    axis.text.y = element_text(size = 10)
  )
sp <- sp +
  labs(subtitle = NULL) +
  theme(plot.subtitle = element_blank())
sp

combined_plots <- combine_plots(
  list(col1, colfi, sp),
  annotation.args = list(
    title = "Domain Property Differences Across Classes",
    caption = "Kruskal-Wallis tests: Col1 χ²=14.52, p=0.00227; COLFI χ²=11.41, p=0.00972; SP χ²=10.37, p=0.00157"
  )
) +
  theme(plot.margin = margin(20, 20, 20, 20))
combined_plots

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
          d1 <- (domain_data %>% filter(class == "Hydrozoans"))$domain_prop
          d2 <- (domain_data %>% filter(class == other_class))$domain_prop

          wtest <- wilcox.test(d2, d1, exact = FALSE)

          med1 <- median(d1)
          med2 <- median(d2)
          pct_diff <- ((med2 - med1) / med1) * 100

          tibble(
            domain = dom,
            compared_class = other_class,
            baseline_class = "Hydrozoans",
            pct_diff = pct_diff,
            abs_pct = abs(pct_diff),
            p.raw = wtest$p.value,
            p.adj = wtest$p.value
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
          pct_diff = numeric(),
          abs_pct = numeric(),
          p.raw = numeric(),
          p.adj = numeric(),
          p.signif = character()
        )
      }
    } else {
      tibble(
        domain = dom,
        compared_class = character(),
        baseline_class = character(),
        pct_diff = numeric(),
        abs_pct = numeric(),
        p.raw = numeric(),
        p.adj = numeric(),
        p.signif = character()
      )
    }
  })

  bind_rows(results) %>%
    mutate(collagen_type = collagen_type)
}

all_vs_hydrozoa <- bind_rows(
  get_pairwise_vs_hydrozoa(hcol1_domain_stats, "Hcol1", c("Col1", "Col2", "COLFI", "Signal_Peptide")),
  get_pairwise_vs_hydrozoa(hcol2_domain_stats, "Hcol2", c("Col1", "Col2", "COLFI", "WAP")),
  get_pairwise_vs_hydrozoa(hcol3_domain_stats, "Hcol3", c("Col1", "Col2", "COLFI", "Signal_Peptide", "WAP")),
  get_pairwise_vs_hydrozoa(hcol4_domain_stats, "Hcol4", c("Col1", "Col2", "C4", "Signal_Peptide")),
  get_pairwise_vs_hydrozoa(hcol5_domain_stats, "Hcol5", c("Col1", "Col2", "COLFI", "Signal_Peptide")),
  get_pairwise_vs_hydrozoa(hcol6_domain_stats, "Hcol6", c("Col1", "Col2", "C4", "Signal_Peptide")),
  get_pairwise_vs_hydrozoa(hcol7_domain_stats, "Hcol7", c("Col1", "Col2", "COLFI", "Signal_Peptide", "TSPN")),
  get_pairwise_vs_hydrozoa(hcol8_domain_stats, "Hcol8", c("Col1", "Col2", "COLFI", "Signal_Peptide"))
)

sig_vs_hydrozoa <- all_vs_hydrozoa %>%
  filter(p.adj < 0.05 & !is.na(p.adj))

heatmap_simple <- sig_vs_hydrozoa %>%
  group_by(collagen_type, domain) %>%
  summarise(
    display_class = compared_class[which.max(log_p)],
    pct_vs_hydro = pct_diff[which.max(log_p)],
    log_p = max(log_p),
    sig_stars = p.signif[which.max(log_p)],
    n_comparisons = n(),
    .groups = "drop"
  )

abbr_map <- c(
  "Cubozoans" = "Cubozoan",
  "Hydrozoans" = "Hydrozoan",
  "Scyphozoans" = "Scyphozoan",
  "Staurozoans" = "Staurozoan"
)

heatmap_simple <- heatmap_simple %>%
  mutate(
    short_class = abbr_map[display_class],
    direction = ifelse(pct_vs_hydro > 0, "↑", "↓"),
    label_short = paste0(short_class, "\n", direction, abs(round(pct_vs_hydro)), "%", "\n")
  )

ggplot(heatmap_simple, aes(x = domain, y = collagen_type, fill = pct_vs_hydro)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = label_short), size = 2.5, lineheight = 0.7) +
  scale_fill_gradient2(
    low = "#ff8888", mid = "#ffffff", high = "#88dd88",
    midpoint = 0, limits = c(-30, 30),
    oob = scales::squish,
    name = "% vs\nHydrozoans"
  ) +
  labs(
    x = "Domain", y = "Collagen Type"
  ) +
  theme_basic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
    axis.text.y = element_text(size = 9, face = "bold"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    legend.position = "right"
  )
