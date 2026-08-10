library(ggtree)
library(tidyverse)
library(ape)
library(pegas)
library(pals)
library(ggtext)
library(tidytree)
library(ggforce)
setwd("Documents/year_5/rp2")

tree <- read.tree("phylogeny/hcol1.nwk")
tree <- root(tree, outgroup = "Radianthus crispa")
tipcategories = read.csv("phylogeny/metadata.csv",
                         sep = ",",
                         header = TRUE,
                         stringsAsFactors = FALSE)
colnames(tipcategories)[colnames(tipcategories) == "species"] <- "label"

tree$tip.label <- gsub("^['\"]|['\"]$", "", tree$tip.label)
tipcategories$label <- gsub("^['\"]|['\"]$", "", tipcategories$label)
tree$tip.label <- gsub("_", " ", tree$tip.label)
tree$tip.label <- trimws(tree$tip.label)
tipcategories <- tipcategories %>% select(label, everything())

node_cubo    <- MRCA(tree, tipcategories$label[tipcategories$class == "Cubozoans"])
node_stauro  <- MRCA(tree, tipcategories$label[tipcategories$class == "Staurozoans"])
node_hydro   <- MRCA(tree, tipcategories$label[tipcategories$class == "Hydrozoans"])
node_scypho1 <- 1
node_scypho2 <- 41
node_scypho3 <- 44
node_scypho4 <- 8
node_scypho5 <- 2

highlight_df <- data.frame(
  node = c(node_stauro, node_cubo, node_hydro, node_scypho1, node_scypho2, node_scypho3, node_scypho4, node_scypho5),
  group = c("Staurozoa", "Cubozoa", "Hydrozoa", "Scyphozoa", "Scyphozoa", "Scyphozoa", "Scyphozoa", "Scyphozoa"),
  color = c("purple", "red", "olivedrab", "deepskyblue", "deepskyblue", "deepskyblue", "deepskyblue", "deepskyblue")
)

color_palette <- setNames(
  c("purple", "red", "olivedrab", "deepskyblue"),
  c("Staurozoa", "Cubozoa", "Hydrozoa", "Scyphozoa")
)

p <- ggtree(tree) %<+% tipcategories +
  geom_treescale(fontsize = 3) +
  geom_highlight(
    data = highlight_df,
    aes(node = node, fill = group),
    alpha = 0.15,
    extendto = 3.0
  ) +
  scale_fill_manual(name = "Taxon Group", values = color_palette) +

  geom_tiplab(size = 2.5, align = FALSE, offset = 0.01, color = "black") +
  geom_nodelab(
    aes(label = label),
    hjust = 1.5,
    vjust = -0.4,
    size = 2.5
  )+
  hexpand(0.1) +
  theme(
    legend.position = "top",
    legend.text = element_markdown(),
    legend.key = element_blank(),
    legend.title = element_text(
      hjust = 0.5,
      margin = margin(b = 5)
    )
  )
p
ggsave("phylogeny/hcol1_tree.pdf", p, width = 12, height = 10, units = "in")

