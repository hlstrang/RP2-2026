library(DiagrammeR)
library(flextable)
library(officer)
library(tidyverse)

setwd("Documents/year_5/rp2")
flowchart_tsa <- grViz("
digraph pipeline {
  graph [layout = dot, rankdir = TB, nodesep = 0.5, ranksep = 0.4]

  node [fontname = 'Arial', fontsize = 10, shape = box, style = 'filled,rounded',
        color = '#1A365D', fillcolor = '#F0F4F8', fontcolor = '#0F172A', penwidth = 1.5]

  edge [color = '#334155', penwidth = 1.2, arrowsize = 0.8]

  input1 [label = 'Medusozoan TSA Datasets\n(NCBI Database)', fillcolor = '#E2E8F0', fontcolor = '#0F172A']
  input2 [label = 'Reference Hcol Sequences\n(Known Collagenome)', fillcolor = '#E2E8F0', fontcolor = '#0F172A']

  p1 [label = 'Local BLAST Database\nConstruction (makeblastdb)']
  p2 [label = 'Homology Search\n(tblastn)']
  p3 [label = 'Nucleotide Sequence Extraction\n(Contig Retrieval)']
  p4 [label = 'Open Reading Frame Translation\n(Longest ORF)']
  p5 [label = 'Domain Architecture Verification\n& Functional Annotation']

  t1 [label = 'tblastn flags:\n-seg no\n-soft_masking false', shape = hexagon, fillcolor = '#FEF3C7', color = '#D97706', fontsize = 8]
  t2 [label = 'Custom Script:\npull_fasta.py', shape = hexagon, fillcolor = '#FEF3C7', color = '#D97706', fontsize = 8]
  t3 [label = 'Tool:\norfipy\n(-procs 2)', shape = hexagon, fillcolor = '#FEF3C7', color = '#D97706', fontsize = 8]
  t4 [label = 'Tool:\nInterProScan\n(-dp -goterms -iprlookup)\nSignalP', shape = hexagon, fillcolor = '#FEF3C7', color = '#D97706', fontsize = 8]

  { rank = same; input1; input2 }
  { rank = same; p2; t1 }
  { rank = same; p3; t2 }
  { rank = same; p4; t3 }
  { rank = same; p5; t4 }

  input1 -> p1
  p1 -> p2
  input2 -> p2
  p2 -> p3
  p3 -> p4
  p4 -> p5

  edge [style = dashed, color = '#D97706', arrowhead = none]
  t1 -> p2
  t2 -> p3
  t3 -> p4
  t4 -> p5
}
")

flowchart_tsa

flowchart_wgs <- grViz("
digraph pipeline {

  graph [layout = dot, rankdir = TB, nodesep = 0.5, ranksep = 0.4]

  node [fontname = 'Arial', fontsize = 10, shape = box, style = 'filled,rounded',
        color = '#1A365D', fillcolor = '#F0F4F8', fontcolor = '#0F172A', penwidth = 1.5]

  edge [color = '#334155', penwidth = 1.2, arrowsize = 0.8]

  input1 [label = 'Medusozoan WGS Datasets\n(NCBI Database)', fillcolor = '#E2E8F0', fontcolor = '#0F172A']
  input2 [label = 'Reference Hcol Sequences\n(Known Collagenome)', fillcolor = '#E2E8F0', fontcolor = '#0F172A']

  p1 [label = 'Local BLAST Database\nConstruction (makeblastdb)']
  p2 [label = 'Homology Search Using N- and C-termini as Queries\n(tblastn)']
  p3 [label = 'Nucleotide Sequence Extraction\n(Contig Retrieval)']
  p4 [label = 'Predicted Protein Sequence']
  p5 [label = 'Domain Architecture Verification\n& Functional Annotation']

  t1 [label = 'tblastn flags:\n-seg no\n-soft_masking false', shape = hexagon, fillcolor = '#FEF3C7', color = '#D97706', fontsize = 8]
  t2 [label = 'Custom Script:\ntblastn_analysis.py', shape = hexagon, fillcolor = '#FEF3C7', color = '#D97706', fontsize = 8]
  t3 [label = 'Tool:\nminiprot\n(--aln, --gff, --trans)', shape = hexagon, fillcolor = '#FEF3C7', color = '#D97706', fontsize = 8]
  t4 [label = 'Tool:\nHMMER Domains\nSignalP', shape = hexagon, fillcolor = '#FEF3C7', color = '#D97706', fontsize = 8]

  { rank = same; input1; input2 }
  { rank = same; p2; t1 }
  { rank = same; p3; t2 }
  { rank = same; p4; t3 }
  { rank = same; p5; t4 }

  input1 -> p1
  p1 -> p2
  input2 -> p2
  p2 -> p3
  p3 -> p4
  p4 -> p5

  edge [style = dashed, color = '#D97706', arrowhead = none]
  t1 -> p2
  t2 -> p3
  t3 -> p4
  t4 -> p5
}
")

flowchart_wgs

combined_flowchart <- grViz("
digraph combined_pipeline {
  graph [layout = dot, rankdir = TB, nodesep = 0.6, ranksep = 0.5, compound = true]

  node [fontname = 'Arial', fontsize = 10, shape = box, style = 'filled,rounded',
        color = '#1A365D', fillcolor = '#F0F4F8', fontcolor = '#0F172A', penwidth = 1.5]

  edge [color = '#334155', penwidth = 1.2, arrowsize = 0.8]

  # --- SHARED INPUT ---
  ref_seq [label = 'Reference Hcol Sequences\\n(Known Collagenome)', fillcolor = '#E2E8F0']

  title_tsa [label = 'Transcriptome Shotgun Assembly (TSA) Pipeline',
           shape = plaintext, style = '', color = none, fillcolor = none,
           fontname = 'Arial-Bold', fontsize = 14, fontcolor = '#1A365D']

  title_wgs [label = 'Whole Genome Shotgun (WGS) Pipeline',
           shape = plaintext, style = '', color = none, fillcolor = none,
           fontname = 'Arial-Bold', fontsize = 14, fontcolor = '#1A365D']

  { rank = same; title_tsa; title_wgs }

  # --- WGS PIPELINE ---
  subgraph cluster_wgs_group {
    labeljust = c
    labelloc = t
    fontname = 'Arial-Bold'
    fontsize = 14
    style = ''
    fontcolor = '#1A365D'
    color = none

    input_wgs [label = 'Medusozoan WGS Datasets\\n(NCBI Database)', fillcolor = '#E2E8F0']
    p1_wgs    [label = 'Local BLAST Database\\nConstruction (makeblastdb)']
    p2_wgs    [label = 'Homology Search (N/C-termini)\\n(tblastn)']
    t1_wgs    [label = 'tblastn flags:\\n-seg no\\n-soft_masking false', shape = hexagon, fillcolor = '#FEF3C7', color = '#D97706', fontsize = 8]
    p3_wgs    [label = 'Nucleotide Sequence Extraction\\n(Contig Retrieval)']
    t2_wgs    [label = 'Custom Script:\\ntblastn_analysis.py', shape = hexagon, fillcolor = '#FEF3C7', color = '#D97706', fontsize = 8]
    p4_wgs    [label = 'Spliced Gene Alignment &\\nProtein Translation']
    t3_wgs    [label = 'Tool:\\nminiprot\\n(--aln, --gff, --trans)', shape = hexagon, fillcolor = '#FEF3C7', color = '#D97706', fontsize = 8]

    input_wgs -> p1_wgs
    p1_wgs -> p2_wgs
    p2_wgs -> t1_wgs [style = invis]
    p2_wgs -> p3_wgs
    p3_wgs -> t2_wgs [style = invis]
    p3_wgs -> p4_wgs
    p4_wgs -> t3_wgs [style = invis]
  }

  # --- TSA PIPELINE ---
  subgraph cluster_tsa_group {
    labeljust = c
    labelloc = t
    fontname = 'Arial-Bold'
    fontsize = 14
    style = ''
    fontcolor = '#1A365D'
    color = none

    input_tsa [label = 'Medusozoan TSA Datasets\\n(NCBI Database)', fillcolor = '#E2E8F0']
    p1_tsa    [label = 'Local BLAST Database\\nConstruction (makeblastdb)']

    t1_tsa    [label = 'tblastn flags:\\n-seg no\\n-soft_masking false', shape = hexagon, fillcolor = '#FEF3C7', color = '#D97706', fontsize = 8]
    p2_tsa    [label = 'Homology Search\\n(tblastn)']

    t2_tsa    [label = 'Custom Script:\\npull_fasta.py', shape = hexagon, fillcolor = '#FEF3C7', color = '#D97706', fontsize = 8]
    p3_tsa    [label = 'Nucleotide Sequence Extraction\\n(Contig Retrieval)']

    t3_tsa    [label = 'Tool:\\norfiipy\\n(-procs 2)', shape = hexagon, fillcolor = '#FEF3C7', color = '#D97706', fontsize = 8]
    p4_tsa    [label = 'Open Reading Frame Translation\\n(Longest ORF)']

    input_tsa -> p1_tsa

    t1_tsa -> p2_tsa [style = invis]
    p1_tsa -> p2_tsa

    t2_tsa -> p3_tsa [style = invis]
    p2_tsa -> p3_tsa

    t3_tsa -> p4_tsa [style = invis]
    p3_tsa -> p4_tsa
  }
  { rank = same; t1_tsa; p2_tsa }
  { rank = same; t2_tsa; p3_tsa }
  { rank = same; t3_tsa; p4_tsa }

  title_tsa -> input_tsa [style = invis, weight = 10]
  title_wgs -> input_wgs [style = invis, weight = 10]

  p_shared [label = 'Domain Architecture Verification\\n& Functional Annotation', fillcolor = '#DBEAFE', penwidth = 2.0]
  t_shared [label = 'Tools:\\nHMMER (Domain Search)\\nSignalP', shape = hexagon, fillcolor = '#FEF3C7', color = '#D97706', fontsize = 8]

  { rank = same; t1_wgs; p2_wgs }
  { rank = same; t2_wgs; p3_wgs }
  { rank = same; t3_wgs; p4_wgs }
  { rank = same; p_shared; t_shared }

  ref_seq -> p2_wgs
  ref_seq -> p2_tsa
  p4_wgs -> p_shared
  p4_tsa -> p_shared

  edge [style = dashed, color = '#D97706', arrowhead = none, constraint = false]
  t1_wgs -> p2_wgs
  t2_wgs -> p3_wgs
  t3_wgs -> p4_wgs
  t1_tsa -> p2_tsa
  t2_tsa -> p3_tsa
  t3_tsa -> p4_tsa
  t_shared -> p_shared
}
")

combined_flowchart

collagen_data <- data.frame(
  `Collagen Type` = c("Hcol1", "Hcol2", "Hcol3", "Hcol4", "Hcol5", "Hcol6", "Hcol7", "Hcol8"),
  `Non-collagenous Domains` = c(
    "No other domains",
    "1 or 2 WAP domains",
    ">1 WAP and >1 VWA domains",
    "No other domains",
    "No other domains",
    ">1 WAP and >1 VWA domains",
    "TSPN domain",
    "No other domains"
  ),
  `Triple Helix Structure` = c(
    "Col1 length = 1010–1030 aa; Col2 length = 50–65 aa; Minor interruption between",
    "Col1 and Col2 present",
    "Multiple segmented sequences",
    "Multiple segmented sequences",
    "Col1 length = 1020–1035 aa; Col2 length = 50–60 aa; Longer interruption between",
    "Multiple segmented sequences",
    "Col1 and Col2 present",
    "Col1 length = ~1000 aa with several short interruptions; No Col2; Longer gap between SP and Col1"
  ),
  `C-terminus` = c("COLFI", "COLFI", "COLFI", "C4", "COLFI", "C4", "COLFI", "COLFI"),
  stringsAsFactors = FALSE
)
colnames(collagen_data) <- c("Collagen Type", "Non-collagenous Domains", "Triple Helix Structure", "C-terminus")
caption_style <- fp_text(color = "#1A365D", font.size = 12)
footer_style <- fp_text(color = "#0F172A")
ft_table <- flextable(collagen_data) %>%
  add_footer_lines("Col1: Major triple-helical domain (>60 amino acids); Col2: Minor triple-helical domain (<60 amino acids). WAP: Whey Acidic Protein domain; VWA: Von Willebrand A domain; TSPN: Thrombospondin N-terminal domain.") %>%
  align(align = "center", part = "all") %>%
  autofit() %>%
  theme_zebra()
ft_table


