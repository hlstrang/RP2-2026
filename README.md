# RP2 2026 Supplementaries
## A Multi-Omics Collagenome Resource for Marine Invertebrate Biotechnology: Expansion of Cnidarian Structural Protein Sequences
# Candidate ID 14226716, University of Manchester
Supplementary materials and source code

## TSA Extraction
To automate this process, Snakemake was used. However, after implementing the Snakefile, it is important to verify the outputs, especially as there is no option to extract scaffolds other than the top for each species. Needs orfipy and pyhmmer dependencies.

Python scripts:
1. pull_fasta.py
2. domain_architecture.py

## WGS Extraction
This process required tblastn using the N- and C-termini as separate reference sequences, which meant knowing the coordinates of the signal peptide. As this was unreliable using HMMER domains, it was done manually using SignalP, and therefore, was not automated.

Instead, it is run using three main Bash scripts:
1. full_extraction.sh
2. miniprot.sh
3. map_domains.sh

Within these, python scripts are used:
1. domain_mapping.py
2. extract_fragments_from_json.py
3. tblastn_analysis.py
4. parse_miniprot_output.py
5. sequence_mapping.py

## Domain Architecture
The domains were mapped against a cutsom HMM database, with accession numbers shown [here](src/added_seq_from_wgs/hmm_domain_accessions.md). The triple helix regex matching expression was:

    r"(G..){3,}"

Which highlights sequences of a Glycine followed by two other amino acids. To account for small interruptions, a max gap of 9 amino acids was allowed where the triple helix was still combined as one. If there was a gap longer than this threshold, this was counted as two triple helix domains (sometimes adjusted depending on collagen type).

## Statistical Analysis
To determine significant variations in median domain lengths, Wilcoxon signed-ranked tests were used with Hydrozoan as the baseline due to its overrepresentation in the dataset. The R script 'statistical_analysis.R' contains individual domain stats for each collagen type including total domains, domain length, and domain proportion. These summaries are combined to conduct the Wilcoxon signed-ranked tests and p-values are adjusted using the Holm method.

## Phylogenetic Analysis
Hcol1 sequences were extracted from the collagenome and once the newick file was obtained, the R script 'phylogeny.R' was used to create the tree topology.