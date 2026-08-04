# RP2 2026 Supplementaries
## A Multi-Omics Collagenome Resource for Marine Invertebrate Biotechnology: Expansion of Cnidarian Structural Protein Sequences
# Candidate ID 14226716, University of Manchester
Supplementary materials and source code

## TSA Extraction
To automate this process, Snakemake was used. However, after implementing the Snakefile, it is important to verify the outputs, especially as there is no option to extract scaffolds other than the top for each species. Needs ofipy and pyhmmer dependencies.

Python scripts:
1. pull_fasta.py
2. domain_architecture.py

## WGS Extraction
This process reuqired tblastn using the N- and C-temrini as separate reference sequences, which meant knowing the coordinates of the signal peptide. As this was unreliable using HMMER domains, it was done manually using SignalP, and therefore, was not automated.

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

