#!/bin/bash
#SBATCH --job-name=hcol_extraction
#SBATCH --output=exon_hcol.out
#SBATCH --error=exon_hcol.err
#SBATCH --time=02:00:00
#SBATCH -p serial

module purge
module load apps/binapps/conda/miniforge3/25.9.1
conda activate /mnt/iusers01/fatpou01/bmh01/msc-bioinf-2025-2026/u78012he/.conda/envs/tblastn

QUERY_FASTA="h_viridissima/ref_seqs/hcol8.fasta"
INPUT_FASTA="h_viridissima/genomic_fragments/hcol8/h_viridissima_hcol8_fragment.fasta"
PREDICTED_PROTEINS="h_viridissima/predicted/hcol8/h_viridissima_hcol8_prediction.txt"
PREDICTED_PROTEINS_FASTA="h_viridissima/predicted/hcol8/h_viridissima_hcol8_prediction.fasta"
EXONS_TSV="h_viridissima/predicted/hcol8/h_viridissima_hcol8_exons.tsv"

echo "=== Starting miniprot extraction workflow ==="
echo "Query: ${QUERY_FASTA}"
echo "Target: ${INPUT_FASTA}"
echo "Running miniprot protein-to-genome alignment..."

miniprot \
  --aln \
  --gff \
  --trans \
  "${INPUT_FASTA}" \
  "${QUERY_FASTA}" > "${PREDICTED_PROTEINS}"

python3 parse_miniprot_output.py "${PREDICTED_PROTEINS}" "${PREDICTED_PROTEINS_FASTA}" "${EXONS_TSV}"

echo "=== Workflow complete ==="
echo "Cleanly spliced proteins saved to: ${PREDICTED_PROTEINS_FASTA}"
echo "Saved exon boundaries to: ${EXONS_TSV}"