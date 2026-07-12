#!/bin/bash
#SBATCH --job-name=protein_prediction
#SBATCH --output=exon_hcol.out
#SBATCH --error=exon_hcol.err
#SBATCH --time=02:00:00
#SBATCH -p serial

module purge
module load apps/binapps/conda/miniforge3/25.9.1
conda activate /mnt/iusers01/fatpou01/bmh01/msc-bioinf-2025-2026/u78012he/.conda/envs/tblastn

HCOL_TYPE="hcol5"
QUERY_FASTA="r_esculentum_${HCOL_TYPE}.fasta"
INPUT_FASTA1="r_luteum_${HCOL_TYPE}_nterm.fasta"
INPUT_FASTA2="r_luteum_${HCOL_TYPE}_cterm.fasta"
PREDICTED_PROTEINS1="r_luteum_${HCOL_TYPE}_prediction1.txt"
PREDICTED_PROTEINS2="r_luteum_${HCOL_TYPE}_prediction2.txt"
PREDICTED_PROTEINS_FASTA1="r_luteum_${HCOL_TYPE}_prediction1.fasta"
PREDICTED_PROTEINS_FASTA2="r_luteum_${HCOL_TYPE}_prediction2.fasta"
FINAL_OUTPUT="r_luteum_${HCOL_TYPE}_prediction.fasta"

echo "=== Starting miniprot extraction workflow ==="
echo "Query: ${QUERY_FASTA}"
echo "Target 1: ${INPUT_FASTA1}"
echo "Running miniprot protein-to-genome alignment on part 1..."

miniprot \
  -t 1 \
  --trans \
  "${INPUT_FASTA1}" \
  "${QUERY_FASTA}" > "${PREDICTED_PROTEINS1}"

echo "Target 2: ${INPUT_FASTA2}"
echo "Running miniprot protein-to-genome alignment on part 2..."

miniprot \
  -t 1 \
  --trans \
  "${INPUT_FASTA2}" \
  "${QUERY_FASTA}" > "${PREDICTED_PROTEINS2}"

echo "Parsing Target 1..."
python3 parse_miniprot_output.py "${PREDICTED_PROTEINS1}" "${PREDICTED_PROTEINS_FASTA1}"
echo "Parsing Target 2..."
python3 parse_miniprot_output.py "${PREDICTED_PROTEINS2}" "${PREDICTED_PROTEINS_FASTA2}"

echo "Combining fragments..."
/mnt/iusers01/fatpou01/bmh01/msc-bioinf-2025-2026/u78012he/.conda/envs/tblastn/bin/python3 stitch_proteins.py "${QUERY_FASTA}" "${PREDICTED_PROTEINS_FASTA1}" "${PREDICTED_PROTEINS_FASTA2}" "${FINAL_OUTPUT}"

echo "=== Workflow complete ==="