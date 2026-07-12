#!/bin/bash
#SBATCH --job-name=hcol1_extraction
#SBATCH --output=exon_hcol1.out
#SBATCH --error=exon_hcol1.err
#SBATCH --time=02:00:00
#SBATCH -p serial

module purge
module load apps/binapps/conda/miniforge3/25.9.1
conda activate /mnt/iusers01/fatpou01/bmh01/msc-bioinf-2025-2026/u78012he/.conda/envs/tblastn

HCOL_TYPE="hcol1"
QUERY_FASTA="r_esculentum_${HCOL_TYPE}.fasta"
INPUT_FASTA="r_luteum_${HCOL_TYPE}_fragment.fasta"
PREDICTED_PROTEINS="r_luteum_${HCOL_TYPE}_prediction.txt"
PREDICTED_PROTEINS_FASTA="r_luteum_${HCOL_TYPE}_prediction.fasta"

echo "=== Starting miniprot extraction workflow ==="
echo "Query: ${QUERY_FASTA}"
echo "Target: ${INPUT_FASTA}"
echo "Running miniprot protein-to-genome alignment..."

miniprot \
  -t 1 \
  --trans \
  "${INPUT_FASTA}" \
  "${QUERY_FASTA}" > "${PREDICTED_PROTEINS}"

if [ $? -eq 0 ] && [ -s "${PREDICTED_PROTEINS}" ]; then
    RESULT_COUNT=$(grep -c "##STA" "${PREDICTED_PROTEINS}")
    echo "SUCCESS: Generated ${RESULT_COUNT} predicted protein sequences."
else
    echo "[ERROR] Failed to extract proteins or output file is empty."
    exit 1
fi

python3 parse_miniprot_output.py "${PREDICTED_PROTEINS}" "${PREDICTED_PROTEINS_FASTA}"

echo "=== Workflow complete ==="
echo "Cleanly spliced proteins saved to: ${PREDICTED_PROTEINS_FASTA}"