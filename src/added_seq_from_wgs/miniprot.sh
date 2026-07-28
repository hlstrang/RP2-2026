#!/bin/bash
#SBATCH --job-name=protein_prediction
#SBATCH --output=h_viridissima/logs/miniprot_%A_%a.out
#SBATCH --error=h_viridissima/logs/miniprot_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH -p multicore
#SBATCH -c 8
#SBATCH --array=0-8

module purge
module load apps/binapps/conda/miniforge3/25.9.1
conda activate /mnt/iusers01/fatpou01/bmh01/msc-bioinf-2025-2026/u78012he/.conda/envs/tblastn

## Change species
SPECIES="h_viridissima"

REF_DIR="${SPECIES}/ref_seqs"
FRAGMENT_DIR="${SPECIES}/genomic_fragments"
PREDICTED_DIR="${SPECIES}/predicted"
mkdir -p "$PREDICTED_DIR" "${SPECIES}"

readarray -t FILES < <(printf '%s\n' "${REF_DIR}"/*.fasta | sort)
current_file=${FILES[$SLURM_ARRAY_TASK_ID]}
HCOL_TYPE=$(basename "$current_file" .fasta)

SAMPLE_OUT_DIR="${PREDICTED_DIR}/${HCOL_TYPE}"
mkdir -p "$SAMPLE_OUT_DIR"

QUERY_FASTA="${REF_DIR}/${HCOL_TYPE}.fasta"
INPUT_FASTA="${FRAGMENT_DIR}/${HCOL_TYPE}/${SPECIES}_${HCOL_TYPE}_fragment.fasta"
PREDICTED_PROTEINS="${SAMPLE_OUT_DIR}/${SPECIES}_${HCOL_TYPE}_prediction.txt"
PREDICTED_PROTEINS_FASTA="${SAMPLE_OUT_DIR}/${SPECIES}_${HCOL_TYPE}_prediction.fasta"
EXONS_TSV="${SAMPLE_OUT_DIR}/${SPECIES}_${HCOL_TYPE}_exons.tsv"

echo "=== Starting miniprot extraction workflow ==="
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Gene Target:   ${HCOL_TYPE}"
echo "Query:         ${QUERY_FASTA}"
echo "Target:        ${INPUT_FASTA}"

if [ ! -f "${INPUT_FASTA}" ]; then
    echo "ERROR: Missing genomic fragment fasta: ${INPUT_FASTA}" >&2
    exit 1
fi

echo "Running miniprot protein-to-genome alignment..."

miniprot \
  --aln \
  --gff \
  --trans \
  "${INPUT_FASTA}" \
  "${QUERY_FASTA}" > "${PREDICTED_PROTEINS}"

echo "Parsing miniprot output..."
python3 parse_miniprot_output.py "${PREDICTED_PROTEINS}" "${PREDICTED_PROTEINS_FASTA}" "${EXONS_TSV}"

echo "=== Workflow complete ==="
echo "Cleanly spliced proteins saved to: ${PREDICTED_PROTEINS_FASTA}"
echo "Saved exon boundaries to: ${EXONS_TSV}"