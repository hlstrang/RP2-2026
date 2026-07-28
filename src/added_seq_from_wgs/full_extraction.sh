#!/bin/bash --login
#SBATCH -p multicore
#SBATCH -t 2-0
#SBATCH -o h_viridissima/logs/extraction_%A_%a.out
#SBATCH -e h_viridissima/logs/extraction_%A_%a.err
#SBATCH -c 8
#SBATCH --array=0-7

module purge
module load apps/binapps/blast/2.17.0

## Change species
SPECIES="h_viridissima"

REF_DIR="${SPECIES}/ref_seqs"
BASE_OUT_DIR="${SPECIES}/genomic_fragments"
mkdir -p "$BASE_OUT_DIR"

FILES=(${REF_DIR}/*.fasta)
current_file=${FILES[$SLURM_ARRAY_TASK_ID]}
HCOL_TYPE=$(basename "$current_file" .fasta)

SAMPLE_OUT_DIR="${BASE_OUT_DIR}/${HCOL_TYPE}"
mkdir -p "$SAMPLE_OUT_DIR"

echo "Task $SLURM_ARRAY_TASK_ID processing sample: $HCOL_TYPE"
echo "Output directory: $SAMPLE_OUT_DIR"

SP_INFO="${SPECIES}/hcol_config.txt"

read SP_START SP_END < <(awk -v target="$HCOL_TYPE" '$1 == target {print $2, $3}' "$SP_INFO")

if [ -z "$SP_START" ] || [ -z "$SP_END" ]; then
    echo "WARNING: No config entry found for ${HCOL_TYPE}, using defaults (0, 0)"
    SP_START=0
    SP_END=0
else
    echo "Found config for ${HCOL_TYPE}: sp_start=${SP_START}, sp_end=${SP_END}"
fi

## Update species in script
python3 domain_mapping.py "${HCOL_TYPE}" --sp_start "${SP_START}" --sp_end "${SP_END}"
python3 extract_fragments_from_json.py "${HCOL_TYPE}" --species "${SPECIES}"

## Change db name
tblastn \
    -query blast_fragments/"${HCOL_TYPE}"_consensus_anchors_combined.fasta \
    -db db/hydra_wgs_db \
    -outfmt "6 qseqid sseqid pident length qcovs qstart qend sstart send evalue bitscore" \
    -evalue 1 \
    -seg no \
    -soft_masking false \
    -num_threads 8 \
    > "${SAMPLE_OUT_DIR}/${SPECIES}_${HCOL_TYPE}_query.tsv"

sed -i '1i qseqid\tsseqid\tpident\tlength\tqcovs\tqstart\tqend\tsstart\tsend\tevalue\tbitscore' "${SAMPLE_OUT_DIR}/${SPECIES}_${HCOL_TYPE}_query.tsv"

## Update species in script
python3 tblastn_analysis.py "${HCOL_TYPE}" --output-fasta "${SAMPLE_OUT_DIR}/${SPECIES}_${HCOL_TYPE}_fragment.fasta"