#!/bin/bash --login
#SBATCH -p multicore_small
#SBATCH -t 2-0
#SBATCH -o extraction.out
#SBATCH -e extraction.err
#SBATCH -c 8

module purge
module load apps/binapps/blast/2.17.0

HCOL_TYPE="hcol8"

echo "=========================================="
echo "STEP 1: Map reference protein domains"
echo "=========================================="
python3 domain_mapping.py "${HCOL_TYPE}" --sp_start 1 --sp_end 21

echo ""
echo "=========================================="
echo "STEP 2: Generate terminal anchor fragments"
echo "=========================================="
python3 extract_fragments_from_json.py "${HCOL_TYPE}"

echo ""
echo "=========================================="
echo "STEP 3: Run consensus TBLASTN"
echo "=========================================="
tblastn \
  -query blast_fragments/"${HCOL_TYPE}"_consensus_anchors_combined.fasta \
  -db db/r_luteum_wgs_db \
  -outfmt "6 qseqid sseqid pident length qcovs qstart qend sstart send evalue bitscore" \
  -evalue 1 \
  -seg no \
  -soft_masking false \
  -num_threads 8 \
> luteum_"${HCOL_TYPE}"_query.tsv

sed -i '1i qseqid\tsseqid\tpident\tlength\tqcovs\tqstart\tqend\tsstart\tsend\tevalue\tbitscore' luteum_"${HCOL_TYPE}"_query.tsv

echo ""
echo "=========================================="
echo "STEP 4: Analyze consensus scoring"
echo "=========================================="
python3 tblastn_analysis.py "${HCOL_TYPE}" --cushion 50000

echo ""
echo "Extraction complete."