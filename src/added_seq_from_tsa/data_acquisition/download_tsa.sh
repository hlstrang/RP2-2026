#!/bin/bash --login
#SBATCH -p multicore
#SBATCH --job-name=blast-db-build
#SBATCH --output=blast_build.out
#SBATCH --time=2-0
#SBATCH --cpus-per-task=4

module purge
module load apps/binapps/blast/2.17.0

DOWNLOAD_DIR="tsa_fasta_files"
mkdir -p "$DOWNLOAD_DIR"

URL_FILE="verified_urls.txt"

echo "=== Starting Download Phase ==="

while read -r acc url; do
    if [[ "$acc" == "Accession" || -z "$acc" ]]; then continue; fi

    echo "Processing $acc..."

    url=$(echo "$url" | tr -d '\r')
    acc=$(echo "$acc" | tr -d '\r')

    base_name=$(basename "$url")

    nt_file="${base_name}.1.fsa_nt.gz"
    pep_file="${base_name}.1.fsa_pep.gz"

    wget -q -O "${DOWNLOAD_DIR}/${nt_file}" "${url}${nt_file}"

    if [ ! -s "${DOWNLOAD_DIR}/${nt_file}" ]; then
        rm -f "${DOWNLOAD_DIR}/${nt_file}"
        echo "   No nucleotide file found. Trying protein sequence..."
        wget -q -O "${DOWNLOAD_DIR}/${pep_file}" "${url}${pep_file}"
    fi

    if [ -s "${DOWNLOAD_DIR}/${nt_file}" ]; then
        echo "   Successfully downloaded Nucleotide: $nt_file"
    elif [ -s "${DOWNLOAD_DIR}/${pep_file}" ]; then
        echo "   Successfully downloaded Protein: $pep_file"
    else
        echo "   WARNING: Could not download sequences for $acc"
    fi

done < "$URL_FILE"

echo "=== Starting BLAST Database Creation ==="

cd "$DOWNLOAD_DIR" || exit

echo "Decompressing files..."
gunzip -f *.gz

if ls *.fsa_nt 1> /dev/null 2>&1; then
    echo "Combining nucleotide sequences..."
    cat *.fsa_nt > combined_tsa_nucl.fasta

    echo "Building Nucleotide BLAST Database..."
    makeblastdb \
        -in combined_tsa_nucl.fasta \
        -dbtype nucl \
        -title "Combined_TSA_Nucl_DB" \
        -out combined_tsa_nucl_db \
        -parse_seqids
fi

if ls *.fsa_pep 1> /dev/null 2>&1; then
    echo "Combining protein sequences..."
    cat *.fsa_pep > combined_tsa_prot.fasta

    echo "Building Protein BLAST Database..."
    makeblastdb \
        -in combined_tsa_prot.fasta \
        -dbtype prot \
        -title "Combined_TSA_Prot_DB" \
        -out combined_tsa_prot_db \
        -parse_seqids
fi

echo "=== Process Complete ==="