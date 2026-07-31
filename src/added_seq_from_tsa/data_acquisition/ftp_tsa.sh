#!/bin/bash --login
#SBATCH -p serial
#SBATCH --job-name=ftp-url-list
#SBATCH --output=ftp_check.out
#SBATCH --time=2-0

OUTPUT_FILE="verified_urls.txt"

echo -e "Accession\tURL" > "$OUTPUT_FILE"

while read -r acc; do
    acc=$(echo "$acc" | tr -d '\r')

    first2="${acc:0:2}"
    next2="${acc:2:2}"
    downloaded=false

    for num in 01 03 04 05; do
        url="https://ftp.ncbi.nlm.nih.gov/sra/wgs${num}/wgs_aux/${first2}/${next2}/${first2}${next2}01/"

        if curl --output /dev/null --silent --head --fail "$url"; then
            echo -e "${acc}\t${url}" >> "$OUTPUT_FILE"
            downloaded=true
            break
        fi
    done

    if [ "$downloaded" = false ]; then
        echo "ERROR: Could not find ${acc} in any directory." >> ftp_check.out
    fi
done < accessions.txt