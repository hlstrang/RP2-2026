import csv
import re
from Bio import SeqIO
import sys
import os

accession_file = snakemake.input.species_list
fasta_file = snakemake.input.tsa_fasta
output_file = snakemake.output[0]

if os.path.getsize(accession_file) == 0:
    with open(output_file, "w") as out:
        pass
    sys.exit(0)

ACCESSION_RE = re.compile(r'^[A-Za-z]{2,6}\d{2,}\.\d+')

def get_accession(text):
    text = str(text).strip().lstrip(">")
    match = ACCESSION_RE.match(text)
    if match:
        return match.group(0)
    token = text.split()[0] if text.split() else text
    if "|" in token:
        parts = [p for p in token.split("|") if p.strip()]
        if parts:
            token = parts[0]
    return token

accession_set = set()

with open(accession_file, mode="r", newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        if "accession" in row and row["accession"].strip():
            accession_set.add(get_accession(row["accession"]))

print(f"Parsed {len(accession_set)} unique target accession IDs to hunt for.", file=sys.stderr)

if not accession_set:
    with open(output_file, "w") as out:
        pass
    sys.exit(0)

records_to_write = []
matched_accessions = set()

for rec in SeqIO.parse(fasta_file, "fasta"):
    rec_accession = get_accession(rec.id)
    if rec_accession in accession_set:
        records_to_write.append(rec)
        matched_accessions.add(rec_accession)

missing = accession_set - matched_accessions
if missing:
    print(f"WARNING: {len(missing)} accessions not found in {fasta_file}:", file=sys.stderr)
    for m in sorted(missing):
        print(f"  missing: {m}", file=sys.stderr)

print(f"Found and extracted {len(records_to_write)} matching sequences.", file=sys.stderr)

with open(output_file, "w") as out:
    SeqIO.write(records_to_write, out, "fasta")
print(f"SUCCESS: Wrote extracted sequences to {output_file}", file=sys.stderr)