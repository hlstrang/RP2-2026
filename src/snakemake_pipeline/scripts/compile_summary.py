import pandas as pd
from Bio import SeqIO
import sys
import os
import re

orf_file = snakemake.input.orfs
new_species_file = snakemake.input.new_species

output_file = snakemake.output.summary
hcol_type = snakemake.params.hcol_type
class_type = snakemake.params.class_type

def clean_accession(acc):
    acc = str(acc).strip()
    if "|" in acc:
        parts = [p for p in acc.split("|") if p.strip()]
        if parts:
            acc = parts[0]
    acc = re.sub(r'(_ORF\d+|\.ORF\.\d+|_mRNA.*|\.p\d+|_path\d+)$', '', acc, flags=re.IGNORECASE)
    return acc.strip()

def split_orf_id(acc):
    acc = str(acc).strip()
    return re.sub(r'_ORF\.\d+$', '', acc, flags=re.IGNORECASE)

def expand_species_keys(species_df, join_col):
    """
    Some qseqid values are compound keys joining multiple accessions with '_'
    (e.g. 'GHHG01000294.1_GGKH01030249.1'). Expand each into one row per
    accession instead of truncating to the first.
    """
    rows = []
    for _, row in species_df.iterrows():
        raw_key = str(row[join_col]).strip()
        parts = re.split(r'(?<=\.\d)_(?=[A-Za-z]+\d+\.\d)', raw_key)
        for part in parts:
            new_row = row.copy()
            new_row[join_col] = clean_accession(part)
            rows.append(new_row)
    return pd.DataFrame(rows)

print(f"--- Starting Final Summary Compilation ---", file=sys.stderr)

if os.path.exists(new_species_file) and os.path.getsize(new_species_file) > 0:
    try:
        species_df = pd.read_csv(new_species_file, sep="\t")
        join_col = "qseqid" if "qseqid" in species_df.columns else species_df.columns[0]
        species_df[join_col] = species_df[join_col].apply(clean_accession)

        if species_df.empty:
            print(f"Warning: No new species found for {hcol_type} {class_type}", file=sys.stderr)
            species_df = pd.DataFrame({join_col: [], "species": []})

    except Exception as e:
        raise ValueError(f"ERROR: Failed to parse species file: {e}")
else:
    raise ValueError(f"ERROR: Target species file missing")

bp_lengths = {}
orf_lengths = {}

if os.path.exists(orf_file) and os.path.getsize(orf_file) > 0:
    try:
        for record in SeqIO.parse(orf_file, "fasta"):
            clean_id = split_orf_id(record.id)
            orf_lengths[clean_id] = len(str(record.seq).replace("*", ""))

            header = record.description
            match = re.search(r'length[:=]\s*(\d+)', header, re.IGNORECASE)
            if match:
                bp_lengths[clean_id] = int(match.group(1))
            else:
                print(f"Warning: Could not find length in header for {clean_id}", file=sys.stderr)

    except Exception as e:
        raise ValueError(f"ERROR: Failed to parse ORF file: {e}")
else:
    print(f"Warning: ORF file is empty for {hcol_type} {class_type}. Writing empty summary.", file=sys.stderr)
    empty_df = pd.DataFrame(columns=["accession", "aa_length", "bp_length", "species"])
    empty_df.to_csv(output_file, index=False)
    sys.exit(0)

bp_df = pd.DataFrame(list(bp_lengths.items()), columns=["accession", "bp_length"])
orf_df = pd.DataFrame(list(orf_lengths.items()), columns=["accession", "aa_length"])

final_df = pd.merge(orf_df, bp_df, on='accession', how='inner')
final_df = pd.merge(
    final_df,
    species_df,
    left_on='accession',
    right_on=join_col,
    how='left'
)

if join_col != 'accession' and join_col in final_df.columns:
    final_df = final_df.drop(columns=[join_col])

final_df.to_csv(output_file, index=False)
print(f"--- Summary Compilation Complete ---", file=sys.stderr)