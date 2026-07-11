import pandas as pd
import re

blast_file = snakemake.input.blast
known_file = snakemake.input.known
output_file = snakemake.output.species
hcol = snakemake.wildcards.hcol

blast_df = pd.read_csv(blast_file, sep="\t")
collagenome_df = pd.read_csv(known_file)

blast_df["species"] = blast_df["species"].str.lower()
collagenome_df["species"] = collagenome_df["species"].str.lower()

hcol_df = collagenome_df[
    collagenome_df["hcol_type"].str.contains(hcol, case=False, na=False)
]
master_species = set(hcol_df["species"].unique())

new_species_df = blast_df[~blast_df["species"].isin(master_species)]

if "bitscore" in new_species_df.columns:
    new_species_df = (
        new_species_df.sort_values("bitscore", ascending=False)
        .drop_duplicates("species")
    )
else:
    new_species_df = new_species_df.drop_duplicates("species")

def clean_sseqid(s):
    s = str(s).strip()
    parts = [p for p in s.split("|") if p.strip()]
    return parts[-1] if len(parts) > 1 else parts[0] if parts else s

new_species_df = new_species_df.copy()
new_species_df["accession"] = new_species_df["sseqid"].apply(clean_sseqid)

new_species_df[["accession", "species"]].to_csv(
    output_file, sep="\t", index=False, header=True
)