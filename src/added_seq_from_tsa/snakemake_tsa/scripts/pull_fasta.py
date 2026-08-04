import pandas as pd
from Bio import SeqIO

blast_file = snakemake.input.blast
tsa_fasta = snakemake.input.tsa_fasta
output_fasta = snakemake.output[0]

blast_df = pd.read_csv(blast_file, sep="\t")

if not blast_df.empty and "sseqid" in blast_df.columns:
    target_ids = set(blast_df["sseqid"].dropna().astype(str))
else:
    target_ids = set()

records_to_write = []
if target_ids:
    for record in SeqIO.parse(tsa_fasta, "fasta"):
        header_id = record.id.split()[0]
        if record.id in target_ids or header_id in target_ids:
            records_to_write.append(record)

SeqIO.write(records_to_write, output_fasta, "fasta")