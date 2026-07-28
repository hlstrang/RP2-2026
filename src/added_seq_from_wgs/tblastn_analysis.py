import argparse
import os
import subprocess
import sys
import pandas as pd

## Change species
SPECIES = "h_viridissima"

BLAST_RESULTS_TEMPLATE = f"{SPECIES}/genomic_fragments/{{}}/{SPECIES}_{{}}_query.tsv"
OUTPUT_FASTA_TEMPLATE = f"{SPECIES}/genomic_fragments/{{}}/{SPECIES}_{{}}_fragment.fasta"

## Update db name
WGS_DB_NAME = "db/hydra_wgs_db"
CUSHION = 5000

def extract_scaffold_name(sseqid):
    """Extract clean scaffold name from various BLAST accession formats."""
    import re
    match = re.search(r'\|([^|]+)\|', str(sseqid))
    if match:
        return match.group(1)
    elif 'contig:' in str(sseqid) or 'SCAFFOLD' in str(sseqid):
        return str(sseqid)
    else:
        return str(sseqid)

def calculate_and_extract(blast_results, output_fasta, cushion_size, rank=1):
    try:
        df = pd.read_csv(blast_results, sep='\t', header=0)
    except FileNotFoundError:
        print(f"ERROR: Could not find BLAST results file: {blast_results}")
        sys.exit(1)

    if df.empty:
        print("ERROR: No hits found, check database or query.")
        sys.exit(1)

    df['sseqid_clean'] = df['sseqid'].apply(extract_scaffold_name)
    print("=" * 60)
    print(f"ORTHOLOG ANALYSIS & SELECTION (RANK: {rank})")
    print("=" * 60)

    identity_threshold = 70.0
    df_filtered = df[df['pident'] >= identity_threshold]

    if df_filtered.empty:
        print(f"WARNING: No hits found above {identity_threshold}% identity. Falling back to all hits.")
        df_filtered = df

    scaffold_scores = df_filtered.groupby('sseqid_clean').agg(
        HSP_count=('bitscore', 'count'),
        Max_Bitscore=('bitscore', 'max'),
        Mean_Identity=('pident', 'mean'),
        Total_Score=('bitscore', 'sum')
    )

    scaffold_scores = scaffold_scores.sort_values(
        by=['Total_Score', 'HSP_count', 'Mean_Identity'],
        ascending=[False, False, False]
    )

    print("\n[Scaffold-Level Aggregation (Filtered for Orthology)]")
    print(scaffold_scores.head(10).to_string())

    print("\n" + "=" * 60)
    print(f"SELECTION STRATEGY: Target Rank {rank} Scaffold")
    print("=" * 60)

    target_index = rank - 1

    if target_index >= len(scaffold_scores):
        print(f"ERROR: Requested rank {rank}, but only {len(scaffold_scores)} distinct scaffolds were found.")
        sys.exit(1)

    top_scaffold = scaffold_scores.index[target_index]
    print(f"Selected Rank {rank} target scaffold: {top_scaffold}")
    print(f"  Mean Identity: {scaffold_scores.loc[top_scaffold, 'Mean_Identity']:.1f}%")
    print(f"  Max Bitscore:  {scaffold_scores.loc[top_scaffold, 'Max_Bitscore']:.1f}")

    print("\nCalculating genomic coordinates...")
    scaffold_hits = df[df['sseqid_clean'] == top_scaffold].copy()

    scaffold_hits['is_forward'] = scaffold_hits['sstart'] < scaffold_hits['send']
    scaffold_hits['hsp_start'] = scaffold_hits[['sstart', 'send']].min(axis=1)
    scaffold_hits['hsp_end'] = scaffold_hits[['sstart', 'send']].max(axis=1)

    scaffold_hits = scaffold_hits.sort_values(by='hsp_start')

    forward_hits = scaffold_hits[scaffold_hits['is_forward']]
    reverse_hits = scaffold_hits[~scaffold_hits['is_forward']]

    f_score = forward_hits['bitscore'].sum()
    r_score = reverse_hits['bitscore'].sum()

    print(f"Hit distribution:")
    print(f"  Forward strand HSPs: {len(forward_hits)} (Total Bitscore: {f_score:,.1f})")
    print(f"  Reverse strand HSPs: {len(reverse_hits)} (Total Bitscore: {r_score:,.1f})")

    if f_score >= r_score:
        dominant_strand_label = "FORWARD"
        winning_hits = forward_hits
        dropped_count = len(reverse_hits)
    else:
        dominant_strand_label = "REVERSE"
        winning_hits = reverse_hits
        dropped_count = len(forward_hits)

    print(f"\nSelected {dominant_strand_label} strand based on higher total bitscore")
    if dropped_count > 0:
        print(f"   Dropped {dropped_count} minority-strand noisy HSPs from coordinate calculations")

    all_coords = pd.concat([winning_hits['sstart'], winning_hits['send']])
    min_coord = int(all_coords.min())
    max_coord = int(all_coords.max())

    print("\n" + "-" * 75)
    print(f"DETAILED HIT INSPECTION FOR SCAFFOLD: {top_scaffold}")
    print("-" * 75)
    print(f"{'Strand':<8} {'Query Region':<15} {'Genomic Span':<22} {'Ident%':<8} {'Bitscore':<8}")
    print("-" * 75)

    for _, row in scaffold_hits.iterrows():
        strand_label = "FORWARD" if row['is_forward'] else "REVERSE"
        query_span = f"{row['qstart']}-{row['qend']}"
        genomic_span = f"{row['sstart']:,}-{row['send']:,}"
        active_marker = "" if row['is_forward'] == (dominant_strand_label == "FORWARD") else "[DROPPED]"

        if not row['is_forward']:
            strand_label = f"*{strand_label}"

        print(f"{strand_label:<8} {query_span:<15} {genomic_span:<22} {row['pident']:<8.1f} {row['bitscore']:<8.1f} {active_marker}")
    print("-" * 75)

    extract_start = max(1, min_coord - cushion_size)
    extract_end = max_coord + cushion_size

    coord_range = f"{extract_start}-{extract_end}"

    extracted_length = extract_end - extract_start + 1

    print(f"\nGene body spans {min_coord:,} to {max_coord:,} on scaffold {top_scaffold}")
    print(f"Applying {cushion_size:,}bp padding -> Extraction Range: {coord_range} ({extracted_length:,} bp)")

    out_dir = os.path.dirname(output_fasta)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    print(f"\nPulling genomic fragment from database...")
    try:
        extract_cmd = [
            "blastdbcmd",
            "-db", WGS_DB_NAME,
            "-entry", top_scaffold,
            "-range", coord_range
        ]
        raw_sequence = subprocess.run(extract_cmd, check=True, capture_output=True, text=True)

        print(f"Applying dustmasker to softmask low-complexity regions...")
        mask_cmd = [
            "dustmasker",
            "-infmt", "fasta",
            "-outfmt", "fasta",
            "-out", output_fasta
        ]
        subprocess.run(mask_cmd, input=raw_sequence.stdout, check=True, text=True)

        print(f"\n{'='*60}")
        print(f"SUCCESS")
        print(f"{'='*60}")
        print(f"  Output file: {output_fasta}")
        print(f"  Region size: {extracted_length:,} bp")
        print(f"  Hits included: {len(scaffold_hits)} HSPs")
        print(f"  Target mean identity: {scaffold_scores.loc[top_scaffold, 'Mean_Identity']:.1f}%")
    except subprocess.CalledProcessError as e:
        print(f"ERROR during extraction: {e}")
        if e.stderr:
            print(f"stderr: {e.stderr}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Extract and mask genomic scaffolds based on aggregated BLAST bitscores.")

    parser.add_argument(
        "hcol_type",
        type=str,
        help="The type of collagen family to process (e.g., hcol1, hcol2a, etc...)"
    )
    parser.add_argument("--blast-results", help="Path to input BLAST TSV file.")
    parser.add_argument("--output-fasta", help="Path to save output softmasked FASTA.")
    parser.add_argument("--cushion", type=int, default=5000, help="BP padding added to flanking regions")
    parser.add_argument("--rank", type=int, default=1, help="The rank of the target scaffold to extract (1 = highest, 2 = second highest, etc.)")

    args = parser.parse_args()
    hcol_type = args.hcol_type

    blast_results = (
        args.blast_results
        if args.blast_results
        else BLAST_RESULTS_TEMPLATE.format(hcol_type, hcol_type)
    )
    output_fasta = (
        args.output_fasta
        if args.output_fasta
        else OUTPUT_FASTA_TEMPLATE.format(hcol_type, hcol_type)
    )

    calculate_and_extract(blast_results, output_fasta, args.cushion, args.rank)

if __name__ == "__main__":
    main()