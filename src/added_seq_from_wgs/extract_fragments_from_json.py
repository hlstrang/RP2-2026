import argparse
import json
from pathlib import Path
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description='Generate TBLASTN anchor query fragments')

    parser.add_argument(
        "hcol_type",
        type=str,
        help="The type of collagen family to process (e.g., hcol1, hcol2, col4a1)"
    )
    parser.add_argument("--coords", help="Path to coordinates JSON (defaults to <hcol_type>_domain_coords.json)")
    parser.add_argument("--outdir", help="Output directory (defaults to blast)")
    args = parser.parse_args()

    hcol_type = args.hcol_type

    coord_file = args.coords if args.coords else f"blast_fragments/{hcol_type}_domain_coords.json"
    out_dir_path = Path(args.outdir) if args.outdir else Path(f"blast_fragments")

    out_dir_path.mkdir(exist_ok=True)

    print(f"Loading domain map from {coord_file}")
    try:
        with open(coord_file) as f:
            data = json.load(f)
    except FileNotFoundError:
        raise FileNotFoundError(f"Could not find coordinates file: {coord_file}. Did you run the mapping pipeline script first?")

    ref_fasta = f"r_esculentum_{hcol_type}.fasta"

    try:
        records = list(SeqIO.parse(ref_fasta, "fasta"))
        full_rec = records[0]
    except FileNotFoundError:
        raise FileNotFoundError(f"Could not find reference fasta file: {ref_fasta}")

    print(f"Generating {len(data['anchor_query_definitions'])} anchor fragments...\n")
    for idx, anchor in enumerate(data['anchor_query_definitions']):
        start_idx = anchor['protein_start'] - 1
        end_idx = anchor['protein_end']
        frag_rec = full_rec.__class__(
            id=anchor['query_label'],
            description=anchor['description'],
            seq=full_rec.seq[start_idx:end_idx]
        )

        out_fasta = out_dir_path / f"{idx+1}_{anchor['query_label']}.fasta"
        SeqIO.write([frag_rec], out_fasta, "fasta")

        evidence_type = anchor.get('primary_evidence', 'unknown')
        domain_note = f"[{anchor.get('underlying_domain', 'none')}]" if 'underlying_domain' in anchor else ""

        print(f"Created: {out_fasta.name:<40} ({len(frag_rec.seq):>4} aa) [{evidence_type}{domain_note}]")

    combined_fasta = out_dir_path / f"{hcol_type}_consensus_anchors_combined.fasta"
    combined_recs = []

    for idx, anchor in enumerate(data['anchor_query_definitions']):
        start_idx = anchor['protein_start'] - 1
        end_idx = anchor['protein_end']

        combined_recs.append(
            full_rec.__class__(
                id=f"ANCHOR_N:{idx+1}",
                description=f"{anchor['query_label']} {anchor.get('primary_evidence','')}",
                seq=full_rec.seq[start_idx:end_idx]
            )
        )

    SeqIO.write(combined_recs, combined_fasta, "fasta")
    print(f"\nCombined file: {combined_fasta.name} ({len(combined_recs)} fragments)")
    print(f"\nReady for TBLASTN search: ")
    print(f" tblastn -query {combined_fasta} -db r_luteum_wgs_db ...")

if __name__ == "__main__":
    main()