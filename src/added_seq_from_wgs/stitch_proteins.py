import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def reference_guided_stitch(ref_path, frag1_path, frag2_path, output_path):
    ref_record = SeqIO.read(ref_path, "fasta")
    frags = []
    for p in [frag1_path, frag2_path]:
        try:
            rec = SeqIO.read(p, "fasta")
            if len(rec.seq) > 0:
                frags.append(rec)
        except Exception:
            print(f"[WARNING] Could not read valid sequence from {p}")

    if not frags:
        print(f"[ERROR] No valid protein predictions found in either fragment file")
        sys.exit(1)
    if len(frags) == 1:
        print(f"[ERROR] Only one fragment contained prediction. Saving only one...")
        SeqIO.write(frags[0], output_path, "fasta")
        return

    mapping = []
    for f in frags:
        seed = str(f.seq[:15])
        start_idx = str(ref_record.seq).find(seed)

        if start_idx == -1:
            seed_tail = str(f.seq[-15:])
            end_idx = str(ref_record.seq).find(seed_tail)
            if end_idx != -1:
                start_idx = end_idx - len(f.seq) + 15

        if start_idx == -1:
            print(f"[ERROR] Could not anchor fragment {f.id} to reference protein")
            sys.exit(1)

        mapping.append({
            'record' : f,
            'start' : start_idx,
            'end' : start_idx + len(f.seq)
        })

    mapping.sort(key=lambda x: x['start'])
    piece1 = mapping[0]
    piece2 = mapping[1]

    print(f"Piece 1 ({piece1['record'].id}) maps to ref amino acids: {piece1['start']}-{piece1['end']}")
    print(f"Piece 2 ({piece2['record'].id}) maps to ref amino acids: {piece2['start']}-{piece2['end']}")

    final_seq_str = str(piece1['record'].seq)
    if piece1['end'] < piece2['start']:
        gap_len = piece2['start'] - piece1['end']
        print(f"[WARNING] Missing {gap_len} amino acids in the assembly. Inserting 'X' gap mask...")
        final_seq_str += ("X" * gap_len) + str(piece2['record'].seq)

    elif piece1['end'] >= piece2['start']:
        overlap_len = piece1['end'] - piece2['start']
        print(f"[WARNING] Trimming {overlap_len} duplicate boundary amino acids")
        final_seq_str += str(piece2['record'].seq)[overlap_len:]

    stitched_rec = SeqRecord(
        Seq(final_seq_str),
        id=f"{frags[0].id}_{frags[1].id}_stitched",
        description="Unified collagen model stitched via reference template guidance"
    )

    SeqIO.write(stitched_rec, output_path, "fasta")
    print(f"[SUCCESS] Stitched protein saved to {output_path} (Length: {len(final_seq_str)} aa)")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python3 stitch_proteins.py <ref.fasta> <frag1.fasta> <frag2.fasta> <output.fasta>")
        sys.exit(1)
    reference_guided_stitch(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])