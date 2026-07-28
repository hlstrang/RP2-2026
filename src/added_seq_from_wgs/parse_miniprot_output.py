import sys
import re

def parse_miniprot(input_file, aa_output_file, exon_output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    paf_line = None
    sta_lines = []
    cds_features = []

    for line in lines:
        line_strip = line.strip()
        if not line_strip:
            continue

        if '\t' in line_strip and line_strip.startswith('##PAF'):
            paf_line = line_strip.replace('##PAF', '').strip()
            continue

        if line_strip.startswith('##STA'):
            seq_part = line_strip.replace('##STA', '').strip()
            if seq_part:
                sta_lines.append(seq_part)
            continue

        if not line_strip.startswith('#') and '\t' in line_strip:
            parts = line_strip.split('\t')
            if len(parts) >= 9 and parts[2] == 'CDS':
                chrom = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                cds_features.append((chrom, start, end, strand))

    if not paf_line:
        print("Error: Could not find ##PAF summary line.")
        sys.exit(1)

    if not cds_features:
        print("Error: Could not find any CDS features in GFF output.")
        sys.exit(1)

    paf_fields = paf_line.split('\t')
    query_id = paf_fields[0]
    raw_target_id = paf_fields[5]
    paf_strand = paf_fields[4]

    coord_match = re.match(r"([^:]+):(\d+)-(\d+)", raw_target_id)
    if coord_match:
        scaffold_name = coord_match.group(1)
        slice_offset = int(coord_match.group(2)) - 1
    else:
        scaffold_name = raw_target_id
        slice_offset = 0

    gene_start = cds_features[0][1] + slice_offset
    gene_end = cds_features[-1][2] + slice_offset

    fasta_header = f">{query_id}_{scaffold_name}_{gene_start}_{gene_end}_{paf_strand}"

    full_sequence = "".join(sta_lines).replace('*', '').replace('-', '').replace(' ', '')

    with open(aa_output_file, 'w') as aa_out:
        aa_out.write(f"{fasta_header}\n")
        for i in range(0, len(full_sequence), 60):
            aa_out.write(f"{full_sequence[i:i+60]}\n")

    with open(exon_output_file, 'w') as exon_out:
        exon_out.write("exon_num\tscaffold\tgenomic_start\tgenomic_end\tstrand\texon_length_bp\n")

        for idx, (chrom, c_start, c_end, c_strand) in enumerate(cds_features, 1):
            abs_start = c_start + slice_offset
            abs_end = c_end + slice_offset
            length = (abs_end - abs_start) + 1

            exon_out.write(f"Exon_{idx}\t{scaffold_name}\t{abs_start}\t{abs_end}\t{c_strand}\t{length}\n")

    print(f"Success!")
    print(f"  -> Amino acid sequence saved to: {aa_output_file}")
    print(f"  -> Exon boundaries saved to:    {exon_output_file}")
    print(f"  -> Locus: {scaffold_name}:{gene_start}-{gene_end} ({paf_strand}) | Total Exons: {len(cds_features)}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 parse_miniprot_output.py <miniprot_output.txt> <output.fasta> <exons.tsv>")
        sys.exit(1)

    parse_miniprot(sys.argv[1], sys.argv[2], sys.argv[3])