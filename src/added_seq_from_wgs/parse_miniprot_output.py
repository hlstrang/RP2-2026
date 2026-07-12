import sys
import re

def parse_miniprot(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    paf_line = None
    sequence_lines = []
    in_sequence_block = False

    for line in lines:
        line_strip = line.strip()
        if not line_strip:
            continue

        if '\t' in line_strip and not line_strip.startswith('#'):
            paf_line = line_strip
            continue

        if line_strip.startswith('##STA'):
            in_sequence_block = True
            seq_start = line_strip.replace('##STA', '').strip()
            if seq_start:
                sequence_lines.append(seq_start)
            continue

        if in_sequence_block:
            sequence_lines.append(line_strip)

    if not paf_line:
        print("Error: Could not find the summary mapping data line.")
        sys.exit(1)

    if not sequence_lines:
        print("Error: Could not find the translated protein sequence block.")
        sys.exit(1)

    fields = paf_line.split('\t')
    query_id = fields[0]
    strand = "plus" if fields[4] == "+" else "minus"
    target_id = fields[5]
    t_start = fields[7]
    t_end = fields[8]

    full_sequence = "".join(sequence_lines).replace('*', '')
    fasta_header = f">{query_id}_{target_id}_{t_start}_{t_end}_{strand}"

    with open(output_file, 'w') as out:
        out.write(f"{fasta_header}\n")
        for i in range(0, len(full_sequence), 60):
            out.write(f"{full_sequence[i:i+60]}\n")

    print(f"Success: Spliced protein successfully parsed into {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 parse_miniprot_output.py <miniprot_raw_output.txt> <output.fasta>")
        sys.exit(1)

    parse_miniprot(sys.argv[1], sys.argv[2])