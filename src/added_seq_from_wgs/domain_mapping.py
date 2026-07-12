import argparse
import json
import pyhmmer
import re
from pathlib import Path
from Bio import SeqIO

REFERENCE_FASTA_TEMPLATE = "r_esculentum_{}.fasta"
HMM_DB = "domains.hmm"
OUTPUT_JSON_TEMPLATE = "blast_fragments/{}_domain_coords.json"
C_TERMINAL_TARGET_DOMAINS = ['COLFI', 'C4']
N_TERMINAL_TARGET_DOMAINS = ['SP']
HEAVY_TRIPLE_HELIX_LENGTH = 60

def classify_domain_category(domain_name):
    domain_upper = domain_name.upper()
    if domain_upper == "SP":
        return "signal_peptide"
    elif domain_upper in C_TERMINAL_TARGET_DOMAINS:
        return "c_terminal_anchor"
    elif domain_upper in ["WAP", "VWA", "TSPN"]:
        return "other_modular"
    else:
        return "undetermined"

def run_hmmer_domains(sequence_str, hmm_db_path=HMM_DB, evalue_threshold=0.05):
    alphabet = pyhmmer.easel.Alphabet.amino()
    text_seq = pyhmmer.easel.TextSequence(
        name=b"hmlr_target",
        sequence=sequence_str.strip().upper()
    )
    digital_seq = text_seq.digitize(alphabet)
    with pyhmmer.plan7.HMMFile(hmm_db_path) as hmm_file:
        hmms = list(hmm_file)
    domain_hits = []

    for top_hits in pyhmmer.hmmer.hmmscan([digital_seq], hmms):
        for hit in top_hits:
            family = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
            for domain in hit.domains:
                if domain.i_evalue < evalue_threshold:
                    start = domain.alignment.target_from + 1
                    end = domain.alignment.target_to
                    hit_dict = {
                        "name" : family,
                        "start" : int(start),
                        "end" : int(end),
                        "length" : int(end-start+1),
                        "evalue": domain.i_evalue,
                        "category" : classify_domain_category(family)
                    }
                    domain_hits.append(hit_dict)

                    is_target = family in C_TERMINAL_TARGET_DOMAINS
                    symbol = "[C-TERM]" if is_target else "[OTHER]"
                    print(f"{symbol} {family}: aa {start}-{end} ({hit_dict['length']} aa), E={domain.i_evalue:.2e}")

    return domain_hits

def detect_triple_helix(sequence_str, min_length=30):
    clean_seq = sequence_str.upper().replace('-', "")
    matches = []
    pattern = r"(G..){3,}"

    for m in re.finditer(pattern, clean_seq):
        start = m.start() + 1
        end = m.end()
        matches.append({
            "name": "TripleHelix",
            "start" : start,
            "end" : end,
            "length" : end-start,
            "evalue" : None
        })
        if end - start >= HEAVY_TRIPLE_HELIX_LENGTH:
            print(f"Regex [TripleHelix_>{HEAVY_TRIPLE_HELIX_LENGTH}aa]: aa {start}-{end}")
        else:
            print(f"Regex [TripleHelix_<{HEAVY_TRIPLE_HELIX_LENGTH}aa]: aa {start}-{end}")

    return matches

def select_c_terminal_anchor(all_domains):
    preferred_candidates = [d for d in all_domains if d['name'].upper() in [dt.upper() for dt in C_TERMINAL_TARGET_DOMAINS]]
    if preferred_candidates:
        preferred_candidates.sort(key=lambda x: x['start'])
        names = ", ".join([d['name'] for d in preferred_candidates])
        positions = ", ".join([f"{d['start']}-{d['end']}" for d in preferred_candidates])
        print(f"\n[PREFERRED C-TERM] Selected {len(preferred_candidates)} domain(s): {names} at aa {positions}")
        return preferred_candidates

    protein_length = max(d['end'] for d in all_domains) if all_domains else 0
    distal_domains = [d for d in all_domains if d['end'] > protein_length * 0.8]
    if distal_domains:
        c_terminal = max(distal_domains, key=lambda x: x['end'])
        print(f"\n[FALLBACK C-TERM] Using {c_terminal['name']} at aa {c_terminal['start']}-{c_terminal['end']} (no COLFI/C4 detected)")
        return c_terminal

    if protein_length:
        fallback_start = int(protein_length * 0.66)
        print(f"\n[EMBEDDED FALLBACK] Using region aa {fallback_start}-end (no suitable domains)")
        return [{
            'name': 'C_terminal_estimated',
            'start' : fallback_start,
            'end' : protein_length,
            'length' : protein_length - fallback_start,
            'method' : 'positional_estimate'
        }]

    return []

def build_tblastn_anchors(sp_coords, all_domains, sequence_length):
    print("\n" + "=" * 70)
    print("BUILDING CONSENSUS ANCHOR FRAGMENTS FOR TBLASTN...")
    print("=" * 70)

    ref_id = list(sp_coords.keys())[0] if sp_coords else "unknown"
    sp_info = sp_coords.get(ref_id)

    modular_domains = [d for d in all_domains if d['category'] == "other_modular"]

    if sp_info:
        n_term_end = min(sp_info['end'] + 100, sequence_length)

        if modular_domains:
            max_modular_end = max(d['end'] for d in modular_domains)
            if n_term_end < (max_modular_end + 20):
                n_term_end = min(max_modular_end + 20, sequence_length)
                print(f"[DYNAMIC] Extended N-term anchor to aa {n_term_end} to safeguard all N-terminal modular domains ({', '.join(set(d['name'] for d in modular_domains))}).")

        n_term_anchor = {
            'query_label':f'{ref_id}_Nterm_SP_start',
            'protein_start':sp_info['start'],
            'protein_end' : n_term_end,
            'description' : f'Signal peptide through early triple helix (aa 1-{n_term_end})',
            'primary_evidence': 'SignalP + dynamic structural extension',
            'span_length' : n_term_end - sp_info['start'] + 1
        }
        print(f"\n[N-TERM ANCHOR]")
        print(f" Coordinates: aa {n_term_anchor['protein_start']}-{n_term_anchor['protein_end']}")
        print(f" Length: {n_term_anchor['span_length']} aa")

    else:
        fallback_end = min(150, sequence_length)
        if modular_domains:
            fallback_end = min(max(d['end'] for d in modular_domains) + 50, sequence_length)

        n_term_anchor = {
            'query_label' : f'{ref_id}_Nterm_AutoExt',
            'protein_start' : 1,
            'protein_end': fallback_end,
            'description' : 'Auto-generated N-terminal fragment (no SignalP)',
            'primary_evidence' : 'position_heuristic',
            'span_length' : fallback_end
        }
        print(f"\n[N_TERM ANCHOR - FALLBACK]")
        print(f" Warning: No signal peptide detected, using structural heuristic up to aa {fallback_end}.")

    c_terminal_domains = select_c_terminal_anchor(all_domains)
    if c_terminal_domains:
        first_domain = c_terminal_domains[0]
        last_domain = c_terminal_domains[-1]
        combined_name = "_".join([d['name'] for d in c_terminal_domains])

        c_term_start = max(first_domain['start'] - 50, n_term_anchor['protein_end'] + 10)
        c_term_end = last_domain['end']

        c_term_anchor = {
            'query_label' : f"{ref_id}_Cterm_{combined_name}",
            'protein_start': c_term_start,
            'protein_end': c_term_end,
            'description' : f"Triple helix -> {combined_name} domain(s)",
            'primary_evidence' : first_domain.get('method', 'HMMER'),
            'span_length' : c_term_end - c_term_start + 1,
            'underlying_domain' : combined_name
        }
        print(f"\n[C_TERM ANCHOR]")
        print(f" Coordinates: aa {c_term_anchor['protein_start']}-{c_term_anchor['protein_end']}")
        print(f" Underlying domain: {combined_name} (aa {first_domain['start']}-{last_domain['end']})")
        print(f" Length: {c_term_anchor['span_length']} aa")

        overlap_check = c_term_anchor['protein_start'] <= n_term_anchor['protein_end'] + 10
        if overlap_check and n_term_anchor['protein_end'] < 200:
            print(f"\n[NOTE] Limited gap between anchors ({c_term_anchor['protein_start']-n_term_anchor['protein_end']}). Consider increasing cushion in downstream extraction.")

    else:
        c_term_anchor = {
            'query_label' : f'{ref_id}_Cterm_Fallback',
            'protein_start' : max(int(sequence_length * 0.66), n_term_anchor['protein_end'] + 50),
            'protein_end' : sequence_length,
            'description' : 'Last third of protein (no COLFI/C4 detected)',
            'primary_evidence' : 'positional_estimate',
            'span_length' : sequence_length - int(sequence_length * 0.66),
            'underlying_domain' : 'none_estimated'
        }
        print(f"\n[C_TERM ANCHOR - FALLBACK]")
        print(f" Using estimated C-terminal region: aa {c_term_anchor['protein_start']}-{c_term_anchor['protein_end']}")

    if c_term_anchor['protein_start'] <= n_term_anchor['protein_end']:
        print(f"\n[WARNING] Overlapping anchor regions detected")
        print(f" Adjustment: Shifting C-term start forward by 10 aa")
        c_term_anchor['protein_start'] = n_term_anchor['protein_end'] + 10

    return [n_term_anchor, c_term_anchor]

def main():
    parser = argparse.ArgumentParser(description="Map collagen domains dynamically based on hcol type.")
    parser.add_argument(
        "hcol_type",
        type=str,
        help="The type of collagen family to process (e.g. hcol1, hcol2a, etc...)"
    )
    parser.add_argument(
        "--sp_start",
        type = int,
        default=1,
        help="Start coordinate for the Signal Peptide from SignalP (default = 1)"
    )
    parser.add_argument(
        "--sp_end",
        type = int,
        default=25,
        help="End coordinate for the Signal Peptide from SignalP (default = 25)"
    )
    args = parser.parse_args()
    hcol_type = args.hcol_type
    reference_fasta = REFERENCE_FASTA_TEMPLATE.format(hcol_type)
    output_json = OUTPUT_JSON_TEMPLATE.format(hcol_type)

    print("=" * 70)
    print("INITIALISING COLLAGEN DOMAIN MAPPING PIPELINE...")
    print("=" * 70)

    records = list(SeqIO.parse(reference_fasta, "fasta"))
    if not records:
        raise FileNotFoundError(f"No sequences found in {reference_fasta}")
    seq_rec = records[0]
    seq_str = str(seq_rec.seq)
    seq_length = len(seq_str)

    print(f"\nReference protein: {seq_rec.id}")
    print(f"Length: {seq_length} amino acids\n")

    sp_coordinates = {
        seq_rec.id: {
            "start": args.sp_start,
            "end": args.sp_end
        }
    }

    all_domains = run_hmmer_domains(seq_str)
    triplex_regions = detect_triple_helix(seq_str)
    anchor_fragments = build_tblastn_anchors(sp_coordinates, all_domains, seq_length)

    final_report = {
        'reference_id' : seq_rec.id,
        'sequence_length' : seq_length,
        'anchor_query_definitions' : anchor_fragments,
        'all_detected_domains' : all_domains,
        'triplex_regions' : triplex_regions,
        'signal_peptide_details' : sp_coordinates,
        'config' : {
            'c_terminal_targets' : C_TERMINAL_TARGET_DOMAINS,
            'triplex_min_length' : HEAVY_TRIPLE_HELIX_LENGTH
        }
    }

    with open(output_json, 'w') as f:
        json.dump(final_report, f, indent=2)

    print("\n" + "=" * 70)
    print("PIPELINE COMPLETE")
    print("=" * 70)
    print(f"\nOutput saved to: {output_json}")

    print("\n" + "-" * 70)
    print("SUMMARY OF ANCHOR FRAGMENTS:")
    print("-" * 70)
    print(f"{'Anchor':<20} {'Coordinates':<25} {'Length':>10}")
    print("-" * 70)
    for frag in anchor_fragments:
        coords = f"{frag['protein_start']}-{frag['protein_end']}"
        length = frag['span_length']
        label = frag['query_label'].replace('_', ' ')
        print(f"{label:<20} {coords:>25} {length:>10}")
    print("-" * 70)

if __name__ == "__main__":
    main()