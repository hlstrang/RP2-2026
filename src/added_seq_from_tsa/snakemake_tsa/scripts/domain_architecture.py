import pyhmmer
import re
import os
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.patches as patches

input_fasta = snakemake.input.pep
hmm_db_path = snakemake.input.hmm_db
output_csv = snakemake.output.csv
output_plot_dir = snakemake.output.plot_dir

os.makedirs(output_plot_dir, exist_ok=True)

HCOL_TYPE = snakemake.wildcards.hcol

DOMAIN_COLOURS = {
    "WAP": "#4C72B0",
    "Triple Helix": "#55A868",
    "COLFI": "#C44E52",
    "TSPN": "#8172B3",
    "VWA": "#CCB974",
    "C4": "#64B5CD",
    "SP" : "#DCEF05"
}
DEFAULT_COLOUR = "#8C8C8C"

def scan_with_hmmer(sequence_str, hmm_db_path="domains.hmm", evalue_threshold=0.05):
    alphabet = pyhmmer.easel.Alphabet.amino()
    text_seq = pyhmmer.easel.TextSequence(
        name=b"target_seq",
        sequence=sequence_str.strip().upper()
    )
    digital_seq = text_seq.digitize(alphabet)

    with pyhmmer.plan7.HMMFile(hmm_db_path) as hmm_file:
        hmms = list(hmm_file)

    domain_hits = []

    print("="* 60)
    print(f"HMMER Domain Feature Search Results (E-value < {evalue_threshold}):")
    print("=" * 60)

    for top_hits in pyhmmer.hmmer.hmmscan([digital_seq], hmms):
        for hit in top_hits:
            if hit.evalue >= evalue_threshold:
                continue

            family = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
            colour = DOMAIN_COLOURS.get(family, DEFAULT_COLOUR)

            valid_domains = []
            for domain in hit.domains:
                if domain.i_evalue < evalue_threshold:
                    valid_domains.append(domain)

            if valid_domains:
                print(f"\n[+] Found Domain: {family}")
                print(f" - Hit E-value: {hit.evalue:.2e}")
                print(f" - Bit Score: {hit.score:.1f}")

                for domain in valid_domains:
                    start = domain.alignment.target_from
                    end = domain.alignment.target_to
                    print(f"   -> Position {start}-{end} (Domain E-value: {domain.i_evalue:.2e})")

                    domain_hits.append({
                        "name": family,
                        "start": start,
                        "end": end,
                        "colour": colour
                    })
    print("=" * 60)
    return domain_hits

def find_triple_helix(sequence, max_gap_allowed=7):
    DOMAIN = r"(G..){3,}"
    clean_seq = "".join(sequence.split()).upper()

    raw_matches = []
    for m in re.finditer(DOMAIN, clean_seq):
        raw_matches.append({"start": m.start(), "end": m.end()})

    if not raw_matches:
        print("\nRegex Triple-Helix Search Results:\nNo triple-helix regions found.\n" + "="*60)
        return []

    merged_matches = [raw_matches[0]]
    for current in raw_matches[1:]:
        last = merged_matches[-1]

        if current["start"] - last["end"] <= max_gap_allowed:
            last["end"] = current["end"]
        else:
            merged_matches.append(current)

    final_domains = []
    print("\nRegex Triple-Helix Search Results:")
    for dom in merged_matches:
        length = dom["end"] - dom["start"]
        name = "Col1" if length >= 60 else "Col2"

        final_domains.append({
            "name": name,
            "start": dom["start"],
            "end": dom["end"],
            "colour": DOMAIN_COLOURS["Triple Helix"]
        })
        print(f"[+] Found Triple Helix: Position {dom['start']}-{dom['end']} (Length: {length} aa)")
    print("=" * 60)

    return final_domains

def add_manual_domain(name, start, end, colour=None):
    return {
        "name": name,
        "start": start,
        "end": end,
        "colour": colour if colour else DOMAIN_COLOURS.get(name, DEFAULT_COLOUR)
    }

def classify_and_validate(domains, expected_family):
    expected_family = expected_family.lower()

    col1 = [d for d in domains if d["name"] == "Col1"]
    col2 = [d for d in domains if d["name"] == "Col2"]
    wap = [d for d in domains if d["name"] == "WAP"]
    vwa = [d for d in domains if d["name"] == "VWA"]
    tspn = [d for d in domains if d["name"] == "TSPN"]
    c4 = [d for d in domains if d["name"] == "C4"]
    colfi = [d for d in domains if d["name"] == "COLFI"]

    col1_len = (col1[0]["end"] - col1[0]["start"] if col1 else 0)
    col2_len = (col2[0]["end"] - col2[0]["start"] if col2 else 0)
    tspn_len = (tspn[0]['end'] - tspn[0]['start'] if tspn else 0)

    has_other_nc = len(wap) > 0 or len(vwa) > 0 or tspn_len > 0

    print(f"\n[VALIDATING] Checking characteristics against expected type: {expected_family.upper()}")

    if expected_family == "hcol1":
        if 1010 <= col1_len <= 1030 and 50 <= col2_len <= 65 and colfi and not has_other_nc:
            return "Fibrillar Hcol1", "PASSED"
        return "Hcol1 Candidate", "WARNING: Coordinates vary slightly from archetype"

    elif expected_family in ["hcol2a", "hcol2b", "hcol2"]:
        wap_count = len(wap)
        if expected_family == "hcol2a" and wap_count == 1:
            return "Fibrillar Hcol2a", "PASSED"
        elif expected_family == "hcol2b" and wap_count == 2:
            return "Fibrillar Hcol2b", "PASSED"
        elif expected_family == "hcol2" and wap_count in [1, 2]:
            return "Fibrillar Hcol2", "PASSED"

        return f"Hcol2 Candidate", "WARNING: Missing required domains (Col1, Col2, WAP, or COLFI)"

    elif expected_family == "hcol3":
        wap_count = len(wap)
        vwa_count = len(vwa)
        if wap_count > 1 and vwa_count > 1 and colfi:
            return "Fibrillar Hcol3", "PASSED"
        return "Hcol3 Candidate", "CHECK: Verify domain structures (Requires >1 WAP, >1 VWA, COLFI)"

    elif expected_family == "hcol4":
        if len(c4) >= 2 and not has_other_nc:
            return "Network-forming Hcol4", "PASSED"
        return "Hcol4 Candidate", "CHECK: Verify C4 presence or unexpected non-collagenous domains found"

    elif expected_family == "hcol5":
        gap = (col2[0]["start"] - col1[0]["end"]) if (col1 and col2) else 0
        if 1020 <= col1_len <= 1035 and 50 <= col2_len <= 60 and gap > 20 and colfi and not has_other_nc:
            return "Fibrillar Hcol5", "PASSED"
        return "Hcol5 Candidate", "CHECK: Inter-domain spacing or domain lengths differs"

    elif expected_family == "hcol6":
        if len(wap) > 1 and len(vwa) > 1 and len(c4) >= 2:
            return "Network-forming Hcol6", "PASSED"
        return "Hcol6 Candidate", "CHECK: Verify domain structures (Requires >1 WAP, >1 VWA, C4)"

    elif expected_family == "hcol7":
        if tspn_len > 0 and colfi:
            return "Fibrillar Hcol7", "PASSED"
        return "Hcol7 Candidate", "CHECK: Verify domain structures (Requires TSPN, COLFI)"

    elif expected_family == "hcol8":
        if colfi and not col2 and 950 <= len(col1) <= 1050 and not has_other_nc:
            return "Fibrillar Hcol8", "PASSED"
        return "Hcol8 Candidate", "WARNING: Domain structure mismatch"

    return "Unidentified Collagen", "UNVERIFIED"

def classify_unknown_collagen(domains, seq_length):
    REFERENCE_PROFILES = {
        "HCOL1": {
            "required": {"col1": 1, "col2": 1, "colfi": True},
            "forbidden": {"wap":0, "vwa":0, "tspn": 0},
            "constraints": {
                "col1_len": (1010, 1030),
                "col2_len": (50, 65)
            }
        },
        "HCOL2": {
            "required": {"col1": 1, "col2": 1, "colfi": True, "wap": [1,2]},
            "optionals": {},
            "constraints": {}
        },
        "HCOL3": {
            "required": {"wap": (">", 1), "vwa": (">", 1), "colfi": True},
            "forbidden": {},
            "constraints": {}
        },
        "HCOL4": {
            "required": {"c4": (">=", 2)},
            "forbidden": {"wap": 0,"vwa": 0, "tspn": 0},
            "constraints": {}
        },
        "HCOL5": {
            "required": {"col1": 1, "col2": 1, "colfi": True},
            "forbidden": {"wap": 0, "vwa": 0, "tspn": 0},
            "constraints": {
                "col1_len": (1020, 1035),
                "col2_len" : (50, 60),
                "col2_start_gap": (">", 20)
            }
        },
        "HCOL6": {
            "required": {"wap": (">", 1), "vwa": (">", 1), "c4": (">=", 2)},
            "forbidden": {},
            "constraints": {}
        },
        "HCOL7": {
            "required": {"tspn": (">=", 1), "colfi": True},
            "forbidden": {},
            "constraints": {}
        },
        "HCOL8": {
            "required": {"col1": 1, "colfi": True},
            "forbidden": {"col2": 0, "wap": 0, "vwa": 0, "tspn": 0},
            "constraints": {
                "col1_len": (950, 1050)
            }
        }
    }

    domain_counts = {
        "col1": len([d for d in domains if d["name"] == "Col1"]),
        "col2": len([d for d in domains if d["name"] == "Col2"]),
        "wap": len([d for d in domains if d["name"] == "WAP"]),
        "vwa": len([d for d in domains if d["name"] == "VWA"]),
        "c4": len([d for d in domains if d["name"] == "C4"]),
        "tspn": len([d for d in domains if d["name"] == "TSPN" and (d["end"] - d["start"]) > 0]),
        "colfi": len([d for d in domains if d["name"] == "COLFI"]) > 0,
    }

    domain_lengths = {}
    for domain_type in ["col1", "col2"]:
        matches = [d for d in domains if d["name"] == domain_type.title()]
        if matches:
            domain_lengths[f"{domain_type}_len"] = matches[0]["end"] - matches[0]["start"]

    domain_gaps = {}
    col1_domains = [d for d in domains if d["name"] == "Col1"]
    col2_domains = [d for d in domains if d["name"] == "Col2"]
    if col1_domains and col2_domains:
        domain_gaps["col2_start_gap"] = col2_domains[0]["start"] - col1_domains[0]["end"]

    scores = {}
    for hcol_type, profile in REFERENCE_PROFILES.items():
        score = 0
        max_possible = 0
        violations = []

        for domain_key, condition in profile["required"].items():
            max_possible += 1

            if isinstance(condition, dict):
                pass
            elif isinstance(condition, list):
                if domain_counts.get(domain_key, 0) in condition:
                    score += 1
                else:
                    violations.append(f"{domain_key}: need one of {condition}, got {domain_counts.get(domain_key, 0)}")
            elif isinstance(condition, tuple) and condition[0] in [">", ">="]:
                if condition[0] == ">":
                    if domain_counts.get(domain_key, 0) > condition[1]:
                        score += 1
                    else:
                        violations.append(f"{domain_key}: need >{condition[1]}, got {domain_counts.get(domain_key, 0)}")
                else:
                    if domain_counts.get(domain_key, 0) >= condition[1]:
                        score += 1
                    else:
                        violations.append(f"{domain_key}: need >={condition[1]}, got {domain_counts.get(domain_key, 0)}")
            else:
                if domain_counts.get(domain_key, 0) == condition:
                    score += 1
                else:
                    violations.append(f"{domain_key}: need {condition}, got {domain_counts.get(domain_key, 0)}")

        forbidden_violation_count = 0
        if "forbidden" in profile:
            for domain_key, value in profile["forbidden"].items():
                if value == 0 and domain_counts.get(domain_key, 0) > 0:
                    violations.append(f"Forbidden {domain_key} found: {domain_counts.get(domain_key, 0)}")
                    forbidden_violation_count += 1

                elif domain_key == "col2_len_lt_30" and value is True:
                    col2_domains_check = [d for d in domains if d["name"] == "Col2"]
                    if col2_domains_check:
                        col2_length = col2_domains_check[0]["end"] - col2_domains_check[0]["start"]
                        if col2_length >= 30:
                            violations.append(
                                f"Forbidden col2 with length >= 30 found: {col2_length} aa"
                            )
                            forbidden_violation_count += 1

        if forbidden_violation_count > 0:
            scores[hcol_type] = {
                "score": 0,
                "matches": "DISQUALIFIED",
                "violations": violations,
                "confidence": "none"
            }
            continue

        if "constraints" in profile:
            for constraint_key, constraint_val in profile["constraints"].items():
                if constraint_key.endswith("_len") and isinstance(constraint_val, tuple):
                    max_possible += 0.5
                    if constraint_key in domain_lengths:
                        length = domain_lengths[constraint_key]
                        min_len, max_len = constraint_val
                        if min_len <= length <= max_len:
                            score += 0.5
                        else:
                            violations.append(f"{constraint_key}: {length} outside [{min_len}, {max_len}]")

                elif constraint_key == "col2_start_gap":
                    max_possible += 0.5
                    if constraint_key in domain_gaps:
                        gap = domain_gaps[constraint_key]
                        if isinstance(constraint_val, tuple) and constraint_val[0] == ">":
                            if gap > constraint_val[1]:
                                score += 0.5
                            else:
                                violations.append(f"col2_start_gap: {gap} (need >{constraint_val[1]})")

        if max_possible > 0:
            scores[hcol_type] = {
                "score": score / max_possible,
                "matches": f"{score}/{max_possible}",
                "violations": violations,
                "confidence": "high" if score / max_possible >= 0.8 else "medium" if score / max_possible >= 0.5 else "low"
            }

    sorted_scores = sorted(scores.items(), key=lambda x: x[1]["score"], reverse=True)

    if sorted_scores:
        best_match = sorted_scores[0]

        print(f"\n{'='*60}")
        print(f"UNKNOWN COLLAGEN CLASSIFICATION RESULTS")
        print(f"{'='*60}")
        print(f"\nBest Match: {best_match[0]}")
        print(f"Match Score: {best_match[1]['matches']}")
        print(f"Confidence: {best_match[1]['confidence'].upper()}")

        if best_match[1]['violations']:
            print(f"\nWarnings/Constraints Violated:")
            for v in best_match[1]['violations']:
                print(f"  - {v}")

        print(f"\nTop Candidates:")
        for hcol_type, data in sorted_scores[:3]:
            marker = ">>> " if hcol_type == best_match[0] else "    "
            print(f"  {marker}{hcol_type}: {data['matches']} ({data['confidence']})")

        print(f"{'='*60}\n")

        return best_match[0], best_match[1]['confidence'], best_match[1]['violations']
    else:
        return "UNKNOWN", "none", ["No reference profiles matched"]

def plot_domain_arcitecture(seq_length, domains, output, title = "Domain Architecture"):
    default_colours = [
        "#4C72B0", "#55A868", "#C44E52", "#8172B3",
        "#CCB974", "#64B5CD", "#8C8C8C"
        ]
    fig, ax = plt.subplots(figsize=(12,2))

    backbone = patches.Rectangle(
        (0, 0.25), seq_length, 0.5,
        linewidth=0,
        facecolor="#939090FF"
        )
    ax.add_patch(backbone)

    for i, dom in enumerate(domains):
        colour = dom.get("colour", default_colours[i % len(default_colours)])
        start, end = dom["start"], dom["end"]
        width = end - start

        if width <= 0:
            continue

        rect = patches.Rectangle(
            (start, 0.25), width, 0.5,
            linewidth = 1, edgecolor = "black", facecolor=colour
        )
        ax.add_patch(rect)

        label = f"{dom['name']} ({width} aa)"

        if i % 2 == 0:
            y = 0.85
            va = "bottom"
        else:
            y = 0.15
            va = "top"

        if width >= 200:
            ax.text(
                start + width/2, 0.5, label,
                ha="center", va="center", fontsize=8, color="black"
            )
        else:
            ax.text(
                start + width/2, y, label,
                ha="center", va=va, fontsize=7, color="black"
            )

    ax.set_xlim(0, seq_length)
    ax.set_ylim(0,1)
    ax.set_yticks([])
    ax.set_xlabel("Amino Acid Position")
    ax.set_title(title)
    plt.tight_layout()
    fig.savefig(f"{output}.png", dpi=300)
    plt.close(fig)

def classify_sequence(domains, seq_length, expected_family=None):
    unknown_type, confidence, violations = classify_unknown_collagen(domains, seq_length)

    if expected_family:
        known_identity, validation_status = classify_and_validate(domains, expected_family)
    else:
        known_identity = "No expected type specified"
        validation_status = "N/A"

    return {
        "known_identity": known_identity,
        "validation_status": validation_status,
        "best_unknown_match": unknown_type,
        "confidence": confidence,
        "violations": violations
    }

records_data = []

if os.path.exists(input_fasta) and os.path.getsize(input_fasta) > 0:
    for record in SeqIO.parse(input_fasta, "fasta"):
        seq_id = record.id
        sequence_str = str(record.seq)
        seq_length = len(sequence_str)

        hmm_domains = scan_with_hmmer(sequence_str, hmm_db_path=hmm_db_path)

        helix_domains = find_triple_helix(sequence_str)

        all_domains = sorted(
            hmm_domains + helix_domains,
            key=lambda d: d["start"]
        )

        results = classify_sequence(all_domains, seq_length, expected_family=HCOL_TYPE)

        plot_path = os.path.join(output_plot_dir, f"{seq_id}")
        plot_title = (
            f"{seq_id} Architecture Map\n"
            f"Validation: {results['known_identity']} ({results['validation_status']}) | "
            f"Best Match: {results['best_unknown_match']}"
        )
        plot_domain_arcitecture(seq_length, all_domains, plot_path, title=plot_title)

        records_data.append({
            "sequence_id": seq_id,
            "seq_length": seq_length,
            "expected_family": HCOL_TYPE,
            "known_identity": results["known_identity"],
            "validation_status": results["validation_status"],
            "best_unknown_match": results["best_unknown_match"],
            "confidence": results["confidence"],
            "violations": "; ".join(results["violations"])
        })

df = pd.DataFrame(records_data)
df.to_csv(output_csv, index=False)