import pyhmmer
import re
import matplotlib.pyplot as plt
import matplotlib.patches as patches

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
    """
    Manually specify a domain region with known coordinates,
    bypassing HMMER/regex detection entirely.
    """
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

    print(f"\n[VALIDATING] Checking characteristics against expected type: {expected_family.upper()}")

    if expected_family == "hcol1":
        if 1010 <= col1_len <= 1030 and 50 <= col2_len <= 65 and colfi:
            return "Fibrillar Hcol1", "PASSED"
        return "Hcol1 Candidate", "WARNING: Coordinates vary slightly from archetype"

    elif expected_family in ["hcol2a", "hcol2b"]:
        wap_count = len(wap)
        if expected_family == "hcol2a" and wap_count == 1:
            return "Fibrillar Hcol2a", "PASSED"
        elif expected_family == "hcol2b" and wap_count == 2:
            return "Fibrillar Hcol2b", "PASSED"

        if wap_count == 1:
            return "Fibrillar Hcol2a structure", "WARNING: Expected hcol2b but found 1x WAP"
        elif wap_count == 2:
            return "Fibrillar Hcol2b structure", "WARNING: Expected hcol2a but found 2x WAP"

        return f"Hcol2 Variant", "WARNING: No WAP domains identified"

    elif expected_family == "hcol3":
        wap_count = len(wap)
        vwa_count = len(vwa)
        if wap_count > 1 and vwa_count > 1 and colfi:
            return "Fibrillar Hcol3", "PASSED"
        return "Hcol3 Candidate", "CHECK: Verify domain structures"

    elif expected_family == "hcol4":
        if len(c4) >= 2 and len(vwa) == 0:
            return "Network-forming Hcol4", "PASSED"
        return "Hcol4 Candidate", "CHECK: Verify domain structures"

    elif expected_family == "hcol5":
        gap = (col2[0]["start"] - col1[0]["end"]) if (col1 and col2) else 0
        if 1020 <= col1_len <= 1035 and gap > 20 and colfi:
            return "Fibrillar Hcol5", "PASSED"
        return "Hcol5 Candidate", "CHECK: Inter-domain spacing differs"

    elif expected_family == "hcol6":
        if len(wap) > 1 and len(vwa) > 1 and len(c4) >= 2:
            return "Network-forming Hcol6", "PASSED"
        return "Hcol6 Candidate", "CHECK: Verify domain structures"

    elif expected_family == "hcol7":
        if tspn_len > 0 and colfi:
            return "Fibrillar Hcol7", "PASSED"
        return "Hcol7 Candidate", "CHECK: Verify domain structures"

    elif expected_family == "hcol8":
        if not col2 or col2_len < 30 and len(col1) >= 3:
            return "Fibrillar Hcol8", "PASSED"
        return "Hcol8 Candidate", "WARNING: Unexpected secondary minor domain found"

    return "Unidentified Collagen", "UNVERIFIED"

def plot_domain_arcitecture(seq_length, domains, output, title = "Domain Architecture"):
    """
    domains: list of dicts with keys:
        - name: str
        - start: int
        - end: int
        - colour: str (optional)"""
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

        if width >= 60:
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

HCOL_TYPE = "hcol8"
with open(f"r_luteum_{HCOL_TYPE}_prediction.fasta", "r") as file:
    target_lines = file.readlines()
USER_SEQUENCE = "".join([line.strip() for line in target_lines if not line.startswith(">")])
TARGET_NAME = f"r_luteum_{HCOL_TYPE}"
SP_START = 0
SP_END = 0

if __name__ == "__main__":
    hmm_domains = scan_with_hmmer(USER_SEQUENCE)
    helix_domains = find_triple_helix(USER_SEQUENCE)
    tspn_domain = add_manual_domain("TSPN", start=0, end=0)
    signal_peptide = add_manual_domain("SP", SP_START, SP_END)

    all_domains = hmm_domains + helix_domains + [signal_peptide, tspn_domain]
    all_domains = sorted(all_domains, key=lambda d: d["start"])
    collagen_identity, validation_verdict = classify_and_validate(all_domains, HCOL_TYPE)

    print("=" * 60)
    print(f"Identity Result: {collagen_identity}")
    print(f"Pipeline Status: {validation_verdict}")
    print("=" * 60)

    seq_length = len(USER_SEQUENCE)
    plot_title = f"{TARGET_NAME.upper()} Architecture Map\nClassification: {collagen_identity}"
    plot_domain_arcitecture(seq_length, all_domains, TARGET_NAME, title=plot_title)
    print(f"Saved structural architecture map image to: {TARGET_NAME}.png")