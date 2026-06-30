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

def find_triple_helix(sequence, max_gap_allowed=1):
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
    Manually specify a domain region with known coordinates -- for SP and TSPN
    """
    return {
        "name": name,
        "start": start,
        "end": end,
        "colour": colour if colour else DOMAIN_COLOURS.get(name, DEFAULT_COLOUR)
    }

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


with open("target_seq.txt", "r") as file:
    target_lines = file.readlines()
USER_SEQUENCE = "".join([line.strip() for line in target_lines if not line.startswith(">")])
TARGET_NAME = "l_quadricornis_hcol8"
SP_START = 0
SP_END = 0

if __name__ == "__main__":
    hmm_domains = scan_with_hmmer(USER_SEQUENCE)
    helix_domains = find_triple_helix(USER_SEQUENCE)
    tspn_domain = add_manual_domain("TSPN", start=0, end=0)
    signal_peptide = add_manual_domain("SP", SP_START, SP_END)

    all_domains = hmm_domains + helix_domains + [signal_peptide, tspn_domain]
    all_domains = sorted(all_domains, key=lambda d: d["start"])
    seq_length = len(USER_SEQUENCE)
    plot_domain_arcitecture(seq_length, all_domains, TARGET_NAME)