import os
import sys
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import pandas as pd

tsv_file = snakemake.input.tsv
output_csv = snakemake.output.csv
hcol = snakemake.wildcards.hcol
plot_dir = f"results/interpro/plots/{hcol}"
os.makedirs(plot_dir, exist_ok=True)

if os.path.getsize(tsv_file) == 0:
    with open(output_csv, "w") as f:
        pass
    sys.exit(0)

DOMAIN_MAP = {
    "triple_helix": {"PF01391", "IPR008160"},
    "colfi": {"PF01410", "IPR000885"},
    "wap": {"PF00095", "IPR008197", "IPR036645"},
    "vwf": {"PF00092", "IPR002035", "IPR036465"},
    "tspn": {"PF18487", "IPR036383"},
    "c4": {"PF01413", "IPR001442"},
}

DOMAIN_COLOURS = {
    "triple_helix": "#4daf4a",
    "colfi": "#377eb8",
    "wap": "#e41a1c",
    "vwf": "#984ea3",
    "tspn": "#ff7f00",
    "c4": "#a65628",
}

CLASS_RULES = {
    "hcol1/5/8": {
        "triple_helix": 1,
        "colfi": 1,
        "wap": 0,
        "vwf": 0,
        "tspn": 0,
        "c4": 0,
    },
    "hcol2a/b": {
        "triple_helix": 1,
        "colfi": 1,
        "wap": 1,
        "vwf": 0,
        "tspn": 0,
        "c4": 0,
    },
    "hcol3": {
        "triple_helix": 1,
        "colfi": 1,
        "wap": 1,
        "vwf": 1,
        "tspn": 0,
        "c4": 0,
    },
    "hcol4": {
        "triple_helix": 1,
        "colfi": 0,
        "wap": 0,
        "vwf": 0,
        "tspn": 0,
        "c4": 1,
    },
    "hcol6": {
        "triple_helix": 1,
        "colfi": 0,
        "wap": 1,
        "vwf": 1,
        "tspn": 0,
        "c4": 1,
    },
    "hcol7": {
        "triple_helix": 1,
        "colfi": 1,
        "wap": 0,
        "vwf": 0,
        "tspn": 1,
        "c4": 0,
    },
}


def calculate_true_length(intervals):
    if not intervals:
        return 0
    intervals.sort(key=lambda x: x[0])
    merged = [intervals[0]]
    for current in intervals[1:]:
        prev_start, prev_end = merged[-1]
        curr_start, curr_end = current
        if curr_start <= prev_end:
            merged[-1] = (prev_start, max(prev_end, curr_end))
        else:
            merged.append(current)
    return sum(end - start + 1 for start, end in merged)


def detect_domain(row, domain_set):
    return (
        str(row.get("signature_accession", "")).upper() in domain_set
        or str(row.get("interpro_accession", "")).upper() in domain_set
    )


def classify_from_rules(features):
    for hcol_name, rule in CLASS_RULES.items():
        if all(features[k] == rule[k] for k in rule):
            return hcol_name
    return "unknown - check SMART"


def classify_protein(group):
    current_id = str(group.name)
    intervals = {k: [] for k in DOMAIN_MAP}

    for _, row in group.iterrows():
        try:
            coords = (int(row["start"]), int(row["end"]))
        except:
            continue
        for domain, accessions in DOMAIN_MAP.items():
            if detect_domain(row, accessions):
                intervals[domain].append(coords)

    merged_intervals = {}
    for domain, coords in intervals.items():
        if coords:
            coords.sort(key=lambda x: x[0])
            merged = [coords[0]]
            for start, end in coords[1:]:
                prev_start, prev_end = merged[-1]
                if start <= prev_end:
                    merged[-1] = (prev_start, max(prev_end, end))
                else:
                    merged.append((start, end))
            merged_intervals[domain] = merged
        else:
            merged_intervals[domain] = []

    lengths = {k: calculate_true_length(v) for k, v in merged_intervals.items()}
    presence = {k: int(lengths[k] > 0) for k in lengths}

    hcol_type = classify_from_rules(presence)

    return pd.Series(
        {
            "Sequence": current_id,
            "Total_hits": len(group),
            **{f"{k}_present": bool(presence[k]) for k in presence},
            **{f"{k}_length": lengths[k] for k in lengths},
            "potential_hcol_type": hcol_type,
            "intervals": merged_intervals,
        }
    )


def plot_domain_architecture(protein_id, intervals, outpath):
    protein_length = max(
        (end for coords in intervals.values() for (_, end) in coords), default=100
    )
    fig, ax = plt.subplots(figsize=(12, 2))
    # Base protein line
    ax.add_patch(
        patches.Rectangle(
            (0, 0.4), protein_length, 0.2, edgecolor="black", facecolor="#dddddd"
        )
    )

    for domain, coords in intervals.items():
        colour = DOMAIN_COLOURS.get(domain, "#999999")
        for start, end in coords:
            ax.add_patch(
                patches.Rectangle(
                    (start, 0.4),
                    end - start + 1,
                    0.2,
                    edgecolor="black",
                    facecolor=colour,
                )
            )

    ax.set_xlim(0, protein_length)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xlabel("Amino acid position")
    ax.set_title(f"Domain architecture: {protein_id}")

    handles = []
    labels = set()
    for domain, colour in DOMAIN_COLOURS.items():
        if intervals.get(domain):
            if domain not in labels:
                handles.append(patches.Patch(color=colour, label=domain))
                labels.add(domain)

    if handles:
        ax.legend(handles=handles, loc="upper right", bbox_to_anchor=(1.15, 1))
    plt.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


df = pd.read_csv(tsv_file, sep="\t")
df["protein_accession"] = df["protein_accession"].astype(str).str.strip()

results = (
    df.groupby("protein_accession", group_keys=False)
    .apply(classify_protein)
    .reset_index(drop=True)
)
results.to_csv(output_csv, index=False)

for _, row in results.iterrows():
    protein_id = row["Sequence"]
    plot_domain_architecture(
        protein_id,
        row["intervals"],
        f"{plot_dir}/{protein_id}_architecture.png",
    )