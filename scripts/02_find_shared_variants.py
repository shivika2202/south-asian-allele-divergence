from pathlib import Path

import pandas as pd

FREQ_DIR = Path("./plink_output")
SAMPLE_LIST_DIR = Path("../panel")
OUT_DIR = Path("./shared_variants_output")

GENES = ["CYP2C19", "HBB", "G6PD"]
POPS = ["GIH", "PJL", "BEB", "ITU"]


def load_freq_file(gene, pop, freq_dir=FREQ_DIR):
    path = freq_dir / f"{gene}_{pop}_filtered_freq.afreq"
    if not path.exists():
        print(f"  WARNING: missing {path}")
        return None

    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={"#CHROM": "CHROM", "ID": "VARIANT_ID"})
    df["POS"] = df["VARIANT_ID"].str.split(":").str[1].astype(int)
    return df[["CHROM", "POS"]]


def intersect_sites(pop_dfs):
    shared = pop_dfs[0]
    for df in pop_dfs[1:]:
        shared = shared.merge(df, on=["CHROM", "POS"], how="inner")
    return shared.drop_duplicates()


def find_shared_variants():
    print("=== Finding variants shared across all 4 populations ===\n")
    OUT_DIR.mkdir(exist_ok=True)
    summary = []

    for gene in GENES:
        print(f"Gene: {gene}")
        pop_dfs = {pop: df for pop in POPS if (df := load_freq_file(gene, pop)) is not None}

        if len(pop_dfs) < len(POPS):
            print(f"  Skipping {gene}: not all {len(POPS)} population files found.\n")
            continue

        shared = intersect_sites([pop_dfs[p] for p in POPS])
        out_path = OUT_DIR / f"{gene}_shared_variants.csv"
        shared.to_csv(out_path, index=False)

        print(f"  {len(shared)} variants shared across {'/'.join(POPS)}")
        print(f"  Saved to {out_path}\n")
        summary.append({"gene": gene, "n_shared_variants": len(shared)})

    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(OUT_DIR / "summary_shared_variant_counts.csv", index=False)
    print("Summary:")
    print(summary_df)


def sample_overlaps(pop_samples):
    """Return {(pop1, pop2): overlapping IDs} for every pair of populations that share a sample."""
    pops = list(pop_samples)
    return {(a, b): pop_samples[a] & pop_samples[b]
            for i, a in enumerate(pops) for b in pops[i + 1:]
            if pop_samples[a] & pop_samples[b]}


def check_sample_list_overlap():
    print("\n=== Sample list QC: checking for cross-population overlap ===\n")
    pop_samples = {}
    for pop in POPS:
        path = SAMPLE_LIST_DIR / f"{pop}_samples.txt"
        if not path.exists():
            print(f"  WARNING: missing {path}")
            continue
        with open(path) as f:
            pop_samples[pop] = {line.strip() for line in f if line.strip()}
        print(f"  {pop}: {len(pop_samples[pop])} samples")

    overlaps = sample_overlaps(pop_samples)
    for (a, b), ids in overlaps.items():
        print(f"  ISSUE: {len(ids)} sample(s) appear in both {a} and {b}: {sorted(ids)}")
    if not overlaps:
        print("  No overlapping sample IDs: each sample belongs to exactly one population.")


if __name__ == "__main__":
    find_shared_variants()
    check_sample_list_overlap()
