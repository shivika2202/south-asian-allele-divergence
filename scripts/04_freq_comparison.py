"""Step 4: compare ALT allele frequencies across populations and rank by spread.

For each variant, max_freq_spread = highest minus lowest population ALT allele
frequency. Variants are sorted by spread so the most divergent sites come first.

Inputs  (cwd = data/vcf): plink_output/<GENE>_<POP>_filtered_freq.afreq
Outputs: freq_comparison/<GENE>_freq_comparison.csv
"""

from pathlib import Path

import pandas as pd

GENES = ["CYP2C19", "HBB", "G6PD"]
POPS = ["GIH", "PJL", "BEB", "ITU"]

PLINK_DIR = Path("plink_output")
OUT_DIR = Path("freq_comparison")


def load_afreq(path: Path, pop: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={"#CHROM": "CHROM"})

    if "POS" not in df.columns:
        df["POS"] = df["ID"].apply(lambda x: int(str(x).split(":")[1]))

    freq_col = "ALT_FREQS" if "ALT_FREQS" in df.columns else "ALT_FREQ"
    keep = df[["CHROM", "POS", "ID", "REF", "ALT", freq_col]].copy()
    keep[freq_col] = pd.to_numeric(keep[freq_col], errors="coerce")
    return keep.rename(columns={freq_col: f"{pop}_ALT_FREQ"})


def add_spread(merged: pd.DataFrame) -> pd.DataFrame:
    """Add max_freq_spread across whichever population columns are present, sorted descending."""
    freq_cols = [f"{pop}_ALT_FREQ" for pop in POPS if f"{pop}_ALT_FREQ" in merged.columns]
    merged = merged.copy()
    merged["max_freq_spread"] = merged[freq_cols].max(axis=1) - merged[freq_cols].min(axis=1)
    return merged.sort_values("max_freq_spread", ascending=False)


def build_comparison_for_gene(gene: str):
    merged = None
    for pop in POPS:
        path = PLINK_DIR / f"{gene}_{pop}_filtered_freq.afreq"
        if not path.exists():
            print(f"  WARNING: {path} not found, skipping {pop} for {gene}")
            continue

        df = load_afreq(path, pop)
        if merged is None:
            merged = df
        else:
            # Outer join: a variant filtered out in one population keeps an empty
            # frequency there rather than being dropped.
            merged = merged.merge(df.drop(columns=["REF", "ALT"]), on=["CHROM", "POS", "ID"], how="outer")

    if merged is None:
        print(f"  WARNING: no frequency files found for {gene}")
        return

    merged = add_spread(merged)
    out_path = OUT_DIR / f"{gene}_freq_comparison.csv"
    merged.to_csv(out_path, index=False)
    print(f"  -> {out_path} ({len(merged)} variants)")
    print(f"     Largest cross-population frequency spread: {merged['max_freq_spread'].iloc[0]:.3f} "
          f"(variant {merged['ID'].iloc[0]})")


def main():
    OUT_DIR.mkdir(exist_ok=True)
    for gene in GENES:
        print(f"=== {gene} ===")
        build_comparison_for_gene(gene)


if __name__ == "__main__":
    main()
