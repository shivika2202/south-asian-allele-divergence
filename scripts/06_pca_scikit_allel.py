"""Step 6: PCA for all three genes with scikit-allel, using X-aware genotype encoding.

PLINK2's PCA excludes chrX, so G6PD needs its own implementation. Running all
three genes through the same method keeps the results comparable. For G6PD,
male samples are hemizygous at non-PAR chrX loci and are encoded on the same
0/1/2 scale as diploid calls; for the autosomal genes the encoding is the
ordinary diploid alt-allele count.

Inputs  (cwd = data/vcf): pca_output/<GENE>_merged.vcf.gz, pca_output/<GENE>_merged.psam,
                          pca_output/population_labels.txt
Outputs: pca_output/<GENE>_pca_scikit.eigenvec.csv, .eigenval.csv, _plot.png
"""

from pathlib import Path

import allel
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

GENES = ["CYP2C19", "HBB", "G6PD"]
X_LINKED_GENES = {"G6PD"}

OUT_DIR = Path("pca_output")
LABELS_PATH = OUT_DIR / "population_labels.txt"

POP_COLORS = {
    "GIH": "#4C72B0",
    "PJL": "#55A868",
    "BEB": "#C44E52",
    "ITU": "#DD8452",
}


def load_sex_map(gene: str) -> dict:
    """Sample ID -> PLINK sex code (1 = male, 2 = female), from the merged .psam.

    The .psam is used because it matches the merged VCF's sample set exactly.
    """
    psam_path = OUT_DIR / f"{gene}_merged.psam"
    sex_df = pd.read_csv(psam_path, sep=r"\s+", engine="python")
    sex_df.columns = [c.lstrip("#") for c in sex_df.columns]
    return dict(zip(sex_df["IID"], sex_df["SEX"]))


def to_allele_counts(gt, is_x_linked: bool, is_male: np.ndarray = None) -> np.ndarray:
    """Convert (variants x samples x ploidy) genotypes to a (variants x samples) alt-allele count matrix.

    Missing calls become -1. scikit-allel pads both true missing calls and the
    second slot of haploid calls with -1, so ploidy cannot be read from the
    array. For X-linked genes, known sample sex decides: males use the first
    allele only, doubled onto the 0/1/2 scale; females use both alleles.
    """
    gt = np.asarray(gt)
    n_variants, n_samples = gt.shape[0], gt.shape[1]
    counts = np.full((n_variants, n_samples), -1, dtype=np.int16)

    for s in range(n_samples):
        sample_gt = gt[:, s, :]
        if is_x_linked and is_male is not None and is_male[s]:
            first_allele = sample_gt[:, 0]
            counts[:, s] = np.where(first_allele >= 0, first_allele * 2, -1)
        else:
            counts[:, s] = np.where(np.all(sample_gt >= 0, axis=1), sample_gt.sum(axis=1), -1)

    return counts


def drop_unusable_variants(ac: np.ndarray):
    """Remove variants with any missing call, then variants monomorphic in the pooled sample.

    Monomorphic variants must go because Patterson scaling divides by
    sqrt(p(1-p)), which is zero when p is 0 or 1. A variant can pass the MAF
    filter within each population and still be (near) monomorphic once pooled.
    Returns (filtered matrix, n_dropped_missing, n_dropped_monomorphic).
    """
    complete = np.all(ac != -1, axis=1)
    ac_complete = ac[complete, :]
    af = ac_complete.mean(axis=1) / 2.0
    polymorphic = (af > 1e-9) & (af < 1 - 1e-9)
    return ac_complete[polymorphic, :], int((~complete).sum()), int((~polymorphic).sum())


def run_pca_for_gene(gene: str, labels_df: pd.DataFrame):
    vcf_path = OUT_DIR / f"{gene}_merged.vcf.gz"
    if not vcf_path.exists():
        print(f"  WARNING: {vcf_path} not found, skipping {gene}")
        return

    print(f"=== {gene} ===")
    callset = allel.read_vcf(str(vcf_path), fields=["samples", "calldata/GT"])
    samples = callset["samples"]
    gt = allel.GenotypeArray(callset["calldata/GT"])
    print(f"  {gt.shape[0]} variants, {gt.shape[1]} samples loaded")

    is_x = gene in X_LINKED_GENES
    is_male = None
    if is_x:
        sex_map = load_sex_map(gene)
        is_male = np.array([sex_map.get(s) == 1 for s in samples])
        n_male = is_male.sum()
        n_unknown = sum(1 for s in samples if s not in sex_map)
        print(f"  {n_male} male, {len(samples) - n_male - n_unknown} female, {n_unknown} unknown sex")

    ac = to_allele_counts(gt, is_x_linked=is_x, is_male=is_male)
    ac_poly, n_missing, n_mono = drop_unusable_variants(ac)
    print(f"  {n_missing} variant(s) dropped for missingness, {n_mono} for being monomorphic after merging; "
          f"{ac_poly.shape[0]} used for PCA")

    if ac_poly.shape[0] < 2:
        print(f"  WARNING: fewer than 2 usable variants for {gene}, skipping PCA")
        return

    n_components = min(5, ac_poly.shape[0] - 1, ac_poly.shape[1] - 1)
    coords, model = allel.pca(ac_poly, n_components=n_components, scaler="patterson")

    scores = pd.DataFrame(coords, columns=[f"PC{i+1}" for i in range(coords.shape[1])])
    scores.insert(0, "sample_id", samples)
    scores = scores.merge(labels_df, on="sample_id", how="left")
    if scores["population"].isna().any():
        print(f"  WARNING: {scores['population'].isna().sum()} sample(s) had no population label")

    scores.to_csv(OUT_DIR / f"{gene}_pca_scikit.eigenvec.csv", index=False)
    pd.DataFrame({
        "PC": [f"PC{i+1}" for i in range(len(model.explained_variance_ratio_))],
        "variance_explained": model.explained_variance_ratio_,
    }).to_csv(OUT_DIR / f"{gene}_pca_scikit.eigenval.csv", index=False)

    print(f"  PC1 variance explained: {model.explained_variance_ratio_[0]:.3f}")
    if n_components > 1:
        print(f"  PC2 variance explained: {model.explained_variance_ratio_[1]:.3f}")

    fig, ax = plt.subplots(figsize=(7, 6))
    for pop, color in POP_COLORS.items():
        subset = scores[scores["population"] == pop]
        ax.scatter(subset["PC1"], subset.get("PC2", 0), label=pop, color=color, alpha=0.75, s=40)
    ax.set_xlabel(f"PC1 ({model.explained_variance_ratio_[0]*100:.1f}%)")
    if n_components > 1:
        ax.set_ylabel(f"PC2 ({model.explained_variance_ratio_[1]*100:.1f}%)")
    ax.set_title(f"{gene} PCA (scikit-allel)")
    ax.legend(title="Population")
    fig.tight_layout()
    fig.savefig(OUT_DIR / f"{gene}_pca_scikit_plot.png", dpi=150)
    plt.close(fig)
    print(f"  -> {OUT_DIR / f'{gene}_pca_scikit_plot.png'}\n")


def main():
    labels_df = pd.read_csv(LABELS_PATH, sep="\t").rename(columns={"IID": "sample_id", "POP": "population"})
    for gene in GENES:
        run_pca_for_gene(gene, labels_df)


if __name__ == "__main__":
    main()
