"""Step 3: label known clinical variants (rsID, star allele) in a frequency file.

Matches on CHROM + POS against data/variant_reference_lookup.csv. Accepts both
vcftools --freq output (.frq, header starts with "CHROM") and PLINK2 --freq
output (.afreq, header starts with "#CHROM"), so the two QC tools' results can
be annotated with the same table and compared.

Usage:
    python 03_annotate_variants.py <freq_file> <lookup_csv> <output_tsv>
"""

import sys

import pandas as pd


def detect_filetype(path):
    with open(path) as f:
        header = f.readline()
    if header.startswith("#CHROM"):
        return "plink"
    if header.startswith("CHROM"):
        return "vcftools"
    raise ValueError(
        f"Unrecognized file format, first line was:\n{header}\n"
        "Expected a vcftools .frq (starts with 'CHROM') or a PLINK2 .afreq (starts with '#CHROM') file."
    )


def load_lookup(lookup_path):
    lookup = pd.read_csv(lookup_path, dtype={"CHROM": str})
    lookup["POS"] = lookup["POS"].astype(int)
    return lookup


def load_vcftools_freq(path):
    """Parse a vcftools .frq file.

    The header has a single {ALLELE:FREQ} placeholder, but each row carries one
    ALLELE:FREQ field per allele (two for biallelic sites, more for multiallelic
    ones), so a fixed-column read would misalign rows. Rows are parsed manually;
    the first two alleles become REF/ALT and the raw fields are kept.
    """
    rows = []
    with open(path) as f:
        f.readline()
        for line in f:
            fields = line.rstrip("\n").split("\t")
            chrom, pos, n_alleles, n_chr = fields[0:4]
            allele_freq_fields = fields[4:]

            alleles, freqs = [], []
            for af in allele_freq_fields:
                allele, freq = af.split(":")
                alleles.append(allele)
                freqs.append(float(freq))

            rows.append({
                "CHROM": chrom,
                "POS": int(pos),
                "N_ALLELES": int(n_alleles),
                "N_CHR": int(n_chr),
                "REF_ALLELE": alleles[0] if len(alleles) > 0 else None,
                "REF_FREQ": freqs[0] if len(freqs) > 0 else None,
                "ALT_ALLELE": alleles[1] if len(alleles) > 1 else None,
                "ALT_FREQ": freqs[1] if len(freqs) > 1 else None,
                "ALL_ALLELES_RAW": ";".join(allele_freq_fields),
            })
    return pd.DataFrame(rows)


def load_plink_afreq(path):
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={"#CHROM": "CHROM"})
    df["CHROM"] = df["CHROM"].astype(str)
    df["POS"] = df["ID"].apply(lambda x: int(str(x).split(":")[1]))
    return df


def annotate(freq_path, lookup_path, output_path):
    filetype = detect_filetype(freq_path)
    lookup = load_lookup(lookup_path)
    df = load_plink_afreq(freq_path) if filetype == "plink" else load_vcftools_freq(freq_path)

    merged = df.merge(lookup, on=["CHROM", "POS"], how="left")

    n_annotated = merged["rsID"].notna().sum()
    print(f"[{filetype}] {n_annotated} of {len(merged)} variants matched a known entry in {lookup_path}")

    merged.to_csv(output_path, sep="\t", index=False)
    print(f"Annotated file written to: {output_path}")
    return merged


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python 03_annotate_variants.py <freq_file> <lookup_csv> <output_tsv>")
        sys.exit(1)
    annotate(sys.argv[1], sys.argv[2], sys.argv[3])
