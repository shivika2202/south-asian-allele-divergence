#!/bin/bash
# Step 1: per-population variant QC with PLINK2.
#
# For each gene x population VCF, applies a minor allele frequency filter and a
# per-variant missingness filter, computes allele frequencies on the retained
# variants, and records exactly which variants were removed.
#
# Inputs  (cwd = data/vcf):  <GENE>_<POP>.vcf.gz, ../sex_info.txt
# Outputs: plink_output/<GENE>_<POP>_filtered.{pgen,pvar,psam}
#          plink_output/<GENE>_<POP>_filtered_freq.afreq
#          plink_output/<GENE>_<POP>_discarded_variants.tsv
#          plink_output/filtration_summary.csv  (before/after counts)

set -e

VCF_DIR="."
OUT_DIR="./plink_output"
MAF_THRESHOLD=0.05
GENO_THRESHOLD=0.05
NEW_ID_MAX_LEN=200            # long structural-variant alleles otherwise break chr:pos:ref:alt IDs
SEX_FILE="../sex_info.txt"

mkdir -p "$OUT_DIR"

GENES=(CYP2C19 HBB G6PD)
POPS=(GIH PJL BEB ITU)

SUMMARY_LOG="${OUT_DIR}/filtration_summary.csv"
echo "gene,population,total_variants_before,retained_variants_after,discarded_variants" > "$SUMMARY_LOG"

for GENE in "${GENES[@]}"; do
  for POP in "${POPS[@]}"; do
    VCF="${VCF_DIR}/${GENE}_${POP}.vcf.gz"
    PREFIX="${OUT_DIR}/${GENE}_${POP}"

    if [ ! -f "$VCF" ]; then
      echo "WARNING: missing $VCF, skipping"
      continue
    fi

    echo "Processing ${GENE} / ${POP} with PLINK2 ..."

    # Unfiltered variant list, kept so the filtered-out set can be reported.
    plink2 --vcf "$VCF" \
      --set-all-var-ids '@:#:$r:$a' \
      --new-id-max-allele-len "$NEW_ID_MAX_LEN" \
      --update-sex "$SEX_FILE" \
      --make-just-pvar \
      --out "${PREFIX}_original"

    # Filtered set: missingness and MAF thresholds.
    plink2 --vcf "$VCF" \
      --set-all-var-ids '@:#:$r:$a' \
      --new-id-max-allele-len "$NEW_ID_MAX_LEN" \
      --update-sex "$SEX_FILE" \
      --geno "$GENO_THRESHOLD" \
      --maf "$MAF_THRESHOLD" \
      --make-pgen \
      --out "${PREFIX}_filtered"

    # Allele frequencies on the filtered set.
    plink2 --pfile "${PREFIX}_filtered" \
      --freq \
      --out "${PREFIX}_filtered_freq"

    # Per-site and per-sample missingness before filtering, for the QC record.
    plink2 --vcf "$VCF" \
      --set-all-var-ids '@:#:$r:$a' \
      --new-id-max-allele-len "$NEW_ID_MAX_LEN" \
      --update-sex "$SEX_FILE" \
      --missing \
      --out "${PREFIX}_missingness"

    # Discarded variants = IDs in the original set but not in the filtered set.
    grep -v '^##' "${PREFIX}_original.pvar"  | tail -n +2 > "${PREFIX}_original_body.tsv"
    grep -v '^##' "${PREFIX}_filtered.pvar"  | tail -n +2 > "${PREFIX}_filtered_body.tsv"

    cut -f3 "${PREFIX}_original_body.tsv" | sort > "${PREFIX}_original_ids.sorted"
    cut -f3 "${PREFIX}_filtered_body.tsv" | sort > "${PREFIX}_filtered_ids.sorted"

    comm -23 "${PREFIX}_original_ids.sorted" "${PREFIX}_filtered_ids.sorted" > "${PREFIX}_discarded_ids.txt"

    {
      grep -v '^##' "${PREFIX}_original.pvar" | head -1
      grep -F -f "${PREFIX}_discarded_ids.txt" "${PREFIX}_original_body.tsv" || true
    } > "${PREFIX}_discarded_variants.tsv"

    TOTAL=$(wc -l < "${PREFIX}_original_ids.sorted" | tr -d ' ')
    RETAINED=$(wc -l < "${PREFIX}_filtered_ids.sorted" | tr -d ' ')
    DISCARDED=$(wc -l < "${PREFIX}_discarded_ids.txt" | tr -d ' ')
    echo "${GENE},${POP},${TOTAL},${RETAINED},${DISCARDED}" >> "$SUMMARY_LOG"

    rm -f "${PREFIX}_original_body.tsv" "${PREFIX}_filtered_body.tsv" \
          "${PREFIX}_original_ids.sorted" "${PREFIX}_filtered_ids.sorted"
  done
done

echo ""
echo "Done. Filtered PLINK2 filesets and frequency tables are in $OUT_DIR"
echo "Before/after variant counts: $SUMMARY_LOG"
