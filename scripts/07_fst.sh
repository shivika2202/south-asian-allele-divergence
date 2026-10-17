#!/bin/bash
# Step 7: pairwise FST between all four populations, per gene, with PLINK2.
#
# Runs both the Hudson and the Weir-Cockerham estimators as a cross-check.
# Hudson is the primary estimate (less biased for small, unequal panels).
# On chrX (G6PD), PLINK2's Weir-Cockerham run completes but writes no
# estimates, so G6PD has a Hudson result only.
#
# Inputs  (cwd = data/vcf): pca_output/<GENE>_merged.{pgen,pvar,psam}, pca_output/population_labels.txt
# Outputs: fst_output/<GENE>_fst_<method>[.x].fst.summary (+ per-variant files)

set -e

GENES=("CYP2C19" "HBB" "G6PD")
PCA_DIR="pca_output"
OUT_DIR="fst_output"

mkdir -p "$OUT_DIR"

for GENE in "${GENES[@]}"; do
    echo "=== $GENE ==="
    for METHOD in hudson wc; do
        plink2 --pfile "${PCA_DIR}/${GENE}_merged" \
            --pheno "${PCA_DIR}/population_labels.txt" \
            --pheno-name POP \
            --fst POP method=${METHOD} report-variants \
            --out "${OUT_DIR}/${GENE}_fst_${METHOD}" \
          || echo "  ${METHOD} not available for ${GENE}"
    done
    echo ""
done

echo "Done. Pairwise FST summaries are in $OUT_DIR"
