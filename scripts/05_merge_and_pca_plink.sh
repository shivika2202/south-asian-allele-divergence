#!/bin/bash
# Step 5: merge populations on the shared sites, then run PLINK2 PCA
# (unpruned and LD-pruned).
#
# Merging uses bcftools merge rather than plink2 --pmerge: plink2's merge is
# built for concatenation (same samples, different variant sets), while here
# the four populations are different samples typed at the same sites.
#
# Restricting to the shared site list from step 2 means every sample has a
# call at every site, so the merge introduces no artificial missingness.
#
# Note: PLINK2's PCA uses autosomal variants only, so G6PD (chrX) produces no
# PLINK2 PCA. Step 6 runs PCA for all three genes with X-aware encoding.
#
# Requires: bcftools, tabix, plink2
# Inputs  (cwd = data/vcf): <GENE>_<POP>.vcf.gz, shared_variants_output/, ../sex_info.txt, ../panel/
# Outputs: pca_output/<GENE>_merged.{vcf.gz,pgen,pvar,psam}
#          pca_output/<GENE>_pca[_pruned].{eigenvec,eigenval}
#          pca_output/population_labels.txt

set -e

VCF_DIR="."
SHARED_DIR="./shared_variants_output"
OUT_DIR="./pca_output"
SEX_FILE="../sex_info.txt"
NEW_ID_MAX_LEN=200
PCA_COMPONENTS=10

# LD pruning: 50-variant window, step 5, drop one of each pair with r^2 > 0.2.
WINDOW=50
STEP=5
R2_THRESHOLD=0.2

mkdir -p "$OUT_DIR"

GENES=(CYP2C19 HBB G6PD)
POPS=(GIH PJL BEB ITU)

for GENE in "${GENES[@]}"; do
  echo "=== ${GENE} ==="
  SHARED_CSV="${SHARED_DIR}/${GENE}_shared_variants.csv"
  if [ ! -f "$SHARED_CSV" ]; then
    echo "WARNING: missing $SHARED_CSV, skipping ${GENE}"
    continue
  fi

  # CHROM <tab> POS file for bcftools -R
  POS_FILE="${OUT_DIR}/${GENE}_shared_positions.txt"
  tail -n +2 "$SHARED_CSV" | awk -F',' '{print $1"\t"$2}' > "$POS_FILE"

  # Subset each population to the shared positions.
  MERGE_INPUTS=()
  for POP in "${POPS[@]}"; do
    VCF="${VCF_DIR}/${GENE}_${POP}.vcf.gz"
    SUBSET="${OUT_DIR}/${GENE}_${POP}_shared.vcf.gz"

    if [ ! -f "$VCF" ]; then
      echo "WARNING: missing $VCF, skipping ${GENE}/${POP}"
      continue
    fi

    bcftools view -R "$POS_FILE" "$VCF" -Oz -o "$SUBSET"
    bcftools index -t -f "$SUBSET"
    MERGE_INPUTS+=("$SUBSET")
  done

  if [ "${#MERGE_INPUTS[@]}" -ne 4 ]; then
    echo "WARNING: expected 4 population files for ${GENE}, got ${#MERGE_INPUTS[@]}, skipping merge"
    continue
  fi

  # One VCF per gene: same sites, all four populations' samples.
  MERGED_VCF="${OUT_DIR}/${GENE}_merged.vcf.gz"
  bcftools merge -Oz -o "$MERGED_VCF" "${MERGE_INPUTS[@]}"
  bcftools index -t -f "$MERGED_VCF"
  echo "Merged VCF: $(bcftools query -l "$MERGED_VCF" | wc -l | tr -d ' ') samples"

  MERGED_PREFIX="${OUT_DIR}/${GENE}_merged"
  plink2 --vcf "$MERGED_VCF" \
    --set-all-var-ids '@:#:$r:$a' \
    --new-id-max-allele-len "$NEW_ID_MAX_LEN" \
    --update-sex "$SEX_FILE" \
    --make-pgen \
    --out "$MERGED_PREFIX"

  # PLINK2 builds its PCA relationship matrix from autosomes only.
  if [ "$(grep -v '^#' "${MERGED_PREFIX}.pvar" | cut -f1 | sort -u)" = "X" ]; then
    echo "Skipping PLINK2 PCA for ${GENE} (chrX); see step 6"
    echo ""
    continue
  fi

  # PCA on the unpruned merged set, kept for comparison with the pruned run.
  plink2 --pfile "$MERGED_PREFIX" \
    --pca "$PCA_COMPONENTS" \
    --out "${OUT_DIR}/${GENE}_pca"
  echo "PCA (unpruned) done for ${GENE}"

  plink2 --pfile "$MERGED_PREFIX" \
    --indep-pairwise "$WINDOW" "$STEP" "$R2_THRESHOLD" \
    --out "${OUT_DIR}/${GENE}_merged_ld"

  N_BEFORE=$(grep -vc '^#' "${MERGED_PREFIX}.pvar" | tr -d ' ')
  N_KEPT=$(wc -l < "${OUT_DIR}/${GENE}_merged_ld.prune.in" | tr -d ' ')
  echo "LD pruning for ${GENE}: kept ${N_KEPT} of ${N_BEFORE} variants (r2 < ${R2_THRESHOLD})"

  PRUNED_PREFIX="${MERGED_PREFIX}_ldpruned"
  plink2 --pfile "$MERGED_PREFIX" \
    --extract "${OUT_DIR}/${GENE}_merged_ld.prune.in" \
    --make-pgen \
    --out "$PRUNED_PREFIX"

  plink2 --pfile "$PRUNED_PREFIX" \
    --pca "$PCA_COMPONENTS" \
    --out "${OUT_DIR}/${GENE}_pca_pruned"
  echo "PCA (LD-pruned) done for ${GENE}"
  echo ""
done

# Sample ID -> population, used for plot colors (step 6, R) and as the FST grouping (step 7).
LABELS_FILE="${OUT_DIR}/population_labels.txt"
echo -e "IID\tPOP" > "$LABELS_FILE"
for POP in "${POPS[@]}"; do
  awk -v pop="$POP" '{print $1"\t"pop}' "../panel/${POP}_samples.txt" >> "$LABELS_FILE"
done

echo ""
echo "Done. PCA results are in $OUT_DIR"
