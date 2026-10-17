#!/bin/bash
# Runs the full analysis in order. Expects the 12 per-population VCFs
# (<GENE>_<POP>.vcf.gz + .tbi) in data/vcf/; see "Input data" in the README.
#
#   bash run_pipeline.sh            # analysis only
#   bash run_pipeline.sh --plots    # also render the R figures

set -e

ROOT="$(cd "$(dirname "$0")" && pwd)"
SCRIPTS="$ROOT/scripts"
PYTHON="${PYTHON:-python3}"

cd "$ROOT/data/vcf"

bash "$SCRIPTS/01_plink_qc.sh"
"$PYTHON" "$SCRIPTS/02_find_shared_variants.py"

mkdir -p annotated
for GENE in CYP2C19 HBB G6PD; do
  for POP in GIH PJL BEB ITU; do
    "$PYTHON" "$SCRIPTS/03_annotate_variants.py" \
      "plink_output/${GENE}_${POP}_filtered_freq.afreq" \
      ../variant_reference_lookup.csv \
      "annotated/${GENE}_${POP}_annotated_plink.tsv"
  done
done

"$PYTHON" "$SCRIPTS/04_freq_comparison.py"
bash "$SCRIPTS/05_merge_and_pca_plink.sh"
"$PYTHON" "$SCRIPTS/06_pca_scikit_allel.py"
bash "$SCRIPTS/07_fst.sh"

if [ "$1" = "--plots" ]; then
  Rscript "$SCRIPTS/plots/plot_pca.R"
  Rscript "$SCRIPTS/plots/plot_freq_comparison.R"
fi

echo "Pipeline complete. Outputs are in data/vcf/."
