# Plot: PLINK2 PCA (PC1 vs PC2, colored by population) and scree plots,
# for the unpruned and LD-pruned runs from step 5.
#
# Run from data/vcf:  Rscript ../../scripts/plots/plot_pca.R
# Reads:    pca_output/<gene>_pca[_pruned].eigenvec, .eigenval, population_labels.txt
# Writes:   pca_output/<gene>_pca[_pruned]_plot.png and _scree.png

library(ggplot2)

genes <- c("CYP2C19", "HBB")  # PLINK2 PCA is autosomal only; G6PD is plotted by step 6
pca_dir <- "pca_output"

# population labels (sample ID -> population), built by pca_analysis.sh
labels <- read.table(file.path(pca_dir, "population_labels.txt"),
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# label_suffix: "" for unpruned, "_pruned" for the LD-pruned version
# plot_note: text describing this version, added to the plot title
plot_pca <- function(gene, label_suffix, plot_note) {

  file_prefix <- paste0(gene, "_pca", label_suffix)
  eigenvec_path <- file.path(pca_dir, paste0(file_prefix, ".eigenvec"))
  eigenval_path <- file.path(pca_dir, paste0(file_prefix, ".eigenval"))

  if (!file.exists(eigenvec_path)) {
    message("Skipping ", gene, label_suffix, " -- missing ", eigenvec_path)
    return(invisible(NULL))
  }

  # --- Load PCA coordinates ---
  # plink2's .eigenvec has a header row starting with '#' (e.g. "#IID PC1 PC2..."
  # or "#FID IID PC1 PC2..." depending on whether family IDs are present).
  # check.names=FALSE stops R from mangling "#IID" before we can clean it up.
  pcs <- read.table(eigenvec_path, header = TRUE, comment.char = "", check.names = FALSE)
  colnames(pcs) <- sub("^#", "", colnames(pcs))

  # --- Load variance explained per PC ---
  eigenvals <- scan(eigenval_path, quiet = TRUE)
  total_var <- sum(eigenvals)
  pc1_pct <- round(100 * eigenvals[1] / total_var, 1)
  pc2_pct <- round(100 * eigenvals[2] / total_var, 1)

  # --- Attach population labels ---
  pcs <- merge(pcs, labels, by = "IID")

  # --- PC1 vs PC2 scatter plot, colored by population ---
  p <- ggplot(pcs, aes(x = PC1, y = PC2, color = POP)) +
    geom_point(size = 2.5, alpha = 0.8) +
    labs(
      title = paste0(gene, " (", plot_note, "): PCA across GIH/PJL/BEB/ITU"),
      x = paste0("PC1 (", pc1_pct, "% variance explained)"),
      y = paste0("PC2 (", pc2_pct, "% variance explained)"),
      color = "Population"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(pca_dir, paste0(file_prefix, "_plot.png")), p, width = 7, height = 5.5, dpi = 300)

  # --- Scree plot: how much variance each PC explains ---
  scree_df <- data.frame(
    PC = factor(paste0("PC", seq_along(eigenvals)), levels = paste0("PC", seq_along(eigenvals))),
    pct_variance = 100 * eigenvals / total_var
  )

  scree <- ggplot(scree_df, aes(x = PC, y = pct_variance)) +
    geom_col(fill = "steelblue") +
    labs(
      title = paste0(gene, " (", plot_note, "): Variance explained per PC"),
      x = "Principal Component", y = "% variance explained"
    ) +
    theme_minimal(base_size = 13)

  ggsave(file.path(pca_dir, paste0(file_prefix, "_scree.png")), scree, width = 6, height = 4.5, dpi = 300)

  cat("\n", gene, "(", plot_note, "): PC1 =", pc1_pct, "%, PC2 =", pc2_pct,
      "%, combined =", pc1_pct + pc2_pct, "%\n")
}

for (gene in genes) {
  plot_pca(gene, "", "unpruned")
  plot_pca(gene, "_pruned", "LD-pruned")
}

cat("\nDone. Plots saved as PNGs in", pca_dir, "\n")
