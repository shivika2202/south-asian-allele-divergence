# Plot: top variants by cross-population ALT allele frequency spread, per gene (step 4 output).
#
# Run from data/vcf:  Rscript ../../scripts/plots/plot_freq_comparison.R
# Writes: freq_comparison/<gene>_freq_comparison_barplot.png

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)

freq_dir <- "freq_comparison"
pops <- c("GIH", "PJL", "BEB", "ITU")
pop_colors <- c(GIH = "#4C72B0", PJL = "#DD8452", BEB = "#55A868", ITU = "#C44E52")

plot_gene_freqs <- function(gene, file, top_n = 8) {
  df <- read_csv(file, show_col_types = FALSE)
  
  freq_cols <- paste0(pops, "_ALT_FREQ")
  
  df <- df %>%
    mutate(across(all_of(freq_cols), as.numeric)) %>%
    filter(if_any(all_of(freq_cols), ~ !is.na(.))) %>%
    arrange(desc(max_freq_spread)) %>%
    # Variants in the same LD block share a frequency profile; keep one per profile.
    distinct(across(all_of(freq_cols)), .keep_all = TRUE) %>%
    slice_head(n = top_n) %>%
    mutate(variant_label = paste0(POS, "\n", REF, ">", ALT))

  df_long <- df %>%
    select(variant_label, all_of(freq_cols)) %>%
    pivot_longer(cols = all_of(freq_cols), names_to = "pop", values_to = "freq") %>%
    mutate(pop = str_remove(pop, "_ALT_FREQ"),
           pop = factor(pop, levels = pops),
           variant_label = factor(variant_label, levels = unique(df$variant_label)))

  p <- ggplot(df_long, aes(x = variant_label, y = freq, fill = pop)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    scale_fill_manual(values = pop_colors, name = "Population") +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(
      title = paste0(gene, ": Top ", nrow(df), " variants by cross-population frequency spread"),
      x = "Variant (POS, REF>ALT)",
      y = "ALT allele frequency"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 8),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )

  ggsave(file.path(freq_dir, paste0(gene, "_freq_comparison_barplot.png")), p, width = 10, height = 6, dpi = 300)
  p
}

plot_gene_freqs("CYP2C19", file.path(freq_dir, "CYP2C19_freq_comparison.csv"))
plot_gene_freqs("HBB", file.path(freq_dir, "HBB_freq_comparison.csv"))
plot_gene_freqs("G6PD", file.path(freq_dir, "G6PD_freq_comparison.csv"))