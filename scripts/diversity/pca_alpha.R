#!/usr/bin/env Rscript

required <- c("phyloseq", "ggplot2", "RColorBrewer")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing required packages: ", paste(missing, collapse = ", "))
}

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(RColorBrewer)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript pca_alpha.R <abundance_table.csv> <output_dir>")
}

input_file <- args[1]
output_dir <- args[2]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data <- read.csv(input_file, row.names = 1, check.names = FALSE)
otu_matrix <- round(as.matrix(data))
storage.mode(otu_matrix) <- "numeric"

otu <- otu_table(otu_matrix, taxa_are_rows = TRUE)

sample_names <- colnames(otu_matrix)
sample_metadata <- data.frame(
  Sample = sample_names,
  Category = factor(sample_names),
  row.names = sample_names
)

physeq <- phyloseq(otu, sample_data(sample_metadata))
physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)

shape_map <- c(
  TVB1 = 24, TVB2 = 24, TVB3 = 24,
  MSTB1 = 22, MSTB2 = 22, MSTB3 = 22,
  MEPB1 = 21, MEPB2 = 21, MEPB3 = 21
)

palette_values <- colorRampPalette(brewer.pal(9, "Set1"))(length(sample_names))
names(palette_values) <- sample_names

p <- plot_richness(
  physeq,
  x = "Category",
  color = "Sample",
  measures = c("Chao1", "Shannon")
) +
  geom_point(aes(shape = Sample, fill = Sample), size = 4) +
  scale_shape_manual(values = shape_map) +
  scale_fill_manual(values = palette_values) +
  scale_color_manual(values = palette_values) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  guides(
    shape = guide_legend(title = "Category"),
    color = guide_legend(title = "Category"),
    fill = guide_legend(title = "Category")
  )

ggsave(
  file.path(output_dir, "alpha_diversity_plot.png"),
  plot = p, width = 8, height = 6, dpi = 300
)
