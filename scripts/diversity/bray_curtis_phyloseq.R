#!/usr/bin/env Rscript

required <- c("phyloseq", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing required packages: ", paste(missing, collapse = ", "))
}

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript bray_curtis_phyloseq.R <abundance_table.csv> <output_dir>")
}

input_file <- args[1]
output_dir <- args[2]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data <- read.csv(input_file, check.names = FALSE)
otu_matrix <- as.matrix(data[, -1, drop = FALSE])
rownames(otu_matrix) <- data[, 1]
storage.mode(otu_matrix) <- "numeric"

otu <- otu_table(otu_matrix, taxa_are_rows = TRUE)

sample_names <- colnames(otu_matrix)
sample_metadata <- data.frame(
  Sample = sample_names,
  Category = factor(sample_names),
  row.names = sample_names
)

physeq <- phyloseq(otu, sample_data(sample_metadata))

bray_curtis_dist <- distance(physeq, method = "bray")
ordination <- ordinate(physeq, method = "PCoA", distance = bray_curtis_dist)

shape_map <- c(
  TVB1 = 24, TVB2 = 24, TVB3 = 24,
  MSTB1 = 22, MSTB2 = 22, MSTB3 = 22,
  MEPB1 = 21, MEPB2 = 21, MEPB3 = 21
)

palette_values <- c(
  TVB1 = "#E41A1C", TVB2 = "#377EB8", TVB3 = "#4DAF4A",
  MSTB1 = "#984EA3", MSTB2 = "#FF7F00", MSTB3 = "#DAA520",
  MEPB1 = "#A65628", MEPB2 = "#F781BF", MEPB3 = "#999999"
)

p <- plot_ordination(
  physeq, ordination,
  color = "Category", shape = "Category"
) +
  geom_point(aes(fill = Category), size = 4, alpha = 0.7) +
  scale_color_manual(values = palette_values) +
  scale_fill_manual(values = palette_values) +
  scale_shape_manual(values = shape_map) +
  theme_bw() +
  guides(
    shape = guide_legend(title = "Category"),
    color = guide_legend(title = "Category"),
    fill = guide_legend(title = "Category")
  )

ggsave(
  file.path(output_dir, "bray_curtis_pcoa.png"),
  plot = p, width = 8, height = 6, dpi = 300
)

write.csv(
  as.matrix(bray_curtis_dist),
  file.path(output_dir, "bray_curtis_distance_matrix.csv"),
  row.names = TRUE
)
