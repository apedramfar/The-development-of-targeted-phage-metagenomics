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
  stop("Usage: Rscript shannon_phyloseq.R <abundance_table.csv> <output_dir> [metadata.csv]")
}

input_file <- args[1]
output_dir <- args[2]
metadata_file <- if (length(args) >= 3) args[3] else NA_character_

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data <- read.csv(input_file, check.names = FALSE)
data <- data[, !grepl("^Unnamed|^X$", colnames(data)), drop = FALSE]

otu_matrix <- as.matrix(data[, -1, drop = FALSE])
rownames(otu_matrix) <- data[, 1]
storage.mode(otu_matrix) <- "numeric"

if (all(rowSums(otu_matrix, na.rm = TRUE) == 0)) {
  stop("The OTU matrix contains only zero values.")
}

otu_matrix[is.na(otu_matrix)] <- 0
otu_matrix <- round(otu_matrix)

otu <- otu_table(otu_matrix, taxa_are_rows = TRUE)

if (!is.na(metadata_file) && file.exists(metadata_file)) {
  sample_metadata <- read.csv(metadata_file, check.names = FALSE)
  if (!"Sample" %in% colnames(sample_metadata)) {
    stop("Metadata file must contain a 'Sample' column.")
  }
  rownames(sample_metadata) <- sample_metadata$Sample
} else {
  sample_names <- colnames(otu_matrix)
  sample_metadata <- data.frame(
    Sample = sample_names,
    Category = factor(sample_names),
    row.names = sample_names
  )
}

physeq <- phyloseq(otu, sample_data(sample_metadata))

shannon_div <- estimate_richness(physeq, measures = "Shannon")
shannon_data <- data.frame(
  Sample = rownames(shannon_div),
  Shannon_Index = shannon_div$Shannon
)

write.csv(
  shannon_data,
  file.path(output_dir, "shannon_diversity_output.csv"),
  row.names = FALSE
)

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

p <- plot_richness(
  physeq,
  x = "Category",
  color = "Sample",
  measures = "Shannon"
) +
  geom_point(aes(shape = Sample, fill = Sample), size = 4) +
  scale_fill_manual(values = palette_values) +
  scale_color_manual(values = palette_values) +
  scale_shape_manual(values = shape_map) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  xlab("") +
  guides(
    shape = guide_legend(title = "Category"),
    color = guide_legend(title = "Category"),
    fill = guide_legend(title = "Category")
  )

ggsave(
  file.path(output_dir, "shannon_diversity_plot.png"),
  plot = p, width = 8, height = 6, dpi = 300
)
