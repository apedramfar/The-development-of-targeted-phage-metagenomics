#!/usr/bin/env Rscript

required <- c("FSA", "dplyr", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing required packages: ", paste(missing, collapse = ", "))
}

suppressPackageStartupMessages({
  library(FSA)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript dunn_test_chao1.R <chao1_csv> <output_dir>")
}

input_file <- args[1]
output_dir <- args[2]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data <- read.csv(input_file, check.names = FALSE)

group_mapping <- c(
  TVB1 = "TVB", TVB2 = "TVB", TVB3 = "TVB",
  MEPB1 = "MEPB", MEPB2 = "MEPB", MEPB3 = "MEPB"
)

data <- data %>%
  mutate(Group = unname(group_mapping[Sample]))

subset_data <- data %>%
  filter(Group %in% c("TVB", "MEPB"))

if (nrow(subset_data) == 0) {
  stop("No TVB/MEPB samples were found in the input data.")
}

dunn_test <- dunnTest(
  Chao1_Index ~ Group,
  data = subset_data,
  method = "bonferroni"
)

write.csv(
  dunn_test$res,
  file.path(output_dir, "dunn_posthoc_chao1.csv"),
  row.names = FALSE
)

p_value <- dunn_test$res$P.adj[1]

p <- ggplot(subset_data, aes(x = Group, y = Chao1_Index, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "Chao1 diversity: TVB vs MEPB",
    x = "Group",
    y = "Chao1 Index"
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(subset_data$Chao1_Index, na.rm = TRUE),
    label = paste0("adjusted p = ", signif(p_value, 3)),
    size = 4
  )

ggsave(
  file.path(output_dir, "chao1_dunn_boxplot.png"),
  plot = p, width = 8, height = 6, dpi = 300
)
