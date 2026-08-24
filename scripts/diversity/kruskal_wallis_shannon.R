#!/usr/bin/env Rscript

required <- c("dplyr", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing required packages: ", paste(missing, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript kruskal_wallis_shannon.R <shannon_csv> <output_dir>")
}

input_file <- args[1]
output_dir <- args[2]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

shannon_data <- read.csv(input_file, check.names = FALSE)

shannon_data <- shannon_data %>%
  mutate(
    Group = case_when(
      grepl("TVB", Sample) ~ "TVB",
      grepl("MEPB", Sample) ~ "MEPB",
      grepl("MSTB", Sample) ~ "MSTB",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group))

run_pair <- function(df, g1, g2) {
  x <- df %>% filter(Group %in% c(g1, g2))
  kruskal.test(Shannon_Index ~ Group, data = x)
}

kruskal_tvb_mepb <- run_pair(shannon_data, "TVB", "MEPB")
kruskal_tvb_mstb <- run_pair(shannon_data, "TVB", "MSTB")
kruskal_mepb_mstb <- run_pair(shannon_data, "MEPB", "MSTB")
kruskal_all_groups <- kruskal.test(Shannon_Index ~ Group, data = shannon_data)

results_df <- data.frame(
  Comparison = c("TVB vs MEPB", "TVB vs MSTB", "MEPB vs MSTB", "All Groups"),
  Statistic = c(
    unname(kruskal_tvb_mepb$statistic),
    unname(kruskal_tvb_mstb$statistic),
    unname(kruskal_mepb_mstb$statistic),
    unname(kruskal_all_groups$statistic)
  ),
  p_value = c(
    kruskal_tvb_mepb$p.value,
    kruskal_tvb_mstb$p.value,
    kruskal_mepb_mstb$p.value,
    kruskal_all_groups$p.value
  )
)

write.csv(
  results_df,
  file.path(output_dir, "kruskal_wallis_shannon_results.csv"),
  row.names = FALSE
)

capture.output(
  kruskal_tvb_mepb,
  kruskal_tvb_mstb,
  kruskal_mepb_mstb,
  kruskal_all_groups,
  file = file.path(output_dir, "kruskal_wallis_shannon_results.txt")
)

p <- ggplot(shannon_data, aes(x = Group, y = Shannon_Index, color = Group)) +
  geom_boxplot(fill = NA, linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Shannon diversity across groups",
    x = "Group",
    y = "Shannon Index"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  file.path(output_dir, "shannon_group_boxplot.png"),
  plot = p, width = 8, height = 6, dpi = 300
)
