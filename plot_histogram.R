
#!/usr/bin/env Rscript

suppressMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript plot_histogram.R <input.tsv> <percent_for_tree>")
}

input_file <- args[1]
selected_percent <- as.numeric(args[2])
output_prefix <- tools::file_path_sans_ext(input_file)

# Read histogram data
data <- read_tsv(input_file, col_types = cols(
  Percent_VCFs = col_integer(),
  Num_Sites = col_integer()
))

# Filter out invalid or out-of-bounds rows
data <- data %>%
  filter(
    !is.na(Percent_VCFs),
    !is.na(Num_Sites),
    Percent_VCFs >= 0,
    Percent_VCFs <= 100
  )

# Compute Y-axis padding
y_max <- max(data$Num_Sites, na.rm = TRUE)
y_pad <- round(y_max * 0.05)
y_limit <- y_max + y_pad

# Get y value for red dot, if available
dot_y <- data %>%
  filter(Percent_VCFs == selected_percent) %>%
  pull(Num_Sites)

# Add padding so the dot is slightly above the bar
dot_y_pos <- dot_y + y_max * 0.02

# Main plot
p <- ggplot(data, aes(x = Percent_VCFs, y = Num_Sites)) +
  geom_col(fill = "steelblue") +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 10)
  ) +
  scale_y_continuous(
    limits = c(0, y_limit),
    expand = c(0, 0)
  ) +
  labs(
    title = "Phylogenetically informative sites",
    x = "Percent of VCFs",
    y = "Number of Sites"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 0)
  )

# Add the red dot above the corresponding bar (if found)
if (length(dot_y) == 1 && !is.na(dot_y)) {
  p <- p + annotate("point", x = selected_percent, y = dot_y_pos, color = "red", size = 3)
}

# Save plots
ggsave(paste0(output_prefix, ".png"), p, width = 10, height = 6, dpi = 300)
ggsave(paste0(output_prefix, ".pdf"), p, width = 10, height = 6)

message("Saved plots to: ", output_prefix, ".png and .pdf")
