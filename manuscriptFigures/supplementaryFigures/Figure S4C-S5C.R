library(ggplot2)
library(dplyr)
library(scales)  # For formatting y-axis
library(svglite)
library(patchwork)

# Paths
temp_csv_forplot <- "zt_temp/PrepPlots_20250429"
save_path <- "Modified Figures"

# Load data
# hg38_variant_counts_per_chromosome.csv
# chm13_variant_counts_per_chromosome.csv
variant_df <- read.csv(file.path(temp_csv_forplot, "hg38_variant_counts_per_chromosome.csv"))

# Ensure chromosome order
chromosome_order <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
variant_df$Chromosome <- factor(variant_df$Chromosome, levels = chromosome_order)

# Custom colors
variant_colors <- c(SNV = "#0077BB", Indel = "#EE7733", SV = "#009988")

## 2025/05/05 new ver plot SVs together with SNVs and indels
variant_df$Variant_Type <- factor(variant_df$Variant_Type, levels = c("SV", "Indel", "SNV"))

## 2026/04/01 plot break y axis
# Upper panel: y from 200 to max
plot_upper <- ggplot(variant_df, aes(x = Chromosome, y = Variant_Count, fill = Variant_Type)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_hline(yintercept = 200, color = "red", linetype = "dashed", linewidth = 1) +
  scale_fill_manual(values = variant_colors) +
  coord_cartesian(ylim = c(200, NA)) +
  scale_y_continuous(labels = comma, expand = c(0, 0)) +
  labs(title = "Variants by Chromosome", x = NULL, y = "Variant Count") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "darkgray"))

# Lower panel: y from 0 to 150
plot_lower <- ggplot(variant_df, aes(x = Chromosome, y = Variant_Count, fill = Variant_Type)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_hline(yintercept = 150, color = "red", linetype = "dashed", linewidth = 1) +
  scale_fill_manual(values = variant_colors) +
  coord_cartesian(ylim = c(0, 150)) +
  scale_y_continuous(labels = comma, expand = c(0, 0), breaks = c(0, 50, 100, 150)) +
  labs(x = "Chromosome", y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        legend.position = "none",
        panel.grid = element_blank(),
        axis.line = element_line(color = "darkgray"))

plot_all_variants <- (plot_upper / plot_lower) +
  plot_layout(heights = c(9.5, 1), guides = "collect")

plot_all_variants

# Hg38: Figure S4C.svg
# CHM13: Figure S5C.svg

ggsave(file.path(save_path, "Figure S4C.svg"), 
       plot = plot_all_variants, width = 8, height = 5.8, units = "in", dpi = 300, bg = "transparent")


