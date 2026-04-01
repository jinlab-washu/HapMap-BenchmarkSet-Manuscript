
library(ggplot2)
library(cowplot)
library(dplyr)
library(svglite)

temp_csv_forplot = "zt_temp/PrepPlots_20250429"
save_path = "V1.6_20250429"

# hg38_SNVcountsForPlot.csv
# chm13_SNVcountsForPlot.csv
snvs_df <- read.csv(file.path(temp_csv_forplot, "chm13_SNVcountsForPlot.csv"))
# hg38_IndelcountsForPlot.csv
# chm13_IndelcountsForPlot.csv
indels_df <- read.csv(file.path(temp_csv_forplot, "chm13_IndelcountsForPlot.csv"))

#######################################################################

# Combine all variants into one dataframe
all_variants_df <- bind_rows(
  snvs_df %>% mutate(type = "SNV"),
  indels_df %>% mutate(type = "Indel")
)

# Summarize total counts per VAF and type
aggregated_vaf <- all_variants_df %>%
  group_by(VAF, type) %>%
  summarize(Total_Count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  mutate(VAF = ifelse(is.na(VAF), "NA", as.character(VAF)))

aggregated_vaf <- aggregated_vaf %>%
  arrange(as.numeric(ifelse(VAF == "NA", NA, VAF))) %>%  # Sort numerically, NAs last
  mutate(VAF = factor(VAF, levels = unique(VAF)))

# Stacked bar plot
stacked_plot <- ggplot(aggregated_vaf, aes(x = VAF, y = Total_Count, fill = type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("SNV" = "#0077BB", "Indel" = "#EE7733")) +
  scale_y_continuous(labels = scales::label_number(accuracy = 1e5),
                     expand = c(0,0)) +
  theme_minimal() +
  labs(
    title = "SNV/Indel VAF Distribution",
    x = "VAF",
    y = "Count",
    fill = "Variant Type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    panel.grid = element_blank(),
    axis.line = element_line(color = "darkgray")
  )

# Show the plot
stacked_plot

ggsave(file.path(save_path, "Figure S5A.svg"), plot = stacked_plot, width = 8, height = 5.8, units = "in", dpi = 300, bg = "transparent")


