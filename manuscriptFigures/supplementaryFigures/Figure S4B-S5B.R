library(ggplot2)
library(svglite)

temp_csv_forplot = "zt_temp/PrepPlots_20250429"
save_path = "V1.6_20250429"

# hg38_mutation_counts.csv
# chm13_mutation_counts.csv
melted_df <- read.csv(file.path(temp_csv_forplot, "chm13_mutation_counts.csv"))

#######################################################################

mutation_type_order <- unique(melted_df$Mutation_Type)
melted_df$Mutation_Type <- factor(melted_df$Mutation_Type, levels = mutation_type_order)

ggplot(melted_df, aes(x = Mutation_Type, y = Count, fill = Tissue_type)) +
  geom_bar(stat = "identity", position = position_dodge(), fill = "#0077BB", 
           width = 0.7) +
  labs(title = "SNV Variant Counts by Mutation Type",
       x = "Mutation Type", y = "Count",
       fill = "Tissue Type") +
  scale_y_continuous(labels = scales::label_number(accuracy = 1e4),
                     expand = c(0,0)) +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.background = element_rect(fill = "transparent", color = NA),  
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_blank(),
        axis.line = element_line(color = "darkgray")
  )


ggsave(file.path(save_path, "Figure S5B.svg"), plot = last_plot(), width = 8, height = 5.8, units = "in", dpi = 300, bg = "transparent")

