library(ggplot2)
library(dplyr)
library(scales)


########################################## Ext 1A2A ##########################################

chromosome_order <- c(paste0("chr", 1:22), "chrX", "chrY")
variant_df$Chromosome <- factor(variant_df$Chromosome, levels = chromosome_order)

variant_colors <- c(SNV = "#0077BB", Indel = "#EE7733", SV = "#009988")

variant_df$Variant_Type <- factor(variant_df$Variant_Type, levels = c("SV", "Indel", "SNV"))

plot_all_variants <- ggplot(variant_df, aes(x = Chromosome, y = Variant_Count, fill = Variant_Type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = variant_colors) +
  scale_y_continuous(labels = comma) +  # Remove scientific notation
  labs(title = "Variants by Chromosome",
       x = "Chromosome", y = "Variant Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        panel.grid = element_blank())

########################################## Ext 1B2B ##########################################

mutation_type_order <- unique(melted_df$Mutation_Type)
melted_df$Mutation_Type <- factor(melted_df$Mutation_Type, levels = mutation_type_order)

ggplot(melted_df, aes(x = Mutation_Type, y = Count, fill = Tissue_type)) +
  geom_bar(stat = "identity", position = position_dodge(), fill = "#0077BB", 
           width = 0.7) +
  labs(title = "SNV Variant Counts by Mutation Type",
       x = "Mutation Type", y = "Count",
       fill = "Tissue Type") +
  scale_y_continuous(labels = scales::label_number(accuracy = 1e4)) +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.background = element_rect(fill = "transparent", color = NA),  
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_blank()

########################################## Ext 1C2C ##########################################





########################################## Ext 1D2D ##########################################






########################################## Ext 1E2E ##########################################

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

stacked_plot <- ggplot(aggregated_vaf, aes(x = VAF, y = Total_Count, fill = type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("SNV" = "#0077BB", "Indel" = "#EE7733")) +
  scale_y_continuous(labels = scales::label_number(accuracy = 1e5)) +
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
    panel.grid = element_blank()  # Remove grid lines
  )

########################################## Ext 1F2F ##########################################






########################################## Ext 1G2G ##########################################







########################################## Ext 1H2H ##########################################







########################################## Ext 1I2I ##########################################








########################################## Ext 1J2J ##########################################
