library(ggplot2)
library(dplyr)
library(scales)
library(stringr)
install.packages("patchwork")
library(patchwork)
library(VariantAnnotation)

########################################## Fig 4A5A ##########################################

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

########################################## Fig 4B5B ##########################################

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

########################################## Fig 4C5C ##########################################



########################################## Fig 4D5D ##########################################

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
########################################## Fig 4E5E ##########################################

# Specify the VCF file path
vcf_file <- "{path to sv vcf}"

# Read the VCF file
vcf <- readVcf(vcf_file)

# Extract REF, ALT, and FORMAT columns
vcf_data <- data.frame(
  REF = as.character(ref(vcf)),
  ALT = unlist(CharacterList(alt(vcf))),
  AF = as.numeric(geno(vcf)$AF),
  Length_Diff = nchar(unlist(CharacterList(alt(vcf)))) - nchar(as.character(ref(vcf)))
)

# Categorize Length_Diff into INS and DEL, keep signed size
vcf_data <- vcf_data %>%
  mutate(
    Type = case_when(
      Length_Diff > 0 ~ "INS",
      Length_Diff < 0 ~ "DEL",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Type)) %>%
  filter(Length_Diff <= -50 | Length_Diff >= 50)

# Equal-width signed bins: 100 bp
bin_width <- 100
max_size <- 10000

# Breaks from -10 kb to +10 kb
breaks <- seq(-max_size, max_size, by = bin_width)

# Bin data numerically, no regex parsing
vcf_data <- vcf_data %>%
  mutate(
    bin_id = cut(
      Length_Diff,
      breaks = breaks,
      include.lowest = TRUE,
      right = FALSE,
      labels = FALSE
    )
  ) %>%
  filter(!is.na(bin_id)) %>%
  mutate(
    Bin_start = breaks[bin_id],
    Bin_end = breaks[bin_id + 1],
    Midpoint = (Bin_start + Bin_end) / 2,
    Bin = paste0("[", Bin_start, ",", Bin_end, ")")
  )

# Count per bin and type
bin_summary <- vcf_data %>%
  group_by(Type, Bin, Midpoint) %>%
  summarise(Count = n(), .groups = "drop")

# Plot line plot: deletions left, insertions right
ggplot(bin_summary, aes(x = Midpoint, y = Count, color = Type, group = Type)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 1.0) +
  scale_x_continuous(
    breaks = c(-10000, -8000, -6000, -4000, -2000, -1000, -500,
               500, 1000, 2000, 4000, 6000, 8000, 10000),
    labels = c("10kb", "8kb", "6kb", "4kb", "2kb", "1kb", "500bp",
               "500bp", "1kb", "2kb", "4kb", "6kb", "8kb", "10kb")
  ) +
  scale_y_log10() +
  scale_color_manual(values = c("INS" = "#44BB99", "DEL" = "#EE99AA")) +
  labs(
    title = "SV size distribution",
    x = "Deletion size                                Insertion size",
    y = "Count",
    color = "SV Type"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")
        
########################################## Fig 4F5F ##########################################
read_vcf_extract_repeats <- function(vcf_path){
  vcf <- read.table(vcf_path, comment.char = "#", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(vcf) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","Sample")
  
  # Extract MEI subtype: second value inside MEI=(class,subtype)
  vcf$repeat_subtype <- str_extract(vcf$INFO, "MEI=\\([^,]*,([^\\)]+)\\)") %>%
    str_replace("MEI=\\([^,]*,", "") %>%
    str_replace("\\)", "")
  
  return(vcf)
}

# Donut plot function
create_donut_chart <- function(data, class_label, threshold = 0) {
  repeat_counts <- data %>%
    group_by(repeat_subtype) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(percentage = round((count / sum(count)) * 100, 1))
  
  other_label <- paste0("other_", class_label)
  repeat_counts <- repeat_counts %>%
    mutate(repeat_grouped = ifelse(percentage < threshold, other_label, repeat_subtype)) %>%
    group_by(repeat_grouped) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    mutate(percentage = round((count / sum(count)) * 100, 1))
  
  rounding_error <- 100 - sum(repeat_counts$percentage)
  repeat_counts <- repeat_counts %>%
    mutate(percentage = ifelse(repeat_grouped == other_label, percentage + rounding_error, percentage))
  
  repeat_counts$repeat_grouped <- factor(repeat_counts$repeat_grouped, levels = unique(repeat_counts$repeat_grouped))
  
  p <- ggplot(repeat_counts, aes(x = "", y = count, fill = repeat_grouped)) +
    geom_bar(stat = "identity", width = 1, color = "black") +
    coord_polar("y", start = 0) +
    theme_void() +
    theme(
      legend.title = element_blank(), 
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(title = paste(class_label, "Insertions Classification")) +
    geom_text(aes(label = ifelse(percentage > threshold, paste0(repeat_grouped, " ", percentage, "%"), "")), 
              position = position_stack(vjust = 0.5), size = 5, color = "black") +
    annotate("text", x = 0, y = 0, 
             label = paste(sum(repeat_counts$count), class_label, "insertions\nare classified"), 
             size = 6, fontface = "bold")
  
  return(p)
}

# Load data
l1_vcf_data <- read_vcf_extract_repeats("{path to l1 vcf}")
alu_vcf_data <- read_vcf_extract_repeats("{path to alu vcf}")
sva_vcf_data <- read_vcf_extract_repeats("{path to sva vcf}")

# Filter and plot
l1_data <- l1_vcf_data %>% filter(str_detect(INFO, "LINE/L1"))
alu_data <- alu_vcf_data %>% filter(str_detect(INFO, "SINE/Alu"))
sva_data <- sva_vcf_data %>% filter(str_detect(INFO, "Retroposon/SVA,SVA_F"))

l1_plot <- create_donut_chart(l1_data, "L1")
alu_plot <- create_donut_chart(alu_data, "Alu")
sva_plot <- create_donut_chart(alu_data, "SVA")

# Extract shared legend
final_plot <- patchwork::wrap_plots(alu_plot, l1_plot, sva_plot,guides = "collect") +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Final plot with shared legend
print(final_plot)


########################################## Fig 4G5G ##########################################
#install.packages("UpSetR")
library(UpSetR)
#reference: https://www.youtube.com/watch?v=n9MRCZxJOfk&t=333s
data = read.delim("SMHTHAPMAP6_GRCh38_v1.2_somatic_benchmark_snvs_filtered.txt", header = TRUE, sep = "\t")
dt = data
        
#make list
x = list(
  HG002_M = dt[grep("HG002_M", dt[,2]),1],
  HG002_P = dt[grep("HG002_P", dt[,2]),1],
  HG00438_M = dt[grep("HG00438_M", dt[,2]),1],
  HG00438_P = dt[grep("HG00438_P", dt[,2]),1],
  HG005_M = dt[grep("HG005_M", dt[,2]),1],
  HG005_P = dt[grep("HG005_P", dt[,2]),1],
  HG02257_M = dt[grep("HG02257_M", dt[,2]),1],
  HG02257_P = dt[grep("HG02257_P", dt[,2]),1],
  HG02486_M = dt[grep("HG02486_M", dt[,2]),1],
  HG02486_P = dt[grep("HG02486_P", dt[,2]),1],
  HG02622_M = dt[grep("HG02622_M", dt[,2]),1],
  HG02622_P = dt[grep("HG02622_P", dt[,2]),1]
)
        
# Merge _M and _P values for each sample
y = list(
  HG002 = c(x$HG002_M, x$HG002_P),
  HG00438 = c(x$HG00438_M, x$HG00438_P),
  HG005 = c(x$HG005_M, x$HG005_P),
  HG02257 = c(x$HG02257_M, x$HG02257_P),
  HG02486 = c(x$HG02486_M, x$HG02486_P),
  HG02622 = c(x$HG02622_M, x$HG02622_P)
)
        
#draw plot
upset(fromList(x), order.by = "freq",
      mainbar.y.label = "Counts", sets.x.label = "Sample",
      sets= c("HG002_M", "HG002_P", "HG00438_M", "HG00438_P", "HG005_M", "HG005_P", "HG02257_M", "HG02257_P", "HG02486_M", "HG02486_P", "HG02622_M", "HG02622_P"))
upset(fromList(y), order.by = "freq",
      mainbar.y.label = "Counts", sets.x.label = "Sample",
      sets= c("HG002", "HG00438", "HG005", "HG02257", "HG02486", "HG02622"))

#save figure
svg("SMHTHAPMAP6_GRCh38_v1.2_somatic_benchmark_snvs_filtered.svg", width=16, height=6)  # Open a new SVG device
upset(fromList(x), order.by = "freq",
      mainbar.y.label = "Counts", sets.x.label = "Sample",
      sets= c("HG002_M", "HG002_P", "HG00438_M", "HG00438_P", "HG005_M", "HG005_P", "HG02257_M", "HG02257_P", "HG02486_M", "HG02486_P", "HG02622_M", "HG02622_P"))
dev.off()

########################################## Fig 4H5H ##########################################
library(UpSetR)
data = read.delim("SMHTHAPMAP6_GRCh38_v1.2_somatic_benchmark_indels_filtered.txt", header = TRUE, sep = "\t")
dt = data
        
#make list
x = list(
  HG002_M = dt[grep("HG002_M", dt[,2]),1],
  HG002_P = dt[grep("HG002_P", dt[,2]),1],
  HG00438_M = dt[grep("HG00438_M", dt[,2]),1],
  HG00438_P = dt[grep("HG00438_P", dt[,2]),1],
  HG005_M = dt[grep("HG005_M", dt[,2]),1],
  HG005_P = dt[grep("HG005_P", dt[,2]),1],
  HG02257_M = dt[grep("HG02257_M", dt[,2]),1],
  HG02257_P = dt[grep("HG02257_P", dt[,2]),1],
  HG02486_M = dt[grep("HG02486_M", dt[,2]),1],
  HG02486_P = dt[grep("HG02486_P", dt[,2]),1],
  HG02622_M = dt[grep("HG02622_M", dt[,2]),1],
  HG02622_P = dt[grep("HG02622_P", dt[,2]),1]
)
        
# Merge _M and _P values for each sample
y = list(
  HG002 = c(x$HG002_M, x$HG002_P),
  HG00438 = c(x$HG00438_M, x$HG00438_P),
  HG005 = c(x$HG005_M, x$HG005_P),
  HG02257 = c(x$HG02257_M, x$HG02257_P),
  HG02486 = c(x$HG02486_M, x$HG02486_P),
  HG02622 = c(x$HG02622_M, x$HG02622_P)
)
        
#draw plot
upset(fromList(x), order.by = "freq",
      mainbar.y.label = "Counts", sets.x.label = "Sample",
      sets= c("HG002_M", "HG002_P", "HG00438_M", "HG00438_P", "HG005_M", "HG005_P", "HG02257_M", "HG02257_P", "HG02486_M", "HG02486_P", "HG02622_M", "HG02622_P"))
upset(fromList(y), order.by = "freq",
      mainbar.y.label = "Counts", sets.x.label = "Sample",
      sets= c("HG002", "HG00438", "HG005", "HG02257", "HG02486", "HG02622"))

#save figure
svg("SMHTHAPMAP6_GRCh38_v1.2_somatic_benchmark_indels_filtered.svg", width=16, height=6)  # Open a new SVG device
upset(fromList(x), order.by = "freq",
      mainbar.y.label = "Counts", sets.x.label = "Sample",
      sets= c("HG002_M", "HG002_P", "HG00438_M", "HG00438_P", "HG005_M", "HG005_P", "HG02257_M", "HG02257_P", "HG02486_M", "HG02486_P", "HG02622_M", "HG02622_P"))
dev.off()

########################################## Fig 4I5I ##########################################
library(UpSetR)
data = read.delim("SMHTHAPMAP6_GRCh38_v1.2_somatic_benchmark_svs_filtered.txt", header = TRUE, sep = "\t")
dt = data
        
#make list
x = list(
  HG002_M = dt[grep("HG002_M", dt[,2]),1],
  HG002_P = dt[grep("HG002_P", dt[,2]),1],
  HG00438_M = dt[grep("HG00438_M", dt[,2]),1],
  HG00438_P = dt[grep("HG00438_P", dt[,2]),1],
  HG005_M = dt[grep("HG005_M", dt[,2]),1],
  HG005_P = dt[grep("HG005_P", dt[,2]),1],
  HG02257_M = dt[grep("HG02257_M", dt[,2]),1],
  HG02257_P = dt[grep("HG02257_P", dt[,2]),1],
  HG02486_M = dt[grep("HG02486_M", dt[,2]),1],
  HG02486_P = dt[grep("HG02486_P", dt[,2]),1],
  HG02622_M = dt[grep("HG02622_M", dt[,2]),1],
  HG02622_P = dt[grep("HG02622_P", dt[,2]),1]
)
        
# Merge _M and _P values for each sample
y = list(
  HG002 = c(x$HG002_M, x$HG002_P),
  HG00438 = c(x$HG00438_M, x$HG00438_P),
  HG005 = c(x$HG005_M, x$HG005_P),
  HG02257 = c(x$HG02257_M, x$HG02257_P),
  HG02486 = c(x$HG02486_M, x$HG02486_P),
  HG02622 = c(x$HG02622_M, x$HG02622_P)
)
        
#draw plot
upset(fromList(x), order.by = "freq",
      mainbar.y.label = "Counts", sets.x.label = "Sample",
      sets= c("HG002_M", "HG002_P", "HG00438_M", "HG00438_P", "HG005_M", "HG005_P", "HG02257_M", "HG02257_P", "HG02486_M", "HG02486_P", "HG02622_M", "HG02622_P"))
upset(fromList(y), order.by = "freq",
      mainbar.y.label = "Counts", sets.x.label = "Sample",
      sets= c("HG002", "HG00438", "HG005", "HG02257", "HG02486", "HG02622"))

#save figure
svg("SMHTHAPMAP6_GRCh38_v1.2_somatic_benchmark_svs_filtered.svg", width=16, height=6)  # Open a new SVG device
upset(fromList(x), order.by = "freq",
      mainbar.y.label = "Counts", sets.x.label = "Sample",
      sets= c("HG002_M", "HG002_P", "HG00438_M", "HG00438_P", "HG005_M", "HG005_P", "HG02257_M", "HG02257_P", "HG02486_M", "HG02486_P", "HG02622_M", "HG02622_P"))
dev.off()

########################################## Fig 4J5J ##########################################
library(UpSetR)
data = read.delim("SMHTHAPMAP6_GRCh38_v1.2_somatic_benchmark_svs_mei_filtered.txt", header = TRUE, sep = "\t")
dt = data
        
#make list
x = list(
  HG002_M = dt[grep("HG002_M", dt[,2]),1],
  HG002_P = dt[grep("HG002_P", dt[,2]),1],
  HG00438_M = dt[grep("HG00438_M", dt[,2]),1],
  HG00438_P = dt[grep("HG00438_P", dt[,2]),1],
  HG005_M = dt[grep("HG005_M", dt[,2]),1],
  HG005_P = dt[grep("HG005_P", dt[,2]),1],
  HG02257_M = dt[grep("HG02257_M", dt[,2]),1],
  HG02257_P = dt[grep("HG02257_P", dt[,2]),1],
  HG02486_M = dt[grep("HG02486_M", dt[,2]),1],
  HG02486_P = dt[grep("HG02486_P", dt[,2]),1],
  HG02622_M = dt[grep("HG02622_M", dt[,2]),1],
  HG02622_P = dt[grep("HG02622_P", dt[,2]),1]
)
        
# Merge _M and _P values for each sample
y = list(
  HG002 = c(x$HG002_M, x$HG002_P),
  HG00438 = c(x$HG00438_M, x$HG00438_P),
  HG005 = c(x$HG005_M, x$HG005_P),
  HG02257 = c(x$HG02257_M, x$HG02257_P),
  HG02486 = c(x$HG02486_M, x$HG02486_P),
  HG02622 = c(x$HG02622_M, x$HG02622_P)
)
        
#draw plot
upset(fromList(x), order.by = "freq",
      mainbar.y.label = "Counts", sets.x.label = "Sample",
      sets= c("HG002_M", "HG002_P", "HG00438_M", "HG00438_P", "HG005_M", "HG005_P", "HG02257_M", "HG02257_P", "HG02486_M", "HG02486_P", "HG02622_M", "HG02622_P"))
upset(fromList(y), order.by = "freq",
      mainbar.y.label = "Counts", sets.x.label = "Sample",
      sets= c("HG002", "HG00438", "HG005", "HG02257", "HG02486", "HG02622"))

#save figure
svg("SMHTHAPMAP6_GRCh38_v1.2_somatic_benchmark_svs_mei_filtered.svg", width=16, height=6)  # Open a new SVG device
upset(fromList(x), order.by = "freq",
      mainbar.y.label = "Counts", sets.x.label = "Sample",
      sets= c("HG002_M", "HG002_P", "HG00438_M", "HG00438_P", "HG005_M", "HG005_P", "HG02257_M", "HG02257_P", "HG02486_M", "HG02486_P", "HG02622_M", "HG02622_P"))
dev.off()
