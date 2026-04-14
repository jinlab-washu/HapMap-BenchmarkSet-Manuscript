############################################ Fig4A ####################################################
# Load necessary libraries
install.packages("ggpubr", repos = "http://cran.us.r-project.org")
install.packages("jsonlite", repos = "http://cran.us.r-project.org")
install.packages("ggplot2", repos = "http://cran.us.r-project.org")

library(ggpubr)
library(jsonlite)
library(ggplot2)

# Function to read SNV & Indel summary.txt files
read_summary_file <- function(directory, caller, genome, vartype) {
  file_path <- file.path(directory, "summary.txt")
  
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  lines <- readLines(file_path)
  lines <- lines[-2]  # Remove problematic line
  
  column_names <- unlist(strsplit(lines[1], "\\s+"))
  
  data_lines <- lines[2:length(lines)]
  data_split <- strsplit(data_lines, "\\s+")
  
  data_cleaned <- lapply(data_split, function(x) {
    x_numeric <- as.numeric(x[-1])
    x_numeric[-1]
  })

  
  data <- do.call(rbind, data_cleaned)
  data <- as.data.frame(data)
  colnames(data) <- c("True_pos_baseline","True_pos_call","False_pos","False_neg","Precision","Sensitivity","Fmeasure")

  data$Caller <- caller
  data$Genome <- genome
  data$VarType <- vartype
  print(paste("Debugging data_cleaned for:", caller, genome, vartype))
  print(data)  
  return(data)
}

# Function to read SV summary.json files
read_summary_json <- function(directory, caller, genome, vartype) {
  file_path <- file.path(directory, "summary.json")
  
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  json_data <- fromJSON(file_path)
  # Extract required fields and store them in the correct format
  data <- data.frame(
    True_pos_baseline = ifelse(!is.null(json_data$`TP-base`), json_data$`TP-base`, NaN),
    True_pos_call = ifelse(!is.null(json_data$`TP-comp`), json_data$`TP-comp`, NaN),
    False_pos = ifelse(!is.null(json_data$`FP`), json_data$`FP`, NaN),
    False_neg = ifelse(!is.null(json_data$`FN`), json_data$`FN`,NaN),
    Precision = ifelse(!is.null(json_data$precision), as.numeric(json_data$precision), NaN),
    Sensitivity = ifelse(!is.null(json_data$recall), as.numeric(json_data$recall), NaN),  # recall = sensitivity
    Fmeasure = ifelse(!is.null(json_data$f1), as.numeric(json_data$f1), NaN),
    Caller = caller,
    Genome = genome,
    VarType = vartype
  )
  
  print(paste("Debugging data_cleaned for:", caller, genome, vartype))
  print(data)
  return(data)
}


# Define base paths separately for SNV, Indel, and SV
base_paths <- list(
  SNV = "{path to SNV comparison result dir}",
  Indel = "{path to indel comparison result dir}",
  SV = "{path to SV comparison result dir}"
)

# Function to construct file paths
construct_paths <- function(caller, genome, metric, vartype) {
  base <- base_paths[[vartype]]
  file.path(base, paste0(toupper(vartype), "_", caller, "_", genome, "_af_", 1:9, "_", metric))
}

# Define variant callers for each type
snvindel_callers <- c("mutect2", "strelka2", "varscan2")
sv_callers <- c("sniffles_tn", "severus")
genomes <- c("hg38", "chm13")
vartypes <- c("SNV", "Indel", "SV")

precision_data_list <- list()
sensitivity_data_list <- list()

# Loop through all variant types
for (vartype in vartypes) {
  base_path <- base_paths[[vartype]]
  
  construct_paths <- function(caller, genome, metric) {
    file.path(base_path, paste0(vartype, "_", caller, "_", genome, "_af_", 1:9, "_", metric))
  }
  
  if (vartype == "SV") {
    callers <- sv_callers  # Use the full vector
  } else {
    callers <- snvindel_callers
  }
  
  for (caller in callers) {
    for (genome in genomes) {
      if (vartype == "SV") {
        print(paste(caller, genomes, vartype))
        precision_paths <- construct_paths(caller, genome, "precision")  # No separate "precision"
        sensitivity_paths <- construct_paths(caller, genome, "sensitivity")
        precision_data <- lapply(precision_paths, read_summary_json, caller=caller, genome=genome, vartype=vartype)
        sensitivity_data <- lapply(sensitivity_paths, read_summary_json, caller=caller, genome=genome, vartype=vartype)
        precision_combined <- do.call(rbind, Map(cbind, precision_data, file = 1:9))
        sensitivity_combined <- do.call(rbind, Map(cbind, sensitivity_data, file = 1:9))  
      } else {
        precision_paths <- construct_paths(caller, genome, "precision")
        sensitivity_paths <- construct_paths(caller, genome, "sensitivity")
        
        precision_data <- lapply(precision_paths, read_summary_file, caller=caller, genome=genome, vartype=vartype)
        sensitivity_data <- lapply(sensitivity_paths, read_summary_file, caller=caller, genome=genome, vartype=vartype)
        
        precision_combined <- do.call(rbind, Map(cbind, precision_data, file = 1:9))
        sensitivity_combined <- do.call(rbind, Map(cbind, sensitivity_data, file = 1:9))
      }
      
      precision_data_list[[paste0(caller, "_", genome, "_", vartype)]] <- precision_combined
      sensitivity_data_list[[paste0(caller, "_", genome, "_", vartype)]] <- sensitivity_combined
    }
  }
}

# Combine all data
precision_combined_data <- do.call(rbind, precision_data_list)
sensitivity_combined_data <- do.call(rbind, sensitivity_data_list)

# Define VAF labels
vaf_labels <- c("0-0.005", "0.005-0.01", "0.01-0.015", "0.015-0.02", "0.02-0.03",
                "0.03-0.04", "0.04-0.05", "0.05-0.1", "0.1-0.165")

# Assign VAF labels
precision_combined_data$file <- factor(precision_combined_data$file, labels=vaf_labels)
sensitivity_combined_data$file <- factor(sensitivity_combined_data$file, labels=vaf_labels)

# Define VarType ordering for facet_grid
precision_combined_data$VarType <- factor(precision_combined_data$VarType, levels = c("SNV", "Indel", "SV"))
sensitivity_combined_data$VarType <- factor(sensitivity_combined_data$VarType, levels = c("SNV", "Indel", "SV"))

# Precision Plot
precision_plot <- ggplot(precision_combined_data, 
                         aes(x=file, y=Precision, color=Caller, linetype=Genome, group=interaction(Caller, Genome))) +
  geom_point(size=3) +
  geom_line(linewidth=1) +
  labs(title="Precision Comparison (SNV, Indel & SV)", x="VAF", y="Precision") +
  ylim(0, 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  facet_grid(. ~ VarType)+
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_blank()
  )

# Sensitivity Plot
sensitivity_plot <- ggplot(sensitivity_combined_data, 
                           aes(x = file, y = Sensitivity, color = Caller, 
                               group = interaction(Caller, Genome), linetype=Genome)) +
  geom_point(size=3) +
  geom_line(linewidth=1) +
  labs(title="Sensitivity Comparison (SNV, Indel & SV)", x="VAF", y="Sensitivity") +
  ylim(0, 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  facet_grid(. ~ VarType)+
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_blank()
  )

# Combine plots
final_combined_plot <- ggarrange( sensitivity_plot, precision_plot, 
                                 ncol = 1, nrow = 2,
                                 common.legend = TRUE, legend="bottom")

# Display the plot
print(final_combined_plot)
############################################ Fig4B ####################################################
library(jsonlite)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggpubr)


# Set working directory
setwd("{path to MEI caller comparison dir}")  # <-- update this

# Get only summary.json files in directories that match MEI_* pattern
folders <- list.dirs(path = ".", recursive = FALSE)
target_dirs <- grep("MEI_.*_(ALU|L1|SVA)_(chm13|hg38)$", folders, value = TRUE)
files <- file.path(target_dirs, "summary.json")
files <- files[file.exists(files)]

# Load and parse JSONs
data <- lapply(files, function(file) {
  json <- fromJSON(file)
  label <- basename(dirname(file))
  parts <- str_split(label, "_", simplify = TRUE)
  data.frame(
    tool = parts[2],
    element = parts[3],
    genome = parts[4],
    precision = json$precision,
    recall = json$recall,
    f1 = json$f1,
    label = label
  )
}) %>% bind_rows()

# Plot Precision
precision <- ggplot(data, aes(x = tool, y = precision, fill = genome)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  facet_wrap(~element, scales = "free_x") +
  theme_minimal() +
  labs(
       y = "Precision", x = "Tool") +
  scale_fill_manual(values = c("chm13" = "steelblue", "hg38" = "tomato")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot Recall
recall <- ggplot(data, aes(x = tool, y = recall, fill = genome)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  facet_wrap(~element, scales = "free_x") +
  theme_minimal() +
  labs(
       y = "Recall", x = "Tool") +
  scale_fill_manual(values = c("chm13" = "steelblue", "hg38" = "tomato")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# Combine plots
final_combined_plot <- ggarrange(recall, precision,
                                 ncol = 1, nrow = 2,
                                 common.legend = TRUE, legend = "bottom")

# Show plot
print(final_combined_plot)

############################################ Fig4C ####################################################
# Load necessary libraries
install.packages("ggpubr", repos = "http://cran.us.r-project.org")
install.packages("jsonlite", repos = "http://cran.us.r-project.org")
install.packages("ggplot2", repos = "http://cran.us.r-project.org")
install.packages("svglite", repos = "http://cran.us.r-project.org")
library(svglite)
library(ggpubr)
library(jsonlite)
library(ggplot2)


read_summary_file <- function(directory, caller, genome, vartype) {
  file_path <- file.path(directory, "summary.txt")
  
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  lines <- readLines(file_path)
  lines <- lines[-2]  # Remove problematic line
  
  column_names <- unlist(strsplit(lines[1], "\\s+"))
  
  data_lines <- lines[2:length(lines)]
  data_split <- strsplit(data_lines, "\\s+")
  
  data_cleaned <- lapply(data_split, function(x) {
    x_numeric <- as.numeric(x[-1])
    x_numeric[-1]
  })
  
  
  data <- do.call(rbind, data_cleaned)
  data <- as.data.frame(data)
  colnames(data) <- c("True_pos_baseline","True_pos_call","False_pos","False_neg","Precision","Sensitivity","Fmeasure")
  
  data$Caller <- caller
  data$Genome <- genome
  data$VarType <- vartype
  print(paste("Debugging data_cleaned for:", caller, genome, vartype))
  print(data)  
  return(data)
}

# Function to read SV summary.json files
read_summary_json <- function(directory, caller, genome, vartype) {
  file_path <- file.path(directory, "summary.json")
  
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  json_data <- fromJSON(file_path)
  # Extract required fields and store them in the correct format
  data <- data.frame(
    True_pos_baseline = ifelse(!is.null(json_data$`TP-base`), json_data$`TP-base`, NaN),
    True_pos_call = ifelse(!is.null(json_data$`TP-comp`), json_data$`TP-comp`, NaN),
    False_pos = ifelse(!is.null(json_data$`FP`), json_data$`FP`, NaN),
    False_neg = ifelse(!is.null(json_data$`FN`), json_data$`FN`,NaN),
    Precision = ifelse(!is.null(json_data$precision), as.numeric(json_data$precision), NaN),
    Sensitivity = ifelse(!is.null(json_data$recall), as.numeric(json_data$recall), NaN),  # recall = sensitivity
    Fmeasure = ifelse(!is.null(json_data$f1), as.numeric(json_data$f1), NaN),
    Caller = caller,
    Genome = genome,
    VarType = vartype
  )
  
  print(paste("Debugging data_cleaned for:", caller, genome, vartype))
  print(data)
  return(data)
}


# Set base paths for SNV/Indel/SV
base_paths <- list(
  SNV = "{path to snv caller comparison dir}",
  Indel = "{path to indel caller comparison dir}",
  SV = "{path to sv caller comparison dir}"
)

# Define callers and genomes
snvindel_callers <- c("mutect2", "strelka2", "varscan2")
sv_callers <- c("sniffles", "severus")
genomes <- c("hg38", "chm13")
vartypes <- c("SNV","Indel","SV")
#vartypes <- c("SV")
kmer_labels <- c("kmer_k24.umap", "kmer_k36.umap", "kmer_k50.umap", "kmer_k100.umap", "kmer_all")

# Function to build directory paths
construct_kmer_dirs <- function(base_path, caller, genome, vartype) {
  suffixes <- c("kmer_k24.umap", "kmer_k36.umap", "kmer_k50.umap", "kmer_k100.umap", "kmer_all")
  file.path(base_path, paste0(vartype, "_", caller, "_", genome, "_", suffixes))
}

# Collect data
precision_data_list <- list()
sensitivity_data_list <- list()

for (vartype in vartypes) {
  base_path <- base_paths[[vartype]]
  callers <- if (vartype == "SV") sv_callers else snvindel_callers
  
  for (caller in callers) {
    for (genome in genomes) {
      dirs <- construct_kmer_dirs(base_path, caller, genome, vartype)
      
      if (vartype == "SV") {
        precision_data <- lapply(dirs, read_summary_json, caller = caller, genome = genome, vartype = vartype)
        sensitivity_data <- lapply(dirs, read_summary_json, caller = caller, genome = genome, vartype = vartype)
      } else {
        precision_data <- lapply(dirs, read_summary_file, caller = caller, genome = genome, vartype = vartype)
        sensitivity_data <- lapply(dirs, read_summary_file, caller = caller, genome = genome, vartype = vartype)
      }
      
      precision_combined <- do.call(rbind, Map(cbind, precision_data, file = 1:5))
      sensitivity_combined <- do.call(rbind, Map(cbind, sensitivity_data, file = 1:5))
      
      precision_data_list[[paste0(caller, "_", genome, "_", vartype)]] <- precision_combined
      sensitivity_data_list[[paste0(caller, "_", genome, "_", vartype)]] <- sensitivity_combined
    }
  }
}

# Combine and label k-mer bins
precision_combined_data <- do.call(rbind, precision_data_list)
sensitivity_combined_data <- do.call(rbind, sensitivity_data_list)

precision_combined_data$file <- factor(precision_combined_data$file, labels = kmer_labels)
sensitivity_combined_data$file <- factor(sensitivity_combined_data$file, labels = kmer_labels)

# Set order for facet plots
precision_combined_data$VarType <- factor(precision_combined_data$VarType, levels = c("SNV", "Indel", "SV"))
sensitivity_combined_data$VarType <- factor(sensitivity_combined_data$VarType, levels = c("SNV", "Indel", "SV"))
# Update labels for plotting
precision_combined_data$file <- factor(precision_combined_data$file,
                                       levels = c("kmer_k24.umap", "kmer_k36.umap", 
                                                  "kmer_k50.umap", "kmer_k100.umap", "kmer_all"),
                                       labels = c("k24", "k36", "k50", "k100", "remaining"))

sensitivity_combined_data$file <- factor(sensitivity_combined_data$file,
                                         levels = c("kmer_k24.umap", "kmer_k36.umap", 
                                                    "kmer_k50.umap", "kmer_k100.umap", "kmer_all"),
                                         labels = c("k24", "k36", "k50", "k100", "remaining"))


# Plot Precision
precision_plot <- ggplot(precision_combined_data,
                         aes(x = file, y = Precision, color = Caller, linetype = Genome, group = interaction(Caller, Genome))) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  labs(title = "", x = "k-mer Bin", y = "Precision") +
  ylim(0, 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(. ~ VarType) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.background = element_blank())

# Plot Sensitivity
sensitivity_plot <- ggplot(sensitivity_combined_data,
                           aes(x = file, y = Sensitivity, color = Caller, linetype = Genome, group = interaction(Caller, Genome))) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  labs(title = "Sensitivity and Precision Comparison by mappabiltiy (SNV, Indel & SV)", x = "", y = "Sensitivity") +
  ylim(0, 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(. ~ VarType) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.background = element_blank())

# Combine plots
final_combined_plot <- ggarrange(sensitivity_plot, precision_plot,
                                 ncol = 1, nrow = 2,
                                 common.legend = TRUE, legend = "bottom")

# Show plot
print(final_combined_plot)


############################################ Fig4D ####################################################
# Load the data
data <- read.table("summary_umap_composition.txt", header = TRUE, sep = "\t")

# Calculate exclusive proportions for each kmer
data <- data %>%
  mutate(entire_genome_excl = c(entire_genome_bp[1], diff(entire_genome_bp)),
         omim_gene_excl = c(omim_gene_bp[1], diff(omim_gene_bp)),
         truthset_snv_excl = c(truthset_snv_bp[1], diff(truthset_snv_bp))) %>%
  select(kmer, entire_genome_excl, omim_gene_excl, truthset_snv_excl)

# Convert the data to long format for plotting
data_long <- data %>%
  pivot_longer(cols = starts_with("entire_genome_excl"):starts_with("truthset_snv_excl"),
               names_to = "category",
               values_to = "bp") %>%
  mutate(category = gsub("_excl", "", category)) %>%
  mutate(category = factor(category, levels = c("entire_genome", "omim_gene", "truthset_snv"),
                           labels = c("Entire genome\n(n=3.1G)", "Omim gene \n(n=4,978)", "Truthset snv\n(n=5.6M)")))

# Define the order for kmer
data_long <- data_long %>%
  mutate(kmer = factor(kmer, 
                       levels = c("hg38_all_chr", "k100.umap", "k50.umap", "k36.umap", "k24.umap"),
                       labels = c("Remaining", "k100 unique", "k50 unique", "k36 unique", "k24")))

# Calculate percentages
data_long <- data_long %>%
  group_by(category) %>%
  mutate(percent = bp / sum(bp) * 100)

# Plot the bar chart
plot<- ggplot(data_long, aes(x = category, y = percent, fill = kmer)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Composition in mappability",
       x = "Base pairs",
       y = "Percentage",
       fill = "UMAP mappability") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # X-axis text
    axis.text.y = element_text(size = 14),  # Y-axis text
    axis.title.x = element_text(size = 16),  # X-axis title
    axis.title.y = element_text(size = 16),  # Y-axis title
    legend.title = element_text(size = 14),  # Legend title
    legend.text = element_text(size = 12),  # Legend labels
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)  # Plot title
  )
# Show plot
print(plot)

