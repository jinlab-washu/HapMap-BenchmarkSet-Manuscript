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
  
  data_lines <- lines[length(lines):length(lines)]
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
  SNV = "/Volumes/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/snv",
  Indel = "/Volumes/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/indels",
  SV = "/Volumes/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/sv"
)

# Function to construct file paths
construct_paths <- function(caller, genome, metric, vartype) {
  base <- base_paths[[vartype]]
  file.path(base, paste0(toupper(vartype), "_", caller, "_", genome, "_af_", 1:9, "_", metric))
}

# Define variant callers for each type
snvindel_callers <- c("mutect2", "strelka2", "varscan2", "neusomatic", "deepsomatic")
sv_callers <- c("sniffles", "severus", "pbsv", "svision", "savana")
genomes <- c("GRCh38", "CHM13")
vartypes <- c("SNV","Indel", "SV")

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

snv_indel_colors <- c(
  "clairS"      = "#E41A1C",
  "deepsomatic" = "#FF7F00",
  "mutect2"     = "#984EA3",
  "neusomatic"  = "#4DAF4A",
  "strelka2"    = "#377EB8",
  "varscan2"    = "#F781BF"
)

sv_colors <- c(
  "pbsv"     = "#1B9E77",
  "savana"   = "#D95F02",
  "severus"  = "#7570B3",
  "sniffles" = "#E7298A",
  "svision"  = "yellow"
)
color_pal <- c(snv_indel_colors, sv_colors)
# Precision Plot
precision_plot <- ggplot(precision_combined_data, 
                         aes(x=file, y=Precision, color=Caller, linetype=Genome, group=interaction(Caller, Genome))) +
  geom_point(size=3) +
  geom_line(linewidth=1) +
  labs(title="Precision Comparison (SNV, Indel & SV)", x="VAF", y="Precision") +
  ylim(0, 1) +
  theme_minimal() +
  scale_color_manual(values = color_pal) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  facet_grid(. ~ VarType)+
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.background = element_blank()
  )

# Sensitivity Plot
sensitivity_plot <- ggplot(sensitivity_combined_data, 
                           aes(x = file, y = Sensitivity, color = Caller, 
                               group = interaction(Caller, Genome), linetype=Genome)) +
  geom_point(size=3) +
  geom_line(linewidth=1) +
  scale_color_manual(values = color_pal) +
  labs(title="Recall Comparison (SNV, Indel & SV)", x="VAF", y="Recall") +
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

