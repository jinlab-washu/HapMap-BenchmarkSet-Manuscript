# Load necessary libraries
install.packages("ggpubr", repos = "http://cran.us.r-project.org")
install.packages("jsonlite", repos = "http://cran.us.r-project.org")
install.packages("ggplot2", repos = "http://cran.us.r-project.org")
install.packages("svglite", repos = "http://cran.us.r-project.org")
install.packages("dplyr", repos = "http://cran.us.r-project.org")
install.packages("tibble", repos = "http://cran.us.r-project.org")
install.packages("purrr", repos = "http://cran.us.r-project.org")

library(svglite)
library(ggpubr)
library(jsonlite)
library(ggplot2)
library(dplyr)
library(tibble)
library(purrr)

read_summary_file <- function(directory, caller, genome, vartype) {
  file_path <- file.path(directory, "summary.txt")
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  lines <- readLines(file_path)
  if (length(lines) < 3) {
    warning(paste("summary.txt malformed:", file_path))
    return(NULL)
  }
  # Remove problematic line (as in your original)
  if (length(lines) >= 2) lines <- lines[-2]
  
  column_names <- unlist(strsplit(lines[1], "\\s+"))
  data_lines <- lines[2:length(lines)]
  data_split <- strsplit(data_lines, "\\s+")
  
  data_cleaned <- lapply(data_split, function(x) {
    x_numeric <- suppressWarnings(as.numeric(x[-1]))
    x_numeric[-1]
  })
  
  data <- do.call(rbind, data_cleaned)
  data <- as.data.frame(data)
  colnames(data) <- c("True_pos_baseline","True_pos_call","False_pos","False_neg","Precision","Sensitivity","Fmeasure")
  
  data$Caller <- caller
  data$Genome <- genome
  data$VarType <- vartype
  
  message(paste("Parsed:", file_path, "->", caller, genome, vartype))
  return(data)
}

# Function to read SV summary.json files
read_summary_json <- function(directory, caller, genome, vartype) {
  file_path <- file.path(directory, "summary.json")
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  json_data <- tryCatch(fromJSON(file_path), error = function(e) NULL)
  if (is.null(json_data)) {
    warning(paste("JSON parse failed:", file_path))
    return(NULL)
  }
  
  data <- data.frame(
    True_pos_baseline = ifelse(!is.null(json_data$`TP-base`), json_data$`TP-base`, NaN),
    True_pos_call     = ifelse(!is.null(json_data$`TP-comp`), json_data$`TP-comp`, NaN),
    False_pos         = ifelse(!is.null(json_data$`FP`), json_data$`FP`, NaN),
    False_neg         = ifelse(!is.null(json_data$`FN`), json_data$`FN`, NaN),
    Precision         = ifelse(!is.null(json_data$precision), as.numeric(json_data$precision), NaN),
    Sensitivity       = ifelse(!is.null(json_data$recall), as.numeric(json_data$recall), NaN),
    Fmeasure          = ifelse(!is.null(json_data$f1), as.numeric(json_data$f1), NaN),
    Caller = caller,
    Genome = genome,
    VarType = vartype
  )
  
  message(paste("Parsed:", file_path, "->", caller, genome, vartype))
  return(data)
}

# Set base paths for SNV/Indel/SV
base_paths <- list(
  SNV   = "/Volumes/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/snv/kmer",
  Indel = "/Volumes/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/indels/kmer",
  SV    = "/Volumes/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/sv/kmer"
)

# Define callers and genomes
snvindel_callers <- c("mutect2", "strelka2", "varscan2", "neusomatic", "deepsomatic")
sv_callers <- c("sniffles", "severus", "pbsv", "svision", "savana")
genomes <- c("GRCh38", "CHM13")
vartypes <- c("SNV", "Indel", "SV")


# --- NEW: Per-genome suffix/label specs  ---
kmer_specs <- function(genome) {
  if (genome == "GRCh38") {
    tibble(
      suffix = c("kmer_k24.umap", "kmer_k36.umap", "kmer_k50.umap", "kmer_k100.umap", "kmer_all"),
      label  = c("k24", "k36", "k50", "k100", "remaining")
    )
  } else if (genome == "CHM13") {
    tibble(
      suffix = c("kmer_k24.umap", "kmer_k36.umap", "kmer_k50.umap", "kmer_k100.umap",
                 "kmer_all_wo_chm13only", "kmer_all_w_chm13only"),
      label  = c("k24", "k36", "k50", "k100", "remaining", "remaining+chm13only")
    )
  } else {
    stop("Unknown genome: ", genome)
  }
}

# Collect data
precision_data_list   <- list()
sensitivity_data_list <- list()

for (vartype in vartypes) {
  base_path <- base_paths[[vartype]]
  callers   <- if (vartype == "SV") sv_callers else snvindel_callers
  
  for (caller in callers) {
    for (genome in genomes) {
      specs <- kmer_specs(genome)
      dirs  <- file.path(base_path, paste0(vartype, "_", caller, "_", genome, "_", specs$suffix))
      
      # Read data (SNV/Indel use summary.txt; SV uses summary.json)
      data_list <- if (vartype == "SV") {
        lapply(dirs, read_summary_json, caller = caller, genome = genome, vartype = vartype)
      } else {
        lapply(dirs, read_summary_file,  caller = caller, genome = genome, vartype = vartype)
      }
      
      # Drop NULLs and align labels accordingly
      keep_idx <- which(!vapply(data_list, is.null, logical(1)))
      if (length(keep_idx) == 0) next
      
      data_list_kept <- data_list[keep_idx]
      labels_kept    <- specs$label[keep_idx]
      
      # Bind and attach file label
      combined <- do.call(rbind, Map(cbind, data_list_kept, file = labels_kept))
      
      # In your original code precision/sensitivity came from the same file
      precision_data_list[[paste0(caller, "_", genome, "_", vartype)]]   <- combined
      sensitivity_data_list[[paste0(caller, "_", genome, "_", vartype)]] <- combined
    }
  }
}

# Combine and label k-mer bins
precision_combined_data   <- do.call(rbind, precision_data_list)
sensitivity_combined_data <- do.call(rbind, sensitivity_data_list)

# Unify the full set of possible labels across genomes (CHM13 has an extra)
all_levels <- c("k24", "k36", "k50", "k100", "remaining", "remaining+chm13only")
precision_combined_data$file   <- factor(precision_combined_data$file,   levels = all_levels)
sensitivity_combined_data$file <- factor(sensitivity_combined_data$file, levels = all_levels)

# Set order for facet plots
precision_combined_data$VarType   <- factor(precision_combined_data$VarType,   levels = c("SNV", "Indel", "SV"))
sensitivity_combined_data$VarType <- factor(sensitivity_combined_data$VarType, levels = c("SNV", "Indel", "SV"))
#color
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
# Plot Precision
precision_plot <- ggplot(precision_combined_data,
                         aes(x = file, y = Precision,
                             color = Caller, linetype = Genome, shape = Genome,
                             group = interaction(Caller, Genome))) +
  geom_point(size = 7) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = color_pal) +
  labs(title = "", x = "k-mer Bin", y = "Precision") +
  ylim(0, 1) +
  scale_x_discrete(drop = FALSE) +
  scale_shape_manual(values = c(chm13 = 16, hg38 = 17)) +  # 16 = filled circle, 17 = triangle
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(. ~ VarType) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.background = element_blank())

# Plot Sensitivity
sensitivity_plot <- ggplot(sensitivity_combined_data,
                           aes(x = file, y = Sensitivity,
                               color = Caller, linetype = Genome, shape = Genome,
                               group = interaction(Caller, Genome))) +
  geom_point(size = 7) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = color_pal) +
  labs(title = "Recall and Precision Comparison by mappability (SNV, Indel & SV)",
       x = "", y = "Recall") +
  ylim(0, 1) +
  scale_x_discrete(drop = FALSE) +
  scale_shape_manual(values = c(chm13 = 16, hg38 = 17)) +  # match above
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

# Save
ggsave("/Users/jinlab/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/Washu/Rotation/PeterJin/Data/Rscript/SMaHT_truthset/Code_for_figures/kmer_precision_sensitivity_plot2.svg",
       final_combined_plot, width = 12, height = 10)
#personal mac
ggsave("/Users/nanakong/Library/CloudStorage/Box-Box/Genetics_Writing_Group/SMaHT_Truth_Set/Figures/Source/Figure4_source/kmer_precision_sensitivity_plot2.svg",
       final_combined_plot, width = 12, height = 10)

write.csv(sensitivity_combined_data, "/Users/nanakong/Library/CloudStorage/Box-Box/Genetics_Writing_Group/SMaHT_Truth_Set/Figures/Source/Figure4_source/4B_kmer_caller_data.csv", row.names = FALSE)
