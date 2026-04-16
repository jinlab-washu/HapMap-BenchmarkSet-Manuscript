library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(viridis)
library(readxl)
library(svglite)

# Set save path
save_path <- "20250908_additional/"

# Load SV hotspots once (global variable)
load_sv_hotspots <- function() {
  sv_hotspots_raw <- read_excel("abf7117_ebert_tables-s1-s56.xlsx", 
                                sheet = "14", 
                                skip = 2,  # Skip first 2 rows
                                col_names = TRUE)  # Use row 3 (now row 1) as column names
  
  # Remove last 2 rows and clean the data
  sv_hotspots_clean <- sv_hotspots_raw %>%
    # Remove last 2 rows
    slice(1:(n()-2)) %>%
    # Select the relevant columns
    select(Chr, Start, End) %>%
    # Clean up the data
    filter(!is.na(Chr), !is.na(Start), !is.na(End)) %>%
    mutate(
      # Convert position to Mb
      start = as.numeric(Start),
      end = as.numeric(End),
      chr = Chr,
      start_mb = start / 1e6,
      end_mb = end / 1e6,
      midpoint_mb = (start + end) / 2 / 1e6,
      # Create chromosome number for filtering
      chr_num = case_when(
        chr == "chrX" ~ 23,
        chr == "chrY" ~ 24, 
        chr == "chrM" ~ 25,
        TRUE ~ as.numeric(gsub("chr", "", chr))
      )
    ) %>%
    filter(!is.na(chr_num), chr_num <= 24) %>%  # Remove chrM, keep only chr1-22, X, Y
    arrange(chr_num, start)
  
  return(sv_hotspots_clean)
}

# Load hotspots once at the beginning
sv_hotspots_clean <- load_sv_hotspots()

print("SV hotspots loaded:")
print(paste("Total hotspots:", nrow(sv_hotspots_clean)))
print("Chromosomes with hotspots:")
print(table(sv_hotspots_clean$chr))

# Main plotting function
create_variant_density_plot <- function(csv_file_path, plot_title, include_sv_hotspots = FALSE, save_plot = TRUE, output_filename = NULL) {
  
  # Read the data
  data <- read_csv(csv_file_path)
  
  # Prepare data for plotting
  plot_data <- data %>%
    mutate(
      chr_num = case_when(
        chr == "chrX" ~ 23,
        chr == "chrY" ~ 24, 
        chr == "chrM" ~ 25,
        TRUE ~ as.numeric(gsub("chr", "", chr))
      ),
      pos_mb = (start + end) / 2 / 1e6,
      # Create proper chromosome labels
      chr_label = case_when(
        chr_num <= 22 ~ as.character(chr_num),
        chr == "chrX" ~ "X",
        chr == "chrY" ~ "Y",
        TRUE ~ as.character(chr_num)
      ),
      chr_label = factor(chr_label, levels = c(as.character(1:22), "X", "Y"))
    ) %>%
    filter(!is.na(chr_num), chr_num <= 24) %>%  # Remove chrM (25)
    arrange(chr_num, start)
  
  # Print count ranges for each variant type and category (if columns exist)
  if(all(c("variant_type", "category") %in% colnames(plot_data))) {
    count_ranges <- plot_data %>%
      group_by(variant_type, category) %>%
      summarise(
        min_count = min(count),
        max_count = max(count),
        mean_count = round(mean(count), 2),
        windows_with_variants = sum(count > 0),
        .groups = "drop"
      ) %>%
      arrange(variant_type, category)
    
    print("Count ranges by variant type and category:")
    print(count_ranges)
  }
  
  # Prepare hotspots data if requested
  hotspots_for_plot <- data.frame()  # Empty by default
  if(include_sv_hotspots) {
    hotspots_for_plot <- sv_hotspots_clean %>%
      filter(chr %in% plot_data$chr) %>%
      mutate(chr_label = case_when(
        chr_num <= 22 ~ as.character(chr_num),
        chr == "chrX" ~ "X",
        chr == "chrY" ~ "Y",
        TRUE ~ as.character(chr_num)
      ),
      chr_label = factor(chr_label, levels = c(as.character(1:22), "X", "Y")))
  }
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = pos_mb, y = count)) +
    {if(include_sv_hotspots && nrow(hotspots_for_plot) > 0) 
      geom_rect(
        data = hotspots_for_plot,
        aes(xmin = start_mb, xmax = end_mb, ymin = -Inf, ymax = Inf),
        fill = "red", alpha = 0.4,
        inherit.aes = FALSE
      )
    } +
    geom_point(aes(color = count), size = 0.6, alpha = 0.8) +
    scale_color_viridis_c(
      name = "Variant\nCount",
      option = "viridis",  # Dark = low, bright = high
      direction = 1,       # More variants = brighter color
      trans = "sqrt",
      labels = function(x) round(x, 0)
    ) +
    facet_wrap(~chr_label, scales = "free_x", ncol = 6, nrow = 4,
               labeller = labeller(chr_label = function(x) paste("Chr", x))) +
    labs(
      title = plot_title,
      x = "Position (Mb)",
      y = "Variants per 1Mb Window"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 9, face = "bold"),
      axis.text.x = element_text(size = 6),
      axis.text.y = element_text(size = 7),
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  # Save plot if requested
  if(save_plot) {
    if(is.null(output_filename)) {
      # Generate filename from plot title if not provided
      clean_title <- gsub("[^A-Za-z0-9_-]", "_", plot_title)
      output_filename <- paste0(clean_title, "_density.svg")
    }
    
    ggsave(file.path(save_path, output_filename),
           plot = p, width = 8, height = 6, dpi = 300, bg = "white")
    
    print(paste("Plot saved as:", file.path(save_path, output_filename)))
  }
  
  print(p)
  return(p)
}

## Generate plots for truth set paper 20250910

WashU_nest_wHotspot <- create_variant_density_plot(
  csv_file_path = "SMHTHAPMAP6_GRCh38_v1.0.0_somatic_benchmark_svs_nested_varCount.csv",
  plot_title = "Nested SV Variant Density",
  include_sv_hotspots = TRUE,
  save_plot = TRUE
)

WashU_UNnest_wHotspot <- create_variant_density_plot(
  csv_file_path = "SMHTHAPMAP6_GRCh38_v1.0.0_somatic_benchmark_svs_unnested_varCount.csv",
  plot_title = "Unnested SV Variant Density",
  include_sv_hotspots = TRUE,
  save_plot = TRUE
)

SNV_nested <- create_variant_density_plot(
  csv_file_path = "SMHTHAPMAP6_GRCh38_v1.0.0_somatic_benchmark_snvs_nested_varCount.csv",
  plot_title = "Nested SNV Variant Density",
  include_sv_hotspots = FALSE,
  save_plot = TRUE
)

SNV_unnested <- create_variant_density_plot(
  csv_file_path = "SMHTHAPMAP6_GRCh38_v1.0.0_somatic_benchmark_snvs_unnested_varCount.csv",
  plot_title = "Unnested SNV Variant Density",
  include_sv_hotspots = FALSE,
  save_plot = TRUE
)

indel_nested <- create_variant_density_plot(
  csv_file_path = "SMHTHAPMAP6_GRCh38_v1.0.0_somatic_benchmark_indels_nested_varCount.csv",
  plot_title = "Nested Indel Variant Density",
  include_sv_hotspots = FALSE,
  save_plot = TRUE
)

indel_unnested <- create_variant_density_plot(
  csv_file_path = "SMHTHAPMAP6_GRCh38_v1.0.0_somatic_benchmark_indels_unnested_varCount.csv",
  plot_title = "Unnested Indel Variant Density",
  include_sv_hotspots = FALSE,
  save_plot = TRUE
)
