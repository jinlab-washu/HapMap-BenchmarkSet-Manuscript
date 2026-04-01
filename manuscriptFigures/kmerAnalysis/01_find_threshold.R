# 01_find_threshold.R
# Find optimal N: SV count > N windows maximize bp-overlap with SV hotspot regions.
#
# Inputs:
#   sv_windows_csv  : the *_varCount.csv files (cols: chr, start, end, count)
#   hotspot_xlsx    : abf7117_ebert_tables-s1-s56.xlsx  (sheet "14", skip 2 rows)
#
# Outputs (written to output_dir):
#   threshold_optimization.pdf   – precision / recall / F1 curve
#   threshold_sweep.tsv          – full metrics table
#   high_sv_windows.bed          – windows with count > optimal N
#   hotspot_high_sv_overlap.bed  – hotspot regions that overlap high-SV windows

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(valr)
})

# ── loaders ───────────────────────────────────────────────────────────────────

load_sv_windows <- function(csv_path) {
  read_csv(csv_path, show_col_types = FALSE) %>%
    filter(!is.na(chr), !is.na(start), !is.na(end), !is.na(count)) %>%
    mutate(
      chr_num = case_when(
        chr == "chrX" ~ 23L,
        chr == "chrY" ~ 24L,
        TRUE ~ as.integer(gsub("chr", "", chr))
      )
    ) %>%
    filter(!is.na(chr_num), chr_num <= 24) %>%   # drop chrM
    select(chrom = chr, start, end, sv_count = count)
}

load_sv_hotspots <- function(path,
                             sheet = "14",
                             skip  = 2) {
  ext <- tolower(tools::file_ext(path))
  
  if (ext %in% c("xlsx", "xls")) {
    df <- read_excel(path, sheet = sheet, skip = skip, col_names = TRUE) %>%
      slice(1:(n() - 2)) %>%
      select(Chr, Start, End) %>%
      filter(!is.na(Chr), !is.na(Start), !is.na(End)) %>%
      mutate(chrom = Chr,
             start = as.integer(Start),
             end   = as.integer(End))
    
  } else if (ext == "bed" || ext == "txt" || ext == "") {
    df <- read_tsv(path, col_names = c("chrom", "start", "end"),
                   col_types = "cii", show_col_types = FALSE) %>%
      filter(!is.na(chrom), !is.na(start), !is.na(end))
    
  }
  
  df %>%
    mutate(
      chr_num = case_when(
        chrom == "chrX" ~ 23L,
        chrom == "chrY" ~ 24L,
        TRUE ~ as.integer(gsub("chr", "", chrom))
      )
    ) %>%
    filter(!is.na(chr_num), chr_num <= 24) %>%
    select(chrom, start, end) %>%
    arrange(chrom, start)
}

# ── interval intersection helper ──────────────────────────────────────────────

# Replaces bed_intersect() to avoid valr column naming issues (.start.x etc.)
intersect_bp <- function(a, b) {
  joined <- inner_join(
    a %>% rename(a_start = start, a_end = end),
    b %>% rename(b_start = start, b_end = end),
    by = "chrom", relationship = "many-to-many"
  ) %>%
    filter(a_start < b_end, b_start < a_end) %>%
    mutate(start = pmax(a_start, b_start),
           end   = pmin(a_end,   b_end)) %>%
    select(chrom, start, end)
  
  if (nrow(joined) == 0) return(0L)
  
  joined %>%
    bed_merge() %>%
    summarise(bp = sum(end - start)) %>%
    pull(bp)
}

# ── threshold sweep ───────────────────────────────────────────────────────────

sweep_thresholds <- function(windows, hotspots) {
  
  hotspots_m     <- bed_merge(hotspots)
  bp_hotspot_tot <- sum(hotspots_m$end - hotspots_m$start)
  
  # Only test thresholds at actual observed count values
  thresholds <- sort(unique(windows$sv_count))
  
  map_dfr(thresholds, function(n) {
    
    high_sv <- windows %>%
      filter(sv_count > n) %>%
      select(chrom, start, end)
    
    if (nrow(high_sv) == 0) return(NULL)
    
    high_sv_m  <- bed_merge(high_sv)
    bp_high_sv <- sum(high_sv_m$end - high_sv_m$start)
    bp_ix      <- intersect_bp(high_sv_m, hotspots_m)
    
    precision <- if (bp_high_sv   > 0) bp_ix / bp_high_sv     else 0
    recall    <- if (bp_hotspot_tot > 0) bp_ix / bp_hotspot_tot else 0
    f1        <- if ((precision + recall) > 0)
      2 * precision * recall / (precision + recall) else 0
    
    tibble(n         = n,
           n_windows = nrow(high_sv),
           bp_high   = bp_high_sv,
           bp_ix     = bp_ix,
           precision = precision,
           recall    = recall,
           f1        = f1)
  })
}

# ── plot ──────────────────────────────────────────────────────────────────────

plot_threshold_sweep <- function(results, optimal_n, out_pdf) {
  
  best <- filter(results, n == optimal_n)
  
  p <- results %>%
    pivot_longer(c(precision, recall, f1),
                 names_to = "metric", values_to = "score") %>%
    mutate(metric = factor(metric,
                           levels = c("f1", "precision", "recall"),
                           labels = c("F1", "Precision", "Recall"))) %>%
    ggplot(aes(x = n, y = score, color = metric)) +
    geom_line(linewidth = 0.9) +
    geom_vline(xintercept = optimal_n,
               linetype = "dashed", color = "black", alpha = 0.7) +
    annotate("text",
             x     = optimal_n,
             y     = min(results$f1, na.rm = TRUE),
             label = sprintf("N = %d\nF1 = %.3f", optimal_n, best$f1),
             hjust = -0.08, vjust = 0, size = 3.4, color = "black") +
    scale_color_manual(
      values = c(F1 = "#4DAF4A", Precision = "#E41A1C", Recall = "#377EB8"),
      name   = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    labs(
      title = "Threshold Optimization: High-SV Windows vs SV Hotspot Overlap",
      x     = "SV Count Threshold (N)",
      y     = "Score (base-pair level)"
    ) +
    theme_classic(base_size = 13) +
    theme(legend.position   = "top",
          plot.background   = element_rect(fill = "white", color = NA),
          panel.background  = element_rect(fill = "white", color = NA))
  
  ggsave(out_pdf, p, width = 7, height = 5, bg = "white")
  message("Saved: ", out_pdf)
}

# ── main ──────────────────────────────────────────────────────────────────────

main <- function(
    sv_windows_csv,
    hotspot_path,
    output_dir    = "output",
    hotspot_sheet = "14",
    hotspot_skip  = 2,
    metric        = "f1",   # optimize by: "f1", "precision", or "recall"
    use_threshold = TRUE
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  message("Loading SV windows from: ", basename(sv_windows_csv))
  windows <- load_sv_windows(sv_windows_csv)
  message(sprintf("  %d windows  |  count range: %d - %d",
                  nrow(windows), min(windows$sv_count), max(windows$sv_count)))
  
  message("Loading SV hotspots from: ", basename(hotspot_path))
  hotspots <- load_sv_hotspots(hotspot_path, hotspot_sheet, hotspot_skip)
  message(sprintf("  %d hotspot regions across %d chromosomes",
                  nrow(hotspots), n_distinct(hotspots$chrom)))
  
  # Diagnostic: check chromosome name consistency between the two datasets
  shared <- intersect(unique(windows$chrom), unique(hotspots$chrom))
  only_w <- setdiff(unique(windows$chrom),  unique(hotspots$chrom))
  only_h <- setdiff(unique(hotspots$chrom), unique(windows$chrom))
  message(sprintf("  Shared chroms: %d  |  Only in windows: %s  |  Only in hotspots: %s",
                  length(shared),
                  if (length(only_w) > 0) paste(head(only_w, 5), collapse = ", ") else "none",
                  if (length(only_h) > 0) paste(head(only_h, 5), collapse = ", ") else "none"))
  if (length(shared) == 0)
    stop("No shared chromosome names. Check chr naming (e.g. 'chr1' vs '1').")
  
  # Restrict both to shared chromosomes only
  windows  <- filter(windows,  chrom %in% shared)
  hotspots <- filter(hotspots, chrom %in% shared)
  
  message("Sweeping thresholds (this may take a moment) ...")
  results <- sweep_thresholds(windows, hotspots)
  
  # Pick optimal N by chosen metric
  optimal_n <- results %>%
    slice_max(.data[[metric]], n = 1, with_ties = FALSE) %>%
    pull(n)
  
  best <- filter(results, n == optimal_n)
  message(sprintf(
    "\nOptimal N = %d\n  F1        = %.4f\n  Precision = %.4f\n  Recall    = %.4f\n  Windows   = %d  |  Overlap bp = %d",
    optimal_n, best$f1, best$precision, best$recall, best$n_windows, best$bp_ix))
  
  # Save sweep table
  write_tsv(results, file.path(output_dir, "threshold_sweep.tsv"))
  message("Saved: threshold_sweep.tsv")
  
  # Plot
  plot_threshold_sweep(results, optimal_n,
                       file.path(output_dir, "threshold_optimization.pdf"))
  
  # Write high-SV windows BED (count > optimal N)
  high_sv_bed_path <- file.path(output_dir, "high_sv_windows.bed")
  windows %>%
    filter(sv_count > optimal_n) %>%
    arrange(chrom, start) %>%
    write_tsv(high_sv_bed_path, col_names = FALSE)
  message("Saved: high_sv_windows.bed")
  
  # Write overlap BED: regions that are in BOTH hotspots AND high-SV windows
  hotspots_m <- bed_merge(hotspots)
  high_sv_m  <- windows %>%
    filter(sv_count > optimal_n) %>%
    select(chrom, start, end) %>%
    bed_merge()
  
  overlap_bed_path <- file.path(output_dir, "hotspot_high_sv_overlap.bed")
  if (use_threshold) {
    # Intersect hotspots with high-SV windows
    inner_join(
      hotspots_m %>% rename(a_start = start, a_end = end),
      high_sv_m  %>% rename(b_start = start, b_end = end),
      by = "chrom", relationship = "many-to-many"
    ) %>%
      filter(a_start < b_end, b_start < a_end) %>%
      mutate(start = pmax(a_start, b_start),
             end   = pmin(a_end,   b_end)) %>%
      select(chrom, start, end) %>%
      bed_merge() %>%
      write_tsv(overlap_bed_path, col_names = FALSE)
    message("Saved: hotspot_high_sv_overlap.bed (threshold applied: N > ", optimal_n, ")")
  } else {
    # No threshold: output all hotspot regions directly
    hotspots_m %>%
      write_tsv(overlap_bed_path, col_names = FALSE)
    message("Saved: hotspot_high_sv_overlap.bed (threshold NOT applied — all hotspot regions used)")
  }
}

# ── run ───────────────────────────────────────────────────────────────────────

SV_WINDOWS_CSV <- "/Volumes/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/zt_temp/CompareSVs_20250814/generate_plots/SMHTHAPMAP6_GRCh38_v1.0.0_somatic_benchmark_svs_nested_varCount.csv"


# ### With using high SV threshold ###
# 
# ## SV Hotspot vs. Entire genome ##
# result <- main(
#   sv_windows_csv = SV_WINDOWS_CSV,
#   hotspot_path   = "/Users/ztang/Downloads/abf7117_ebert_tables-s1-s56.xlsx",
#   output_dir     = "/Volumes/jin810/Active/testing/ztang/code/SMaHT/true_mut_set/202602_resub_kmer/output/HotSpot_vs_AllGen_N28",
#   metric         = "f1",
#   use_threshold  = TRUE
# )
# 
# ## SV Harboring vs. Entire genome ##
# result <- main(
#   sv_windows_csv = SV_WINDOWS_CSV,
#   hotspot_path   = "/Volumes/jin810/Active/testing/ztang/code/SMaHT/truth_set_comparison/testDifficultRegions_20250828/sv_harbor_AllVar.bed",
#   output_dir     = "/Volumes/jin810/Active/testing/ztang/code/SMaHT/true_mut_set/202602_resub_kmer/output/SVHarb_vs_AllGen_N",
#   metric         = "f1",
#   use_threshold  = TRUE 
# )
# 
# ## noTP vs. Entire genome ##
# result <- main(
#   sv_windows_csv = SV_WINDOWS_CSV,
#   hotspot_path   = "/Volumes/jin810/Active/testing/ztang/code/SMaHT/truth_set_comparison/testDifficultRegions_20250828/noTP_1_AllVar_inSVHarbor.bed",
#   output_dir     = "/Volumes/jin810/Active/testing/ztang/code/SMaHT/true_mut_set/202602_resub_kmer/output/noTP_vs_AllGen_N",
#   metric         = "f1",
#   use_threshold  = TRUE 
# )



### Without using high SV threshold
## SV Hotspot vs. Entire genome ##
result <- main(
  sv_windows_csv = SV_WINDOWS_CSV,
  hotspot_path   = "/Users/ztang/Downloads/abf7117_ebert_tables-s1-s56.xlsx",
  output_dir     = "/Volumes/jin810/Active/testing/ztang/code/SMaHT/true_mut_set/202602_resub_kmer/output/HotSpot_vs_AllGen_NoFilter",
  metric         = "f1",
  use_threshold  = FALSE
)

## SV Harboring vs. Entire genome ##
result <- main(
  sv_windows_csv = SV_WINDOWS_CSV,
  hotspot_path   = "/Volumes/jin810/Active/testing/ztang/code/SMaHT/truth_set_comparison/testDifficultRegions_20250828/sv_harbor_AllVar.bed",
  output_dir     = "/Volumes/jin810/Active/testing/ztang/code/SMaHT/true_mut_set/202602_resub_kmer/output/SVHarb_vs_AllGen_NoFilter",
  metric         = "f1",
  use_threshold  = FALSE 
)

## noTP vs. Entire genome ##
result <- main(
  sv_windows_csv = SV_WINDOWS_CSV,
  hotspot_path   = "/Volumes/jin810/Active/testing/ztang/code/SMaHT/truth_set_comparison/testDifficultRegions_20250828/noTP_1_AllVar_inSVHarbor.bed",
  output_dir     = "/Volumes/jin810/Active/testing/ztang/code/SMaHT/true_mut_set/202602_resub_kmer/output/noTP_vs_AllGen_NoFilter",
  metric         = "f1",
  use_threshold  = FALSE 
)




# result <- main(
#   sv_windows_csv = SV_WINDOWS_CSV,
#   hotspot_path   = HOTSPOT_XLSX,
#   output_dir     = OUTPUT_DIR,
#   metric         = "f1",
#   use_threshold  = FALSE 
# )