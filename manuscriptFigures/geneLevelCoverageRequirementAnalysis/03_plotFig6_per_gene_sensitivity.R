#!/usr/bin/env Rscript

#################################################################
# Fit per-gene sensitivity vs coverage curves and generate Figure 6 panels and Supplementary table 7-8.
#
# This script:
#   1. Fits an exponential saturation model per gene:
#        Sensitivity = 1 - exp(-b * Coverage)
#   2. Estimates coverage thresholds for target sensitivities
#   3. Writes an output TSV with fitted parameters
#   4. Generates:
#        - gene-level fitted curve plot
#        - histogram of predicted sensitivity at 30X
#        - histogram of coverage required for 80% sensitivity
#
# Usage:
#   Rscript 03_plotFig6_per_gene_sensitivity.R
#
# Optional path edits are needed below if running in a new environment.
#
# Required packages:
#   readr, dplyr, tidyr, purrr, ggplot2, ggbreak
#
# Author: Nahyun Kong
# Contact: nahyun@wustl.edu
#################################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(ggbreak)
})

############################
# Input / output paths
############################

infile <- "/Volumes/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/snv/per_gene/per_gene_vcfeval_runs_final/per_gene_sensitivity.tsv"

outfile_tsv <- "/Volumes/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/snv/per_gene/per_gene_vcfeval_runs_final/per_gene_sensitivity_with_cov_targets_exponential2.tsv"

outfile_svg_curve <- "/Users/nanakong/Library/CloudStorage/Box-Box/Genetics_Writing_Group/SMaHT_Truth_Set/Figures/Source/Figure6/6C_sensitivity_coverage_gene_plot.svg"

outfile_hist_sens30 <- "/Users/nanakong/Library/CloudStorage/Box-Box/Genetics_Writing_Group/SMaHT_Truth_Set/Figures/Source/Figure6/revision/6A_sens_cov30.svg"

outfile_hist_cov80 <- "/Users/nanakong/Library/CloudStorage/Box-Box/Genetics_Writing_Group/SMaHT_Truth_Set/Figures/Source/Figure6/revision/6B_hist_cov.svg"

############################
# Load and clean input table
############################

df_raw <- read_tsv(infile, show_col_types = FALSE)

df <- df_raw %>% distinct()
message("Removed ", nrow(df_raw) - nrow(df), " duplicate rows (exact matches).")

############################
# Coverage per lab
############################

coverage_map <- c(
  washu = 510,
  BCM   = 463,
  NYGC  = 271,
  BROAD = 175
)

############################
# Reshape to long format
############################

df_long <- df %>%
  filter(Truthset_number >= 10) %>%
  pivot_longer(
    cols = starts_with("Sensitivity"),
    names_to = "Lab",
    values_to = "Sensitivity"
  ) %>%
  mutate(
    Coverage = case_when(
      Lab == "Sensitivity_washu" ~ coverage_map["washu"],
      Lab == "Sensitivity_BCM"   ~ coverage_map["BCM"],
      Lab == "Sensitivity_NYGC"  ~ coverage_map["NYGC"],
      Lab == "Sensitivity_BROAD" ~ coverage_map["BROAD"],
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(Sensitivity), !is.na(Coverage))

############################
# Fit exponential model per gene
############################

fit_one_gene <- function(dat) {
  if (nrow(dat) < 2 || dplyr::n_distinct(dat$Coverage) < 2) {
    return(tibble(
      cov_75 = NA_real_,
      cov_80 = NA_real_,
      cov_90 = NA_real_,
      cov_95 = NA_real_,
      b = NA_real_,
      sens_cov30 = NA_real_
    ))
  }

  b_start <- 0.005

  fit <- try(
    nls(
      Sensitivity ~ 1 - exp(-b * Coverage),
      data = dat,
      start = list(b = b_start),
      algorithm = "port",
      lower = c(b = 1e-6),
      upper = c(b = 1.0)
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {
    return(tibble(
      cov_75 = NA_real_,
      cov_80 = NA_real_,
      cov_90 = NA_real_,
      cov_95 = NA_real_,
      b = NA_real_,
      sens_cov30 = NA_real_
    ))
  }

  b <- unname(coef(fit)[["b"]])

  cov_fun <- function(target) {
    if (is.na(b) || b <= 0 || target <= 0 || target >= 1) {
      return(NA_real_)
    }
    x <- -log(1 - target) / b
    if (!is.finite(x) || x < 0) NA_real_ else x
  }

  sens_at <- function(cov) {
    if (is.na(b) || b <= 0) {
      return(NA_real_)
    }
    val <- 1 - exp(-b * cov)
    pmin(pmax(val, 0), 1)
  }

  tibble(
    cov_75 = cov_fun(0.75),
    cov_80 = cov_fun(0.80),
    cov_90 = cov_fun(0.90),
    cov_95 = cov_fun(0.95),
    b = b,
    sens_cov30 = sens_at(30)
  )
}

predictions <- df_long %>%
  group_by(gene) %>%
  group_modify(~ fit_one_gene(.x)) %>%
  ungroup()

############################
# Save augmented TSV
############################

orig_cols <- names(df)[1:min(7, ncol(df))]

out <- df %>%
  select(all_of(orig_cols)) %>%
  left_join(predictions, by = "gene")

write_tsv(out, outfile_tsv)
message("Wrote TSV file: ", outfile_tsv)

############################
# Add "all" genomic region row
############################

df_long_plot <- df_long %>%
  bind_rows(
    tibble(
      gene = "all",
      Truthset_number = 2000160L,
      Lab = c(
        "Sensitivity_washu",
        "Sensitivity_BCM",
        "Sensitivity_NYGC",
        "Sensitivity_BROAD"
      ),
      Sensitivity = c(0.1227, 0.1255, 0.1040, 0.0521),
      Coverage = c(510, 463, 271, 175)
    )
  )

############################
# Plot 1: selected genes with fitted curves
############################

top_rows <- df %>%
  drop_na(Truthset_number) %>%
  slice_max(order_by = Truthset_number, n = 5, with_ties = FALSE)

genes_to_plot <- top_rows %>% pull(gene)
genes_to_plot <- c("BCL2L1-AS1", "MUC1", "PMS2", "SBDS", "UNC93A", "all")

xmax_curve <- 600

curve_grid <- predictions %>%
  filter(gene %in% genes_to_plot) %>%
  select(gene, b) %>%
  tidyr::expand_grid(Coverage = seq(0, xmax_curve, by = 1)) %>%
  mutate(
    Sens_pred = ifelse(is.finite(b), 1 - exp(-b * Coverage), NA_real_)
  )

label_pos_x <- 0.88 * xmax_curve

label_df <- predictions %>%
  filter(gene %in% genes_to_plot) %>%
  mutate(Coverage = label_pos_x) %>%
  left_join(
    curve_grid %>% select(gene, Coverage, Sens_pred),
    by = c("gene", "Coverage")
  ) %>%
  mutate(
    label = ifelse(
      is.na(b),
      "fit failed",
      sprintf("y = 1 - e^{-b·x}\nb = %.4f", b)
    ),
    Sens_pred = pmin(pmax(Sens_pred, 0), 1)
  )

p_curve <- ggplot() +
  geom_point(
    data = df_long_plot %>% filter(gene %in% genes_to_plot),
    aes(x = Coverage, y = Sensitivity, color = gene),
    size = 3
  ) +
  geom_line(
    data = curve_grid,
    aes(x = Coverage, y = Sens_pred, color = gene),
    linewidth = 1
  ) +
  geom_text(
    data = label_df,
    aes(x = Coverage, y = Sens_pred, label = label, color = gene),
    hjust = 1,
    vjust = -0.6,
    size = 3.5,
    show.legend = FALSE,
    lineheight = 1.05
  ) +
  scale_x_continuous(limits = c(0, xmax_curve)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Variant Detection Rate (MuTect2) vs Coverage (SR)",
    subtitle = paste("Genes:", paste(genes_to_plot, collapse = ", ")),
    x = "Coverage (X)",
    y = "Variant Detection Rate"
  ) +
  theme_minimal(base_size = 12)

print(p_curve)
ggsave(outfile_svg_curve, p_curve, width = 7, height = 5, dpi = 300)
message("Wrote SVG file: ", outfile_svg_curve)

############################
# Plot 2: histogram of sensitivity at 30X
############################

predictions_filtered <- out %>% filter(Truthset_number >= 10)

df_sens30 <- predictions_filtered %>%
  select(gene, sens_cov30) %>%
  filter(!is.na(sens_cov30))

mean_sens30 <- mean(df_sens30$sens_cov30, na.rm = TRUE)

p_hist_sens30 <- ggplot(df_sens30, aes(x = sens_cov30)) +
  geom_histogram(
    bins = 30,
    boundary = 0,
    closed = "left",
    fill = "steelblue"
  ) +
  geom_vline(
    xintercept = mean_sens30,
    linetype = "dashed",
    linewidth = 0.9
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "SNV Detection Rate at 30× for Clinically Relevent Gene (VAF 1%)",
    subtitle = sprintf("Mean = %.3f", mean_sens30),
    x = "Estimated MuTect2 SNV Detection Rate at 30×",
    y = "Gene count"
  ) +
  theme_minimal(base_size = 12)

print(p_hist_sens30)
ggsave(outfile_hist_sens30, p_hist_sens30, width = 6.5, height = 4.5, dpi = 300)
message("Wrote SVG file: ", outfile_hist_sens30)

############################
# Plot 3: histogram of coverage needed for 80% sensitivity
############################

n_total <- sum(!is.na(predictions$cov_80))
n_over <- sum(predictions$cov_80 > 500, na.rm = TRUE)
pct_over <- 100 * n_over / n_total

message(sprintf(
  "Genes requiring >500X: %d / %d (%.2f%%)",
  n_over, n_total, pct_over
))

xmax_hist <- max(predictions$cov_80, na.rm = TRUE)

p_hist_cov80 <- predictions %>%
  select(gene, cov_80) %>%
  filter(!is.na(cov_80)) %>%
  ggplot(aes(x = cov_80)) +
  geom_histogram(
    binwidth = 20,
    boundary = 0,
    closed = "left",
    fill = "steelblue"
  ) +
  geom_vline(xintercept = 500, linetype = "dashed") +
  labs(
    title = "Distribution of Coverage Needed for 80% Sensitivity",
    subtitle = sprintf(
      "Vertical line at 500×; %.2f%% of genes require >500× (n=%d/%d)",
      pct_over, n_over, n_total
    ),
    x = "Coverage for 80% detection (cov_80, X)",
    y = "Gene count"
  ) +
  theme_minimal(base_size = 12)

if (is.finite(xmax_hist) && xmax_hist > 51000) {
  p_hist_cov80 <- p_hist_cov80 +
    scale_x_break(c(2000, 32500)) +
    scale_x_break(c(32700, 51600)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.02)))
} else {
  p_hist_cov80 <- p_hist_cov80 +
    coord_cartesian(xlim = c(0, xmax_hist)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.02)))
}

print(p_hist_cov80)
ggsave(outfile_hist_cov80, p_hist_cov80, width = 6.5, height = 4.5, dpi = 300)
message("Wrote SVG file: ", outfile_hist_cov80)
