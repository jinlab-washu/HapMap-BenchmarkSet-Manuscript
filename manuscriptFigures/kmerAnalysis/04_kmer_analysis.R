# 04_kmer_plots.R
#
# Reads the pre-computed per-k differential k-mer TSV files produced by
# 04_kmer_differential.py and generates publication-ready plots.
#
# For each requested k the script produces:
#   1. Paired bar plot  — top N SIGNIFICANT k-mers (Bonferroni < 0.05) by
#                         p-value, showing log2(hotspot/genome) for both
#                         conditions with significance annotations.
#                         If fewer than top_n_bars pass the threshold, only
#                         those are plotted (no non-significant fill).
#   2. Delta strip      — flanking panel showing delta_log2 (noTP − SVHarb).
#   3. Scatter plot     — log2_ratio_svharb vs log2_ratio_notp for all k-mers
#                         (only for k <= k_scatter_max; labeled by significance).
#
# Significance tiers (Bonferroni-corrected p-value):
#   ***  p < 0.001
#   **   p < 0.01
#   *    p < 0.05
#
# NOTE: This script does NO computation — it only reads and plots.
#       Run 04_kmer_differential.py first to generate the TSV files.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(svglite)
  library(patchwork)
})


# ── helpers ───────────────────────────────────────────────────────────────────

sig_label <- function(p) {
  case_when(
    p <  0.001 ~ "***",
    p <  0.01  ~ "**",
    p <  0.05  ~ "*",
    TRUE       ~ NA_character_   # non-significant → not plotted
  )
}

load_differential_k <- function(input_dir, k) {
  f <- file.path(input_dir, sprintf("kmer_differential_k%d.tsv", k))
  if (!file.exists(f)) stop("File not found: ", f,
                            "\n  Run 04_kmer_differential.py first.")
  read_tsv(f, col_types = cols(), show_col_types = FALSE) %>%
    mutate(k         = as.integer(k),
           delta_log2 = -delta_log2)   # Python wrote SVHarb-noTP; flip to noTP-SVHarb
}


# ── Plot 1: Paired bar + delta strip ─────────────────────────────────────────
#
# Only k-mers with Bonferroni p < 0.05 are included, up to top_n.
# Bars show log2(hotspot/genome) for SVHarb (red) and noTP (blue).
# The right strip shows delta_log2 (noTP − SVHarb) on a diverging scale.
# Significance stars are drawn above each bar pair.

plot_bars_k <- function(diff_tbl, k_val, top_n = 30) {
  
  pd <- diff_tbl %>%
    filter(k == k_val) %>%
    # Keep only significant k-mers, then take top N by p-value
    filter(!is.na(pvalue_bonferroni), pvalue_bonferroni < 0.05) %>%
    arrange(pvalue_bonferroni, desc(abs(delta_log2))) %>%
    slice_head(n = top_n) %>%
    mutate(
      sig  = sig_label(pvalue_bonferroni),
      kmer = fct_reorder(kmer, delta_log2)   # order by delta (noTP - SVHarb)
    )
  
  if (nrow(pd) == 0) {
    message(sprintf("  k = %d: no significant k-mers (Bonferroni < 0.05) — skipping bar plot.",
                    k_val))
    return(NULL)
  }
  
  message(sprintf("  k = %d: plotting %d significant k-mers (of %d total).",
                  k_val, nrow(pd),
                  diff_tbl %>% filter(k == k_val) %>% nrow()))
  
  # Long form for grouped bars
  pd_long <- pd %>%
    select(kmer, delta_log2, sig,
           `SV-harboring` = log2_ratio_svharb,
           `noTP`         = log2_ratio_notp) %>%
    pivot_longer(cols      = c(`SV-harboring`, `noTP`),
                 names_to  = "condition",
                 values_to = "log2_ratio") %>%
    mutate(condition = factor(condition, levels = c("SV-harboring", "noTP")))
  
  # Per-kmer annotation y position (just above the taller bar)
  y_range <- max(abs(pd_long$log2_ratio), na.rm = TRUE)
  sig_pos <- pd_long %>%
    group_by(kmer, sig) %>%
    summarise(y_max = max(log2_ratio, na.rm = TRUE), .groups = "drop") %>%
    mutate(y_label = y_max + 0.04 * y_range)
  
  # ── Panel A: grouped bars ──────────────────────────────────────────────────
  p_bars <- ggplot(pd_long,
                   aes(x = kmer, y = log2_ratio, fill = condition)) +
    geom_col(position = position_dodge(width = 0.7),
             width = 0.65, alpha = 0.85) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "gray30") +
    geom_text(data = sig_pos,
              aes(x = kmer, y = y_label, label = sig,
                  fill = NULL, color = NULL),
              size = 3.2, vjust = 0.5, hjust = -0.1,
              inherit.aes = FALSE) +
    scale_x_discrete(position = "top") +
    scale_fill_manual(
      values = c("SV-harboring" = "#2980B9", "noTP" = "#C0392B"),
      name   = "Region"
    ) +
    labs(
      x     = "k-mer",
      y     = expression(log[2]~"(region freq / genome freq)"),
      title = sprintf(
        "k = %d  |  %d significant k-mers (Bonferroni < 0.05)  |  \u0394 = noTP \u2212 SVHarb",
        k_val, nrow(pd))
    ) +
    coord_flip() +
    theme_classic(base_size = 11) +
    theme(
      legend.position  = "top",
      plot.background  = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA)
    )
  
  # ── Panel B: delta strip ───────────────────────────────────────────────────
  p_delta <- ggplot(pd,
                    aes(x = kmer, y = delta_log2, fill = delta_log2)) +
    geom_col(width = 0.65) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "gray30") +
    scale_fill_gradient2(
      low      = "#2980B9",   # negative = more SVHarb-enriched
      mid      = "gray90",
      high     = "#C0392B",   # positive = more noTP-enriched
      midpoint = 0,
      name     = expression(Delta~log[2])
    ) +
    labs(x = NULL,
         y = expression(Delta~log[2]~"(noTP \u2212 SVHarb)")) +
    coord_flip() +
    theme_classic(base_size = 11) +
    theme(
      axis.text.y      = element_blank(),
      axis.ticks.y     = element_blank(),
      plot.background  = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA)
    )
  
  p_bars + p_delta + plot_layout(widths = c(3, 1))
}


# ── main workflow ─────────────────────────────────────────────────────────────

# ╔══════════════════════════════════════════════════════════════════════════════
# ║  USER SETTINGS — edit only this block
# ╠══════════════════════════════════════════════════════════════════════════════

# Directory where 04_kmer_differential.py wrote its output
input_dir  <- "/Volumes/jin810/Active/testing/ztang/code/SMaHT/true_mut_set/202602_resub_kmer/output/Differential_SVHarb_vs_noTP"

# Where to save the plots produced by this script
output_dir <- "/Users/ztang/Library/CloudStorage/Box-Box/MacBookAir_WorkFromHome/SMaHTannualMeetingPlots/20260320_forRevision"

# K values to plot — each must have a corresponding kmer_differential_k{k}.tsv
k_plot <- c(31)

# How many top SIGNIFICANT k-mers (Bonferroni < 0.05) to show in the bar plot.
# If fewer pass the threshold, only those are shown.
top_n_bars <- 30

# Save plots as SVG? Set TRUE once results are final.
save_plots <- TRUE

# ╚══════════════════════════════════════════════════════════════════════════════

# dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

for (k_val in k_plot) {
  
  message(sprintf("=== k = %d ===", k_val))
  
  diff_tbl <- tryCatch(
    load_differential_k(input_dir, k_val),
    error = function(e) { message("  ERROR: ", conditionMessage(e)); NULL }
  )
  if (is.null(diff_tbl)) next
  
  n_total <- nrow(diff_tbl)
  n_sig   <- sum(!is.na(diff_tbl$pvalue_bonferroni) &
                   diff_tbl$pvalue_bonferroni < 0.05)
  message(sprintf("  Loaded %d k-mers  |  Bonferroni < 0.05: %d", n_total, n_sig))
  
  top_n <- min(top_n_bars, 4^k_val)
  
  # ── Paired bar + delta strip ──
  p_bar <- plot_bars_k(diff_tbl, k_val, top_n = top_n)
  if (!is.null(p_bar)) {
    print(p_bar)
    if (save_plots) {
      n_shown  <- diff_tbl %>%
        filter(k == k_val, !is.na(pvalue_bonferroni), pvalue_bonferroni < 0.05) %>%
        nrow() %>% min(top_n)
      out_path <- file.path(output_dir, sprintf("diff_k%02d_bars.svg", k_val))
      ggsave(out_path, plot = p_bar,
             width  = 11,
             height = max(4, n_shown * 0.25 + 2),
             bg     = "transparent")
      message("  Saved: ", out_path)
    }
  }
  
}

message("\nDone. Set save_plots <- TRUE to write SVG files.")