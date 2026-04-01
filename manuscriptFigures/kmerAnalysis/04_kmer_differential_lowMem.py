
"""
04_kmer_differential.py

For each k, compares k-mer enrichment between SVHarb and noTP hotspot regions,
both expressed relative to their respective (separate) genome backgrounds.

Statistical test per k-mer
--------------------------
Primary 2x2 contingency table (tests differential proportion between regions):

                | this k-mer | all other k-mers |
  SVHarb hotspot|     a      |        b         |
  noTP hotspot  |     c      |        d         |

  Chi-squared (no Yates correction) is used by default.
  Falls back to Fisher's exact test when any expected cell count < 5.
  Bonferroni correction is applied within each k independently.

Enrichment metadata (reported, not tested separately)
------------------------------------------------------
  For each condition vs its genome background:
    log2_ratio    = log2(freq_hotspot / freq_genome)         -- frequency-based
    enrichment_or = (count_h / rest_h) / (count_g / rest_g) -- odds-ratio-based

  These let you distinguish two cases downstream in R:
    Case 1: enriched in both vs genome, but more so in noTP
    Case 2: not enriched (or depleted) in SVHarb vs genome, enriched in noTP

Memory design
-------------
  Uses plain Python dicts {kmer: count} instead of pandas DataFrames for the
  four raw count files. This avoids pandas column overhead and the temporary
  copies created during pd.merge, reducing peak memory by ~4-6x compared to
  the naive DataFrame approach. The dicts are explicitly deleted before the
  result DataFrame is constructed so peak memory is:
    4 dicts + result rows list  (not 4 dicts + enrichment DFs + merged DF).

Output
------
  One TSV per k in <output_dir>/kmer_differential_k{k}.tsv
  Sorted by pvalue_bonferroni ascending (most significant first).
  If the output file already exists for a given k, that k is skipped.

  Columns:
    kmer
    count_svharb, count_notp                         -- raw hotspot counts
    count_genome_sv, count_genome_notp               -- raw genome counts
    enrichment_or_svharb, log2_ratio_svharb          -- SVHarb vs its genome
    enrichment_or_notp,   log2_ratio_notp            -- noTP   vs its genome
    delta_log2                                       -- log2_ratio_svharb - log2_ratio_notp
    OR_differential                                  -- OR from the 2x2 differential table
    pvalue_raw, pvalue_bonferroni
    test_used                                        -- "chi2", "fisher", or "skipped"

Usage
-----
  python 04_kmer_differential.py \
    --svharb_hotspot_prefix /path/to/SVHarb_vs_AllGen_NoFilter/hotspot \
    --svharb_genome_prefix  /path/to/SVHarb_vs_AllGen_NoFilter/genome  \
    --notp_hotspot_prefix   /path/to/noTP_vs_AllGen_NoFilter/hotspot   \
    --notp_genome_prefix    /path/to/noTP_vs_AllGen_NoFilter/genome    \
    --output_dir            /path/to/Differential_SVHarb_vs_noTP       \
    --k_values 2 3 4 5      # omit to default to 2-31
"""

import math
import os
import sys
import argparse

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact


# ── argument parsing ──────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="Differential k-mer enrichment: SVHarb vs noTP",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument("--svharb_hotspot_prefix", required=True,
                        help="File prefix for SVHarb hotspot k-mer counts "
                             "(files expected as <prefix>_k{k}.tsv)")
    parser.add_argument("--svharb_genome_prefix",  required=True,
                        help="File prefix for SVHarb genome k-mer counts")
    parser.add_argument("--notp_hotspot_prefix",   required=True,
                        help="File prefix for noTP hotspot k-mer counts")
    parser.add_argument("--notp_genome_prefix",    required=True,
                        help="File prefix for noTP genome k-mer counts")
    parser.add_argument("--output_dir", required=True,
                        help="Directory where per-k TSV files will be written")
    parser.add_argument("--k_values", nargs="+", type=int,
                        default=list(range(2, 32)),
                        help="k values to process (default: 2 through 31)")
    return parser.parse_args()


# ── I/O helpers ───────────────────────────────────────────────────────────────

def load_kmer_counts(prefix: str, k: int) -> tuple[dict, int]:
    """
    Load a k-mer count TSV directly into a dict {kmer: count}.
    Returns (counts_dict, total_count).
    Much lower memory than a pandas DataFrame — no column overhead, no copies.
    """
    path = f"{prefix}_k{k}.tsv"
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing file: {path}")
    counts = {}
    total = 0
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            kmer, count_str = line.split("\t")
            c = int(count_str)
            counts[kmer] = c
            total += c
    return counts, total


# ── statistical test ──────────────────────────────────────────────────────────

def run_contingency_test(a: float, b: float,
                         c: float, d: float) -> tuple[float, float, str]:
    """
    Test whether k-mer proportions differ between SVHarb and noTP hotspots.

    2x2 table:
        [[a, b],   <- SVHarb hotspot: (this k-mer count, rest-of-k-mers count)
         [c, d]]   <- noTP   hotspot: (this k-mer count, rest-of-k-mers count)

    Uses chi-squared (no Yates correction) when all expected cell counts >= 5,
    otherwise falls back to Fisher's exact test (two-sided).

    Returns
    -------
    pvalue   : float
    or_value : float  -- odds ratio (a*d)/(b*c)
    test_used: str    -- "chi2" | "fisher" | "skipped"
    """
    table = np.array([[a, b], [c, d]], dtype=float)

    # Guard: any cell <= 0 that would break the test
    if b <= 0 or d <= 0 or a < 0 or c < 0:
        return float("nan"), float("nan"), "skipped"

    # Expected cell counts under H0
    row_sums = table.sum(axis=1, keepdims=True)
    col_sums = table.sum(axis=0, keepdims=True)
    grand    = table.sum()
    if grand == 0:
        return float("nan"), float("nan"), "skipped"
    expected = row_sums * col_sums / grand

    or_value = (a * d) / (b * c) if (b * c) > 0 else float("nan")

    if (expected < 5).any():
        # Fisher's exact (two-sided)
        res = fisher_exact(table.astype(int), alternative="two-sided")
        return float(res.pvalue), float(res.statistic), "fisher"
    else:
        chi2, pval, _, _ = chi2_contingency(table, correction=False)
        return float(pval), float(or_value), "chi2"


# ── per-k processing ──────────────────────────────────────────────────────────

def process_k(k: int,
              svharb_hotspot_prefix: str, svharb_genome_prefix: str,
              notp_hotspot_prefix:   str, notp_genome_prefix:   str) -> pd.DataFrame:
    """
    Memory-efficient pipeline for one k value.
    Uses dicts instead of DataFrames to avoid pandas merge overhead and copies.
    Streams through the kmer intersection in a single pass.
    """
    # ── load all four files as dicts ──────────────────────────────────────────
    sv_hot,   total_sv_hot   = load_kmer_counts(svharb_hotspot_prefix, k)
    sv_gen,   total_sv_gen   = load_kmer_counts(svharb_genome_prefix,  k)
    notp_hot, total_notp_hot = load_kmer_counts(notp_hotspot_prefix,   k)
    notp_gen, total_notp_gen = load_kmer_counts(notp_genome_prefix,    k)

    # ── find kmers present in all four files ──────────────────────────────────
    shared_kmers = (sv_hot.keys()
                    & sv_gen.keys()
                    & notp_hot.keys()
                    & notp_gen.keys())
    n_tests = len(shared_kmers)
    if n_tests == 0:
        return pd.DataFrame()

    # ── stream through shared kmers, compute everything in one pass ───────────
    rows = []
    for kmer in shared_kmers:
        a_sv     = sv_hot[kmer]
        a_sv_g   = sv_gen[kmer]
        a_notp   = notp_hot[kmer]
        a_notp_g = notp_gen[kmer]

        # --- enrichment vs genome (SVHarb) ---
        freq_sv_hot = a_sv   / total_sv_hot
        freq_sv_gen = a_sv_g / total_sv_gen
        if freq_sv_hot > 0 and freq_sv_gen > 0:
            log2_sv = math.log2(freq_sv_hot / freq_sv_gen)
        else:
            log2_sv = float("nan")
        rest_sv_hot = total_sv_hot - a_sv
        rest_sv_gen = total_sv_gen - a_sv_g
        or_sv = ((a_sv / rest_sv_hot) / (a_sv_g / rest_sv_gen)
                 if rest_sv_hot > 0 and rest_sv_gen > 0 and a_sv_g > 0
                 else float("nan"))

        # --- enrichment vs genome (noTP) ---
        freq_notp_hot = a_notp   / total_notp_hot
        freq_notp_gen = a_notp_g / total_notp_gen
        if freq_notp_hot > 0 and freq_notp_gen > 0:
            log2_notp = math.log2(freq_notp_hot / freq_notp_gen)
        else:
            log2_notp = float("nan")
        rest_notp_hot = total_notp_hot - a_notp
        rest_notp_gen = total_notp_gen - a_notp_g
        or_notp = ((a_notp / rest_notp_hot) / (a_notp_g / rest_notp_gen)
                   if rest_notp_hot > 0 and rest_notp_gen > 0 and a_notp_g > 0
                   else float("nan"))

        # --- differential 2x2 test ---
        b = total_sv_hot   - a_sv    # rest of SVHarb hotspot
        d = total_notp_hot - a_notp  # rest of noTP hotspot
        pval, or_diff, test = run_contingency_test(a_sv, b, a_notp, d)

        rows.append({
            "kmer":                  kmer,
            "count_svharb":          a_sv,
            "count_notp":            a_notp,
            "count_genome_sv":       a_sv_g,
            "count_genome_notp":     a_notp_g,
            "enrichment_or_svharb":  or_sv,
            "log2_ratio_svharb":     log2_sv,
            "enrichment_or_notp":    or_notp,
            "log2_ratio_notp":       log2_notp,
            "delta_log2":            (log2_sv - log2_notp
                                      if not (math.isnan(log2_sv) or math.isnan(log2_notp))
                                      else float("nan")),
            "OR_differential":       or_diff,
            "pvalue_raw":            pval,
            "test_used":             test,
        })

    # ── free the dicts immediately — no longer needed ─────────────────────────
    del sv_hot, sv_gen, notp_hot, notp_gen

    # ── build result DataFrame and apply Bonferroni ───────────────────────────
    result = pd.DataFrame(rows)
    valid  = result["pvalue_raw"].notna()
    result["pvalue_bonferroni"] = float("nan")
    result.loc[valid, "pvalue_bonferroni"] = (
        result.loc[valid, "pvalue_raw"] * n_tests
    ).clip(upper=1.0)

    return (result
            .sort_values("pvalue_bonferroni", ascending=True, na_position="last")
            .reset_index(drop=True))


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"Output directory: {args.output_dir}")
    print(f"Processing k values: {args.k_values}\n")

    for k in args.k_values:

        out_path = os.path.join(args.output_dir, f"kmer_differential_k{k}.tsv")

        try:
            if os.path.exists(out_path):
                print(f"  k = {k:>2d} ... SKIPPED — output already exists: {out_path}")
                continue

            print(f"  k = {k:>2d} ... ", end="", flush=True)

            result = process_k(
                k,
                args.svharb_hotspot_prefix, args.svharb_genome_prefix,
                args.notp_hotspot_prefix,   args.notp_genome_prefix,
            )

            result.to_csv(out_path, sep="\t", index=False, float_format="%.6g")

            n_sig  = int((result["pvalue_bonferroni"] < 0.05).sum())
            n_chi2 = int((result["test_used"] == "chi2").sum())
            n_fish = int((result["test_used"] == "fisher").sum())
            n_skip = int((result["test_used"] == "skipped").sum())
            print(f"{len(result):>6} k-mers  |  "
                  f"sig (Bonf<0.05): {n_sig:>5}  |  "
                  f"chi2: {n_chi2}  fisher: {n_fish}  skipped: {n_skip}")

        except FileNotFoundError as e:
            print(f"SKIPPED — {e}", file=sys.stderr)

    print("\nDone.")


if __name__ == "__main__":
    main()