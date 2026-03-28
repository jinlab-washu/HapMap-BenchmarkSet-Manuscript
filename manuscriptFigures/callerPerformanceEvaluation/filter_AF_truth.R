#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

options(scipen = 999)

option_list <- list(
  make_option(
    c("-f", "--vcf_file"),
    type = "character",
    help = "Path to VCF body text file"
  ),
  make_option(
    c("-t", "--threads"),
    type = "integer",
    help = "Number of threads for data.table"
  ),
  make_option(
    c("-r", "--reference"),
    type = "character",
    help = "Genome reference name (e.g. GRCh38, CHM13)"
  ),
  make_option(
    c("-c", "--class"),
    type = "character",
    help = "Variant class (snvs, indels, svs)"
  ),
  make_option(
    c("-v", "--version"),
    type = "character",
    help = "Truth set version"
  ),
  make_option(
    c("-o", "--outdir"),
    type = "character",
    help = "Output directory (do not include trailing slash)"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

data.table::setDTthreads(threads = opt$threads)

# Read VCF body text (no header lines, tab-delimited)
vcf_txt <- fread(opt$vcf_file, header = FALSE, sep = "\t")

vcf_txt_header <- c(
  "chr", "pos", "id", "ref", "alt",
  "qual", "filter", "info", "format", "sample"
)
setnames(vcf_txt, vcf_txt_header)

extract_AF <- function(x) {
  parts <- strsplit(x, ":", fixed = TRUE)[[1]]
  if (length(parts) < 2) return(NA_character_)
  parts[2]
}

# Extract AF from SAMPLE column assuming GT:AF-style format
vcf_txt[, AF := sapply(sample, extract_AF)]
vcf_txt <- vcf_txt[AF != "."]
vcf_txt[, AF := as.numeric(AF)]
vcf_txt <- vcf_txt[!is.na(AF)]

# Bin definitions
vcf_txt_new_1 <- vcf_txt[AF <= 0.005, .SD, .SDcols = 1:10]
vcf_txt_new_2 <- vcf_txt[AF > 0.005 & AF <= 0.01, .SD, .SDcols = 1:10]
vcf_txt_new_3 <- vcf_txt[AF > 0.01 & AF <= 0.015, .SD, .SDcols = 1:10]
vcf_txt_new_4 <- vcf_txt[AF > 0.015 & AF <= 0.02, .SD, .SDcols = 1:10]
vcf_txt_new_5 <- vcf_txt[AF > 0.02 & AF <= 0.03, .SD, .SDcols = 1:10]
vcf_txt_new_6 <- vcf_txt[AF > 0.03 & AF <= 0.04, .SD, .SDcols = 1:10]
vcf_txt_new_7 <- vcf_txt[AF > 0.04 & AF <= 0.05, .SD, .SDcols = 1:10]
vcf_txt_new_8 <- vcf_txt[AF > 0.05 & AF <= 0.1, .SD, .SDcols = 1:10]
vcf_txt_new_9 <- vcf_txt[AF > 0.1 & AF <= 0.165, .SD, .SDcols = 1:10]

vcf_txt_new_list <- list(
  vcf_txt_new_1, vcf_txt_new_2, vcf_txt_new_3,
  vcf_txt_new_4, vcf_txt_new_5, vcf_txt_new_6,
  vcf_txt_new_7, vcf_txt_new_8, vcf_txt_new_9
)

out_dir <- opt$outdir
ref <- opt$reference
variant_class <- opt$class
version <- opt$version

for (i in 1:9) {
  out_path <- paste0(
    out_dir, "/truth_set_", version, "_",
    variant_class, "_", ref, "_", i, ".txt"
  )
  fwrite(
    vcf_txt_new_list[[i]],
    out_path,
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE
  )
}
