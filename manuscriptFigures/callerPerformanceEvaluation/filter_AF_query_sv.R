library(optparse)
library(data.table)

options(scipen = 999)

option_list <- list(
  make_option(c("-f", "--vcf_file"),
              type = "character",
              help = "Path to vcf txt file"),
  make_option(c("-t", "--threads"),
              type = "integer",
              help = "number of threads"),
  make_option(c("-n", "--name"),
              type = "character",
              help = "caller name"),
  make_option(c("-c", "--class"),
              type = "character",
              help = "variants class"),
  make_option(c("-o", "--outdir"),
              type = "character",
              help = "output file directory, don't plus /!!!")
)

opt <- parse_args(OptionParser(option_list = option_list))

out_dir <- opt$outdir
name <- opt$name
class <- opt$class

data.table::setDTthreads(threads = opt$threads)

vcf_txt_path <- opt$vcf_file
vcf_txt <- fread(vcf_txt_path, header = FALSE, sep = "\t")

ncol_vcf <- ncol(vcf_txt)

if (ncol_vcf == 10) {
  # one-sample VCF
  setnames(vcf_txt, c("chr", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", "sample"))
  keep_cols <- 1:10
} else if (ncol_vcf == 11) {
  # two-sample VCF; for svision, use first sample column
  setnames(vcf_txt, c("chr", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", "sample1", "sample2"))
  vcf_txt[, sample := sample1]
  keep_cols <- 1:11
} else {
  stop(paste("Unexpected number of columns:", ncol_vcf, "Expected 10 or 11"))
}

vcf_txt_new <- copy(vcf_txt)

extract_format_value <- function(format_str, sample_str, keys = c("AF", "VAF")) {
  format_parts <- strsplit(format_str, ":", fixed = TRUE)[[1]]
  sample_parts <- strsplit(sample_str, ":", fixed = TRUE)[[1]]
  
  idx <- which(format_parts %in% keys)
  if (length(idx) == 0) return(NA_real_)
  
  idx <- idx[1]
  if (idx > length(sample_parts)) return(NA_real_)
  
  val <- sample_parts[idx]
  if (is.na(val) || val %in% c(".", "", "NA")) return(NA_real_)
  
  val <- strsplit(val, ",", fixed = TRUE)[[1]][1]
  suppressWarnings(as.numeric(val))
}

extract_info_value <- function(info_str, key) {
  parts <- strsplit(info_str, ";", fixed = TRUE)[[1]]
  hit <- parts[startsWith(parts, paste0(key, "="))]
  if (length(hit) == 0) return(NA_real_)
  
  val <- sub(paste0("^", key, "="), "", hit[1])
  if (is.na(val) || val %in% c(".", "", "NA")) return(NA_real_)
  
  val <- strsplit(val, ",", fixed = TRUE)[[1]][1]
  suppressWarnings(as.numeric(val))
}

extract_pbsv_AF <- function(format_str, sample_str) {
  format_parts <- strsplit(format_str, ":", fixed = TRUE)[[1]]
  sample_parts <- strsplit(sample_str, ":", fixed = TRUE)[[1]]
  
  ad_idx <- which(format_parts == "AD")
  dp_idx <- which(format_parts == "DP")
  
  if (length(ad_idx) == 0 || length(dp_idx) == 0) return(NA_real_)
  ad_idx <- ad_idx[1]
  dp_idx <- dp_idx[1]
  
  if (ad_idx > length(sample_parts) || dp_idx > length(sample_parts)) return(NA_real_)
  
  ad_val <- sample_parts[ad_idx]
  dp_val <- sample_parts[dp_idx]
  
  if (ad_val %in% c(".", "", "NA") || dp_val %in% c(".", "", "NA")) return(NA_real_)
  
  ad_parts <- strsplit(ad_val, ",", fixed = TRUE)[[1]]
  if (length(ad_parts) < 2) return(NA_real_)
  
  alt_count <- suppressWarnings(as.numeric(ad_parts[2]))
  dp <- suppressWarnings(as.numeric(dp_val))
  
  if (is.na(alt_count) || is.na(dp) || dp == 0) return(NA_real_)
  
  alt_count / dp
}

extract_dr_dv_AF <- function(format_str, sample_str) {
  format_parts <- strsplit(format_str, ":", fixed = TRUE)[[1]]
  sample_parts <- strsplit(sample_str, ":", fixed = TRUE)[[1]]
  
  dr_idx <- which(format_parts == "DR")
  dv_idx <- which(format_parts == "DV")
  
  if (length(dr_idx) == 0 || length(dv_idx) == 0) return(NA_real_)
  dr_idx <- dr_idx[1]
  dv_idx <- dv_idx[1]
  
  if (dr_idx > length(sample_parts) || dv_idx > length(sample_parts)) return(NA_real_)
  
  dr <- suppressWarnings(as.numeric(sample_parts[dr_idx]))
  dv <- suppressWarnings(as.numeric(sample_parts[dv_idx]))
  
  if (is.na(dr) || is.na(dv) || (dr + dv) == 0) return(NA_real_)
  
  dv / (dr + dv)
}

if (grepl("^sniffles", name, ignore.case = TRUE)) {
  message("Computing AF for Sniffles from FORMAT DR/DV")
  vcf_txt_new[, AF := mapply(extract_dr_dv_AF, format, sample)]
  
} else if (grepl("^pbsv", name, ignore.case = TRUE)) {
  message("Computing AF for pbsv from FORMAT AD/DP")
  vcf_txt_new[, AF := mapply(extract_pbsv_AF, format, sample)]
  
} else if (grepl("^severus", name, ignore.case = TRUE)) {
  message("Computing AF for Severus from FORMAT VAF")
  vcf_txt_new[, AF := mapply(function(f, s) extract_format_value(f, s, keys = c("VAF")), format, sample)]
  
} else if (grepl("^savana", name, ignore.case = TRUE)) {
  message("Computing AF for Savana from INFO TUMOUR_AF")
  vcf_txt_new[, AF := sapply(info, extract_info_value, key = "TUMOUR_AF")]
  
} else if (grepl("^svision", name, ignore.case = TRUE)) {
  message("Computing AF for svision from first-sample FORMAT DR/DV")
  vcf_txt_new[, AF := mapply(extract_dr_dv_AF, format, sample)]
  
} else {
  stop(paste("Unsupported caller name:", name))
}

if (!"AF" %in% colnames(vcf_txt_new)) {
  stop("AF column was not created.")
}

vcf_txt_new <- vcf_txt_new[!is.na(AF)]

vcf_txt_new_1 <- vcf_txt_new[AF <= 0.005, .SD, .SDcols = keep_cols]
vcf_txt_new_2 <- vcf_txt_new[AF > 0.005 & AF <= 0.01, .SD, .SDcols = keep_cols]
vcf_txt_new_3 <- vcf_txt_new[AF > 0.01 & AF <= 0.015, .SD, .SDcols = keep_cols]
vcf_txt_new_4 <- vcf_txt_new[AF > 0.015 & AF <= 0.02, .SD, .SDcols = keep_cols]
vcf_txt_new_5 <- vcf_txt_new[AF > 0.02 & AF <= 0.03, .SD, .SDcols = keep_cols]
vcf_txt_new_6 <- vcf_txt_new[AF > 0.03 & AF <= 0.04, .SD, .SDcols = keep_cols]
vcf_txt_new_7 <- vcf_txt_new[AF > 0.04 & AF <= 0.05, .SD, .SDcols = keep_cols]
vcf_txt_new_8 <- vcf_txt_new[AF > 0.05 & AF <= 0.1, .SD, .SDcols = keep_cols]
vcf_txt_new_9 <- vcf_txt_new[AF > 0.1 & AF <= 0.165, .SD, .SDcols = keep_cols]

vcf_txt_new_list <- list(
  vcf_txt_new_1, vcf_txt_new_2, vcf_txt_new_3, vcf_txt_new_4,
  vcf_txt_new_5, vcf_txt_new_6, vcf_txt_new_7, vcf_txt_new_8, vcf_txt_new_9
)

vcf_txt_new_output_path <- vector("character", length = 9)
for (i in 1:9) {
  vcf_txt_new_output_path[i] <- paste0(out_dir, "/", name, "_", class, "_", i, ".txt")
}

for (i in 1:9) {
  fwrite(vcf_txt_new_list[[i]],
         vcf_txt_new_output_path[i],
         sep = "\t",
         col.names = FALSE,
         row.names = FALSE)
}