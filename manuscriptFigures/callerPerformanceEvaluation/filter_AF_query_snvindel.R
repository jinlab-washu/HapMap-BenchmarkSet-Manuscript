library(optparse)

options(scipen = 999)

option_list <- list(
  make_option(c("-f", "--vcf_file"),
              type = "character",
              help = "Path to vcf txt file"),
  make_option(c("-t", "--threads"),
              type = "integer",
              help = "number of threads"),
  make_option(c("-c", "--class"),
              type = "character",
              help = "variants class"),
  make_option(c("-n", "--name"),
              type = "character",
              help = "variants name"),
  make_option(c("-o", "--outdir"),
              type = "character",
              help = "output file directory, don't plus /!!!")
)

opt <- parse_args(OptionParser(option_list = option_list))

library(data.table)

data.table::setDTthreads(threads = opt$threads)

vcf_txt_path <- opt$vcf_file
vcf_txt <- fread(vcf_txt_path, header = FALSE, sep = "\t")

ncol_vcf <- ncol(vcf_txt)

if (ncol_vcf == 10) {
  # one-sample VCF
  vcf_txt_header <- c("chr", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", "sample")
  setnames(vcf_txt, vcf_txt_header)
  keep_cols <- 1:10
} else if (ncol_vcf == 11) {
  # two-sample VCF
  vcf_txt_header <- c("chr", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", "normal", "sample")
  setnames(vcf_txt, vcf_txt_header)
  keep_cols <- 1:11
} else {
  stop(paste("Unexpected number of columns in VCF txt:", ncol_vcf,
             ". Expected 10 (1 sample) or 11 (2 samples)."))
}

vcf_txt_new <- copy(vcf_txt)

extract_AF <- function(format_str, sample_str) {
  format_parts <- strsplit(format_str, ":", fixed = TRUE)[[1]]
  sample_parts <- strsplit(sample_str, ":", fixed = TRUE)[[1]]
  
  af_index <- which(format_parts %in% c("AF", "VAF"))
  
  if (length(af_index) == 0) {
    return(NA_real_)
  }
  
  af_index <- af_index[1]
  
  if (af_index > length(sample_parts)) {
    return(NA_real_)
  }
  
  af_value <- sample_parts[af_index]
  
  if (af_value %in% c(".", "", NA)) {
    return(NA_real_)
  }
  
  # if multiallelic like 0.12,0.03, take first value
  af_value <- strsplit(af_value, ",", fixed = TRUE)[[1]][1]
  
  suppressWarnings(as.numeric(af_value))
}

vcf_txt_new[, AF := mapply(extract_AF, format, sample)]
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

out_dir <- opt$outdir
name <- opt$name
class <- opt$class

vcf_txt_new_output_path <- vector("character", length = 9)
for (i in 1:9) {
  vcf_txt_new_output_path[i] <- paste(out_dir, "/", name, "_", class, "_", i, ".txt", sep = "")
}

for (i in 1:9) {
  fwrite(vcf_txt_new_list[[i]],
         vcf_txt_new_output_path[i],
         sep = "\t",
         col.names = FALSE,
         row.names = FALSE)
}