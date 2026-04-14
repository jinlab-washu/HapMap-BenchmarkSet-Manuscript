library(tidyverse)
library(jsonlite)

base_dir <- "/Volumes/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation"

# --------------------------------------------------
# Get ONLY k-mer summary files
# --------------------------------------------------
files <- Sys.glob(file.path(base_dir, "*/kmer/*/summary*"))

# --------------------------------------------------
# Filter unwanted files
# --------------------------------------------------
files <- files %>%
  keep(~ !str_detect(.x, "/mei/kmer/")) %>%
  keep(~ !str_detect(.x, "wseq")) %>%
  keep(~ !str_detect(.x, "_Test/summary"))

# --------------------------------------------------
# Parse TXT summary (SNV / INDEL)
# --------------------------------------------------
parse_txt <- function(fp) {
  lines <- readLines(fp)
  
  data_lines <- lines[str_detect(lines, "^(\\s*None|\\s*[0-9])")]
  if (length(data_lines) == 0) return(NULL)
  
  row <- data_lines[length(data_lines)]
  parts <- str_split(str_squish(row), "\\s+")[[1]]
  
  tibble(
    Precision = as.numeric(parts[6]),
    Recall    = as.numeric(parts[7]),
    F1        = as.numeric(parts[8])
  )
}

# --------------------------------------------------
# Parse JSON summary (SV / MEI)
# --------------------------------------------------
parse_json <- function(fp) {
  d <- fromJSON(fp)
  
  tibble(
    Precision = d$precision,
    Recall    = d$recall,
    F1        = d$f1
  )
}

# --------------------------------------------------
# Extract metadata
# --------------------------------------------------
parse_meta <- function(fp) {
  folder <- basename(dirname(fp))
  
  vartype <- case_when(
    str_detect(fp, "/snv/")    ~ "SNV",
    str_detect(fp, "/indels/") ~ "Indel",
    str_detect(fp, "/sv/")     ~ "SV",
    TRUE ~ NA_character_
  )
  
  genome <- case_when(
    str_detect(folder, "GRCh38|hg38") ~ "GRCh38",
    str_detect(folder, "CHM13|chm13") ~ "CHM13",
    TRUE ~ NA_character_
  )
  
  region <- case_when(
    str_detect(folder, "_total$")  ~ "total",
    str_detect(folder, "kmer_k24\\.umap$")  ~ "k24",
    str_detect(folder, "kmer_k36\\.umap$")  ~ "k36",
    str_detect(folder, "kmer_k50\\.umap$")  ~ "k50",
    str_detect(folder, "kmer_k100\\.umap$") ~ "k100",
    genome == "GRCh38" & str_detect(folder, "kmer_all$") ~ "remaining",
    genome == "CHM13"  & str_detect(folder, "kmer_all_wo_chm13only$") ~ "remaining",
    genome == "CHM13"  & str_detect(folder, "kmer_all_w_chm13only$") ~ "remaining+chm13only",
    TRUE ~ NA_character_
  )
  
  tokens <- str_split(folder, "_")[[1]]
  
  caller <- case_when(
    length(tokens) >= 2 ~ tokens[2],
    TRUE ~ NA_character_
  )
  
  tibble(
    VarType = vartype,
    Genome  = genome,
    Region  = region,
    Caller  = caller,
    File    = fp
  )
}

# --------------------------------------------------
# Parse all files
# --------------------------------------------------
df <- map_dfr(files, function(fp) {
  meta <- parse_meta(fp)
  
  metrics <- if (str_detect(fp, "\\.txt$")) {
    parse_txt(fp)
  } else {
    parse_json(fp)
  }
  
  if (is.null(metrics)) return(NULL)
  bind_cols(meta, metrics)
})

df <- df %>%
  filter(!is.na(VarType), !is.na(Genome), !is.na(Region))

# --------------------------------------------------
# Final mean F1 table
# --------------------------------------------------
region_levels <- c("total","k24", "k36", "k50", "k100", "remaining", "remaining+chm13only")

paper_table <- df %>%
  filter(VarType %in% c("SNV", "Indel", "SV")) %>%
  group_by(Region, Genome, VarType) %>%
  summarise(F1 = mean(F1, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = VarType, values_from = F1) %>%
  mutate(Region = factor(Region, levels = region_levels)) %>%
  arrange(Region, Genome)

print(paper_table)

# --------------------------------------------------
# Plot mean F1 table
# --------------------------------------------------
plot_df <- paper_table %>%
  pivot_longer(
    cols = c(SNV, Indel, SV),
    names_to = "VarType",
    values_to = "MeanF1"
  ) %>%
  mutate(
    Region = factor(Region, levels = region_levels),
    Region_label = recode(as.character(Region),
                          "total" = "Total",
                          "k24" = "k24",
                          "k36" = "k36",
                          "k50" = "k50",
                          "k100" = "k100",
                          "remaining" = "Remaining",
                          "remaining+chm13only" = "Remaining + CHM13-only"
    ),
    Region_label = factor(
      Region_label,
      levels = c("Total", "k24", "k36", "k50", "k100", "Remaining", "Remaining + CHM13-only")
    ),
    Genome = factor(Genome, levels = c("GRCh38", "CHM13")),
    VarType = factor(VarType, levels = c("SNV", "Indel", "SV"))
  )

p <- ggplot(plot_df, aes(x = Genome, y = MeanF1, fill = Genome)) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = sprintf("%.3f", MeanF1)),
    vjust = -0.25,
    size = 3,
    family = "Helvetica"
  ) +
  facet_grid(Region_label ~ VarType, switch = "y", drop = FALSE) +
  scale_fill_manual(values = c("GRCh38" = "#F8766D", "CHM13" = "#00BFC4")) +
  scale_y_continuous(limits = c(0, 1.05), expand = c(0, 0)) +
  labs(
    title = "Average F1 Across Callers by Minimum Unique K-mer",
    x = NULL,
    y = "Mean F1"
  ) +
  theme_bw(base_family = "Helvetica") +
  theme(
    text = element_text(family = "Helvetica"),
    strip.background = element_rect(fill = "white"),
    strip.placement = "outside",
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor.y = element_line(color = "grey92"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p)

ggsave(
   "/Users/nanakong/Library/CloudStorage/Box-Box/Genetics_Writing_Group/SMaHT_Truth_Set/2026_Rebuttal/Revision/ModifiedFigures/FigureS13.pdf",
   p
 )
