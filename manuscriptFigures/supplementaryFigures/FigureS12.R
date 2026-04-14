library(jsonlite)
library(purrr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(readxl)

snv_indel_colors <- c(
  "Deepsomatic" = "#FF7F00",
  "Mutect2"     = "#984EA3",
  "Neusomatic"  = "#4DAF4A",
  "Strelka2"    = "#377EB8",
  "Varscan2"    = "#F781BF"
)

sv_colors <- c(
  "pbsv"        = "#1B9E77",
  "Savana"      = "#D95F02",
  "Severus"     = "#7570B3",
  "Sniffles2"    = "#E7298A",
  "SVision-pro" = "#FFFF00"
)


df <- read_excel("Table S1-6_0403.xlsx",sheet = "Table S4",skip=1) %>%
  fill(Variant_Type, Caller, Version, Command, .direction = "down") %>%
  rename(Truvari.rtg.tool.parameter = `Truvari/rtg-tool parameter`) %>%
  select(Variant_Type, Caller, Truvari.rtg.tool.parameter, Recall, Precision, F1) %>%
  filter(
    !is.na(Recall), !is.na(Precision), !is.na(F1),
    !Variant_Type %in% c("MEI", ""),
    Variant_Type != "Variant_Type"  
  ) %>%
  mutate(
    Caller = trimws(Caller),
    Recall = as.numeric(Recall),
    Precision = as.numeric(Precision),
    F1 = as.numeric(F1),
    Row_Group = case_when(
      Variant_Type == "SNV" ~ "SNV",
      Variant_Type == "INDEL" ~ "INDEL",
      Variant_Type == "SV" & grepl("no sequence", Truvari.rtg.tool.parameter) ~ "SV (w/o seq comparison)",
      Variant_Type == "SV" & grepl("sequence comparison", Truvari.rtg.tool.parameter) & !grepl("no sequence", Truvari.rtg.tool.parameter) ~ "SV (w/ seq comparison)"
    ) 
  ) %>%
  filter(!is.na(Row_Group)) %>%
  pivot_longer(cols = c(Recall, Precision, F1), names_to = "Metric", values_to = "Value") %>%
  mutate(
    Row_Group = factor(Row_Group, levels = c("SNV", "INDEL", "SV (w/o seq comparison)", "SV (w/ seq comparison)")),
    Metric = factor(Metric, levels = c("Recall", "Precision", "F1"))
  )


make_plots <- function(sub, color_pal, row_label) {
  metrics <- c("Recall", "Precision", "F1")
  plots <- map(metrics, function(m) {
    d <- sub %>% filter(Metric == m) %>% arrange(desc(Value)) %>%
      mutate(Caller = factor(Caller, levels = Caller))
    
    p <- ggplot(d, aes(x = Caller, y = Value, fill = Caller)) +
      geom_col() +
      scale_fill_manual(values = color_pal) +
      scale_y_continuous(limits = c(0, 1)) +
      theme_bw() +
      theme(text = element_text(family = "Helvetica"))+
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "none",
        strip.text = element_text(size = 9),
        plot.title = element_text(size = 9)
      ) +
      labs(x = NULL, y = NULL, title = m)
    
    
    if (m == "Recall") p <- p + labs(y = row_label)
    p
  })
  wrap_plots(plots, nrow = 1)
}

# SNV
snv_sub <- df %>% filter(Variant_Type == "SNV")
p_snv <- make_plots(snv_sub, snv_indel_colors, "SNV")

# INDEL
indel_sub <- df %>% filter(Variant_Type == "INDEL")
p_indel <- make_plots(indel_sub, snv_indel_colors, "INDEL")

# SV w/ seq
sv_seq_sub <- df %>%
  filter(Variant_Type == "SV", grepl("^sequence comparison", Truvari.rtg.tool.parameter)) %>%
  group_by(Caller, Metric) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")
p_sv_seq <- make_plots(sv_seq_sub, sv_colors, "SV (w/ seq)")

# SV w/o seq
sv_noseq_sub <- df %>%
  filter(Variant_Type == "SV", grepl("no sequence comparison", Truvari.rtg.tool.parameter)) %>%
  group_by(Caller, Metric) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")
p_sv_noseq <- make_plots(sv_noseq_sub, sv_colors, "SV (w/o seq)")

p = p_snv / p_indel / p_sv_seq / p_sv_noseq

library(svglite)
ggsave("FigureS12.svg", plot = p, device = svglite, width = 9, height = 10)

