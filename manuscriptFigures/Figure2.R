library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(reshape2)

########################################## Fig2A ##########################################
# Data
coverage <- c(500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500)
validation_rate <- c(0.9555, 0.9711, 0.9778, 0.9799, 0.9787, 0.9813, 0.9817, 0.9819,0.9819)
data <- data.frame(Coverage = coverage, ValidationRate = validation_rate)

# Plot
ggplot(data, aes(x = Coverage, y = ValidationRate)) +
  geom_line(color = "orange", size = 1) +
  geom_point(color = "orange", size = 2) +
  labs(title = "Variant Validation Rate per Coverage", 
       x = "Coverage (x)", 
       y = "Validation Rate") +
  ylim(0, 1) +
  theme_minimal() +
  theme(text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 20, face = "bold"))

########################################## Fig2B ##########################################

data <- data.frame(
  expectedVAFrange = rep(c("VAF<0.01", "0.01≤VAF<0.02", "0.02≤VAF<0.05", "0.05≤VAF<0.10", "VAF≥0.10"), 2),
  coverage = factor(c(
    rep("4000x", 5),
    rep("500x", 5)
  ), levels = c("500x", "4000x")),
  recall_rate = c(
    # 4000x
    0.9742519717, 0.9560608736, 0.9867675554, 0.9762574652, 0.9909701658,
    # 500x
    0.6548925822, 0.9341421427, 0.9808647101, 0.971674132, 0.9882678066
  ),
  variant_count = c(
    # SNVs counts
    372145, 2203321, 679164, 1863491, 424814,
    # Same counts for 500x
    372145, 2203321, 679164, 1863491, 424814
  )
)

# Set the correct factor levels to ensure proper ordering
data$expectedVAFrange <- factor(data$expectedVAFrange,
                                levels = c("VAF<0.01", "0.01≤VAF<0.02", "0.02≤VAF<0.05", "0.05≤VAF<0.10", "VAF≥0.10"))

# Create an x-axis with sample counts included for all points
data$x_label <- paste0(data$expectedVAFrange, "\n(n=", data$variant_count, ")")

data$x_label <- factor(data$x_label, 
                       levels = unique(data$x_label[order(data$expectedVAFrange)]))

# Create the line plot
vaf_plot <- ggplot(data, aes(x = x_label, y = recall_rate, 
                             group = coverage, color = coverage)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  # Add text labels for recall rates
  geom_text(aes(label = sprintf("%.4f", recall_rate), 
                vjust = ifelse(coverage == "4000x", -1.8, 2.5)), 
            hjust = 0.5, 
            size = 4) +
  # Set specific colors for 4000x and 500x
  scale_color_manual(values = c("500x" = "#0077BB", "4000x" = "#EE3377")) +
  # Y-axis settings
  scale_y_continuous(
    limits = c(0.6, 1.0),  # Convert to decimal range
    breaks = seq(0.6, 1.0, 0.1),  # Convert breaks to decimals
    labels = function(x) sprintf("%.1f", x),  # Format as decimal
    expand = expansion(mult = c(0, 0.1))
  ) +
  # Labels
  labs(
    title = "SNV Validation Rate per Benchmark Set VAF",
    x = "Expected VAF",
    y = "Recall Rate",
    color = "Coverage"
  ) +
  # Custom theme settings (maintained from original)
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 13),
    plot.title = element_text(size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )


########################################## Fig2C ##########################################

data <- data.frame(
  Variant_Type = c("SNVs (N=6,037,703)", "Indels (N=1,832,989)", "SVs (N=51,006)"),
  Alignment_Recall = c(0.9325612737, 0.8669380995, 0.9372230718),
  Percent_Rescued = c(0.007846030187, 0.1011609999, 0.02235031173)
)

# IMPORTANT FIX: The factor levels in your code don't match the actual data values
# Set the factor levels to match the actual values in your data
data$Variant_Type <- factor(data$Variant_Type, 
                            levels = c("SNVs (N=6,037,703)", 
                                       "Indels (N=1,832,989)", 
                                       "SVs (N=51,006)"))

# Calculate total values for each variant type
data$Total <- data$Alignment_Recall + data$Percent_Rescued

# Reshape the data to long format
melted_data <- melt(data, id.vars = c("Variant_Type", "Total"), 
                    variable.name = "Metric",
                    value.name = "Value")

# Reorder the levels of Metric to match the stacking order
melted_data$Metric <- factor(melted_data$Metric, 
                             levels = c("Percent_Rescued", 
                                        "Alignment_Recall"))

# Create the plot
validation_plot <- ggplot(melted_data, aes(x = Variant_Type, y = Value)) +
  geom_bar(aes(fill = Metric), stat = "identity", width = 0.7) +
  # Add text annotations for totals
  geom_text(data = data,  # Use original data frame
            aes(y = Total, label = sprintf("%.1f%%", Total * 100)),
            vjust = -0.5,
            size = 4) +
  geom_text(data = data,
            aes(y = Alignment_Recall/2, 
                label = sprintf("%.1f%%", Alignment_Recall * 100)),
            size = 4,
            color = "white") +
  scale_fill_manual(values = c("Percent_Rescued" = "#EE3377",       
                               "Alignment_Recall" = "#0077bb"),    
                    labels = c("Assembly-based Validation",
                               "Alignment-based Validation")) +
  scale_y_continuous(limits = c(0, 1.05),  # Increased upper limit to accommodate labels
                     breaks = seq(0, 1, 0.25),
                     labels = scales::percent,
                     expand = c(0, 0)) +
  labs(title = "Recall Across Variant Types",
       x = "",
       y = "Recall",
       fill = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 13),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.position = "top",
    legend.justification = "center",
    panel.grid = element_blank(),  # Remove all gridlines
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

########################################## Fig2D ##########################################

########################################## Fig2E ##########################################
# In separate Python scripts under the current directory.

