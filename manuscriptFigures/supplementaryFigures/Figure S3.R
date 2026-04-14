library(ggplot2)
suppressMessages(library(tidyverse))
install.packages("svglite")
library(dplyr)
################Reliable and Unreliable Counts per Sample#########
# Read the data from the text file
data <- read.table("reli_unreli_length_final.txt",  header = FALSE, col.names = c("Sample", "Type", "Count"))

# Transform, reorder, and rename the data
data <- data %>%
  mutate(Sample = as.factor(Sample),
         Sample = fct_recode(Sample,
                             "common_chm13" = "6hapmap",
                             "common_liftover_hg38" = "6hapmap_liftover"),
         Sample = fct_relevel(Sample, "common_chm13", "common_liftover_hg38", after = Inf),
         Type = factor(Type, levels = c("reliable", "unreliable"))) %>%
  arrange(desc(Type))

# Create the stacked bar plot
plot <- ggplot(data, aes(x = Sample, y = Count, fill = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("reliable" = "#00AF92", "unreliable" = "#DF536B")) +
  geom_vline(xintercept = which(levels(data$Sample) == "common") - 0.5, linetype = "dashed", color = "black") +
  labs(title = "Size of Reliable and Unreliable Regions per Sample",
       x = "Sample",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # Increase x-axis label size and rotate for visibility
        axis.text.y = element_text(size = 12), # Increase y-axis label size
        axis.title.x = element_text(size = 14), # Increase x-axis title size
        axis.title.y = element_text(size = 14), # Increase y-axis title size
        plot.title = element_text(size = 16, hjust = 0.5), # Increase and center the plot title
        legend.title = element_text(size = 12), # Increase legend title size
        legend.text = element_text(size = 10)) # Increase legend text size

print(plot)
