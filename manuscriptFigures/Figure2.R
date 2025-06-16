library(ggplot2)

#Fig2A
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

#Fig2B

#Fig2C

#Fig2D


