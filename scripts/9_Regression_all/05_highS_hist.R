##################################################################################
## Select highest slopes
## Plot along chromosomes
## 
## Author Daniel Anstett
## 
## 
## Last Modified 07072025
###################################################################################
#Import libraries
library(tidyverse)

#Import top plink clumped 215 SNPs per pop
p1_215 <- read_csv("data/genomic_data/p01_215.csv")
p2_215 <- read_csv("data/genomic_data/p02_215.csv")
p3_215 <- read_csv("data/genomic_data/p03_215.csv")
p4_215 <- read_csv("data/genomic_data/p04_215.csv")
p5_215 <- read_csv("data/genomic_data/p05_215.csv")
p6_215 <- read_csv("data/genomic_data/p06_215.csv")
p7_215 <- read_csv("data/genomic_data/p07_215.csv")
p8_215 <- read_csv("data/genomic_data/p08_215.csv")
p9_215 <- read_csv("data/genomic_data/p09_215.csv")
p10_215 <- read_csv("data/genomic_data/p10_215.csv")
p11_215 <- read_csv("data/genomic_data/p11_215.csv")


#Make plot for pop 1 for top 215 SNPs
ggplot(p1_215, aes(x = Position, y = Slope_abs)) +
  geom_point(shape = 1, size = 2, stroke = 0.9, color = "black", alpha = 0.7) +
  facet_wrap(~ chr, scales = "free_x", nrow = 2, ncol = 4) +
  theme_classic() +
  labs(
    title = "Top 215 SNPs for Site 1 by Chromosome",
    x = "Genomic Position",
    y = "Slope"
  )
#ggsave("Graphs/Regression_all/p1_215_plot.pdf", width = 10,height = 6,units = "in")


# For all pop for 215
for (i in 1:11) {
  # Dynamically get the data frame (e.g. p1_215, p2_215, ...)
  df_name <- paste0("p", i, "_215")
  df <- get(df_name)
  
  # Create plot
  p <- ggplot(df, aes(x = Position, y = Slope_abs)) +
    geom_point(shape = 1, size = 2, stroke = 0.9, color = "black", alpha = 0.7) +
    facet_wrap(~ chr, scales = "free_x", nrow = 2, ncol = 4) +
    theme_classic() +
    labs(
      title = paste("Top 215 SNPs for Site", i, "by Chromosome"),
      x = "Genomic Position",
      y = "Slope"
    )
  
  # Save plot
  ggsave(
    filename = paste0("p", i, "_215_plot.pdf"),
    plot = p,
    path = "Graphs/Regression_all/215",
    width = 10,
    height = 6,
    units = "in"
  )
}


