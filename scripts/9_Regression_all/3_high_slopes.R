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

#Import
snp_all <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_all.csv") %>% filter(SE<5)

#Get top 212 SNPs per pop
p1_212 <- snp_all %>% filter(Site==1) %>% arrange(desc(Slope)) %>% slice_head(n = 212)
p2_212 <- snp_all %>% filter(Site == 2) %>% arrange(desc(Slope)) %>% slice_head(n = 212)
p3_212 <- snp_all %>% filter(Site == 3) %>% arrange(desc(Slope)) %>% slice_head(n = 212)
p4_212 <- snp_all %>% filter(Site == 4) %>% arrange(desc(Slope)) %>% slice_head(n = 212)
p5_212 <- snp_all %>% filter(Site == 5) %>% arrange(desc(Slope)) %>% slice_head(n = 212)
p6_212 <- snp_all %>% filter(Site == 6) %>% arrange(desc(Slope)) %>% slice_head(n = 212)
p7_212 <- snp_all %>% filter(Site == 7) %>% arrange(desc(Slope)) %>% slice_head(n = 212)
p8_212 <- snp_all %>% filter(Site == 8) %>% arrange(desc(Slope)) %>% slice_head(n = 212)
p9_212 <- snp_all %>% filter(Site == 9) %>% arrange(desc(Slope)) %>% slice_head(n = 212)
p10_212 <- snp_all %>% filter(Site == 10) %>% arrange(desc(Slope)) %>% slice_head(n = 212)
p11_212 <- snp_all %>% filter(Site == 11) %>% arrange(desc(Slope)) %>% slice_head(n = 212)

#Get top 2000 SNPs per pop
p1_2000 <- snp_all %>% filter(Site == 1) %>% arrange(desc(Slope)) %>% slice_head(n = 2000)
p2_2000 <- snp_all %>% filter(Site == 2) %>% arrange(desc(Slope)) %>% slice_head(n = 2000)
p3_2000 <- snp_all %>% filter(Site == 3) %>% arrange(desc(Slope)) %>% slice_head(n = 2000)
p4_2000 <- snp_all %>% filter(Site == 4) %>% arrange(desc(Slope)) %>% slice_head(n = 2000)
p5_2000 <- snp_all %>% filter(Site == 5) %>% arrange(desc(Slope)) %>% slice_head(n = 2000)
p6_2000 <- snp_all %>% filter(Site == 6) %>% arrange(desc(Slope)) %>% slice_head(n = 2000)
p7_2000 <- snp_all %>% filter(Site == 7) %>% arrange(desc(Slope)) %>% slice_head(n = 2000)
p8_2000 <- snp_all %>% filter(Site == 8) %>% arrange(desc(Slope)) %>% slice_head(n = 2000)
p9_2000 <- snp_all %>% filter(Site == 9) %>% arrange(desc(Slope)) %>% slice_head(n = 2000)
p10_2000 <- snp_all %>% filter(Site == 10) %>% arrange(desc(Slope)) %>% slice_head(n = 2000)
p11_2000 <- snp_all %>% filter(Site == 11) %>% arrange(desc(Slope)) %>% slice_head(n = 2000)


#Make plot for pop 1 for top 212 SNPs
ggplot(p1_212, aes(x = Position, y = Slope)) +
  geom_point(shape = 1, size = 2, stroke = 0.9, color = "black", alpha = 0.7) +
  facet_wrap(~ chr, scales = "free_x", nrow = 2, ncol = 4) +
  theme_classic() +
  labs(
    title = "Top 212 SNPs for Site 1 by Chromosome",
    x = "Genomic Position",
    y = "Slope"
  )
ggsave("Graphs/Regression_all/p1_212_plot.pdf", width = 10,height = 6,units = "in")


# For all pop for 212
for (i in 1:11) {
  # Dynamically get the data frame (e.g. p1_212, p2_212, ...)
  df_name <- paste0("p", i, "_212")
  df <- get(df_name)
  
  # Create plot
  p <- ggplot(df, aes(x = Position, y = Slope)) +
    geom_point(shape = 1, size = 2, stroke = 0.9, color = "black", alpha = 0.7) +
    facet_wrap(~ chr, scales = "free_x", nrow = 2, ncol = 4) +
    theme_classic() +
    labs(
      title = paste("Top 212 SNPs for Site", i, "by Chromosome"),
      x = "Genomic Position",
      y = "Slope"
    )
  
  # Save plot
  ggsave(
    filename = paste0("p", i, "_212_plot.pdf"),
    plot = p,
    path = "Graphs/Regression_all/212",
    width = 10,
    height = 6,
    units = "in"
  )
}


# For all pop for 2000
for (i in 1:11) {
  # Dynamically get the data frame (e.g. p1_2000, p2_2000, ...)
  df_name <- paste0("p", i, "_2000")
  df <- get(df_name)
  
  # Create plot
  p <- ggplot(df, aes(x = Position, y = Slope)) +
    geom_point(shape = 1, size = 2, stroke = 0.9, color = "black", alpha = 0.7) +
    facet_wrap(~ chr, scales = "free_x", nrow = 2, ncol = 4) +
    theme_classic() +
    labs(
      title = paste("Top 200 SNPs for Site", i, "by Chromosome"),
      x = "Genomic Position",
      y = "Slope"
    )
  
  # Save plot
  ggsave(
    filename = paste0("p", i, "_2000_plot.pdf"),
    plot = p,
    path = "Graphs/Regression_all/2000",
    width = 10,
    height = 6,
    units = "in"
  )
}
