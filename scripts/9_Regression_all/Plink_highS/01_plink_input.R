#############################################################################################################
## Calculate Major and Minor Allele for High S SNP set
## Author Daniel Anstett
## 
## 
## 07112015
#############################################################################################################
#Import libraries
library(tidyverse)

#Import data
snp_all <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_all.csv") %>% filter(SE<5)

#Get top 2000 SNPs per pop
p01_2000 <- snp_all %>% filter(Site == 1) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p02_2000 <- snp_all %>% filter(Site == 2) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p03_2000 <- snp_all %>% filter(Site == 3) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p04_2000 <- snp_all %>% filter(Site == 4) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p05_2000 <- snp_all %>% filter(Site == 5) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p06_2000 <- snp_all %>% filter(Site == 6) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p07_2000 <- snp_all %>% filter(Site == 7) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p08_2000 <- snp_all %>% filter(Site == 8) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p09_2000 <- snp_all %>% filter(Site == 9) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p10_2000 <- snp_all %>% filter(Site == 10) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p11_2000 <- snp_all %>% filter(Site == 11) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
rm(snp_all)

#Import Baseline snp table to get 55 pop abudnance data
pop_order_base<-read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", 
                           header=F, sep="\t")
snp_base<-read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table", 
                     header=F, sep=" ")
loci_base<-read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", 
                      header=F, sep="\t")
colnames(loci_base) <- c("Chromosome","SNP")
loci_united_base <- loci_base %>% unite(chr_snp,"Chromosome","SNP",sep="_")
loci_snp <-cbind(loci_united_base,snp_base) #add snp lables to rows


#############################################################################################################

#Left join to get 55 pop abundance data only for strong SNP set
abundance_p01 <- p01_2000 %>% left_join(loci_snp, by = c("snp_ID" = "chr_snp"))
abundance_p02 <- p02_2000 %>% left_join(loci_snp, by = c("snp_ID" = "chr_snp"))
abundance_p03 <- p03_2000 %>% left_join(loci_snp, by = c("snp_ID" = "chr_snp"))
abundance_p04 <- p04_2000 %>% left_join(loci_snp, by = c("snp_ID" = "chr_snp"))
abundance_p05 <- p05_2000 %>% left_join(loci_snp, by = c("snp_ID" = "chr_snp"))
abundance_p06 <- p06_2000 %>% left_join(loci_snp, by = c("snp_ID" = "chr_snp"))
abundance_p07 <- p07_2000 %>% left_join(loci_snp, by = c("snp_ID" = "chr_snp"))
abundance_p08 <- p08_2000 %>% left_join(loci_snp, by = c("snp_ID" = "chr_snp"))
abundance_p09 <- p09_2000 %>% left_join(loci_snp, by = c("snp_ID" = "chr_snp"))
abundance_p10 <- p10_2000 %>% left_join(loci_snp, by = c("snp_ID" = "chr_snp"))
abundance_p11 <- p11_2000 %>% left_join(loci_snp, by = c("snp_ID" = "chr_snp"))

#Calc empirical p-value based on BF rank of all SNP
abundance_p01 <- abundance_p01 %>% mutate(Empirical_p = 1 - rank(abundance_p01$Slope_abs) / length(abundance_p01$Slope_abs))
abundance_p02 <- abundance_p02 %>% mutate(Empirical_p = 1 - rank(abundance_p02$Slope_abs) / length(abundance_p02$Slope_abs))
abundance_p03 <- abundance_p03 %>% mutate(Empirical_p = 1 - rank(abundance_p03$Slope_abs) / length(abundance_p03$Slope_abs))
abundance_p04 <- abundance_p04 %>% mutate(Empirical_p = 1 - rank(abundance_p04$Slope_abs) / length(abundance_p04$Slope_abs))
abundance_p05 <- abundance_p05 %>% mutate(Empirical_p = 1 - rank(abundance_p05$Slope_abs) / length(abundance_p05$Slope_abs))
abundance_p06 <- abundance_p06 %>% mutate(Empirical_p = 1 - rank(abundance_p06$Slope_abs) / length(abundance_p06$Slope_abs))
abundance_p07 <- abundance_p07 %>% mutate(Empirical_p = 1 - rank(abundance_p07$Slope_abs) / length(abundance_p07$Slope_abs))
abundance_p08 <- abundance_p08 %>% mutate(Empirical_p = 1 - rank(abundance_p08$Slope_abs) / length(abundance_p08$Slope_abs))
abundance_p09 <- abundance_p09 %>% mutate(Empirical_p = 1 - rank(abundance_p09$Slope_abs) / length(abundance_p09$Slope_abs))
abundance_p10 <- abundance_p10 %>% mutate(Empirical_p = 1 - rank(abundance_p10$Slope_abs) / length(abundance_p10$Slope_abs))
abundance_p11 <- abundance_p11 %>% mutate(Empirical_p = 1 - rank(abundance_p11$Slope_abs) / length(abundance_p11$Slope_abs))

#Export as table delimited file
write_tsv(abundance_p01, "data/genomic_data/plink_input_p01.tsv")
write_tsv(abundance_p02, "data/genomic_data/plink_input_p02.tsv")
write_tsv(abundance_p03, "data/genomic_data/plink_input_p03.tsv")
write_tsv(abundance_p04, "data/genomic_data/plink_input_p04.tsv")
write_tsv(abundance_p05, "data/genomic_data/plink_input_p05.tsv")
write_tsv(abundance_p06, "data/genomic_data/plink_input_p06.tsv")
write_tsv(abundance_p07, "data/genomic_data/plink_input_p07.tsv")
write_tsv(abundance_p08, "data/genomic_data/plink_input_p08.tsv")
write_tsv(abundance_p09, "data/genomic_data/plink_input_p09.tsv")
write_tsv(abundance_p10, "data/genomic_data/plink_input_p10.tsv")
write_tsv(abundance_p11, "data/genomic_data/plink_input_p11.tsv")









write_csv(snp_set_abundance_p01, "data/genomic_data/snp_set_abundance_p01.csv")
write_csv(snp_set_abundance_p02, "data/genomic_data/snp_set_abundance_p02.csv")
write_csv(snp_set_abundance_p03, "data/genomic_data/snp_set_abundance_p03.csv")
write_csv(snp_set_abundance_p04, "data/genomic_data/snp_set_abundance_p04.csv")
write_csv(snp_set_abundance_p05, "data/genomic_data/snp_set_abundance_p05.csv")
write_csv(snp_set_abundance_p06, "data/genomic_data/snp_set_abundance_p06.csv")
write_csv(snp_set_abundance_p07, "data/genomic_data/snp_set_abundance_p07.csv")
write_csv(snp_set_abundance_p08, "data/genomic_data/snp_set_abundance_p08.csv")
write_csv(snp_set_abundance_p09, "data/genomic_data/snp_set_abundance_p09.csv")
write_csv(snp_set_abundance_p10, "data/genomic_data/snp_set_abundance_p10.csv")
write_csv(snp_set_abundance_p11, "data/genomic_data/snp_set_abundance_p11.csv")



