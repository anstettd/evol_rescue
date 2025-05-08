#############################################################################################################
## Calculate Major and Minor Allele for SNP set
## Done on full WZA BF filtered SNPs (no removal of non-linear clim SNPs)
## Author Daniel Anstett
## 
## 
## 
#############################################################################################################

library(tidyverse)

#Import SNP set per env variable. These have been filtered for Bayes, WZA cut offs and linear clim association
snp_set_env1 <- read_csv("data/genomic_data/snp_set_env1.csv") %>% select(chr_snp)
snp_set_env2 <- read_csv("data/genomic_data/snp_set_env2.csv") %>% select(chr_snp)
snp_set_env3 <- read_csv("data/genomic_data/snp_set_env3.csv") %>% select(chr_snp)
snp_set_env4 <- read_csv("data/genomic_data/snp_set_env4.csv") %>% select(chr_snp)
snp_set_env5 <- read_csv("data/genomic_data/snp_set_env5.csv") %>% select(chr_snp)
snp_set_env6 <- read_csv("data/genomic_data/snp_set_env6.csv") %>% select(chr_snp)
snp_set_env7 <- read_csv("data/genomic_data/snp_set_env7.csv") %>% select(chr_snp)
snp_set_env8 <- read_csv("data/genomic_data/snp_set_env8.csv") %>% select(chr_snp)
snp_set_env9 <- read_csv("data/genomic_data/snp_set_env9.csv") %>% select(chr_snp)

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

#Left join to get 55 pop abudnance data only for strong SNP set 
snp_set_abundance_env1 <- snp_set_env1 %>% left_join(loci_snp,by="chr_snp")
snp_set_abundance_env2 <- snp_set_env2 %>% left_join(loci_snp,by="chr_snp")
snp_set_abundance_env3 <- snp_set_env3 %>% left_join(loci_snp,by="chr_snp")
snp_set_abundance_env4 <- snp_set_env4 %>% left_join(loci_snp,by="chr_snp")
snp_set_abundance_env5 <- snp_set_env5 %>% left_join(loci_snp,by="chr_snp")
snp_set_abundance_env6 <- snp_set_env6 %>% left_join(loci_snp,by="chr_snp")
snp_set_abundance_env7 <- snp_set_env7 %>% left_join(loci_snp,by="chr_snp")
snp_set_abundance_env8 <- snp_set_env8 %>% left_join(loci_snp,by="chr_snp")
snp_set_abundance_env9 <- snp_set_env9 %>% left_join(loci_snp,by="chr_snp")

write_csv(snp_set_abundance_env1,"data/genomic_data/snp_set_abundance_full_env1.csv")
write_csv(snp_set_abundance_env2,"data/genomic_data/snp_set_abundance_full_env2.csv")  
write_csv(snp_set_abundance_env3,"data/genomic_data/snp_set_abundance_full_env3.csv")  
write_csv(snp_set_abundance_env4,"data/genomic_data/snp_set_abundance_full_env4.csv")  
write_csv(snp_set_abundance_env5,"data/genomic_data/snp_set_abundance_full_env5.csv")  
write_csv(snp_set_abundance_env6,"data/genomic_data/snp_set_abundance_full_env6.csv")  
write_csv(snp_set_abundance_env7,"data/genomic_data/snp_set_abundance_full_env7.csv")  
write_csv(snp_set_abundance_env8,"data/genomic_data/snp_set_abundance_full_env8.csv")  
write_csv(snp_set_abundance_env9,"data/genomic_data/snp_set_abundance_full_env9.csv")  
  


