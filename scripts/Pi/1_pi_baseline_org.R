##################################################################################
## Pi data prep
## Author Daniel Anstett
## 
## 
## Last Modified July 11, 2023
###################################################################################
##Libraries
library(tidyverse)

#Import SNPlist
snp_set <- read_csv("data/genomic_data/snp_set_env.csv")

#Import 2.1 M loci
loci_base<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", 
                      header=F, sep="\t")
colnames(loci_base) <- c("Chromosome","SNP")
loci_united <- loci_base %>% unite(chr_snp,"Chromosome","SNP",sep="_")
rm(loci_base)

###################################################################################

##Calculate snp set and global pi for all 55 baseline populations

pi_df <- as.data.frame(1:55)


for (i in 1:55){
  fname_in<-paste("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/Pi/pop_baseline_", i,".sites.pi", sep="")
  pop_temp<-read.table(fname_in, header=F, sep="\t")
  colnames(pop_temp) <- c("Chromosome","SNP","PI")
  pop_temp <- pop_temp %>% unite(chr_snp,"Chromosome","SNP",sep="_")
  pop_filter <-pop_temp %>% filter (chr_snp %in% as.character(loci_united$chr_snp))
  snp_set <-pop_filter %>% filter (chr_snp %in% as.character(snp_set$chr_snp))
  
  pi_df[i,2] <- mean(as.numeric(snp_set$PI))
  pi_df[i,3] <- mean(as.numeric(pop_filter$PI))
  rm(pop_temp)
  print(i)
}

colnames(pi_df) <- c("Site","pi_snp_set","pi_all_snps")

write_csv(pi_df, "data/genomic_data/baseline_pi.csv")

#demography pi

pi_demo <- pi_df %>% filter(between(Site,1,12) | Site==14 | Site==15| Site==17| 
                              Site==27| Site==28| Site==29| Site==55)


write_csv(pi_demo, "data/genomic_data/raw_pi.csv")




