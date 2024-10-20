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

##Import full Pi data set

#pop1
pop1<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/Pi/1_2010.sites.pi", header=F, sep="\t")
colnames(pop1) <- c("Chromosome","SNP","PI")
pop1_all <- pop1 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
rm(pop1)
##Filter Pi datatset by 2.1 Million SNPs
pop1_filter <-pop1_all %>% filter (chr_snp %in% as.character(loci_united$chr_snp))
rm(pop1_all)

#pop2
pop2<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/Pi/2_2010.sites.pi", header=F, sep="\t")
colnames(pop2) <- c("Chromosome","SNP","PI")
pop2_all <- pop2 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
rm(pop2)
##Filter Pi datatset by 2.1 Million SNPs
pop2_filter <-pop2_all %>% filter (chr_snp %in% as.character(loci_united$chr_snp))
rm(pop2_all)

#pop3
pop3<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/Pi/3_2010.sites.pi", header=F, sep="\t")
colnames(pop3) <- c("Chromosome","SNP","PI")
pop3_all <- pop3 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
rm(pop3)
##Filter Pi datatset by 2.1 Million SNPs
pop3_filter <-pop3_all %>% filter (chr_snp %in% as.character(loci_united$chr_snp))
rm(pop3_all)

#pop4
pop4<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/Pi/4_2010.sites.pi", header=F, sep="\t")
colnames(pop4) <- c("Chromosome","SNP","PI")
pop4_all <- pop4 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
rm(pop4)
##Filter Pi datatset by 2.1 Million SNPs
pop4_filter <-pop4_all %>% filter (chr_snp %in% as.character(loci_united$chr_snp))
rm(pop4_all)

#pop5
pop5<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/Pi/5_2010.sites.pi", header=F, sep="\t")
colnames(pop5) <- c("Chromosome","SNP","PI")
pop5_all <- pop5 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
rm(pop5)
##Filter Pi datatset by 2.1 Million SNPs
pop5_filter <-pop5_all %>% filter (chr_snp %in% as.character(loci_united$chr_snp))
rm(pop5_all)

#pop6
pop6<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/Pi/6_2010.sites.pi", header=F, sep="\t")
colnames(pop6) <- c("Chromosome","SNP","PI")
pop6_all <- pop6 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
rm(pop6)
##Filter Pi datatset by 2.1 Million SNPs
pop6_filter <-pop6_all %>% filter (chr_snp %in% as.character(loci_united$chr_snp))
rm(pop6_all)

#pop7
pop7<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/Pi/7_2010.sites.pi", header=F, sep="\t")
colnames(pop7) <- c("Chromosome","SNP","PI")
pop7_all <- pop7 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
rm(pop7)
##Filter Pi datatset by 2.1 Million SNPs
pop7_filter <-pop7_all %>% filter (chr_snp %in% as.character(loci_united$chr_snp))
rm(pop7_all)

#pop8
pop8<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/Pi/8_2011.sites.pi", header=F, sep="\t")
colnames(pop8) <- c("Chromosome","SNP","PI")
pop8_all <- pop8 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
rm(pop8)
##Filter Pi datatset by 2.1 Million SNPs
pop8_filter <-pop8_all %>% filter (chr_snp %in% as.character(loci_united$chr_snp))
rm(pop8_all)

#pop9
pop9<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/Pi/9_2010.sites.pi", header=F, sep="\t")
colnames(pop9) <- c("Chromosome","SNP","PI")
pop9_all <- pop9 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
rm(pop9)
##Filter Pi datatset by 2.1 Million SNPs
pop9_filter <-pop9_all %>% filter (chr_snp %in% as.character(loci_united$chr_snp))
rm(pop9_all)

#pop10
pop10<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/Pi/10_2011.sites.pi", header=F, sep="\t")
colnames(pop10) <- c("Chromosome","SNP","PI")
pop10_all <- pop10 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
rm(pop10)
##Filter Pi datatset by 2.1 Million SNPs
pop10_filter <-pop10_all %>% filter (chr_snp %in% as.character(loci_united$chr_snp))
rm(pop10_all)

#pop11
pop11<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/Pi/11_2010.sites.pi", header=F, sep="\t")
colnames(pop11) <- c("Chromosome","SNP","PI")
pop11_all <- pop11 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
rm(pop11)
##Filter Pi datatset by 2.1 Million SNPs
pop11_filter <-pop11_all %>% filter (chr_snp %in% as.character(loci_united$chr_snp))
rm(pop11_all)

#pop12
pop12<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/Pi/12_2010.sites.pi", header=F, sep="\t")
colnames(pop12) <- c("Chromosome","SNP","PI")
pop12_all <- pop12 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
rm(pop12)
##Filter Pi datatset by 2.1 Million SNPs
pop12_filter <-pop12_all %>% filter (chr_snp %in% as.character(loci_united$chr_snp))
rm(pop12_all)

#Demography only pops

#pop14 Hauser Creek
pop14<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/Pi/pop_demo_14.sites.pi", header=F, sep="\t")
colnames(pop14) <- c("Chromosome","SNP","PI")
pop14_all <- pop14 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
rm(pop14)
##Filter Pi datatset by 2.1 Million SNPs
pop14_filter <-pop14_all %>% filter (chr_snp %in% as.character(loci_united$chr_snp))
rm(pop14_all)

#pop15 Kitchen Creek
pop15<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/Pi/pop_demo_15.sites.pi", header=F, sep="\t")
colnames(pop15) <- c("Chromosome","SNP","PI")
pop15_all <- pop15 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
rm(pop15)
##Filter Pi datatset by 2.1 Million SNPs
pop15_filter <-pop15_all %>% filter (chr_snp %in% as.character(loci_united$chr_snp))
rm(pop15_all)

#pop17 Whitewater Caynon
pop17<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/Pi/pop_demo_17.sites.pi", header=F, sep="\t")
colnames(pop17) <- c("Chromosome","SNP","PI")
pop17_all <- pop17 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
rm(pop17)
##Filter Pi datatset by 2.1 Million SNPs
pop17_filter <-pop17_all %>% filter (chr_snp %in% as.character(loci_united$chr_snp))
rm(pop17_all)

#pop27 Buck Meadows
pop27<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/Pi/pop_demo_27.sites.pi", header=F, sep="\t")
colnames(pop27) <- c("Chromosome","SNP","PI")
pop27_all <- pop27 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
rm(pop27)
##Filter Pi datatset by 2.1 Million SNPs
pop27_filter <-pop27_all %>% filter (chr_snp %in% as.character(loci_united$chr_snp))
rm(pop27_all)

#pop28 Carlon
pop28<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/Pi/pop_demo_28.sites.pi", header=F, sep="\t")
colnames(pop28) <- c("Chromosome","SNP","PI")
pop28_all <- pop28 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
rm(pop28)
##Filter Pi datatset by 2.1 Million SNPs
pop28_filter <-pop28_all %>% filter (chr_snp %in% as.character(loci_united$chr_snp))
rm(pop28_all)

#pop29 Rainbow Pool
pop29<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/Pi/pop_demo_29.sites.pi", header=F, sep="\t")
colnames(pop29) <- c("Chromosome","SNP","PI")
pop29_all <- pop29 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
rm(pop29)
##Filter Pi datatset by 2.1 Million SNPs
pop29_filter <-pop29_all %>% filter (chr_snp %in% as.character(loci_united$chr_snp))
rm(pop29_all)

#pop55 Coast Fork Willamette 
pop55<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/Pi/pop_demo_55.sites.pi", header=F, sep="\t")
colnames(pop55) <- c("Chromosome","SNP","PI")
pop55_all <- pop55 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
rm(pop55)
##Filter Pi datatset by 2.1 Million SNPs
pop55_filter <-pop55_all %>% filter (chr_snp %in% as.character(loci_united$chr_snp))
rm(pop55_all)



#Filter out snp set
snp_set_p1 <-pop1_filter %>% filter (chr_snp %in% as.character(snp_set$chr_snp))
snp_set_p2 <-pop2_filter %>% filter (chr_snp %in% as.character(snp_set$chr_snp))
snp_set_p3 <-pop3_filter %>% filter (chr_snp %in% as.character(snp_set$chr_snp))
snp_set_p4 <-pop4_filter %>% filter (chr_snp %in% as.character(snp_set$chr_snp))
snp_set_p5 <-pop5_filter %>% filter (chr_snp %in% as.character(snp_set$chr_snp))
snp_set_p6 <-pop6_filter %>% filter (chr_snp %in% as.character(snp_set$chr_snp))
snp_set_p7 <-pop7_filter %>% filter (chr_snp %in% as.character(snp_set$chr_snp))
snp_set_p8 <-pop8_filter %>% filter (chr_snp %in% as.character(snp_set$chr_snp))
snp_set_p9 <-pop9_filter %>% filter (chr_snp %in% as.character(snp_set$chr_snp))
snp_set_p10 <-pop10_filter %>% filter (chr_snp %in% as.character(snp_set$chr_snp))
snp_set_p11 <-pop11_filter %>% filter (chr_snp %in% as.character(snp_set$chr_snp))
snp_set_p12 <-pop12_filter %>% filter (chr_snp %in% as.character(snp_set$chr_snp))

snp_set_p14 <-pop14_filter %>% filter (chr_snp %in% as.character(snp_set$chr_snp))
snp_set_p15 <-pop15_filter %>% filter (chr_snp %in% as.character(snp_set$chr_snp))
snp_set_p17 <-pop17_filter %>% filter (chr_snp %in% as.character(snp_set$chr_snp))
snp_set_p27 <-pop27_filter %>% filter (chr_snp %in% as.character(snp_set$chr_snp))
snp_set_p28 <-pop28_filter %>% filter (chr_snp %in% as.character(snp_set$chr_snp))
snp_set_p29 <-pop29_filter %>% filter (chr_snp %in% as.character(snp_set$chr_snp))
snp_set_p55 <-pop55_filter %>% filter (chr_snp %in% as.character(snp_set$chr_snp))


#Calc pi
pi_df <- 1:12
pi_df[13] <- 14
pi_df[14] <- 15
pi_df[15] <- 17
pi_df[16] <- 27
pi_df[17] <- 28
pi_df[18] <- 29
pi_df[19] <- 55
pi_df <- as.data.frame(pi_df)

#SNP set PI
pi_df[1,2] <- mean(as.numeric(snp_set_p1$PI))
pi_df[2,2] <- mean(as.numeric(snp_set_p2$PI))
pi_df[3,2] <- mean(as.numeric(snp_set_p3$PI))
pi_df[4,2] <- mean(as.numeric(snp_set_p4$PI))
pi_df[5,2] <- mean(as.numeric(snp_set_p5$PI))
pi_df[6,2] <- mean(as.numeric(snp_set_p6$PI))
pi_df[7,2] <- mean(as.numeric(snp_set_p7$PI))
pi_df[8,2] <- mean(as.numeric(snp_set_p8$PI))
pi_df[9,2] <- mean(as.numeric(snp_set_p9$PI))
pi_df[10,2] <- mean(as.numeric(snp_set_p10$PI))
pi_df[11,2] <- mean(as.numeric(snp_set_p11$PI))
pi_df[12,2] <- mean(as.numeric(snp_set_p12$PI))

pi_df[13,2] <- mean(as.numeric(snp_set_p14$PI))
pi_df[14,2] <- mean(as.numeric(snp_set_p15$PI))
pi_df[15,2] <- mean(as.numeric(snp_set_p17$PI))
pi_df[16,2] <- mean(as.numeric(snp_set_p27$PI))
pi_df[17,2] <- mean(as.numeric(snp_set_p28$PI))
pi_df[18,2] <- mean(as.numeric(snp_set_p29$PI))
pi_df[19,2] <- mean(as.numeric(snp_set_p55$PI))


#Full set PI
pi_df[1,3] <- mean(as.numeric(pop1_filter$PI))
pi_df[2,3] <- mean(as.numeric(pop2_filter$PI))
pi_df[3,3] <- mean(as.numeric(pop3_filter$PI))
pi_df[4,3] <- mean(as.numeric(pop4_filter$PI))
pi_df[5,3] <- mean(as.numeric(pop5_filter$PI))
pi_df[6,3] <- mean(as.numeric(pop6_filter$PI))
pi_df[7,3] <- mean(as.numeric(pop7_filter$PI))
pi_df[8,3] <- mean(as.numeric(pop8_filter$PI))
pi_df[9,3] <- mean(as.numeric(pop9_filter$PI))
pi_df[10,3] <- mean(as.numeric(pop10_filter$PI))
pi_df[11,3] <- mean(as.numeric(pop11_filter$PI))
pi_df[12,3] <- mean(as.numeric(pop12_filter$PI))

pi_df[13,3] <- mean(as.numeric(pop14_filter$PI))
pi_df[14,3] <- mean(as.numeric(pop15_filter$PI))
pi_df[15,3] <- mean(as.numeric(pop17_filter$PI))
pi_df[16,3] <- mean(as.numeric(pop27_filter$PI))
pi_df[17,3] <- mean(as.numeric(pop28_filter$PI))
pi_df[18,3] <- mean(as.numeric(pop29_filter$PI))
pi_df[19,3] <- mean(as.numeric(pop55_filter$PI))

colnames(pi_df) <- c("Site","pi_snp_set","pi_all_snps")

write_csv(pi_df, "data/genomic_data/raw_pi.csv")
