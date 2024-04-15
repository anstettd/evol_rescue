#############################################################################################################
## Calc WZA for windows
## Get empirical p-value
## Author Daniel Anstett
## 
## 
## Modified from Tom Booker WZA Vignette
## Last Modified August 3, 2022
#############################################################################################################
#Import libraries

library(tidyverse)

#Import files

snps_env1 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_env1_bf.csv")
snps_env2 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_env2_bf.csv")
snps_env3 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_env3_bf.csv")
snps_env4 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_env4_bf.csv")
snps_env5 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_env5_bf.csv")

snps_env6 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_env6_bf.csv")
snps_env7 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_env7_bf.csv")
snps_env8 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_env8_bf.csv")
snps_env9 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_env9_bf.csv")




###########################################################################################################

#Filter needed variables for WZA
#snps_env1_input <- snps_env1 %>% select(BF,win,MAF)
#snps_env2_input <- snps_env2 %>% select(BF,win,MAF)
#snps_env3_input <- snps_env3 %>% select(BF,win,MAF)
#snps_env4_input <- snps_env4 %>% select(BF,win,MAF)
#snps_env5_input <- snps_env5 %>% select(BF,win,MAF)

#snps_env6_input <- snps_env6 %>% select(BF,win,MAF)
#snps_env7_input <- snps_env7 %>% select(BF,win,MAF)
#snps_env8_input <- snps_env8 %>% select(BF,win,MAF)
#snps_env9_input <- snps_env9 %>% select(BF,win,MAF)


#Too large to store on github. Store locally
#write_csv(snps_env1_input, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_env1_input_bf.csv")
#write_csv(snps_env2_input, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_env2_input_bf.csv")
#write_csv(snps_env3_input, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_env3_input_bf.csv")
#write_csv(snps_env4_input, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_env4_input_bf.csv")
#write_csv(snps_env5_input, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_env5_input_bf.csv")   

#write_csv(snps_env6_input, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_env6_input_bf.csv")
#write_csv(snps_env7_input, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_env7_input_bf.csv")
#write_csv(snps_env8_input, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_env8_input_bf.csv")
#write_csv(snps_env9_input, "/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_env9_input_bf.csv")

#Run python script run_WZA.txt
#Import WZA scores
WZA_df_env1 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/env1_WZA_bf.csv")
WZA_df_env2 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/env2_WZA_bf.csv")
WZA_df_env3 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/env3_WZA_bf.csv")
WZA_df_env4 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/env4_WZA_bf.csv")
WZA_df_env5 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/env5_WZA_bf.csv")

WZA_df_env6 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/env6_WZA_bf.csv")
WZA_df_env7 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/env7_WZA_bf.csv")
WZA_df_env8 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/env8_WZA_bf.csv")
WZA_df_env9 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/env9_WZA_bf.csv")

###########################################################################################################
#Re-name windows
names(WZA_df_env1)[names(WZA_df_env1) == 'gene'] <- 'win'
names(WZA_df_env2)[names(WZA_df_env2) == 'gene'] <- 'win'
names(WZA_df_env3)[names(WZA_df_env3) == 'gene'] <- 'win'
names(WZA_df_env4)[names(WZA_df_env4) == 'gene'] <- 'win'
names(WZA_df_env5)[names(WZA_df_env5) == 'gene'] <- 'win'

names(WZA_df_env6)[names(WZA_df_env6) == 'gene'] <- 'win'
names(WZA_df_env7)[names(WZA_df_env7) == 'gene'] <- 'win'
names(WZA_df_env8)[names(WZA_df_env8) == 'gene'] <- 'win'
names(WZA_df_env9)[names(WZA_df_env9) == 'gene'] <- 'win'

#Make DF that has chr ID for each window
chr_env1 <- snps_env1 %>% select(chr,win) %>% distinct()
chr_env2 <- snps_env2 %>% select(chr,win) %>% distinct()
chr_env3 <- snps_env3 %>% select(chr,win) %>% distinct()
chr_env4 <- snps_env4 %>% select(chr,win) %>% distinct()
chr_env5 <- snps_env5 %>% select(chr,win) %>% distinct()

chr_env6 <- snps_env6 %>% select(chr,win) %>% distinct()
chr_env7 <- snps_env7 %>% select(chr,win) %>% distinct()
chr_env8 <- snps_env8 %>% select(chr,win) %>% distinct()
chr_env9 <- snps_env9 %>% select(chr,win) %>% distinct()


#Get chromosome information into WZA score
WZA_df_env1_chr <- left_join(WZA_df_env1,chr_env1,by="win")
WZA_df_env2_chr <- left_join(WZA_df_env2,chr_env2,by="win")
WZA_df_env3_chr <- left_join(WZA_df_env3,chr_env3,by="win")
WZA_df_env4_chr <- left_join(WZA_df_env4,chr_env4,by="win")
WZA_df_env5_chr <- left_join(WZA_df_env5,chr_env5,by="win")

WZA_df_env6_chr <- left_join(WZA_df_env6,chr_env6,by="win")
WZA_df_env7_chr <- left_join(WZA_df_env7,chr_env7,by="win")
WZA_df_env8_chr <- left_join(WZA_df_env8,chr_env8,by="win")
WZA_df_env9_chr <- left_join(WZA_df_env9,chr_env9,by="win")


#Input position for each window
WZA_df_env1_chr$pos <- (WZA_df_env1_chr$win+0.5)*10000
WZA_df_env2_chr$pos <- (WZA_df_env2_chr$win+0.5)*10000
WZA_df_env3_chr$pos <- (WZA_df_env3_chr$win+0.5)*10000
WZA_df_env4_chr$pos <- (WZA_df_env4_chr$win+0.5)*10000
WZA_df_env5_chr$pos <- (WZA_df_env5_chr$win+0.5)*10000

WZA_df_env6_chr$pos <- (WZA_df_env6_chr$win+0.5)*10000
WZA_df_env7_chr$pos <- (WZA_df_env7_chr$win+0.5)*10000
WZA_df_env8_chr$pos <- (WZA_df_env8_chr$win+0.5)*10000
WZA_df_env9_chr$pos <- (WZA_df_env9_chr$win+0.5)*10000


WZA_df_env1_chr <- na.omit(WZA_df_env1_chr)
WZA_df_env2_chr <- na.omit(WZA_df_env2_chr)
WZA_df_env3_chr <- na.omit(WZA_df_env3_chr)
WZA_df_env4_chr <- na.omit(WZA_df_env4_chr)
WZA_df_env5_chr <- na.omit(WZA_df_env5_chr)

WZA_df_env6_chr <- na.omit(WZA_df_env6_chr)
WZA_df_env7_chr <- na.omit(WZA_df_env7_chr)
WZA_df_env8_chr <- na.omit(WZA_df_env8_chr)
WZA_df_env9_chr <- na.omit(WZA_df_env9_chr)


write_csv(WZA_df_env1_chr, "data/genomic_data/WZA_win_env1_bf.csv")
write_csv(WZA_df_env2_chr, "data/genomic_data/WZA_win_env2_bf.csv")
write_csv(WZA_df_env3_chr, "data/genomic_data/WZA_win_env3_bf.csv")
write_csv(WZA_df_env4_chr, "data/genomic_data/WZA_win_env4_bf.csv")
write_csv(WZA_df_env5_chr, "data/genomic_data/WZA_win_env5_bf.csv")

write_csv(WZA_df_env6_chr, "data/genomic_data/WZA_win_env6_bf.csv")
write_csv(WZA_df_env7_chr, "data/genomic_data/WZA_win_env7_bf.csv")
write_csv(WZA_df_env8_chr, "data/genomic_data/WZA_win_env8_bf.csv")
write_csv(WZA_df_env9_chr, "data/genomic_data/WZA_win_env9_bf.csv")






