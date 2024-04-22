##################################################################################
## Randomize major and minor alleles prior to GLM
## Author Daniel Anstett
## 
## 
## Last Modified Dec 30, 2023
###################################################################################
#Import libraries
library(tidyverse)

###################################################################################
#Functions

#Make Randomized tables
#test_swiss<-swiss_minor_major[1:10,]

rand_50_50<-function(df, j){
  fifty_fifty <- data.frame()
for (i in 1:dim(df)[1]){
  set.seed(i+j)
  rand_int <- sample(0:1, 1)
  tmp_chr_snp<-df %>% filter (chr_snp==as.character(df$chr_snp)[i])
  print(i)
  
  if(rand_int==0){
    fifty_fifty<-rbind(fifty_fifty,tmp_chr_snp)
  }else{
    
    for(j in 0:61){
      A<-tmp_chr_snp[1,3+(2*j)]
      B<-tmp_chr_snp[1,2+(2*j)]
      
      tmp_chr_snp[1,2+(2*j)] <- A
      tmp_chr_snp[1,3+(2*j)] <- B
    }
    fifty_fifty<-rbind(fifty_fifty,tmp_chr_snp)
    
  }
}
  return(fifty_fifty)
}

###################################################################################

#Randomize whether major or minor allele is selected

seed_num<-0
swiss_1 <- read.csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_1.csv")
fifty_fifty_1<-rand_50_50(swiss_1, seed_num)
write_csv(fifty_fifty_1,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_50_50_abund_table_1.csv")
seed_num<-seed_num + dim(swiss_1)[1]
rm(swiss_1)
rm(fifty_fifty_1)


swiss_2 <- read.csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_2.csv")
fifty_fifty_2<-rand_50_50(swiss_2, seed_num)
write_csv(fifty_fifty_2,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_50_50_abund_table_2.csv")
seed_num<-seed_num + dim(swiss_2)[1]
rm(swiss_2)
rm(fifty_fifty_2)


swiss_3 <- read.csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_3.csv")
fifty_fifty_3<-rand_50_50(swiss_3, seed_num)
write_csv(fifty_fifty_3,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_50_50_abund_table_3.csv")
seed_num<-seed_num + dim(swiss_3)[1]
rm(swiss_3)
rm(fifty_fifty_3)


swiss_4 <- read.csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_4.csv")
fifty_fifty_4<-rand_50_50(swiss_4, seed_num)
write_csv(fifty_fifty_4,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_50_50_abund_table_4.csv")
seed_num<-seed_num + dim(swiss_4)[1]
rm(swiss_4)
rm(fifty_fifty_4)


swiss_5 <- read.csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_5.csv")
fifty_fifty_5<-rand_50_50(swiss_5, seed_num)
write_csv(fifty_fifty_5,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_50_50_abund_table_5.csv")
seed_num<-seed_num + dim(swiss_5)[1]
rm(swiss_5)
rm(fifty_fifty_5)


swiss_6 <- read.csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_6.csv")
fifty_fifty_6<-rand_50_50(swiss_6, seed_num)
write_csv(fifty_fifty_6,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_50_50_abund_table_6.csv")
seed_num<-seed_num + dim(swiss_6)[1]
rm(swiss_6)
rm(fifty_fifty_6)


swiss_7 <- read.csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_7.csv")
fifty_fifty_7<-rand_50_50(swiss_7, seed_num)
write_csv(fifty_fifty_7,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_50_50_abund_table_7.csv")
seed_num<-seed_num + dim(swiss_7)[1]
rm(swiss_7)
rm(fifty_fifty_7)


swiss_8 <- read.csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_8.csv")
fifty_fifty_8<-rand_50_50(swiss_8, seed_num)
write_csv(fifty_fifty_8,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_50_50_abund_table_8.csv")
seed_num<-seed_num + dim(swiss_8)[1]
rm(swiss_8)
rm(fifty_fifty_8)


swiss_9 <- read.csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_9.csv")
fifty_fifty_9<-rand_50_50(swiss_9, seed_num)
write_csv(fifty_fifty_9,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_50_50_abund_table_9.csv")
seed_num<-seed_num + dim(swiss_9)[1]
rm(swiss_9)
rm(fifty_fifty_9)


#Reset environment 

seed_num <- 900000
swiss_10 <- read.csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_10.csv")
fifty_fifty_10<-rand_50_50(swiss_10, seed_num)
write_csv(fifty_fifty_10,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_50_50_abund_table_10.csv")
seed_num<-seed_num + dim(swiss_10)[1]
rm(swiss_10)
rm(fifty_fifty_10)


swiss_11 <- read.csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_11.csv")
fifty_fifty_11<-rand_50_50(swiss_11, seed_num)
write_csv(fifty_fifty_11,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_50_50_abund_table_11.csv")
seed_num<-seed_num + dim(swiss_11)[1]
rm(swiss_11)
rm(fifty_fifty_11)


swiss_12 <- read.csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_12.csv")
fifty_fifty_12<-rand_50_50(swiss_12, seed_num)
write_csv(fifty_fifty_12,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_50_50_abund_table_12.csv")
seed_num<-seed_num + dim(swiss_12)[1]
rm(swiss_12)
rm(fifty_fifty_12)


swiss_13 <- read.csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_13.csv")
fifty_fifty_13<-rand_50_50(swiss_13, seed_num)
write_csv(fifty_fifty_13,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_50_50_abund_table_13.csv")
seed_num<-seed_num + dim(swiss_13)[1]
rm(swiss_13)
rm(fifty_fifty_13)


swiss_14 <- read.csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_14.csv")
fifty_fifty_14<-rand_50_50(swiss_14, seed_num)
write_csv(fifty_fifty_14,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_50_50_abund_table_14.csv")
seed_num<-seed_num + dim(swiss_14)[1]
rm(swiss_14)
rm(fifty_fifty_14)


swiss_15 <- read.csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_15.csv")
fifty_fifty_15<-rand_50_50(swiss_15, seed_num)
write_csv(fifty_fifty_15,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_50_50_abund_table_15.csv")
seed_num<-seed_num + dim(swiss_15)[1]
rm(swiss_15)
rm(fifty_fifty_15)


swiss_16 <- read.csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_16.csv")
fifty_fifty_16<-rand_50_50(swiss_16,seed_num)
write_csv(fifty_fifty_16,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_50_50_abund_table_16.csv")
seed_num<-seed_num + dim(swiss_16)[1]
rm(swiss_16)
rm(fifty_fifty_16)


swiss_17 <- read.csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_17.csv")
fifty_fifty_17<-rand_50_50(swiss_17, seed_num)
write_csv(fifty_fifty_17,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_50_50_abund_table_17.csv")
rm(swiss_17)
rm(fifty_fifty_17)










