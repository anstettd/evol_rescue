##################################################################################
## Finding which allele is associated with climate change
## For peak window SNPS (WZA) with BF (BayPass) > 10 and all BF>30
## Done for all alleles positively associated with climate change
## Author Daniel Anstett
## 
## 
## Last Modified May 25, 2023
###################################################################################
#Import libraries
library(tidyverse)
library(boot) #this has the inv.logit function

###################################################################################

#Functions
theme_spaghetti <- function(){ 
  
  theme_classic() %+replace%    #replace elements we want to change
    theme(legend.position = "none",
          axis.text.x = element_text(size = 14, face = "bold", angle = 45,hjust = 1, vjust = 1), 
          axis.title = element_text(size =16, face = "bold"), 
          axis.text.y = element_text(size = 14, face = "bold"),legend.title = element_blank(),
          legend.text = element_text(size=12,face="bold"),
          strip.background = element_blank(), 
          strip.text.x=element_text(size=14,face="bold",hjust=0,vjust=-1.2))
}

###################################################################################


#Read In Data


#Import Baseline climate 
climate <- read_csv("data/climate_pop.csv")

#Import pop names (site_year names)
pop_order<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")

pop_order_base<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", 
                           header=F, sep="\t")

#Read in timeseries SNP abundances
snp2_time <- read_csv("data/snp_set_time_env2.csv")

#Read in baseline SNP abundances
env2_base <- read_csv("data/snp_set_base_env2.csv")

###################################################################################

#Ensure baseline has same SNPs as timeseries and remove any not in the timeseries 
env2_base <- env2_base %>% filter (chr_snp %in% as.character(snp2_time$chr_snp))


#Make abundance table for SNP A
snp_prop_A<- env2_base %>% select (chr_snp)
counter=2
pop_num=1
for (i in seq(1,dim(env2_base)[2]-1,2)) {
  
  snpA<-paste("V",i, sep="") #sets up string for snpA for pop num
  P<-paste("P", pop_num, sep="") #sets up paper_ID (population ID)
  tmp<-env2_base %>% select(snpA) #graps all snps per paper_ID
  
  colnames(tmp)<-"A" #renames column headers
  
  snp_prop_A[,counter]<-tmp$A #assign A
  colnames (snp_prop_A)[counter]<-P 
  
  counter<-counter+1
  pop_num<-pop_num+1
}


#Make abundance table for SNP B

snp_prop_B<- env2_base %>% select (chr_snp)
counter=2
pop_num=1
for (i in seq(1,dim(env2_base)[2]-1,2)) {
  
  snpB<-paste("V",i+1, sep="") #sets up string for snpB for pop num
  P<-paste("P", pop_num, sep="") #sets up paper_ID (population ID)
  tmp<-env2_base %>% select(snpB) #graps all snps per paper_ID
  
  colnames(tmp)<-"B" #renames column headers
  
  snp_prop_B[,counter]<-tmp$B #assign B
  colnames (snp_prop_B)[counter]<-P 
  
  counter<-counter+1
  pop_num<-pop_num+1
}


## Make table with climate change associated SNP abundances for timeseries using baseline climate


freq.temp2 <- data.frame()
for (i in 1:dim(snp_prop_A)[1]) { #for each SNP
  
  env_pop <- climate %>% select(Paper_ID,MAP) #Select population ID and climate data
  snp_prop_A_tmp<-snp_prop_A %>% filter(chr_snp==snp_prop_A$chr_snp[i]) %>% 
    select(-chr_snp) #filter prop A data 
  snp_prop_B_tmp<-snp_prop_B %>% filter(chr_snp==snp_prop_B$chr_snp[i]) %>%
    select(-chr_snp) #filter prop B data 
  
  env_pop <- cbind(env_pop,t(snp_prop_A_tmp), t(snp_prop_B_tmp))
  colnames(env_pop)[3] <-"prop_A"
  colnames(env_pop)[4]<-"prop_B"
  
  df_A <- cbind(env_pop$prop_A,env_pop$prop_B)
  df_B <- cbind(env_pop$prop_B,env_pop$prop_A)
  
  
  #This is the glm where we associate the A and B SNPs with climate, in this case it's MAP
  
  lm.temp_A <- glm(df_A~env_pop[,2],family=binomial) # save glm of climate predicting prop A
  lm.temp_B <- glm(df_B~env_pop[,2],family=binomial) # save glm of climate predicting prop B
  
  #Here, I'm just printing out the invserse logit of the glm 
  print(paste(inv.logit(lm.temp_A$coefficients[2]), inv.logit(lm.temp_B$coefficients[2]), sep=" "))
  
  #This is where we've run into problems, I'm not sure what condition to use to force the selection of the 
  #climate associated SNP. I've coded it this way because the designation of A and B are arbitrary, and I want
  #to force A to always be the climate associated state. This code reorders the abudances of the SNPs so that
  #the number of climate associated SNPs will always be under column "True_SNP_A", and the abundance of non-climate
  #associated SNPs will always be under column "True_SNP_B"
  
  #This should be the case where SNP A is climate associated
  if(inv.logit(lm.temp_A$coefficients[2])<0.5){ 
    print(paste("A, ", chr_snp=snp_prop_A$chr_snp[i]))
    tmp_in<-env_pop %>% select(Paper_ID, prop_A, prop_B) %>% 
      mutate (ENV="MAP", chr_snp=snp_prop_A$chr_snp[i], SNP_Select="A", 
              Coeff_Pick=inv.logit(lm.temp_A$coefficients[2]), 
              Coeff_Unpick=inv.logit(lm.temp_B$coefficients[2]))
    
    colnames(tmp_in)[2]<-"True_SNP_A"
    colnames(tmp_in)[3]<-"True_SNP_B"
    
    freq.temp2<-rbind(freq.temp2,tmp_in)
    #    print(summary(lm.temp_A)) 
    #    print(lm.temp_A) 
    
    
    #This should be the case where SNP A is climate associated
  } else if (inv.logit(lm.temp_B$coefficients[2])<0.5){
    print(paste("B, ", chr_snp=snp_prop_A$chr_snp[i]))
    
    tmp_in<-env_pop %>% select(Paper_ID, prop_B, prop_A) %>% 
      mutate (ENV="MAP", chr_snp=snp_prop_A$chr_snp[i], SNP_Select="B",
              Coeff_Pick=inv.logit(lm.temp_B$coefficients[2]), 
              Coeff_Unpick=inv.logit(lm.temp_A$coefficients[2]))
    
    colnames(tmp_in)[2]<-"True_SNP_A"
    colnames(tmp_in)[3]<-"True_SNP_B"
    
    freq.temp2<-rbind(freq.temp2,tmp_in)
    #    print(summary(lm.temp_B)) 
    #    print(lm.temp_B) 
    
  }
}



#Convert the abundance counts to be number of 1s and 0s

#output df
abund_table2<-data.frame()

for (i in 1:dim(freq.temp2)[1]){
  Binomial_A<-c(rep(1, freq.temp2[i, "True_SNP_A"]), rep(0,freq.temp2[i, "True_SNP_B"]))
  tmp_df<- as.data.frame(Binomial_A) %>% mutate (Paper_ID=freq.temp2$Paper_ID[i], 
                                                 chr_snp=freq.temp2$chr_snp[i], 
                                                 ENV=freq.temp2$ENV[i])
  abund_table2<-rbind(abund_table2, tmp_df)
  
}

#Merge with climate data
abund_clim <- left_join(abund_table2,climate,by="Paper_ID")


#Plot frequency vs climate
plot_env_2<-ggplot(data=abund_clim,aes(MAP,Binomial_A,group=chr_snp)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.2,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="MAP")  + theme_spaghetti()

plot_env_2