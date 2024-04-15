##################################################################################
## Plot BayPass Manhattan Plots
## Author Daniel Anstett & Julia Anstett
## 
##
## Last Modified June 1, 2023
###################################################################################
#Functions
theme_man <- function(){ 
  theme_classic() %+replace%    #replace elements we want to change
    theme(
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(hjust = 0.8),
      strip.background = element_blank(),
      strip.placement = "outside",
      axis.text.x = element_text(size=15, face="bold"),
      axis.text.y = element_text(size=15,face="bold"),
      axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
      axis.title.y = element_text(color="black", size=18,vjust = 1.4, face="bold",hjust=0.5))
}

###################################################################################
#Import libraries
library(tidyverse)
library(qqman)

#Import chromosome size
chr_size <- read_csv("data/genomic_data/chr_size.csv")
chr_size[,3] <- cumsum(chr_size$size) #get cumulative chromosome position
colnames(chr_size)[3] <- "poz"


#Import snp env associations
env1 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_1_trim.tsv",header=F, sep=" ")
env2 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_2_trim.tsv",header=F, sep=" ")
env3 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_3_trim.tsv",header=F, sep=" ")
env4 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_4_trim.tsv",header=F, sep=" ")
env5 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_5_trim.tsv",header=F, sep=" ")
env6 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_6_trim.tsv",header=F, sep=" ")
env7 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_7_trim.tsv",header=F, sep=" ")
env8 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_8_trim.tsv",header=F, sep=" ")
env9 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_9_trim.tsv",header=F, sep=" ")


#Name Columns
colnames(env1) <- c("Chromosome","SNP","Env","BF")
colnames(env2) <- c("Chromosome","SNP","Env","BF")
colnames(env3) <- c("Chromosome","SNP","Env","BF")
colnames(env4) <- c("Chromosome","SNP","Env","BF")
colnames(env5) <- c("Chromosome","SNP","Env","BF")
colnames(env6) <- c("Chromosome","SNP","Env","BF")
colnames(env7) <- c("Chromosome","SNP","Env","BF")
colnames(env8) <- c("Chromosome","SNP","Env","BF")
colnames(env9) <- c("Chromosome","SNP","Env","BF")


########################################################################################################################
#Make Chromosome SNP variable for easier left_joining later
#Make a cumulative basepair value
env1_united <- env1 %>% unite(chr_snp,Chromosome,SNP)
env1_united <- env1_united %>% select(chr_snp)
env1_united <- cbind(env1_united,env1)
env1_united <- env1_united%>% mutate(CHR=as.numeric(gsub("CE10_chr","",Chromosome))) %>%
  mutate(BP = ifelse(CHR == 1, SNP,
                     ifelse(CHR==2, SNP+chr_size$poz[1],
                            ifelse(CHR==3, SNP+chr_size$poz[2],
                                   ifelse(CHR==4, SNP+chr_size$poz[3],
                                          ifelse(CHR==5, SNP+chr_size$poz[4],
                                                 ifelse(CHR==6, SNP+chr_size$poz[5],
                                                        ifelse(CHR==7, SNP+chr_size$poz[6],
                                                               SNP+chr_size$poz[7]))))))))

env2_united <- env2 %>% unite(chr_snp,Chromosome,SNP)
env2_united <- env2_united %>% select(chr_snp)
env2_united <- cbind(env2_united,env2)
env2_united <- env2_united %>% mutate(CHR=as.numeric(gsub("CE10_chr","",Chromosome))) %>%
  mutate(BP = ifelse(CHR == 1, SNP,
                     ifelse(CHR==2, SNP+chr_size$poz[1],
                            ifelse(CHR==3, SNP+chr_size$poz[2],
                                   ifelse(CHR==4, SNP+chr_size$poz[3],
                                          ifelse(CHR==5, SNP+chr_size$poz[4],
                                                 ifelse(CHR==6, SNP+chr_size$poz[5],
                                                        ifelse(CHR==7, SNP+chr_size$poz[6],
                                                               SNP+chr_size$poz[7]))))))))

env3_united <- env3 %>% unite(chr_snp,Chromosome,SNP)
env3_united <- env3_united %>% select(chr_snp)
env3_united <- cbind(env3_united,env3)
env3_united <- env3_united %>% mutate(CHR=as.numeric(gsub("CE10_chr","",Chromosome))) %>%
  mutate(BP = ifelse(CHR == 1, SNP,
                     ifelse(CHR==2, SNP+chr_size$poz[1],
                            ifelse(CHR==3, SNP+chr_size$poz[2],
                                   ifelse(CHR==4, SNP+chr_size$poz[3],
                                          ifelse(CHR==5, SNP+chr_size$poz[4],
                                                 ifelse(CHR==6, SNP+chr_size$poz[5],
                                                        ifelse(CHR==7, SNP+chr_size$poz[6],
                                                               SNP+chr_size$poz[7]))))))))

env4_united <- env4 %>% unite(chr_snp,Chromosome,SNP)
env4_united <- env4_united %>% select(chr_snp)
env4_united <- cbind(env4_united,env4)
env4_united <- env4_united %>% mutate(CHR=as.numeric(gsub("CE10_chr","",Chromosome))) %>%
  mutate(BP = ifelse(CHR == 1, SNP,
                     ifelse(CHR==2, SNP+chr_size$poz[1],
                            ifelse(CHR==3, SNP+chr_size$poz[2],
                                   ifelse(CHR==4, SNP+chr_size$poz[3],
                                          ifelse(CHR==5, SNP+chr_size$poz[4],
                                                 ifelse(CHR==6, SNP+chr_size$poz[5],
                                                        ifelse(CHR==7, SNP+chr_size$poz[6],
                                                               SNP+chr_size$poz[7]))))))))

env5_united <- env5 %>% unite(chr_snp,Chromosome,SNP)
env5_united <- env5_united %>% select(chr_snp)
env5_united <- cbind(env5_united,env5)
env5_united <- env5_united %>% mutate(CHR=as.numeric(gsub("CE10_chr","",Chromosome))) %>%
  mutate(BP = ifelse(CHR == 1, SNP,
                     ifelse(CHR==2, SNP+chr_size$poz[1],
                            ifelse(CHR==3, SNP+chr_size$poz[2],
                                   ifelse(CHR==4, SNP+chr_size$poz[3],
                                          ifelse(CHR==5, SNP+chr_size$poz[4],
                                                 ifelse(CHR==6, SNP+chr_size$poz[5],
                                                        ifelse(CHR==7, SNP+chr_size$poz[6],
                                                               SNP+chr_size$poz[7]))))))))

env6_united <- env6 %>% unite(chr_snp,Chromosome,SNP)
env6_united <- env6_united %>% select(chr_snp)
env6_united <- cbind(env6_united,env6)
env6_united <- env6_united %>% mutate(CHR=as.numeric(gsub("CE10_chr","",Chromosome))) %>%
  mutate(BP = ifelse(CHR == 1, SNP,
                     ifelse(CHR==2, SNP+chr_size$poz[1],
                            ifelse(CHR==3, SNP+chr_size$poz[2],
                                   ifelse(CHR==4, SNP+chr_size$poz[3],
                                          ifelse(CHR==5, SNP+chr_size$poz[4],
                                                 ifelse(CHR==6, SNP+chr_size$poz[5],
                                                        ifelse(CHR==7, SNP+chr_size$poz[6],
                                                               SNP+chr_size$poz[7]))))))))

env7_united <- env7 %>% unite(chr_snp,Chromosome,SNP)
env7_united <- env7_united %>% select(chr_snp)
env7_united <- cbind(env7_united,env7)
env7_united <- env7_united %>% mutate(CHR=as.numeric(gsub("CE10_chr","",Chromosome))) %>%
  mutate(BP = ifelse(CHR == 1, SNP,
                     ifelse(CHR==2, SNP+chr_size$poz[1],
                            ifelse(CHR==3, SNP+chr_size$poz[2],
                                   ifelse(CHR==4, SNP+chr_size$poz[3],
                                          ifelse(CHR==5, SNP+chr_size$poz[4],
                                                 ifelse(CHR==6, SNP+chr_size$poz[5],
                                                        ifelse(CHR==7, SNP+chr_size$poz[6],
                                                               SNP+chr_size$poz[7]))))))))

env8_united <- env8 %>% unite(chr_snp,Chromosome,SNP)
env8_united <- env8_united %>% select(chr_snp)
env8_united <- cbind(env8_united,env8)
env8_united <- env8_united %>% mutate(CHR=as.numeric(gsub("CE10_chr","",Chromosome))) %>%
  mutate(BP = ifelse(CHR == 1, SNP,
                     ifelse(CHR==2, SNP+chr_size$poz[1],
                            ifelse(CHR==3, SNP+chr_size$poz[2],
                                   ifelse(CHR==4, SNP+chr_size$poz[3],
                                          ifelse(CHR==5, SNP+chr_size$poz[4],
                                                 ifelse(CHR==6, SNP+chr_size$poz[5],
                                                        ifelse(CHR==7, SNP+chr_size$poz[6],
                                                               SNP+chr_size$poz[7]))))))))

env9_united <- env9 %>% unite(chr_snp,Chromosome,SNP)
env9_united <- env9_united %>% select(chr_snp)
env9_united <- cbind(env9_united,env9)
env9_united <- env9_united %>% mutate(CHR=as.numeric(gsub("CE10_chr","",Chromosome))) %>%
  mutate(BP = ifelse(CHR == 1, SNP,
                     ifelse(CHR==2, SNP+chr_size$poz[1],
                            ifelse(CHR==3, SNP+chr_size$poz[2],
                                   ifelse(CHR==4, SNP+chr_size$poz[3],
                                          ifelse(CHR==5, SNP+chr_size$poz[4],
                                                 ifelse(CHR==6, SNP+chr_size$poz[5],
                                                        ifelse(CHR==7, SNP+chr_size$poz[6],
                                                               SNP+chr_size$poz[7]))))))))



########################################################################################################################
#Test plot with only two chromosomes
#env1_united_chr1_2 <- env1_united %>% filter(CHR<3)

#axisdf_1 <- env1_united_chr1_2 %>% group_by(CHR) %>% summarize(center=( max(BP) + min(BP) ) / 2 )

#ggplot(env1_united_chr1_2, aes(x=BP, y=BF)) +
  
  # Show all points
#  geom_point(aes(color=as.factor(CHR)), alpha=0.8,size=0.5) +
#  scale_color_manual(values = rep(c("black", "magenta2"), 22 )) +
#  geom_hline(yintercept=10, linetype="dashed",color = "deepskyblue1", size=0.9) +
  # custom X axis:
#  scale_x_continuous(expand = c(0, 0), label = axisdf_1$CHR, breaks= axisdf_1$center) +
#  scale_y_continuous(breaks = c(-20,-10,0,10,20,30,40)) +
#  theme_classic() +
#  labs(
#    y = "Bayes Factor",
#    x = "Position")+
#  theme(
#    legend.position="none",
#    panel.border = element_blank(),
#    panel.grid.major = element_blank(),
#    panel.grid.minor = element_blank(),
#    panel.background = element_blank(),
#    plot.title = element_text(hjust = 0.8),
#    strip.background = element_blank(),
#    strip.placement = "outside",
#    axis.text.x = element_text(size=15, face="bold"),
#    axis.text.y = element_text(size=15,face="bold"),
#    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
#    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5)
#  )




########################################################################################################################
#Make Manhattan Plots

#ENV 1 - MAT 
ggplot(env1_united, aes(x=BP/1e6, y=BF)) +
  geom_point(aes(color=as.factor(CHR)), alpha=0.8,size=1.4) +
  scale_color_manual(values = rep(c("black", "darkgoldenrod"), 22 )) +
  geom_hline(yintercept=10, linetype="dashed",color = "red",lty = 2, lwd = 1) +
  geom_hline(yintercept=30, linetype="dashed",color = "deepskyblue", lty = 2, lwd = 1) +
  scale_y_continuous (breaks = c(-20,-10,0,10,20,30,40,50,60), limits=c(-20,60))+
  labs( y = "-log10(BF)", x = "Position (Mbp)")+
  theme_classic() +
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.8),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.4, face="bold",hjust=0.5))
ggsave("Graphs/Manhattan_BF/bf_manhattan_1.png",width=13, height=4.5)

#ENV 2 - MAP 
ggplot(env2_united, aes(x=BP/1e6, y=BF)) +
  geom_point(aes(color=as.factor(CHR)), alpha=0.8,size=1.4) +
  scale_color_manual(values = rep(c("black", "deepskyblue"), 22 )) +
  geom_hline(yintercept=10, linetype="dashed",color = "red",lty = 2, lwd = 1) +
  geom_hline(yintercept=30, linetype="dashed",color = "deepskyblue", lty = 2, lwd = 1) +
  scale_y_continuous (breaks = c(-20,-10,0,10,20,30,40))+
  labs( y = "-log10(BF)", x = "Position (Mbp)")+
  theme_classic() +
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.8),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.4, face="bold",hjust=0.5))
ggsave("Graphs/Manhattan_BF/bf_manhattan_2.png",width=13, height=4.5)


#ENV 3 - PAS 
ggplot(env3_united, aes(x=BP/1e6, y=BF)) +
  geom_point(aes(color=as.factor(CHR)), alpha=0.8,size=1.4) +
  scale_color_manual(values = rep(c("black", "deepskyblue"), 22 )) +
  geom_hline(yintercept=10, linetype="dashed",color = "red",lty = 2, lwd = 1) +
  geom_hline(yintercept=30, linetype="dashed",color = "deepskyblue", lty = 2, lwd = 1) +
  scale_y_continuous (breaks = c(-20,-10,0,10,20,30,40,50,60), limits=c(-20,60))+
  labs( y = "-log10(BF)", x = "Position (Mbp)")+
  theme_classic() +
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.8),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.4, face="bold",hjust=0.5))
ggsave("Graphs/Manhattan_BF/bf_manhattan_3.png",width=13, height=4.5)


#ENV 4 - EXT 
ggplot(env4_united, aes(x=BP/1e6, y=BF)) +
  geom_point(aes(color=as.factor(CHR)), alpha=0.8,size=1.4) +
  scale_color_manual(values = rep(c("black", "darkgoldenrod"), 22 )) +
  geom_hline(yintercept=10, linetype="dashed",color = "red",lty = 2, lwd = 1) +
  geom_hline(yintercept=30, linetype="dashed",color = "deepskyblue", lty = 2, lwd = 1) +
  scale_y_continuous (breaks = c(-20,-10,0,10,20,30,40,50,60), limits=c(-20,60))+
  labs( y = "-log10(BF)", x = "Position (Mbp)")+
  theme_classic() +
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.8),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.4, face="bold",hjust=0.5))
ggsave("Graphs/Manhattan_BF/bf_manhattan_4.png",width=13, height=4.5)

#ENV 5 - CMD 
ggplot(env5_united, aes(x=BP/1e6, y=BF)) +
  geom_point(aes(color=as.factor(CHR)), alpha=0.8,size=1.4) +
  scale_color_manual(values = rep(c("black", "magenta3"), 22 )) +
  geom_hline(yintercept=10, linetype="dashed",color = "red",lty = 2, lwd = 1) +
  geom_hline(yintercept=30, linetype="dashed",color = "deepskyblue", lty = 2, lwd = 1) +
  scale_y_continuous (breaks = c(-20,-10,0,10,20,30,40,50,60), limits=c(-20,60))+
  labs( y = "-log10(BF)", x = "Position (Mbp)")+
  theme_classic() +
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.8),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.4, face="bold",hjust=0.5))
ggsave("Graphs/Manhattan_BF/bf_manhattan_5.png",width=13, height=4.5)


#ENV 6 - Tave_wt 
ggplot(env6_united, aes(x=BP/1e6, y=BF)) +
  geom_point(aes(color=as.factor(CHR)), alpha=0.8,size=1.4) +
  scale_color_manual(values = rep(c("black", "darkgoldenrod"), 22 )) +
  geom_hline(yintercept=10, linetype="dashed",color = "red",lty = 2, lwd = 1) +
  geom_hline(yintercept=30, linetype="dashed",color = "deepskyblue", lty = 2, lwd = 1) +
  scale_y_continuous (breaks = c(-20,-10,0,10,20,30,40,50,60), limits=c(-20,60))+
  labs( y = "-log10(BF)", x = "Position (Mbp)")+
  theme_classic() +
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.8),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.4, face="bold",hjust=0.5))
ggsave("Graphs/Manhattan_BF/bf_manhattan_6.png",width=13, height=4.5)


#ENV 7 - Tave_sm 
ggplot(env7_united, aes(x=BP/1e6, y=BF)) +
  geom_point(aes(color=as.factor(CHR)), alpha=0.8,size=1.4) +
  scale_color_manual(values = rep(c("black", "darkgoldenrod"), 22 )) +
  geom_hline(yintercept=10, linetype="dashed",color = "red",lty = 2, lwd = 1) +
  geom_hline(yintercept=30, linetype="dashed",color = "deepskyblue", lty = 2, lwd = 1) +
  scale_y_continuous (breaks = c(-20,-10,0,10,20,30,40,50,60), limits=c(-20,60))+
  labs( y = "-log10(BF)", x = "Position (Mbp)")+
  theme_classic() +
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.8),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.4, face="bold",hjust=0.5))
ggsave("Graphs/Manhattan_BF/bf_manhattan_7.png",width=13, height=4.5)


#ENV 8 - PPT_wt 
ggplot(env8_united, aes(x=BP/1e6, y=BF)) +
  geom_point(aes(color=as.factor(CHR)), alpha=0.8,size=1.4) +
  scale_color_manual(values = rep(c("black", "deepskyblue"), 22 )) +
  geom_hline(yintercept=10, linetype="dashed",color = "red",lty = 2, lwd = 1) +
  geom_hline(yintercept=30, linetype="dashed",color = "deepskyblue", lty = 2, lwd = 1) +
  scale_y_continuous (breaks = c(-20,-10,0,10,20,30,40,50,60), limits=c(-20,60))+
  labs( y = "-log10(BF)", x = "Position (Mbp)")+
  theme_classic() +
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.8),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.4, face="bold",hjust=0.5))
ggsave("Graphs/Manhattan_BF/bf_manhattan_8.png",width=13, height=4.5)



#ENV 9 - PPT_sm 
ggplot(env9_united, aes(x=BP/1e6, y=BF)) +
  geom_point(aes(color=as.factor(CHR)), alpha=0.8,size=1.4) +
  scale_color_manual(values = rep(c("black", "deepskyblue"), 22 )) +
  geom_hline(yintercept=10, linetype="dashed",color = "red",lty = 2, lwd = 1) +
  geom_hline(yintercept=30, linetype="dashed",color = "deepskyblue", lty = 2, lwd = 1) +
  scale_y_continuous (breaks = c(-20,-10,0,10,20,30,40,50,60), limits=c(-20,60))+
  labs( y = "-log10(BF)", x = "Position (Mbp)")+
  theme_classic() +
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.8),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.4, face="bold",hjust=0.5))
ggsave("Graphs/Manhattan_BF/bf_manhattan_9.png",width=13, height=4.5)





















#ENV 2 - MAP
ggplot(env2_united, aes(x=BP/1e6, y=BF)) +
  geom_point(aes(color=as.factor(CHR)), alpha=0.8,size=0.5) +
  scale_color_manual(values = rep(c("black", "skyblue"), 22 )) +
  geom_hline(yintercept=10, linetype="dashed",color = "red", size=2) +
  geom_hline(yintercept=30, linetype="dashed",color = "skyblue", size=2) +
  scale_x_continuous("Position (Mbp)")+
  scale_y_continuous("-log10(BF)", breaks = c(-20,-10,0,10,20,30,40)) +
  theme_man() 
ggsave("Graphs/Manhattan_BF/bf_manhattan_2.png",width=9, height=6)






#ENV 2 - MAP 
axisdf_2 <- env2_united %>% group_by(CHR) %>% summarize(center=( max(BP) + min(BP) ) / 2 )

ggplot(env2_united, aes(x=BP, y=BF)) +
  
  # Show all points
  geom_point(aes(color=as.factor(CHR)), alpha=0.8,size=0.5) +
  scale_color_manual(values = rep(c("black", "deepskyblue"), 22 )) +
  geom_hline(yintercept=10, linetype="dashed",color = "red", size=2) +
  # custom X axis:
  scale_x_continuous(expand = c(0, 0), label = axisdf_2$CHR, breaks= axisdf_2$center) +
  scale_y_continuous(breaks = c(-20,-10,0,10,20,30,40)) +
  theme_classic() +
  labs(
    y = "Bayes Factor",
    x = "Position")+
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.8),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size=25, face="bold"),
    axis.text.y = element_text(size=25,face="bold"),
    axis.title.x = element_text(color="black", size=0, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=0,vjust = 2, face="bold",hjust=0.5)
  )

#Export. File is very large and somewhat hard to hande/open. 
ggsave("Graphs/Manhattan_BF/graph_manhattan_2.pdf",width=9, height=6)




#ENV 5 - CMD 
axisdf_5 <- env5_united %>% group_by(CHR) %>% summarize(center=( max(BP) + min(BP) ) / 2 )

ggplot(env5_united, aes(x=BP, y=BF)) +
  
  # Show all points
  geom_point(aes(color=as.factor(CHR)), alpha=0.8,size=0.5) +
  scale_color_manual(values = rep(c("black", "magenta3"), 22 )) +
  geom_hline(yintercept=10, linetype="dashed",color = "blue", size=2) +
  # custom X axis:
  scale_x_continuous(expand = c(0, 0), label = axisdf_5$CHR, breaks= axisdf_5$center) +
  scale_y_continuous(breaks = c(-20,-10,0,10,20,30,40)) +
  theme_classic() +
  labs(
    y = "Bayes Factor",
    x = "Position")+
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.8),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size=25, face="bold"),
    axis.text.y = element_text(size=25,face="bold"),
    axis.title.x = element_text(color="black", size=0, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=0,vjust = 2, face="bold",hjust=0.5)
  )

#Export. File is very large and somewhat hard to hande/open. 
ggsave("/Users/daniel_anstett/Dropbox/AM_Workshop/large_graphs/Manhattan/graph_manhattan_5.pdf",width=9, height=6)
































