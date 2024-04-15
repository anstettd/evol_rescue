#############################################################################################################
## #View empirical p-value plots using plotly
## Author Daniel Anstett
## 
## 
## 
## Last Modified March 14, 2023
#############################################################################################################
#Import libraries

library(tidyverse)
library(Kendall)
library(plotly)

#Import files

#Windows
#wza_win_mat <- read_csv("data/genomic_data/WZA_win_mat_bf.csv")
wza_win_env1 <- read_csv("data/genomic_data/WZA_win_env1_bf.csv")
wza_win_env2 <- read_csv("data/genomic_data/WZA_win_env2_bf.csv")
wza_win_env3 <- read_csv("data/genomic_data/WZA_win_env3_bf.csv")
wza_win_env4 <- read_csv("data/genomic_data/WZA_win_env4_bf.csv")
wza_win_env5 <- read_csv("data/genomic_data/WZA_win_env5_bf.csv")
wza_win_env6 <- read_csv("data/genomic_data/WZA_win_env6_bf.csv")
wza_win_env7 <- read_csv("data/genomic_data/WZA_win_env7_bf.csv")
wza_win_env8 <- read_csv("data/genomic_data/WZA_win_env8_bf.csv")
wza_win_env9 <- read_csv("data/genomic_data/WZA_win_env9_bf.csv")

##############################################################################################################
##############################################################################################################
#View empirical p-value plots using plotly

#Annual

#ENV 1 - MAT
wza_empri_mat <- ggplot(data = wza_win_env1, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  geom_line()+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)",limits=c(64.5,64.7))+ #MODIFY RANGE HERE TO VIEW DATA UP CLOSE
  geom_hline(aes(yintercept = -log10(0.05/dim(wza_win_env1)[1])), col = "red", lty = 2, lwd = 1)+
  scale_color_manual(values = rep(c("black", "darkgoldenrod"), 22 )) +
  theme_classic()+
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
    axis.title.y = element_text(color="black", size=18,vjust = 1.4, face="bold",hjust=0.5)
  ) 
#wza_empri_mat
ggplotly(wza_empri_mat)
##############################################################################################################

#ENV2 - MAP
wza_empri_map <- ggplot(data = wza_win_env2, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  geom_line()+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)",limits=c(38,39))+ #MODIFY RANGE HERE TO VIEW DATA UP CLOSE
  geom_hline(aes(yintercept = -log10(0.05/dim(wza_win_env2)[1])), col = "red", lty = 2, lwd = 1)+
  scale_color_manual(values = rep(c("black", "deepskyblue"), 22 )) +
  theme_classic()+
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
    axis.title.y = element_text(color="black", size=18,vjust = 1.4, face="bold",hjust=0.5)
  ) 
#wza_empri_map
ggplotly(wza_empri_map)
##############################################################################################################

#ENV3
wza_empri_env3 <- ggplot(data = wza_win_env3, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  geom_line()+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)",limits=c(377.0,377.2))+ #MODIFY RANGE HERE TO VIEW DATA UP CLOSE
  geom_hline(aes(yintercept = -log10(0.05/dim(wza_win_env3)[1])), col = "red", lty = 2, lwd = 1)+
  scale_color_manual(values = rep(c("black", "deepskyblue"), 22 )) +
  theme_classic()+
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
    axis.title.y = element_text(color="black", size=18,vjust = 1.4, face="bold",hjust=0.5)
  ) 
#wza_empri_map
ggplotly(wza_empri_env3)
##############################################################################################################

#ENV4
wza_empri_env4 <- ggplot(data = wza_win_env4, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  geom_line()+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)",limits=c(108.5,108.55))+ #MODIFY RANGE HERE TO VIEW DATA UP CLOSE
  geom_hline(aes(yintercept = -log10(0.05/dim(wza_win_env4)[1])), col = "red", lty = 2, lwd = 1)+
  scale_color_manual(values = rep(c("black", "deepskyblue"), 22 )) +
  theme_classic()+
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
    axis.title.y = element_text(color="black", size=18,vjust = 1.4, face="bold",hjust=0.5)
  ) 
#wza_empri_map
ggplotly(wza_empri_env4)
##############################################################################################################

#ENV 5 - CMD
wza_empri_cmd <- ggplot(data = wza_win_env5, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  geom_line()+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)",limits=c(64.5,64.55))+ #MODIFY RANGE HERE TO VIEW DATA UP CLOSE
  geom_hline(aes(yintercept = -log10(0.05/dim(wza_win_env5)[1])), col = "red", lty = 2, lwd = 1)+
  scale_color_manual(values = rep(c("black", "magenta3"), 22 )) +
  theme_classic()+
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
    axis.title.y = element_text(color="black", size=18,vjust = 1.4, face="bold",hjust=0.5)
  ) 
#wza_empri_cmd
ggplotly(wza_empri_cmd)
##############################################################################################################

#Seasonal
#ENV6
wza_empri_env6 <- ggplot(data = wza_win_env6, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  geom_line()+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)",limits=c(99.35,99.41))+ #MODIFY RANGE HERE TO VIEW DATA UP CLOSE
  geom_hline(aes(yintercept = -log10(0.05/dim(wza_win_env6)[1])), col = "red", lty = 2, lwd = 1)+
  scale_color_manual(values = rep(c("black", "deepskyblue"), 22 )) +
  theme_classic()+
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
    axis.title.y = element_text(color="black", size=18,vjust = 1.4, face="bold",hjust=0.5)
  ) 
#wza_empri_map
ggplotly(wza_empri_env6)
##############################################################################################################

#ENV7
wza_empri_env7 <- ggplot(data = wza_win_env7, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  geom_line()+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)",limits=c(64.60,64.70))+ #MODIFY RANGE HERE TO VIEW DATA UP CLOSE
  geom_hline(aes(yintercept = -log10(0.05/dim(wza_win_env7)[1])), col = "red", lty = 2, lwd = 1)+
  scale_color_manual(values = rep(c("black", "deepskyblue"), 22 )) +
  theme_classic()+
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
    axis.title.y = element_text(color="black", size=18,vjust = 1.4, face="bold",hjust=0.5)
  ) 
#wza_empri_map
ggplotly(wza_empri_env7)
##############################################################################################################

#ENV8
wza_empri_env8 <- ggplot(data = wza_win_env8, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  geom_line()+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)",limits=c(64.5,64.6))+ #MODIFY RANGE HERE TO VIEW DATA UP CLOSE
  geom_hline(aes(yintercept = -log10(0.05/dim(wza_win_env8)[1])), col = "red", lty = 2, lwd = 1)+
  scale_color_manual(values = rep(c("black", "deepskyblue"), 22 )) +
  theme_classic()+
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
    axis.title.y = element_text(color="black", size=18,vjust = 1.4, face="bold",hjust=0.5)
  ) 
#wza_empri_map
ggplotly(wza_empri_env8)
##############################################################################################################

#ENV9
wza_empri_env9 <- ggplot(data = wza_win_env9, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  geom_line()+
  scale_y_continuous("-log10(WZA Empirical p-value)", limits=c(0,20))+
  scale_x_continuous("Position (Mbp)",limits=c(38.78,38.85))+ #MODIFY RANGE HERE TO VIEW DATA UP CLOSE
  geom_hline(aes(yintercept = -log10(0.05/dim(wza_win_env9)[1])), col = "red", lty = 2, lwd = 1)+
  scale_color_manual(values = rep(c("black", "deepskyblue"), 22 )) +
  theme_classic()+
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
    axis.title.y = element_text(color="black", size=18,vjust = 1.4, face="bold",hjust=0.5)
  ) 
#wza_empri_map
ggplotly(wza_empri_env9)


