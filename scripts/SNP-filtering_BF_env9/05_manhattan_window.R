#############################################################################################################
## Plot Manhattan for WZA windows
## Author Daniel Anstett
## 
## 
## Modified from Tom Booker WZA Vignette
## Last Modified August 3, 2022
#############################################################################################################
#Import libraries

library(tidyverse)
library(Kendall)

#Import files

#Windows
wza_win_env1 <- read_csv("data/genomic_data/WZA_win_env1_bf.csv")
wza_win_env2 <- read_csv("data/genomic_data/WZA_win_env2_bf.csv")
wza_win_env3 <- read_csv("data/genomic_data/WZA_win_env3_bf.csv")
wza_win_env4 <- read_csv("data/genomic_data/WZA_win_env4_bf.csv")
wza_win_env5 <- read_csv("data/genomic_data/WZA_win_env5_bf.csv")

wza_win_env6 <- read_csv("data/genomic_data/WZA_win_env6_bf.csv")
wza_win_env7 <- read_csv("data/genomic_data/WZA_win_env7_bf.csv")
wza_win_env8 <- read_csv("data/genomic_data/WZA_win_env8_bf.csv")
wza_win_env9 <- read_csv("data/genomic_data/WZA_win_env9_bf.csv")

###########################################################################################################

# WZA Empirical p-value

#ENV1
wza_empri_mat <- ggplot(data = wza_win_env1, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  scale_y_continuous("-log10(p-value)", limits=c(0,10))+
  scale_x_continuous("Position (Mbp)")+
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
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.4, face="bold",hjust=0.5)
  )
wza_empri_mat
ggsave("Graphs/WZA/wza_env1_bf.png", wza_empri_mat, width=13, height = 4.5, units = "in")


#ENV2
wza_empri_map <- ggplot(data = wza_win_env2, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  scale_y_continuous("-log10(p-value)", limits=c(0,10))+
  scale_x_continuous("Position (Mbp)")+
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
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.4, face="bold",hjust=0.5)
  )
wza_empri_map
ggsave("Graphs/WZA/wza_env2_bf.png", wza_empri_map, width=13, height = 4.5, units = "in")

#ENV3
wza_empri_map <- ggplot(data = wza_win_env3, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  scale_y_continuous("-log10(p-value)", limits=c(0,10))+
  scale_x_continuous("Position (Mbp)")+
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
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.4, face="bold",hjust=0.5)
  )
wza_empri_map
ggsave("Graphs/WZA/wza_env3_bf.png", wza_empri_map, width=13, height = 4.5, units = "in")

#ENV4
wza_empri_map <- ggplot(data = wza_win_env4, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  scale_y_continuous("-log10(p-value)", limits=c(0,10))+
  scale_x_continuous("Position (Mbp)")+
  geom_hline(aes(yintercept = -log10(0.05/dim(wza_win_env4)[1])), col = "red", lty = 2, lwd = 1)+
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
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.4, face="bold",hjust=0.5)
  )
wza_empri_map
ggsave("Graphs/WZA/wza_env4_bf.png", wza_empri_map, width=13, height = 4.5, units = "in")


#ENV5
wza_empri_cmd <- ggplot(data = wza_win_env5, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  scale_y_continuous("-log10(p-value)", limits=c(0,10))+
  scale_x_continuous("Position (Mbp)")+
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
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.4, face="bold",hjust=0.5)
  )
wza_empri_cmd
ggsave("Graphs/WZA/wza_env5_bf.png", wza_empri_cmd, width=13, height = 4.5, units = "in")

#ENV6
wza_empri_map <- ggplot(data = wza_win_env6, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  scale_y_continuous("-log10(p-value)", limits=c(0,10))+
  scale_x_continuous("Position (Mbp)")+
  geom_hline(aes(yintercept = -log10(0.05/dim(wza_win_env6)[1])), col = "red", lty = 2, lwd = 1)+
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
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.4, face="bold",hjust=0.5)
  )
wza_empri_map
ggsave("Graphs/WZA/wza_env6_bf.png", wza_empri_map, width=13, height = 4.5, units = "in")

#ENV7
wza_empri_map <- ggplot(data = wza_win_env7, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  scale_y_continuous("-log10(p-value)", limits=c(0,10))+
  scale_x_continuous("Position (Mbp)")+
  geom_hline(aes(yintercept = -log10(0.05/dim(wza_win_env7)[1])), col = "red", lty = 2, lwd = 1)+
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
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.4, face="bold",hjust=0.5)
  )
wza_empri_map
ggsave("Graphs/WZA/wza_env7_bf.png", wza_empri_map, width=13, height = 4.5, units = "in")

#ENV8
wza_empri_map <- ggplot(data = wza_win_env8, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  scale_y_continuous("-log10(p-value)", limits=c(0,10))+
  scale_x_continuous("Position (Mbp)")+
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
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.4, face="bold",hjust=0.5)
  )
wza_empri_map
ggsave("Graphs/WZA/wza_env8_bf.png", wza_empri_map, width=13, height = 4.5, units = "in")

#ENV9
wza_empri_map <- ggplot(data = wza_win_env9, aes( x = pos/1e6, y = -log10(Z_pVal)))+
  geom_point(aes(color=as.factor(chr), alpha=0.9))+
  scale_y_continuous("-log10(p-value)", limits=c(0,10))+
  scale_x_continuous("Position (Mbp)")+
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
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.4, face="bold",hjust=0.5)
  )
wza_empri_map
ggsave("Graphs/WZA/wza_env9_bf.png", wza_empri_map, width=13, height = 4.5, units = "in")

