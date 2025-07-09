##################################################
### Phenotypic trait change over time
## Break down assocaitions between Genomic strength of selection and trait data
##################################################
#Import libraries
library(tidyverse)
library(MASS)
library(sfsmisc)
library(egg)

#Import Data
demo_pop <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv") %>% filter(Paper_ID<12)
trait_pop <- read_csv("data/trait_data/trait_pop.csv")
trait_pop_graph <- cbind(demo_pop,trait_pop)  %>% dplyr::select(Latitude,Site,Paper_ID, Lat.Color,Median,SLA_D,FT_D,A_D,SC_D,WC_D)

###################################################################################
# Stat testing looking at the combined effect of Assimilation and Stomatal Conducance

mod.dry.gasx <- lm(Median ~ A_D + SC_D, data=trait_pop)
summary(mod.dry.gasx) #positive effects of evolution towards greater carbon assimilation, decreased stomatal conductance
rob.mod.dry.gasx <- rlm(Median ~ A_D + SC_D, data=trait_pop)
f.robftest(rob.mod.dry.gasx, var="A_D") #*
f.robftest(rob.mod.dry.gasx, var="SC_D") #***


###################################################################################
# Graphs of all 5 traits

###################################################################################
trait_pop_graph$Lat.Color<-as.factor(trait_pop_graph$Lat.Color)
trait_pop_graph$Lat.Color<-factor(trait_pop_graph$Lat.Color,levels=trait_pop_graph$Lat.Color)

a <- ggplot(trait_pop_graph, aes(x=SLA_D, y=Median)) + 
  geom_smooth(method=lm,color="black", size=1.8, linetype="dashed", fill="gray75")+
  geom_smooth(method=MASS::rlm, color="black", size=1.8, fill="grey50") +
  #geom_smooth(data=trait_pop_graph_cull1, aes(x=Median, y=log(mean.lambda.recovery+0.5)), method=lm, color="black", size=1.8, linetype="longdash", fill="grey35") +
  geom_point(aes(fill=Lat.Color), shape=21, size =6)+
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.7) +
  scale_x_continuous(name="SLA")+#,breaks=c(0.5,1,1.5,2,2.5))+
  scale_y_continuous(name="Median S")+#,breaks=c(-0.1,0,0.1,0.2))+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(name = "Latitude (°N)",labels=round(trait_pop_graph$Latitude,1),
                    values=as.character(trait_pop_graph$Lat.Color)) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=16, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=16,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines"),
    legend.position = "none"#Reduce height
  )

b <- ggplot(trait_pop_graph, aes(x=FT_D, y=Median)) + 
  geom_smooth(method=lm,color="black", size=1.8, linetype="dashed", fill="gray75")+
  geom_smooth(method=MASS::rlm, color="black", size=1.8, fill="grey50") +
  #geom_smooth(data=trait_pop_graph_cull1, aes(x=Median, y=log(mean.lambda.recovery+0.5)), method=lm, color="black", size=1.8, linetype="longdash", fill="grey35") +
  geom_point(aes(fill=Lat.Color), shape=21, size =6)+
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.7) +
  scale_x_continuous(name="Date of First Flower")+#,breaks=c(0.5,1,1.5,2,2.5))+
  scale_y_continuous(name="Median S")+#,breaks=c(-0.1,0,0.1,0.2))+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(name = "Latitude (°N)",labels=round(trait_pop_graph$Latitude,1),
                    values=as.character(trait_pop_graph$Lat.Color)) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=16, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=16,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines"), #Reduce height
    legend.position = "none"#Reduce height
  )

c <- ggplot(trait_pop_graph, aes(x=A_D, y=Median)) + 
  geom_smooth(method=lm,color="black", size=1.8, linetype="dashed", fill="gray75")+
  geom_smooth(method=MASS::rlm, color="black", size=1.8, fill="grey50") +
  #geom_smooth(data=trait_pop_graph_cull1, aes(x=Median, y=log(mean.lambda.recovery+0.5)), method=lm, color="black", size=1.8, linetype="longdash", fill="grey35") +
  geom_point(aes(fill=Lat.Color), shape=21, size =6)+
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.7) +
  scale_x_continuous(name="Carbon Assimilation")+#,breaks=c(0.5,1,1.5,2,2.5))+
  scale_y_continuous(name="Median S")+#,breaks=c(-0.1,0,0.1,0.2))+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(name = "Latitude (°N)",labels=round(trait_pop_graph$Latitude,1),
                    values=as.character(trait_pop_graph$Lat.Color)) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=16, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=16,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines"), #Reduce height
    legend.position = "none"#Reduce height
  )

d <- ggplot(trait_pop_graph, aes(x=SC_D, y=Median)) + 
  geom_smooth(method=lm,color="black", size=1.8, linetype="dashed", fill="gray75")+
  geom_smooth(method=MASS::rlm, color="black", size=1.8, fill="grey50") +
  #geom_smooth(data=trait_pop_graph_cull1, aes(x=Median, y=log(mean.lambda.recovery+0.5)), method=lm, color="black", size=1.8, linetype="longdash", fill="grey35") +
  geom_point(aes(fill=Lat.Color), shape=21, size =6)+
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.7) +
  scale_x_continuous(name="Stomatal Conductance")+#,breaks=c(0.5,1,1.5,2,2.5))+
  scale_y_continuous(name="Median S")+#,breaks=c(-0.1,0,0.1,0.2))+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(name = "Latitude (°N)",labels=round(trait_pop_graph$Latitude,1),
                    values=as.character(trait_pop_graph$Lat.Color)) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=16, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=16,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce height
  )

e <- ggplot(trait_pop_graph, aes(x=Median, y=WC_D)) + 
  geom_smooth(method=lm,color="black", size=1.8, linetype="dashed", fill="gray75")+
  geom_smooth(method=MASS::rlm, color="black", size=1.8, fill="grey50") +
  #geom_smooth(data=trait_pop_graph_cull1, aes(x=Median, y=log(mean.lambda.recovery+0.5)), method=lm, color="black", size=1.8, linetype="longdash", fill="grey35") +
  geom_point(aes(fill=Lat.Color), shape=21, size =6)+
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.7) +
  scale_x_continuous(name="Water Content")+#,breaks=c(0.5,1,1.5,2,2.5))+
  scale_y_continuous(name="Median S")+#,breaks=c(-0.1,0,0.1,0.2))+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(name = "Latitude (°N)",labels=round(trait_pop_graph$Latitude,1),
                    values=as.character(trait_pop_graph$Lat.Color)) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=16, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=16,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines"),
    legend.position = "none"#Reduce height#Reduce height
  )

trait_5 <- ggarrange(a,b,c,d,e)
trait_2 <- ggarrange(c,d,ncol=2)

ggsave("Graphs/Traits/trait_5.pdf",trait_5, width=10, height = 12, units = "in")

ggsave("Graphs/Traits/trait_2.pdf",trait_2,width=10, height = 5, units = "in")












