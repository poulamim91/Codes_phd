rm(list=ls())

library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)

data_new = read.csv(some_file, header = T)

data_new <- data_new %>% mutate(key6 = paste0(key2,key4), key5 = paste0(key2,key3))
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", 
               "#CC79A7")

val <- c("Absence of Toilet (95% CI)", "Absence of Toilet (Estimate)", "Presence of Toilet (95% CI)",
         "Presence of Toilet (Estimate)")
gg <- c("Age of Children is > 12 months", "Age of Children is <= 12 months", "")

data_new1 <- data_new %>% 
  filter(grepl("mu00|mu01", key)) %>% 
  select(R_u, key6, value) %>%
  spread(key6, value) %>%
  as.data.frame()

data_new1 = data_new1[,c(1,4,2,3,7,5,6)]

lt=c(4,4,1,3,3,2)
temp<-data.frame()
p1<-ggplot(temp) + geom_point()+
geom_line(aes(x=data_new1$R_u,y=data_new1$group2notCI,linetype="Presence of Toilet (Estimate)"))+
geom_line(aes(x=data_new1$R_u,y=data_new1$group1notCI,linetype="Absence of Toilet (Estimate)"))+
geom_line(aes(x=data_new1$R_u,y=data_new1$group2CI_LL,linetype="Presence of Toilet (95% CI)"))+
geom_line(aes_string(x=data_new1$R_u,y=data_new1$group2CI_UL),linetype=3)+
geom_line(aes(x=data_new1$R_u,y=data_new1$group1CI_LL,linetype="Absence of Toilet (95% CI)"))+
geom_line(aes_string(x=data_new1$R_u,y=data_new1$group1CI_UL),linetype=4)+
ggtitle("Age of Children is > 12 months")+
scale_x_continuous(name = "Time (days)") +
scale_y_continuous(name = "Mean Function Estimates", limits = c(0, 1.5)) +
theme(panel.background = element_rect(fill = "white",colour = "black",
                                        size = 0.5, linetype = "solid"),
        axis.text=element_text(size=20, face = "bold"), 
        axis.title=element_text(size=20, face = "bold"),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
      legend.text=element_text(size=17.5),
      legend.title=element_blank(), 
      legend.position = c(0.27,0.912))+
  theme(panel.border= element_blank())+
  theme(axis.line.x = element_line(color="black", size = 1.2),
        axis.line.y = element_line(color="black", size = 1.2),
        legend.key.width = unit(5, "line"))+
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
  #scale_color_manual(values=c("Y1"="black","Y2"="black","Y3"="black",
   #                             "Y4"="black","Y5"="black","Y6"="black"))+
  scale_linetype_manual(name = "",
                        values=c("Absence of Toilet (Estimate)"=1,
                                 "Presence of Toilet (Estimate)"=2,
                                 "Absence of Toilet (95% CI)"=4,
                                 "Presence of Toilet (95% CI)"=3),
                        breaks=c("Absence of Toilet (Estimate)",
                                 "Absence of Toilet (95% CI)",
                                 "Presence of Toilet (Estimate)",
                                 "Presence of Toilet (95% CI)"))

data_new2 <- data_new %>% 
  filter(grepl("mu10|mu11", key)) %>% 
  select(R_u, key6, value) %>%
  spread(key6, value) %>%
  as.data.frame()


p2<-ggplot(temp) + geom_point()+
  geom_line(aes(x=data_new2$R_u,y=data_new2$group4notCI,linetype="Presence of Toilet (Estimate)"))+
  geom_line(aes(x=data_new2$R_u,y=data_new2$group3notCI,linetype="Absence of Toilet (Estimate)"))+
  geom_line(aes(x=data_new2$R_u,y=data_new2$group4CI_LL,linetype="Presence of Toilet (95% CI)"))+
  geom_line(aes_string(x=data_new2$R_u,y=data_new2$group4CI_UL),linetype=3)+
  geom_line(aes(x=data_new2$R_u,y=data_new2$group3CI_LL,linetype="Absence of Toilet (95% CI)"))+
  geom_line(aes_string(x=data_new2$R_u,y=data_new2$group3CI_UL),linetype=4)+
  ggtitle("Age of Children is <= 12 months")+
  scale_x_continuous(name = "Time (days)") +
  scale_y_continuous(name = "Mean Function Estimates", limits = c(0, 1.5)) +
  theme(panel.background = element_rect(fill = "white",colour = "black",
                                        size = 0.5, linetype = "solid"),
        axis.text=element_text(size=20, face = "bold"), 
        axis.title=element_text(size=20, face = "bold"),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
        legend.text=element_text(size=17.5),
        legend.title=element_blank(), 
        legend.position = c(0.27,0.912))+
  theme(panel.border= element_blank())+
  theme(axis.line.x = element_line(color="black", size = 1.2),
        axis.line.y = element_line(color="black", size = 1.2),
        legend.key.width = unit(5, "line"))+
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
  #scale_color_manual(values=c("Y1"="black","Y2"="black","Y3"="black",
  #                             "Y4"="black","Y5"="black","Y6"="black"))+
  scale_linetype_manual(name = "",
                        values=c("Absence of Toilet (Estimate)"=1,
                                 "Presence of Toilet (Estimate)"=2,
                                 "Absence of Toilet (95% CI)"=4,
                                 "Presence of Toilet (95% CI)"=3),
                        breaks=c("Absence of Toilet (Estimate)",
                                 "Absence of Toilet (95% CI)",
                                 "Presence of Toilet (Estimate)",
                                 "Presence of Toilet (95% CI)"))

postscript("data_example_ratefn_new.eps", horizontal=F, onefile = F,
           paper = "special", width = 18, height = 8)
print(grid.arrange(p1,p2,nrow=1))
dev.off()

pdf("data_example_ratefn_new.pdf", width = 18, height = 10)
print(grid.arrange(p1,p2,nrow=1))
dev.off()