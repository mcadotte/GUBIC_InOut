###plots Sorensen for native and exotic in cities and buffers
setwd("~/Dropbox/2018WorkingFiles/Marc2018/GUBIC/sGUBIK")
source("In_Out/R/data_import.R")

all<-read.csv("In_Out/data/Sor_buff+urb.csv",row.names=1)
N<-all[all$Taxa=="Native",]
I<-all[all$Taxa=="Intro",]
N_I_both<-rbind(N, I)

#quartz()


NIB<-ggplot(N_I_both, aes(x=Spat,y=Sor,color=Taxa))+
  geom_point(data = subset(N_I_both, Taxa == "Native"),
             aes(x = Spat, y = Sor), color = "darkolivegreen4", alpha = 0.6, size = 2) +
  geom_point(data = subset(N_I_both, Taxa == "Intro"),
             aes(x = Spat, y = Sor), color = "darkorange3", alpha = 0.6, size = 2) +
  geom_smooth(data = subset(N_I_both, Taxa == "Native"),
              aes(x = Spat, y = Sor), color = "#243b29", se = FALSE, linewidth = 1.4) +
  geom_smooth(data = subset(N_I_both, Taxa == "Intro"),
              aes(x = Spat, y = Sor), color = "darkorange4", se = FALSE, linewidth = 1.4) +
  facet_wrap(~Taxa)+
  theme_minimal(base_size = 14) +
  theme(axis.text.x=element_text(size=11,angle=45,hjust=1))+
  labs(x = "Spatial distance", y = "Sorensen dissimilarity",
       title="A) Urban + 5 km buffer")  
#NIB

all.buff<-read.csv("In_Out/data/Sor_buff.csv",row.names=1)
N<-all.buff[all.buff$Taxa=="Native",]
I<-all.buff[all.buff$Taxa=="Intro",]
N_I_buff<-rbind(N, I)

NI5<-ggplot(N_I_buff, aes(x=Spat,y=Sor,color=Taxa))+
  geom_point(data = subset(N_I_buff, Taxa == "Native"),
             aes(x = Spat, y = Sor), color = "darkolivegreen4", alpha = 0.6, size = 2) +
  geom_point(data = subset(N_I_buff, Taxa == "Intro"),
             aes(x = Spat, y = Sor), color = "darkorange3", alpha = 0.6, size = 2) +
  geom_smooth(data = subset(N_I_buff, Taxa == "Native"),
              aes(x = Spat, y = Sor), color = "#243b29", se = FALSE, linewidth = 1.4) +
  geom_smooth(data = subset(N_I_buff, Taxa == "Intro"),
              aes(x = Spat, y = Sor), color = "darkorange4", se = FALSE, linewidth = 1.4) +
  facet_wrap(~Taxa)+
  theme_minimal(base_size = 14) +
  theme(axis.text.x=element_text(size=11,angle=45,hjust=1))+
  labs(x = "Spatial distance", y = "Sorensen dissimilarity",
       title="C) 5 km buffer")  
#NI5

all.urb<-read.csv("In_Out/data/Sor_urb.csv",row.names=1)
N<-all.urb[all.urb$Taxa=="Native",]
I<-all.urb[all.urb$Taxa=="Intro",]
N_I_urb<-rbind(N, I)

NIU<-ggplot(N_I_urb, aes(x=Spat,y=Sor,color=Taxa))+
  geom_point(data = subset(N_I_urb, Taxa == "Native"),
             aes(x = Spat, y = Sor), color = "darkolivegreen4", alpha = 0.6, size = 2) +
  geom_point(data = subset(N_I_urb, Taxa == "Intro"),
             aes(x = Spat, y = Sor), color = "darkorange3", alpha = 0.6, size = 2) +
  geom_smooth(data = subset(N_I_urb, Taxa == "Native"),
              aes(x = Spat, y = Sor), color = "#243b29", se = FALSE, linewidth = 1.4) +
  geom_smooth(data = subset(N_I_urb, Taxa == "Intro"),
              aes(x = Spat, y = Sor), color = "darkorange4", se = FALSE, linewidth = 1.4) +
  facet_wrap(~Taxa)+
  theme_minimal(base_size = 14) +
  theme(axis.text.x=element_text(size=11,angle=45,hjust=1))+
  labs(x = "Spatial distance", y = "Sorensen dissimilarity",
       title="B) Urban ")  
#NIU

quartz()
grid.arrange(NIB,NIU,NI5, nrow=1)
