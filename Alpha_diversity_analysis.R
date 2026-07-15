setwd("~/Dropbox/2018WorkingFiles/Marc2018/GUBIC/sGUBIK")
source("In_Out/R/data_import.R")
library(rethinking)

all.out<-readRDS("In_Out/data/R_Phy_Trait_output.RData")

delta.PD.n<-all.out$Metrics.cities

t.test(delta.PD.n$Nat.D.rich)
t.test(delta.PD.n$Nat.D.MPD)
t.test(delta.PD.n$Inv.D.rich)
t.test(delta.PD.n$Inv.D.MPD)


quartz()
dens(delta.PD.n$Inv.D.rich,lwd=3,col="darkorange3")
dens(delta.PD.n$Nat.D.rich,lwd=3,col="darkolivegreen4",add=TRUE)

dens(delta.PD.n$Inv.D.MPD,lwd=3,col="darkorange3")
dens(delta.PD.n$Nat.D.MPD,lwd=3,col="darkolivegreen4",add=TRUE)

sum.stat<-function(x){
  n<-length(x)
  mean<-mean(x,na.rm=TRUE)
  SE<-sd(x,na.rm=TRUE)/sqrt(length(x))
  CI95.L<-mean - qt(0.975, df = n - 1) * SE
  CI95.U<-mean + qt(0.975, df = n - 1) * SE
  return(data.frame(n,mean,SE,CI95.L,CI95.U))
}

sum.dat<-do.call(rbind,apply(delta.PD.n,2,sum.stat))

metric<-rep(c("Richness","PD","MPD", "FDis", "FDiv", "FEve"),2)
status<-c(rep("Introduced",6),rep("Native",6))

sum.PD<-cbind(metric, status, sum.dat)
#write.csv(sum.PD,"In_Out/data/In_Out_allMetrics.csv")


##test difference from 0
Ps<-apply(sum.PD[,3:5],1,function(x) 2*pt(-abs(x[2]/x[3]),x[1]-1))
sum.PDs<-sum.PD[sum.PD$metric!="PD" & sum.PD$metric!="FDiv",]
sum.PDs<-sum.PDs[sum.PDs$metric!="FDis" & sum.PDs$metric!="FEve",]
#write.csv(sum.PDs,"In_Out/data/In_Out_main4Metrics.csv")


out_metrics<-ggplot(sum.PDs, aes(x = mean, y = status, colour = status)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = CI95.L, xmax = CI95.U), height = 0.2,linewidth = 1) +
  facet_wrap(~ metric, scales="free_x") +
  scale_colour_manual(values = c("Introduced" = "darkorange3", "Native" = "darkolivegreen4")) +
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(size = 16, face = "bold")) +
  labs(x = "Mean ± SE", y = "Non-native/Native",colour = "status") +
  geom_vline(xintercept = 0,linetype = "dashed",colour = "grey45",linewidth = 1.2)+
  
  geom_segment( data = data.frame(metric = "MPD"),
                aes(x = 0, xend = -3.5, y = 2.4, yend = 2.4),
                arrow = arrow(length = unit(0.4, "cm")),linewidth = 1.2,
                inherit.aes = FALSE,colour = "grey45" ) +
  
  geom_text(data = data.frame(metric = "MPD"),
            aes(x = -1.75, y = 2.25),  # midpoint of arrow, slightly below
            label = "Greater in buffer",  inherit.aes = FALSE,
            size = 4, colour = "grey45" ) +
  
  # geom_segment( data = data.frame(metric = "FDis"),
  #   aes(x = 0, xend = 0.06, y = 2.4, yend = 2.4),
  #   arrow = arrow(length = unit(0.4, "cm")),linewidth = 1.2,
  #   inherit.aes = FALSE,colour = "grey45" ) +
  
  #  geom_text(data = data.frame(metric = "FDis"),
  #   aes(x = 0.03, y = 2.25),  # midpoint of arrow, slightly below
  #   label = "Greater in urban",  inherit.aes = FALSE,
  #   size = 4, colour = "grey45" ) +
  
  geom_text(data = data.frame(metric = "Richness"),aes(x=-106,y=2.1),
            label="***",inherit.aes = FALSE,size=8,colour = "grey25") +
  geom_text(data = data.frame(metric = "Richness"), aes(x=15.4,y=1.1),
            label="**",inherit.aes = FALSE,size=8,colour = "grey25")+
  # geom_text(data = data.frame(metric = "PD"),aes(x=1.636,y=2.1),
  #           label="***",inherit.aes = FALSE,size=8,colour = "grey25") +
  # geom_text(data = data.frame(metric = "PD"), aes(x=0.015,y=1.15),
  #           label="ns",inherit.aes = FALSE,size=4,colour = "grey25")+
  geom_text(data = data.frame(metric = "MPD"),aes(x=-3.89,y=2.1),
            label="***",inherit.aes = FALSE,size=8,colour = "grey25") +
  geom_text(data = data.frame(metric = "MPD"), aes(x=-2.77,y=1.1),
            label="***",inherit.aes = FALSE,size=8,colour = "grey25")+
  #  geom_text(data = data.frame(metric = "FEve"),aes(x=0.046,y=2.1),
  #      label="***",inherit.aes = FALSE,size=8,colour = "grey25") +
  # geom_text(data = data.frame(metric = "FEve"), aes(x=0.0125,y=1.1),
  #          label="**",inherit.aes = FALSE,size=8,colour = "grey25")+
  # geom_text(data = data.frame(metric = "FDiv"),aes(x=-0.001,y=2.15),
  #           label="ns",inherit.aes = FALSE,size=4,colour = "grey25") +
  # geom_text(data = data.frame(metric = "FDiv"), aes(x=-0.0026,y=1.1),
  #             label="*",inherit.aes = FALSE,size=8,colour = "grey25")+
  # geom_text(data = data.frame(metric = "FDis"),aes(x=0.0577,y=2.1),
  #           label="***",inherit.aes = FALSE,size=8,colour = "grey25") +
  # geom_text(data = data.frame(metric = "FDis"), aes(x=0.0376,y=1.1),
  #             label="**",inherit.aes = FALSE,size=8,colour = "grey25")
  
  quartz()
out_metrics

#ggsave(path = "In_Out/figures",filename="delta_MPD_PD_rich.pdf", plot=out_metrics, dpi=300, units="cm")

###CWM analyses



library(dplyr)
library(purrr)
library(patchwork)

CWM.list.i<-all.out$CWM.Int
CWM.list.n<-all.out$CWM.Nat

CWM.list.n <- CWM.list.n[sapply(CWM.list.n,length)>1]
delta.CWM.n.L <- lapply(CWM.list.n, function(df) {
  df[1, ] - df[2, ]})

##pos = larger in buffer, neg, larger in urban
delta.CWM.n <- map_dfr(delta.CWM.n.L, ~.x, .id = "list_id")


CWM.list.i <- CWM.list.i[sapply(CWM.list.i,length)>1]
delta.CWM.i.L <- lapply(CWM.list.i, function(df) {
  df[1, ] - df[2, ]})

##pos = larger in buffer, neg, larger in urban
delta.CWM.i <- map_dfr(delta.CWM.i.L, ~.x, .id = "list_id")

##combine
delta.CWM.n<-delta.CWM.n[,2:ncol(delta.CWM.n)]
delta.CWM.i<-delta.CWM.i[,2:ncol(delta.CWM.i)]

plots <- lapply(names(delta.CWM.n), function(i) {
  d <- data.frame(
    value = c(delta.CWM.n[[i]], delta.CWM.i[[i]]),
    group = rep(c("Native", "Intro"), c(nrow(delta.CWM.n), nrow(delta.CWM.i))))

  ggplot(d, aes(value, colour = group, fill = group)) +
    geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed", linewidth = 0.8)+
    geom_histogram(
  aes(y = after_stat(count / sum(count)), group = group),
  position = "identity",
  alpha = 0.4,
  bins = 30
)+
    scale_colour_manual(values = c("darkorange3","darkolivegreen4")) +
    scale_fill_manual(values = c("darkorange3","darkolivegreen4")) +
    labs(title = i,x = NULL,y = "Relative frequency") +
    theme_minimal(base_size = 14)+
    theme(legend.position = "top")
})

quartz()
wrap_plots(plots, ncol = 3)

results0 <- lapply(names(delta.CWM.n), function(i){

  data.frame(
    Variable = i,
    Native.m = t.test(delta.CWM.n[[i]], mu = 0)$estimate,
    Native.p = t.test(delta.CWM.n[[i]], mu = 0)$p.value,
    Intro.m = t.test(delta.CWM.i[[i]], mu = 0)$estimate,
    Intro.p  = t.test(delta.CWM.i[[i]], mu = 0)$p.value
  )
})

do.call(rbind, results0)


results.D <- lapply(names(delta.CWM.n), function(i){

  data.frame(
    Variable = i,
    t.p = t.test(delta.CWM.n[[i]], delta.CWM.i[[i]])$p.value,
    wilcox.p = wilcox.test(delta.CWM.n[[i]], delta.CWM.i[[i]])$p.value,
    ks.p = ks.test(delta.CWM.n[[i]], delta.CWM.i[[i]])$p.value
  )
})

do.call(rbind, results.D)


delta.CWM.n$source <- "Native"
delta.CWM.i$source <- "Intro"

comb.delt.CWM <- bind_rows(delta.CWM.n,delta.CWM.i) %>%
  pivot_longer(cols=-source, names_to = "Var", values_to = "Value")





quartz()
dens(delta.PD.n$Inv.D.rich,lwd=3,col="darkorange3")
dens(delta.PD.n$Nat.D.rich,lwd=3,col="darkolivegreen4",add=TRUE)

dens(delta.PD.n$Inv.D.MPD,lwd=3,col="darkorange3")
dens(delta.PD.n$Nat.D.MPD,lwd=3,col="darkolivegreen4",add=TRUE)



out.i<-CWM.list.i[[1]][1,]

for (i in 1:length(CWM.list.i)){
  tmp<-CWM.list.i[[i]]
  if (class(tmp)=="data.frame"){
    if (ncol(tmp)>0){
      out.i[i,]<-tmp[2,]-tmp[1,]
    }
  }  
}
row.names(out.i)<-NULL
out.i<-na.omit(out.i)

out.n<-CWM.list.n[[1]][1,]

for (i in 1:length(CWM.list.n)){
  tmp<-CWM.list.n[[i]]
  if (class(tmp)=="data.frame"){
    if (ncol(tmp)>0){
      out.n[i,]<-tmp[2,]-tmp[1,]
    }
  }  
}
row.names(out.n)<-NULL
out.n<-na.omit(out.n)

out.i$Taxa <-"I"
out.n$Taxa <- "N"

out<-rbind(out.i,out.n)
out$Taxa<-as.factor(out$Taxa)
labs<-c("Seed mass","SLA","Height","Leaf area","Stem density","Leaf N", "Taxa")
colnames(out)<-labs
#stats
out.t<-list(NULL)
for (i in 1:6){
  out.t[[i]]<-t.test(out[,i]~out$Taxa)
  names(out.t)[i]<-colnames(out)[i]
}
apply(out.i[1:6],2,t.test)
apply(out.n[1:6],2,t.test)


#stats -log transformed
Lout.t<-list(NULL)
for (i in 1:6){
  Lout.t[[i]]<-t.test(log(out[,i])~out$Taxa)
  names(Lout.t)[i]<-colnames(out)[i]
}
apply(log(out.i[1:6]),2,t.test)
apply(log(out.n[1:6]),2,t.test)

quartz()
par(mfrow=c(2,3),cex.lab=1.4,cex.main=1.5)
plot(out$Taxa,out$`Seed mass`,ylab = "City-Buffer",
     xlab="Status",col=c("darkorange3","darkolivegreen4"),
     main="Seed mass")
abline(h=0,col="gray50",lty="dashed",lwd=2)
text(1.5,3700,"*",cex=4)
text(0.75,500,"*",cex=2.5,col="darkorange3")
text(1.75,400,"ns",cex=1.5,col="darkolivegreen4")
text(2,3700, "Greater in buffer", col="gray50")
text(1,-3500, "Greater in urban", col="gray50")

plot(out$Taxa,out$SLA,ylab = "City-Buffer",
     xlab="Status",col=c("darkorange3","darkolivegreen4"),
     main="SLA")
abline(h=0,col="gray50",lty="dashed",lwd=2)
text(1.5,9,"*",cex=4)
text(0.75,3,"*",cex=2.5,col="darkorange3")
text(1.75,3,"ns",cex=1.5,col="darkolivegreen4")

plot(out$Taxa,out$Height,ylab = "City-Buffer",
     xlab="Status",col=c("darkorange3","darkolivegreen4"),
     main="Height")
abline(h=0,col="gray50",lty="dashed",lwd=2)
text(1.5,15,"ns",cex=2)
text(0.75,3.5,"*",cex=2.5,col="darkorange3")
text(1.75,3,"*",cex=2.5,col="darkolivegreen4")

t.la<-out[out$`Leaf area`<40000,]
t.test(t.la$`Leaf area`[t.la$Taxa=="N"])
plot(t.la$Taxa,t.la$`Leaf area`,ylab = "City-Buffer",
     xlab="Status",col=c("darkorange3","darkolivegreen4"),
     main="Leaf area")
abline(h=0,col="gray50",lty="dashed",lwd=2)
text(1.5,15000,"ns",cex=2)
text(0.75,4000,"*",cex=2.5,col="darkorange3")
text(1.75,4000,"*",cex=2.5,col="darkolivegreen4")

plot(out$Taxa,out$`Stem density`,ylab = "City-Buffer",
     xlab="Status",col=c("darkorange3","darkolivegreen4"),
     main="Stem density")
abline(h=0,col="gray50",lty="dashed",lwd=2)
text(1.5,0.35,"*",cex=4)
text(0.75,0.08,"*",cex=2.5,col="darkorange3")
text(1.75,0.05,"*",cex=2.5,col="darkolivegreen4")

plot(out$Taxa,out$`Leaf N`,ylab = "City-Buffer",
     xlab="Status",col=c("darkorange3","darkolivegreen4"),
     main="Leaf N")
abline(h=0,col="gray50",lty="dashed",lwd=2)
text(1.5,9,"*",cex=4)
text(0.75,2,"*",cex=2.5,col="darkorange3")
text(1.75,2,"*",cex=2.5,col="darkolivegreen4")

#do non-parametric tests for difference from zero and each other


###plot ggplot, y-axis pseudo-log scale
sm<-ggplot(data=out,aes(x=Taxa,y=`Seed mass`,fill=Taxa))+
  geom_boxplot()+
  scale_y_continuous(trans = pseudo_log_trans(base = 10)) +
  #geom_hline(yintercept = 0,linetype = "dashed",colour = "grey45",linewidth = 1.2)+
  scale_fill_manual( values=c("darkorange3", "darkolivegreen4"))+
  xlab("Status") + ylab("Urban - buffer") +
  ggtitle("A) Seed mass")+
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(size = 15))+
  guides(fill = "none") +
  annotate("text",x=1.5,y=2000,label="*",color="grey20",
           fontface="bold",size=8)+
  annotate("text",x=1,y=100,label="*",color="grey90",
           fontface="bold",size=7)+
  annotate("text",x=2,y=15,label="ns",color="grey90",
           fontface="bold",size=4)
  
sla<-ggplot(data=out,aes(x=Taxa,y=SLA,fill=Taxa))+
  geom_boxplot()+
  scale_y_continuous(trans = pseudo_log_trans(base = 10)) +
  #geom_hline(yintercept = 0,linetype = "dashed",colour = "grey45",linewidth = 1.2)+
  scale_fill_manual( values=c("darkorange3", "darkolivegreen4"))+
  xlab("Status") + ylab("Urban - buffer") +
  ggtitle("B) SLA")+
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(size = 15))+
  guides(fill = "none") +
  annotate("text",x=1.5,y=8,label="*",color="grey20",
           fontface="bold",size=8)+
  annotate("text",x=1,y=0.3,label="*",color="grey90",
           fontface="bold",size=7)+
  annotate("text",x=2,y=0.34,label="ns",color="grey90",
           fontface="bold",size=4)

ht<-ggplot(data=out,aes(x=Taxa,y=Height,fill=Taxa))+
  geom_boxplot()+
  scale_y_continuous(trans = pseudo_log_trans(base = 10)) +
  #geom_hline(yintercept = 0,linetype = "dashed",colour = "grey45",linewidth = 1.2)+
  scale_fill_manual( values=c("darkorange3", "darkolivegreen4"))+
  xlab("Status") + ylab("Urban - buffer") +
  ggtitle("C) Height")+
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(size = 15))+
  guides(fill = "none") +
  annotate("text",x=1.5,y=14,label="ns",color="grey20",
           fontface="bold",size=5)+
  annotate("text",x=1,y=1,label="*",color="grey90",
           fontface="bold",size=7)+
  annotate("text",x=2,y=0.6,label="*",color="grey90",
           fontface="bold",size=7)

out.t<-out[out$`Leaf area`<40000,]

la<-ggplot(data=out.t,aes(x=Taxa,y=`Leaf area`,fill=Taxa))+
  geom_boxplot()+
  scale_y_continuous(trans = pseudo_log_trans(base = 10)) +
  #geom_hline(yintercept = 0,linetype = "dashed",colour = "grey45",linewidth = 1.2)+
  scale_fill_manual( values=c("darkorange3", "darkolivegreen4"))+
  xlab("Status") + ylab("Urban - buffer") +
  ggtitle("D) Leaf area")+
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(size = 15))+
  guides(fill = "none") +
  annotate("text",x=1.5,y=14000,label="ns",color="grey20",
           fontface="bold",size=5)+
  annotate("text",x=1,y=15,label="*",color="grey90",
           fontface="bold",size=7)+
  annotate("text",x=2,y=15,label="*",color="grey90",
           fontface="bold",size=7)

sd<-ggplot(data=out,aes(x=Taxa,y=`Stem density`,fill=Taxa))+
  geom_boxplot()+
  scale_y_continuous(trans = pseudo_log_trans(base = 10)) +
  #geom_hline(yintercept = 0,linetype = "dashed",colour = "grey45",linewidth = 1.2)+
  scale_fill_manual( values=c("darkorange3", "darkolivegreen4"))+
  xlab("Status") + ylab("Urban - buffer") +
  ggtitle("E) Stem density")+
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(size = 15))+
  guides(fill = "none") +
  annotate("text",x=1.5,y=0.35,label="*",color="grey20",
           fontface="bold",size=8)+
  annotate("text",x=1,y=0.02,label="*",color="grey90",
           fontface="bold",size=7)+
  annotate("text",x=2,y=0,label="*",color="grey90",
           fontface="bold",size=7)

ln<-ggplot(data=out,aes(x=Taxa,y=`Leaf N`,fill=Taxa))+
  geom_boxplot()+
  scale_y_continuous(trans = pseudo_log_trans(base = 10)) +
  #geom_hline(yintercept = 0,linetype = "dashed",colour = "grey45",linewidth = 1.2)+
  scale_fill_manual( values=c("darkorange3", "darkolivegreen4"))+
  xlab("Status") + ylab("Urban - buffer") +
  ggtitle("F) Leaf N")+
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(size = 15))+
  guides(fill = "none") +
  annotate("text",x=1.5,y=8,label="*",color="grey20",
           fontface="bold",size=8)+
  annotate("text",x=1,y=0.1,label="*",color="grey90",
           fontface="bold",size=7)+
  annotate("text",x=2,y=0.34,label="*",color="grey90",
           fontface="bold",size=7)

quartz()
grid.arrange(sm, sla, ht, la, sd, ln, nrow=2)
