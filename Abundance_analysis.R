setwd("~/Dropbox/2018WorkingFiles/Marc2018/GUBIC/sGUBIK")
source("In_Out/R/data_import.R")

#read in rarified data
Rare.dat<-read.csv("In_Out/R/Rare.dat.csv",row.names = 1)
Rare.dat$log.records<-log(Rare.dat$N_total)

IN_hist<-ggplot(dat.no.0, aes(x = per_change, fill = I_N)) +
  geom_histogram(position = "identity", alpha = 0.7, bins = 30) +
  scale_fill_manual(values = c("I" = "darkorange3","N" = "darkolivegreen4")) +
  labs(x = "Percent Change", y = "Count", fill = "") +
  scale_x_log10() +
  theme_minimal(base_size = 14)

quartz()
IN_hist

CI95<-function(x){
  n  <- length(x)
  m  <- mean(x, na.rm = TRUE)
  se <- sd(x, na.rm = TRUE) / sqrt(n)
  m + c(-1, 1) * qt(0.975, df = n - 1) * se
}


data.frame(Taxa = c("N","I"), 
           Mean_per = c(mean(dat.no.0$per_change[dat.no.0$I_N=="N"],na.rm=TRUE),
                        mean(dat.no.0$per_change[dat.no.0$I_N=="I"],na.rm=TRUE)),
           CI95.L = c(CI95(dat.no.0$per_change[dat.no.0$I_N=="N"])[1],
                      CI95(dat.no.0$per_change[dat.no.0$I_N=="I"])[1]),
           CI95.U = c(CI95(dat.no.0$per_change[dat.no.0$I_N=="N"])[2],
                      CI95(dat.no.0$per_change[dat.no.0$I_N=="I"])[2]),
           Ttest = c(t.test(dat.no.0$per_change[dat.no.0$I_N=="N"])$statistic,
                     t.test(dat.no.0$per_change[dat.no.0$I_N=="I"])$statistic),
           Pval = c(t.test(dat.no.0$per_change[dat.no.0$I_N=="N"])$p.value,
                    t.test(dat.no.0$per_change[dat.no.0$I_N=="I"])$p.value)
           )


lm1<-lmer((per_change)~I_N+(1|ID_HDC_G0_2),data=dat.no.0)
summary(lm1)
emm <- emmeans(lm1, specs =  pairwise~ I_N)
emm$emmeans
emm$contrasts

I_PC<-dat.no.0[dat.no.0$I_N=="I",]
N_PC<-dat.no.0[dat.no.0$I_N=="N",]

diff<-data.frame(I.mean.per.dif = tapply(I_PC$per_change,I_PC$ID_HDC_G0_2,mean,na.rm=TRUE),
           N.mean.per.dif = tapply(N_PC$per_change,N_PC$ID_HDC_G0_2,mean,na.rm=TRUE))
diff$pop<-city$POP_2015[match(rownames(diff),city$ID_HDC_G0_2)]
diff$gdp<-city$gdp_md_est_2016[match(rownames(diff),city$ID_HDC_G0_2)]  
diff$continent<-city$continent[match(rownames(diff),city$ID_HDC_G0_2)]
diff$native<-city$nsp_not_in_GloNAF[match(rownames(diff),city$ID_HDC_G0_2)]  
diff$intro<-city$nsp_naturalized_GloNAF[match(rownames(diff),city$ID_HDC_G0_2)] 
diff$tot.spp<-diff$native+diff$intro

dat.no.0$tot_records<-dat.no.0$buffer_5km_records  + dat.no.0$urban_records
  
diff$tot.records <-  tapply(dat.no.0$tot_records,dat.no.0$ID_HDC_G0_2,sum,na.rm=TRUE)




##

mDiff.N<-lm(log(N.mean.per.dif+100)~log(pop)+tot.records,data=diff)
summary(mDiff.N)

mDiff.I<-lm(log(I.mean.per.dif+100)~log(pop)+log(tot.spp),data=diff)
summary(mDiff.I)

cDiff.I<-lm(log(I.mean.per.dif+100)~continent,data=diff)
summary(cDiff.I)
anova(cDiff.I)

TukeyHSD(aov(log(I.mean.per.dif+100)~continent,data=diff))

cDiff.N<-lm(log(N.mean.per.dif+100)~continent,data=diff)
summary(cDiff.N)
anova(cDiff.N)

TukeyHSD(aov(log(N.mean.per.dif+100)~continent,data=diff))

quartz()
par(mfrow=c(2,1))
plot(as.factor(diff$continent),log(diff$N.mean.per.dif),
     xlab = "Continent",ylab = "Log Ave. record difference",main = "Native")
text()

plot(as.factor(diff$continent),log(diff$I.mean.per.dif),
     xlab = "Continent",ylab = "Log Ave. record difference",main = "Introduced")



quartz()
par(mfrow=c(2,2))

plot(log(diff$pop),log(diff$N.mean.per.dif+100),pch=19,cex=1.3,xlab = "Log Population size",
     ylab = "Ave. record difference",main = "Native")
plot(log(diff$tot.spp),log(diff$N.mean.per.dif+100),pch=19,cex=1.3,xlab = "Log Total species",
     ylab = "Ave. record difference",main = "Native")

plot(log(diff$pop),log(diff$I.mean.per.dif+100),pch=19,cex=1.3,xlab = "Log Population size",
     ylab = "Ave. record difference",main = "Introduced")
plot(log(diff$tot.spp),log(diff$I.mean.per.dif+100),pch=19,cex=1.3,xlab = "Log Total species",
     ylab = "Ave. record difference",main = "Introduced")


plot(log(diff$tot.records),log(diff$N.mean.per.dif+100),pch=19,cex=1.3,xlab = "Log Population size",
     ylab = "Ave. record difference",main = "Native")
plot(log(diff$pop),log(diff$tot.records),pch=19,cex=1.3,xlab = "Log Population size",
     ylab = "Ave. record difference",main = "Native")
summary(lm(tot.spp~pop,data=diff))

dat.no.0.long <- dat.no.0 %>%
  pivot_longer(
    cols = c(buffer_5km_records, urban_records),
    names_to = "location",
    values_to = "records"
  )

dat.no.0.long$log.records<-log(dat.no.0.long$records)


summary_dat.no.0.long <- dat.no.0.long %>%
  group_by(location, I_N) %>%
  summarise(
    mean_rec = mean(log.records),
    se_rec = sd(log.records)/sqrt(n()),
    CI95U = qt(0.975, df = n() - 1) * sd(log.records)/sqrt(n()),
    CI95L = qt(0.975, df = n() - 1) * sd(log.records)/sqrt(n()),
    SD = sd(log.records),
    mean.notrans = mean(records),
    CI95U.notrans = qt(0.975, df = n() - 1) * sd(records)/sqrt(n()),
    CI95L.notrans = qt(0.975, df = n() - 1) * sd(records)/sqrt(n()),
    .groups = 'drop'
  )

summary_Rare <- Rare.dat %>%
  group_by(location, I_N) %>%
  summarise(
    mean_rec = mean(log.records),
    se_rec = sd(log.records)/sqrt(n()),
    CI95U = qt(0.975, df = n() - 1) * sd(log.records)/sqrt(n()),
    CI95L = qt(0.975, df = n() - 1) * sd(log.records)/sqrt(n()),
    SD = sd(log.records),
    mean.notrans = mean(N_total),
    CI95U.notrans = qt(0.975, df = n() - 1) * sd(N_total)/sqrt(n()),
    CI95L.notrans = qt(0.975, df = n() - 1) * sd(N_total)/sqrt(n()),
    .groups = 'drop'
  )

##plot boxplots of percent change
quartz()
p<-ggplot(dat.no.0, aes(x = I_N, y = per_change,fill = I_N)) +
  geom_boxplot()+
  scale_y_continuous(trans = pseudo_log_trans(base = 10)) +
  scale_fill_manual( values=c("darkorange3", "darkolivegreen4"))+
  xlab("Introduced or Native") + ylab("Percent change (pseudo-log scale)") +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  annotate("segment", 
           y = 10, yend = 50000, x = 0.5, xend = 0.5,
           arrow=arrow(length=unit(0.5, "cm"),type="closed"),color="orange3",linewidth=2)+
  annotate("text",x=0.55,y=1000,label="Increasing in buffer",color="orange3",
           fontface="bold",size=5,angle=90)

p 


##plot density distributions for record numbers

rec.buf<-ggplot(dat.no.0, aes(x = buffer_5km_records, fill = I_N)) +
  geom_density(alpha = 0.8) +
  scale_x_continuous(trans = pseudo_log_trans(base = 10)) +
  theme_minimal() +
  ggtitle("A) 5 Km buffer records")+
  scale_fill_manual( values=c("darkorange3", "darkolivegreen4"))+
  labs(x = "Number of records",
       y = "Density")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12,face="bold"))

rec.urb<-ggplot(dat.no.0, aes(x = urban_records, fill = I_N)) +
  geom_density(alpha = 0.8) +
  scale_x_continuous(trans = pseudo_log_trans(base = 10)) +
  theme_minimal() +
  ggtitle("B) Urban records")+
  scale_fill_manual( values=c("darkorange3", "darkolivegreen4"))+
  labs(x = "Number of records",
       y = "Density")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12,face="bold"))
  
quartz()
grid.arrange(rec.buf, rec.urb,nrow=2)


ks.test(dat.no.0$buffer_5km_records[dat.no.0$I_N=="I"],dat.no.0$buffer_5km_records[dat.no.0$I_N=="N"])
ks.test(dat.no.0$urban_records[dat.no.0$I_N=="I"],dat.no.0$urban_record[dat.no.0$I_N=="N"])


##plot density distributions for percent difference

dd<-ggplot(dat.no.0, aes(x = per_change, fill = I_N)) +
  geom_density(alpha = 0.8) +
  scale_x_continuous(trans = pseudo_log_trans(base = 10)) +
  theme_minimal() +
  scale_fill_manual( values=c("darkorange3", "darkolivegreen4"))+
  geom_vline(xintercept = 0,linetype = "dashed",colour = "grey45",linewidth = 1.2)+
  ggtitle("B) Species % diff")+
  labs(x = "Percent change (pseudo-log scale)",
       y = "Density")+
  theme(axis.text.y=element_text(size=11), axis.text.x=element_text(size=9,angle=45,hjust=1),
        axis.title=element_text(size=12,face="bold"))+
  annotate("text",x=-30,y=0.72,label="Greater in urban",color="grey45",
           fontface="bold",size=4)+
  annotate("text",x=100,y=0.72,label="Greater in buffer",color="grey45",
           fontface="bold",size=4)

quartz()
dd

ks.test(dat.no.0$per_change[dat.no.0$I_N=="I"],dat.no.0$per_change[dat.no.0$I_N=="N"])

library(quantreg)

quan.reg<-rq(per_change ~ I_N, tau = c(0.1,0.5,0.9), data = dat.no.0)
summary(quan.reg)


## Ratio N to I species for urban and buffer for each city I:N

rat <- function(x){
  buff<-x$buffer_5km_records
  urb<-x$urban_records
  IN<-x$I_N
  rB<-buff[buff>0]
  INB<-IN[buff>0]
  rU<-urb[urb>0]
  INU<-IN[urb>0]
  R_Buff <-length(rB[INB=="I"])/length(rB[INB=="N"])
  R_Urb <- length(rU[INU=="I"])/length(rU[INU=="N"])
  R_Buff_Tot <- length(rB[INB=="I"])/length(buff)
  R_Urb_Tot <- length(rU[INU=="I"])/length(buff)
  return(data.frame(City = x$ID_HDC_G0_2[1],R_Buff,R_Urb,R_Buff_Tot,R_Urb_Tot))
}

Cities<-unique(dat$ID_HDC_G0_2)
out<-list()

for (i in 1:length(Cities)){
  tmp<-dat[dat$ID_HDC_G0_2==Cities[i],]
  out[[i]]<-rat(tmp)  
}

outR<-do.call(rbind, out)
outR.L<-c(outR$R_Buff,outR$R_Urb)
loc<-c(rep("Buff",nrow(outR)),rep("Urb",nrow(outR)))
outR.L<-data.frame(loc,I.N_ratio=outR.L)

INR_hist<-ggplot(outR.L, aes(y = I.N_ratio,x=loc, fill = loc)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("Buff" = "peru","Urb" = "steelblue")) +
  labs(x = "Location", y = "I:N ratio") +
 # scale_y_log10() +
  theme_minimal(base_size = 14)+
 theme(legend.position = "none")
  

mean(outR$R_Urb,na.rm=TRUE)
mean(outR$R_Buff,na.rm=TRUE)
t.test(outR$R_Urb,outR$R_Buff,paired=TRUE)

quartz()
INR_hist

# 2. Create plot with connecting lines (need rarefaction)

p_BU<-ggplot(summary_dat.no.0.long, aes(x = factor(location), y = mean_rec, 
                                        color = I_N, group = I_N)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_rec - CI95L, ymax = mean_rec + CI95U),
                width = 0.2) +
  geom_line(linewidth = 1.5) +
  scale_x_discrete(labels = c("Buffer","Urban"))+
  labs(x = "Location of records",y = "Log(mean number of species records)",color = "Status") +
  theme_bw()+
  ggtitle("A) Number of records")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12,face="bold"))+
  theme(legend.text = element_text(size = 12))+
  theme(legend.title = element_text(size = 14,face="bold"))+
  scale_color_manual(values = c("I" = "darkorange3", "N" = "darkolivegreen4"))

quartz()
p_BU

quartz()
grid.arrange(p_BU,dd,nrow=2)
#ggsave("In_Out/figures/buff_urb.pdf", plot = p_BU, width = 7, height = 7)

##model for above plot
lm.l<-lmer(log.records~I_N+location+I_N:location+(1|ID_HDC_G0_2),data=dat.no.0.long)
summary(lm.l)
anova(lm.l)
emm <- emmeans(lm.l, specs =  pairwise~ I_N)
emm$emmeans
emm$contrasts

#just look at exotics
ex_long<-dat.no.0.long[dat.no.0.long$I_N=="I",]
lm.ex<-lmer(log.records~location+(1|ID_HDC_G0_2),data=ex_long)
summary(lm.ex)
anova(lm.ex)###not significant

