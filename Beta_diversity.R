

setwd("~/Dropbox/2018WorkingFiles/Marc2018/GUBIC/sGUBIK")
source("In_Out/R/data_import.R")


##by city similarity between buffer and urban

all.out<-NULL
N.out<-NULL
I.out<-NULL

cities<-unique(dat.long.sp.remove$ID_HDC_G0_2)

for (i in 1:length(cities)){
  tmp<-dat.long.sp.remove[dat.long.sp.remove$ID_HDC_G0_2==cities[i],]
  all.tmp.mat<-sample2matrix(tmp[,c("location","records","accepted_taxon_name_binomial")])
  all.tmp.out<-betadiver(all.tmp.mat,"sor")
  all.out[i]<-as.vector(all.tmp.out)
  
  tmp.N<-tmp[tmp$I_N=="N",]
  N.tmp.mat<-sample2matrix(tmp.N[,c("location","records","accepted_taxon_name_binomial")])
  N.tmp.out<-betadiver(N.tmp.mat,"sor")
  N.out[i]<-as.vector(N.tmp.out)
  
  tmp.I<-tmp[tmp$I_N=="I",]
  I.tmp.mat<-sample2matrix(tmp.I[,c("location","records","accepted_taxon_name_binomial")])
  I.tmp.out<-betadiver(I.tmp.mat,"sor")
  I.out[i]<-as.vector(I.tmp.out)
}

all.S<-rbind(data.frame(Sor = all.out,Taxa = "All"),data.frame(Sor = N.out,Taxa = "Native"),
             data.frame(Sor = I.out,Taxa = "Intro"))

lm.S<-lm(Sor~Taxa,data=all.S)
summary(lm.S)
anova(lm.S)

cols<-c("Intro"="#ba882b","Native"="#709c81","All"="#9c7089")

plot.S<-ggplot(all.S, aes(x = Taxa, y = Sor, fill = Taxa)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.fill = "white") +
  scale_fill_manual(values=cols) +
  #scale_y_log10() +
  labs(x = "Taxa group", y = "Buffer-urban dissimilarity (Sorenson)") +
  theme_minimal(base_size = 14)

quartz()
plot.S

#mean values
tapply(all.S$Sor,all.S$Taxa,mean,na.rm=TRUE)
mean(all.S$Sor)
sd(all.S$Sor)
kruskal.test(Sor~Taxa,data=all.S)

##similarity between buffers globally, and cities globally, for the taxa groups
S.mat<-readRDS("In_Out/data/urban_center_distance.rds")

#urban + 5km radius combined
both_all_mat<-sample2matrix(dat.long.sp.remove[,c("ID_HDC_G0_2","records","accepted_taxon_name_binomial")])
#both_all.out<-betadiver(both_all_mat,"sor")

both_all.out<-vegdist(both_all_mat,method = "bray", binary = TRUE)
both_all.out<-as.matrix(both_all.out)

comm_rows <- intersect(rownames(S.mat), rownames(both_all.out))
comm_cols <- intersect(colnames(S.mat), colnames(both_all.out))

Spat <- S.mat[comm_rows, comm_cols, drop = FALSE]
Sor_all <- both_all.out[comm_rows, comm_cols, drop = FALSE]

df.tmp<-data.frame(Spat=as.vector(Spat),Sor = as.vector(Sor_all),Taxa = "All",Comp="Urban+5km")

quartz()
pl<-ggplot(df.tmp, aes(x=Spat,y=Sor))+
  geom_point(color = "lightblue",size=2,alpha=0.8)+
  geom_smooth(method = "loess", se = FALSE, color = "dodgerblue3", fill = "lightblue") +
  theme_minimal() +
  labs(x = "Spatial distance", y = "Sorensen dissimilarity")  
pl

both.N<-dat.long.sp.remove[dat.long.sp.remove$I_N=="N",]
both_N_mat<-sample2matrix(both.N[,c("ID_HDC_G0_2","records","accepted_taxon_name_binomial")])
both_N.out<-vegdist(both_N_mat,method = "bray", binary = TRUE)
both_N.out<-as.matrix(both_N.out)

comm_rows <- intersect(rownames(S.mat), rownames(both_N.out))
comm_cols <- intersect(colnames(S.mat), colnames(both_N.out))

Spat <- S.mat[comm_rows, comm_cols, drop = FALSE]
Sor_N <- both_N.out[comm_rows, comm_cols, drop = FALSE]

df.N.tmp<-data.frame(Spat=as.vector(Spat),Sor = as.vector(Sor_N),Taxa = "Native",Comp="Urban+5km")

both.I<-dat.long.sp.remove[dat.long.sp.remove$I_N=="I",]
both_I_mat<-sample2matrix(both.I[,c("ID_HDC_G0_2","records","accepted_taxon_name_binomial")])
both_I.out<-vegdist(both_I_mat,method = "bray", binary = TRUE)
both_I.out<-as.matrix(both_I.out)

comm_rows <- intersect(rownames(S.mat), rownames(both_I.out))
comm_cols <- intersect(colnames(S.mat), colnames(both_I.out))

Spat <- S.mat[comm_rows, comm_cols, drop = FALSE]
Sor_I <- both_I.out[comm_rows, comm_cols, drop = FALSE]

df.I.tmp<-data.frame(Spat=as.vector(Spat),Sor = as.vector(Sor_I),Taxa = "Intro",Comp="Urban+5km")


all<-rbind(df.tmp,df.N.tmp,df.I.tmp)
N_I_both<-rbind(df.N.tmp,df.I.tmp)

#write.csv(all, 'In_Out/data/Sor_buff+urb.csv')
#all<-read.csv("In_Out/data/Sor_buff+urb.csv",row.names=1)
#N<-all[all$Taxa=="Native",]
#I<-all[all$Taxa=="Intro",]
#N_I_both<-rbind(N, I)
#test N_I_both2 <-N_I_both[sample(1:nrow(N_I_both),3000),]

quartz()
NIB<-ggplot(N_I_both, aes(x=Spat,y=Sor,color=Taxa))+
  geom_point(data = subset(N_I_both, Taxa == "Native"),
             aes(x = Spat, y = Sor), color = "#93ab9c", alpha = 0.6, size = 2) +
  geom_point(data = subset(N_I_both, Taxa == "Intro"),
             aes(x = Spat, y = Sor), color = "#bda06c", alpha = 0.6, size = 2) +
  geom_smooth(data = subset(N_I_both, Taxa == "Native"),
              aes(x = Spat, y = Sor), color = "#709c81", se = FALSE, linewidth = 1.4) +
  geom_smooth(data = subset(N_I_both, Taxa == "Intro"),
              aes(x = Spat, y = Sor), color = "#ba882b", se = FALSE, linewidth = 1.4) +
  theme_minimal(base_size = 14) +
  labs(x = "Spatial distance", y = "Sorensen dissimilarity",
       title="A) Urban + 5 km buffer") 
NIB


#buffer only
buff<-dat.long.sp.remove[dat.long.sp.remove$location=="buffer_5km_records",]
buff_all_mat<-sample2matrix(buff[,c("ID_HDC_G0_2","records","accepted_taxon_name_binomial")])
buff_all.out<-vegdist(buff_all_mat,method = "bray", binary = TRUE)
buff_all.out<-as.matrix(buff_all.out)

comm_rows <- intersect(rownames(S.mat), rownames(buff_all.out))
comm_cols <- intersect(colnames(S.mat), colnames(buff_all.out))

Spat <- S.mat[comm_rows, comm_cols, drop = FALSE]
Sor_buff_all <- buff_all.out[comm_rows, comm_cols, drop = FALSE]

df.buff.tmp<-data.frame(Spat=as.vector(Spat),Sor = as.vector(Sor_buff_all),Taxa = "All",Comp="5km buffer")


buff.N<-buff[buff$I_N=="N",]
buff_N_mat<-sample2matrix(buff.N[,c("ID_HDC_G0_2","records","accepted_taxon_name_binomial")])
buff_N.out<-vegdist(buff_N_mat,method = "bray", binary = TRUE)
buff_N.out<-as.matrix(buff_N.out)

comm_rows <- intersect(rownames(S.mat), rownames(buff_N.out))
comm_cols <- intersect(colnames(S.mat), colnames(buff_N.out))

Spat <- S.mat[comm_rows, comm_cols, drop = FALSE]
Sor_buff_N <- buff_N.out[comm_rows, comm_cols, drop = FALSE]

df.buff.N.tmp<-data.frame(Spat=as.vector(Spat),Sor = as.vector(Sor_buff_N),Taxa = "Native",Comp="5km buffer")


buff.I<-buff[buff$I_N=="I",]
buff_I_mat<-sample2matrix(buff.I[,c("ID_HDC_G0_2","records","accepted_taxon_name_binomial")])
buff_I.out<-vegdist(buff_I_mat,method = "bray", binary = TRUE)
buff_I.out<-as.matrix(buff_I.out)

comm_rows <- intersect(rownames(S.mat), rownames(buff_I.out))
comm_cols <- intersect(colnames(S.mat), colnames(buff_I.out))

Spat <- S.mat[comm_rows, comm_cols, drop = FALSE]
Sor_buff_I <- buff_I.out[comm_rows, comm_cols, drop = FALSE]

df.buff.I.tmp<-data.frame(Spat=as.vector(Spat),Sor = as.vector(Sor_buff_I),Taxa = "Intro",Comp="5km buffer")

all.buff<-rbind(df.buff.tmp,df.buff.N.tmp,df.buff.I.tmp)
N_I_buff<-rbind(df.buff.N.tmp,df.buff.I.tmp)

#write.csv(all.buff, 'In_Out/data/Sor_buff.csv')
#all<-read.csv("In_Out/data/Sor_buff.csv",row.names=1)
#N<-all[all.buff$Taxa=="Native",]
#I<-all[all.buff$Taxa=="Intro",]
#N_I_buff<-rbind(N, I)


quartz()

buff<-ggplot(N_I_buff, aes(x=Spat,y=Sor,color=Taxa))+
  geom_point(data = subset(N_I_buff, Taxa == "Native"),
             aes(x = Spat, y = Sor), color = "#93ab9c", alpha = 0.6, size = 2) +
  geom_point(data = subset(N_I_buff, Taxa == "Intro"),
             aes(x = Spat, y = Sor), color = "#bda06c", alpha = 0.6, size = 2) +
  geom_smooth(data = subset(N_I_buff, Taxa == "Native"),
              aes(x = Spat, y = Sor), color = "#709c81", se = FALSE, linewidth = 1.4) +
  geom_smooth(data = subset(N_I_buff, Taxa == "Intro"),
              aes(x = Spat, y = Sor), color = "#ba882b", se = FALSE, linewidth = 1.4) +
  theme_minimal(base_size = 14) +
  labs(x = "Spatial distance", y = "Sorensen dissimilarity",
       title = "C) 5 km buffer")  
buff

#urban only
urb<-dat.long.sp.remove[dat.long.sp.remove$location=="urban_records",]
urb_all_mat<-sample2matrix(urb[,c("ID_HDC_G0_2","records","accepted_taxon_name_binomial")])
urb_all.out<-vegdist(urb_all_mat,method = "bray", binary = TRUE)
urb_all.out<-as.matrix(urb_all.out)

comm_rows <- intersect(rownames(S.mat), rownames(urb_all.out))
comm_cols <- intersect(colnames(S.mat), colnames(urb_all.out))

Spat <- S.mat[comm_rows, comm_cols, drop = FALSE]
Sor_urb_all <- urb_all.out[comm_rows, comm_cols, drop = FALSE]

df.urb.tmp<-data.frame(Spat=as.vector(Spat),Sor = as.vector(Sor_urb_all),Taxa = "All",Comp="Urban")


urb.N<-urb[urb$I_N=="N",]
urb_N_mat<-sample2matrix(urb.N[,c("ID_HDC_G0_2","records","accepted_taxon_name_binomial")])
urb_N.out<-vegdist(urb_N_mat,method = "bray", binary = TRUE)
urb_N.out<-as.matrix(urb_N.out)

comm_rows <- intersect(rownames(S.mat), rownames(urb_N.out))
comm_cols <- intersect(colnames(S.mat), colnames(urb_N.out))

Spat <- S.mat[comm_rows, comm_cols, drop = FALSE]
Sor_urb_N <- urb_N.out[comm_rows, comm_cols, drop = FALSE]

df.urb.N.tmp<-data.frame(Spat=as.vector(Spat),Sor = as.vector(Sor_urb_N),Taxa = "Native",Comp="Urban")

urb.I<-urb[urb$I_N=="I",]
urb_I_mat<-sample2matrix(urb.I[,c("ID_HDC_G0_2","records","accepted_taxon_name_binomial")])
urb_I.out<-vegdist(urb_I_mat,method = "bray", binary = TRUE)
urb_I.out<-as.matrix(urb_I.out)

comm_rows <- intersect(rownames(S.mat), rownames(urb_I.out))
comm_cols <- intersect(colnames(S.mat), colnames(urb_I.out))

Spat <- S.mat[comm_rows, comm_cols, drop = FALSE]
Sor_urb_I <- urb_I.out[comm_rows, comm_cols, drop = FALSE]

df.urb.I.tmp<-data.frame(Spat=as.vector(Spat),Sor = as.vector(Sor_urb_I),Taxa = "Intro",Comp="Urban")

all.urb<-rbind(df.urb.tmp,df.urb.N.tmp,df.urb.I.tmp)
N_I_urb<-rbind(df.urb.N.tmp,df.urb.I.tmp)

#write.csv(all.urb, 'In_Out/data/Sor_urb.csv')
#all<-read.csv("In_Out/data/Sor_urb.csv",row.names=1)
#N<-all[all.urb$Taxa=="Native",]
#I<-all[all.urb$Taxa=="Intro",]
#N_I_urb<-rbind(N, I)
quartz()

urb<-ggplot(N_I_urb, aes(x=Spat,y=Sor,color=Taxa))+
  geom_point(data = subset(N_I_urb, Taxa == "Native"),
             aes(x = Spat, y = Sor), color = "#93ab9c", alpha = 0.6, size = 2) +
  geom_point(data = subset(N_I_urb, Taxa == "Intro"),
             aes(x = Spat, y = Sor), color = "#bda06c", alpha = 0.6, size = 2) +
  geom_smooth(data = subset(N_I_urb, Taxa == "Native"),
              aes(x = Spat, y = Sor), color = "#709c81", se = FALSE, linewidth = 1.4) +
  geom_smooth(data = subset(N_I_urb, Taxa == "Intro"),
              aes(x = Spat, y = Sor), color = "#ba882b", se = FALSE, linewidth = 1.4) +
  theme_minimal(base_size = 14) +
  labs(x = "Spatial distance", y = "Sorensen dissimilarity",
       title = "B) Urban")  
urb

quartz()
grid.arrange(NIB,urb,buff, ncol = 3)


glob.dat<-rbind(data.frame(Sor = as.vector(both_all.out),Taxa = "All",Comp="Urban + buffer"),
                data.frame(Sor = as.vector(both_N.out),Taxa = "Native",Comp="Urban + buffer"),
                data.frame(Sor = as.vector(both_I.out),Taxa = "Intro",Comp="Urban + buffer"),
                data.frame(Sor = as.vector(buff_all.out),Taxa = "All",Comp="5km buffer"),
                data.frame(Sor = as.vector(buff_N.out),Taxa = "Native",Comp="5km buffer"),
                data.frame(Sor = as.vector(buff_I.out),Taxa = "Intro",Comp="5km buffer"),
                data.frame(Sor = as.vector(urb_all.out),Taxa = "All",Comp="Urban"),
                data.frame(Sor = as.vector(urb_N.out),Taxa = "Native",Comp="Urban"),
                data.frame(Sor = as.vector(urb_I.out),Taxa = "Intro",Comp="Urban"))
glob.dat2<-na.omit(glob.dat)

summary_groups <- glob.dat2 %>%
  group_by(Comp, Taxa) %>%
  summarise(mean = mean(Sor),
            median = median(Sor),
            sd = sd(Sor),
            n = n(),
            .groups = "drop")

summary_groups
###plots include urban + buffer, below I redo without it

cols<-c("All"="#9c7089","Native"="darkolivegreen4","Intro"="darkorange3")

plot.both<-ggplot(glob.dat2, aes(x = Taxa, y = Sor, fill = Taxa)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.fill = "white") +
  scale_fill_manual(values=cols) +
  scale_y_log10() +
  facet_wrap(~Comp)+
  labs(x = "Taxa group", y = "Log Dissimilarity between cities (Sorenson)") +
  theme_minimal(base_size = 14)+
  theme(strip.text = element_text(size = 14, face = "bold"))

quartz()
plot.both

#mod.TxUxB <- betareg(Sor ~ as.factor(Taxa)*as.factor(Comp) , data = glob.dat2)
#summary(mod.TxUxB)

##need to adjust Sor values, can't have 0 or 1s
eps <- 1e-6
glob.dat2$Sor_adj <- pmin(pmax(glob.dat2$Sor, eps), 1 - eps)

modTUB <- glmmTMB(Sor_adj ~ Taxa * Comp,data = glob.dat2,family = beta_family())
summary(modTUB)


emm <- emmeans(modTUB, ~ Taxa * Comp, type = "response")
emm
#m.out<-kruskal.test(Sor ~ as.factor(Comp)*as.factor(Comp) , data = glob.dat2)
#summary(m.out)
pairs(emm, by = "Comp", type = "response")

pairs(emm, by = "Taxa", type = "response")

##with errors

emm_df <- as.data.frame(emmeans(modTUB, ~ Taxa * Comp, type = "response"))

quartz()
ggplot(emm_df, aes(x = Comp, y = response, color = Taxa, group = Taxa)) +
  geom_line(linewidth = 1.4) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1, linewidth = 0.8) +
  scale_color_manual(values = c("Native" = "darkolivegreen4","Intro"  = "darkorange3", "All" = "#9c7089")) +
  theme_minimal(base_size = 14) +
  labs(x = "Comparison type",
       y = "Predicted Sørensen dissimilarity",
       color = "Taxa") +
  theme(legend.position = "top",
        legend.title = element_text(face = "bold"))

###now without Urban+buffer
glob.dat3<-droplevels(subset(glob.dat2, Taxa != "All" & Comp != "Urban + buffer"))


cols<-c("All"="#9c7089","Native"="darkolivegreen4","Intro"="darkorange3")

plot.both<-ggplot(glob.dat3, aes(x = Taxa, y = Sor, fill = Taxa)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.fill = "white") +
  scale_fill_manual(values=cols) +
  scale_y_log10() +
  facet_wrap(~Comp)+
  labs(x = "Taxa group", y = "Log Dissimilarity between cities (Sorenson)") +
  theme_minimal(base_size = 14)+
  theme(strip.text = element_text(size = 14, face = "bold"))

quartz()
plot.both

#mod.TxUxB <- betareg(Sor ~ as.factor(Taxa)*as.factor(Comp) , data = glob.dat2)
#summary(mod.TxUxB)

##need to adjust Sor values, can't have 0 or 1s
eps <- 1e-6
glob.dat3$Sor_adj <- pmin(pmax(glob.dat3$Sor, eps), 1 - eps)

modTUB <- glmmTMB(Sor_adj ~ Taxa * Comp,data = glob.dat3,family = beta_family())
summary(modTUB)


emm <- emmeans(modTUB, ~ Taxa * Comp, type = "response")
emm
#m.out<-kruskal.test(Sor ~ as.factor(Comp)*as.factor(Comp) , data = glob.dat2)
#summary(m.out)
pairs(emm, by = "Comp", type = "response")

pairs(emm, by = "Taxa", type = "response")

##with errors

emm_df <- as.data.frame(emmeans(modTUB, ~ Taxa * Comp, type = "response"))

quartz()
ggplot(emm_df, aes(x = Comp, y = response, color = Taxa, group = Taxa)) +
  geom_line(linewidth = 1.4) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1, linewidth = 0.8) +
  scale_color_manual(values = c("Native" = "darkolivegreen4","Intro"  = "darkorange3", "All" = "#9c7089")) +
  theme_minimal(base_size = 14) +
  labs(x = "Comparison type",
       y = "Predicted Sørensen dissimilarity",
       color = "Taxa") +
  theme(legend.position = "top",
        legend.title = element_text(face = "bold"))
