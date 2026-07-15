library(ape)  
library(betapart) 
library(betareg)
library(dplyr)
library(ecodist)
library(emmeans)
library(FD) 
library(hillR)
library(glmmTMB)
library(ggplot2)
library(gridExtra)
library(lme4)
library(lmerTest)
library(mFD)
library(missForest)
library(picante)
#library(rethinking)
library(scales)
library(statmod)
library(tibble)
library(tidyr)
library(vegan)  
library(VennDiagram)

###GUBIC inside/outside prelim analysis
setwd("~/Dropbox/2018WorkingFiles/Marc2018/GUBIC/sGUBIK")

city<-read.csv('GUBIC data/urban_sf_sgubik.csv')
dat<-read.csv('GUBIC data/gbif_gubic_list_per_city_final2_20241015.csv')

dat$ID_HDC_G0_2<-as.factor(dat$ID_HDC_G0_2)

#add column for Introduced and native
dat$I_N<-"N"
dat$I_N[is.na(dat$status_glonaf)==FALSE]<-"I"


#species found in both urban and buffer (getting rid of singletons)
dat.no.0<-dat[dat$buffer_5km_records > 1 & dat$urban_records > 1,]
dat.no.0<-dat.no.0[is.na(dat.no.0$buffer_5km_records)==FALSE,]

dim(dat.no.0) #353,197 species records in both

#percent change from urban to buffer: positive = lower in urban
#negative = higher in urban; 0 = no difference
per_change<-NULL
for (i in 1:nrow(dat.no.0)){
  per_change[i]<-100*(dat.no.0$buffer_5km_records[i]-dat.no.0$urban_records[i])/dat.no.0$urban_records[i]
}
dat.no.0$per_change<-per_change

##read in phylogeny

phy<-read.tree("GUBIC data/phylo_gubic.tre")

Xunder<-function(names) {
  as.vector(sapply(names,function(x) sub("_"," ",x)))
}

phy$tip.label<-Xunder(phy$tip.label)


#traits
traits<-readRDS("GUBIC data/trait_plant_all_standardized_names.rds")
traits<-column_to_rownames(traits,var = "accepted_wcvp_name_binomial")
traits<-traits[,1:6]


###need data with 0s, make new long format data
dat.all.long <- dat %>%
  pivot_longer(
    cols = c(buffer_5km_records, urban_records),
    names_to = "location",
    values_to = "records"
  )

##remove species not in phylogeny
tmp<-match(dat.all.long$accepted_taxon_name_binomial,phy$tip.label)
dat.long.sp.remove<-dat.all.long[!is.na(tmp),]
