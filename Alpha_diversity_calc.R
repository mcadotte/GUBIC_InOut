setwd("~/Dropbox/2018WorkingFiles/Marc2018/GUBIC/sGUBIK")
source("In_Out/R/data_import.R")



subTree<- function(phy,species) {
  
  all.in<-is.na(match(species,phy$tip.label))
  
  if (sum(all.in)>0) 
    paste("species",as.character(species[all.in]),
          "not in phylogeny")
  
  dropme<-phy$tip.label[!phy$tip.label %in% unique(species)]
  sub.tr<-drop.tip(phy,dropme)
  
}

##delta_mpd and delta_richness for each city for exotics and natives (urban - buffer)

delta.PD<-data.frame(Inv.D.rich = 0, Inv.D.PD = 0, Inv.D.MPD = 0, 
                     Inv.D.FDis = 0, Inv.D.FD = 0, Inv.D.FDeve = 0,
                     Nat.D.rich = 0, Nat.D.PD = 0, Nat.D.MPD = 0,
                     Nat.D.FDis = 0, Nat.D.FD = 0, Nat.D.FDeve = 0)

CWM.list.i<-list(NULL)
CWM.list.n<-list(NULL)
FD_spp<-list(NULL)

#I standardized PD by urban richness, so these are proportion increase/decrease for richness
for (i in 1:length(city$ID_HDC_G0_2)){
  sub<-dat.long.sp.remove[dat.long.sp.remove$ID_HDC_G0_2==city$ID_HDC_G0_2[i],]
  if (nrow(sub)>100) {
    sub.i<-sub[sub$I_N=="I",]
    sub.n<-sub[sub$I_N=="N",]
    samp.i<-sample2matrix(sub.i[c("location","records","accepted_taxon_name_binomial")])
    samp.n<-sample2matrix(sub.n[c("location","records","accepted_taxon_name_binomial")])
    #remove 0 species (summed across both)
    samp.i.n0<-samp.i[,apply(samp.i,2,sum) > 0]
    samp.n.n0<-samp.n[,apply(samp.n,2,sum) > 0]
    if (sum(samp.i[1,])>10 & sum(samp.n[1,])>10)
    {
      sub.phy<-subTree(phy,unique(sub$accepted_taxon_name_binomial))
      #pd.tmp.i<-ses.pd(samp.i,phy,runs=100)
      #pd.tmp.n<-ses.pd(samp.n,phy,runs=100)
      pd.tmp.i<-pd(samp.i.n0,sub.phy,include.root = FALSE)
      pd.tmp.n<-pd(samp.n,sub.phy,include.root = FALSE)
      
      mpd.tmp.i<-mpd(samp.i.n0,cophenetic(sub.phy))
      mpd.tmp.n<-mpd(samp.n,cophenetic(sub.phy))
      
      tr.i<-traits[match(colnames(samp.i.n0),rownames(traits)),]
      tr.i<-na.omit(tr.i)
      shared<-sort(intersect(rownames(tr.i), colnames(samp.i.n0)))
      tr.i <- tr.i[shared, ]
      samp.i.n0 <- samp.i.n0[, shared]
      
      tr.n<-traits[match(colnames(samp.n.n0),rownames(traits)),]
      tr.n<-na.omit(tr.n)
      shared.n<-sort(intersect(rownames(tr.n), colnames(samp.n.n0)))
      tr.n <- tr.n[shared.n, ]
      samp.n.n0 <- samp.n.n0[, shared.n]
      
      if(apply(samp.i.n0,1,sum)[2]>0 & apply(samp.i.n0,1,sum)[1]>0){
        fd.tmp.i<-dbFD(tr.i,samp.i.n0,calc.CWM = TRUE,calc.FDiv = TRUE)
      }
      else{
        fd.tmp.i<-list(FDis = c(0,0),FDiv =c(0,0), FEve = c(0,0), CWM=0, nbsp = c(1,1))
        names(fd.tmp.i$nbsp)<-c("buffer_5km_records","urban_records")
      }
      
      if(apply(samp.n.n0,1,sum)[2]>0 & apply(samp.n.n0,1,sum)[1]>0){
        fd.tmp.n<-dbFD(tr.n,samp.n.n0,calc.CWM = TRUE,calc.FDiv = TRUE)
      }
      else{
        fd.tmp.n<-list(FDis = c(0,0),FDiv =c(0,0), FEve = c(0,0), CWM=0, nbsp = c(1,1))
        names(fd.tmp.n$nbsp)<-c("buffer_5km_records","urban_records")
      }
      
      #tr.n<-traits[match(colnames(samp.n.n0),rownames(traits)),]
      # fd.tmp.n<-dbFD(tr.n,samp.n.n0,calc.CWM = FALSE,calc.FDiv = FALSE)
      
      delta.PD[i,]<-c((pd.tmp.i[2,2]-pd.tmp.i[1,2]),
                      (pd.tmp.i[2,1]/pd.tmp.i[2,2])-(pd.tmp.i[1,1]/pd.tmp.i[1,2]),
                      mpd.tmp.i[2]-mpd.tmp.i[1],
                      fd.tmp.i$FDis[2]-fd.tmp.i$FDis[1],
                      fd.tmp.i$FDiv[2]/fd.tmp.i$nbsp[2]-fd.tmp.i$FDiv[1]/fd.tmp.i$nbsp[1],
                      fd.tmp.i$FEve[2]-fd.tmp.i$FEve[1],
                      (pd.tmp.n[2,2]-pd.tmp.n[1,2]),
                      (pd.tmp.n[2,1]/pd.tmp.n[2,2])-(pd.tmp.n[1,1]/pd.tmp.n[1,2]),
                      mpd.tmp.n[2]-mpd.tmp.n[1],
                      fd.tmp.n$FDis[2]-fd.tmp.n$FDis[1],
                      fd.tmp.n$FDiv[2]/fd.tmp.n$nbsp[2]-fd.tmp.n$FDiv[1]/fd.tmp.n$nbsp[1],
                      fd.tmp.n$FEve[2]-fd.tmp.n$FEve[1])
      
      CWM.list.i[[i]]<-fd.tmp.i$CWM
      CWM.list.n[[i]]<-fd.tmp.n$CWM
      FD_spp[[i]]<-data.frame(N_I = c("N","N","I","I"),
                              Num_SP=c(fd.tmp.n$nbsp,fd.tmp.i$nbsp),
                              Loc=c(names(fd.tmp.n$nbsp),names(fd.tmp.i$nbsp)))
    }
    
    
    
  }
  print(paste("iter",i,"of",length(city$ID_HDC_G0_2),"done"))
}

delta.PD.n<-na.omit(delta.PD)
summary(delta.PD.n)
all.out<-list(Metrics.cities=delta.PD.n,CWM.Nat=CWM.list.n,CWM.Int=CWM.list.i,FD.city=FD_spp)
saveRDS(all.out,"In_Out/data/R_Phy_Trait_output.RData")

########


