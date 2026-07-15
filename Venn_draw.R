###draw Venn diagram for number of species by city records in GUBIC and GBIF
setwd("~/Dropbox/2018WorkingFiles/Marc2018/GUBIC/sGUBIK")
source("In_Out/R/data_import.R")


#introduced only
datI<-dat[dat$I_N=="I",]
datI$ID<-1:nrow(datI)
GBIF <- datI$ID[grepl("GBIF",datI$sources)]
GUBIC <- datI$ID[grepl("GUBIC",datI$sources)]



#venn diagram for GUBIC and GBIF
venn.plot <- venn.diagram(
  x = list("GBIF" = GBIF, "GUBIC" = GUBIC),
  filename = NULL,  
  fill = c("slateblue3", "yellowgreen"),
  alpha = 0.5,
  cat.cex = 2,
  cex=2
)
quartz()
grid::grid.draw(venn.plot)