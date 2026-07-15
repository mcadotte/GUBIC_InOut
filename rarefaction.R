###rarefaction
setwd("~/Dropbox/2018WorkingFiles/Marc2018/GUBIC/sGUBIK")
source("In_Out/R/data_import.R")

library(iNEXT)
library(vegan)

#simplify data
dat_sum <- dat.no.0.long %>%
  group_by(ID_HDC_G0_2,location,I_N,accepted_taxon_name_binomial) %>%
  summarise(records = sum(records),.groups = "drop")

#split data matrices
comm_list <- dat_sum %>%
  pivot_wider(
    names_from  = accepted_taxon_name_binomial,
    values_from = records,
    values_fill = 0 ) 

meta <- comm_list %>%
  select(ID_HDC_G0_2,location,I_N)

comm_mat <- comm_list %>%
  select(-ID_HDC_G0_2,-location,-I_N) %>%
  as.matrix()


#check sample sizes, some very low
N <- rowSums(comm_mat)
summary(N)
hist(log(N))
length(N)
length(N[N>400]) #so will use 400 as cut off.
#will use 

keep <- N >= 400
comm_mat <- comm_mat[keep,]
meta    <- meta[keep, ]
N         <- N[keep]

rare_S <- rarefy(comm_mat, sample = min(N))



Rare.dat <- meta %>%
  mutate(
    N_total = N,
    S_rare = rare_S
  )

write.csv(Rare.dat,"In_Out/R/Rare.dat.csv")



###let's redo figure 2
Rare.dat<-read.csv("In_Out/R/Rare.dat.csv",row.names = 1)
Rare.dat$log.records<-log(Rare.dat$N_total)

delta_S <- Rare.dat %>%
  select(ID_HDC_G0_2, I_N, location,  S_rare) %>%
  pivot_wider(
    names_from  = location,
    values_from = S_rare
  ) %>%
  mutate(
    delta_S = buffer_5km_records - urban_records
  )


delta_S %>%
  group_by(I_N) %>%
  summarise(
    mean_delta = mean(delta_S, na.rm = TRUE),
    sd_delta   = sd(delta_S, na.rm = TRUE),
    n          = n(),
    se         = sd_delta / sqrt(n),
    t_crit     = qt(0.975, df = n - 1),
    ci_lower   = mean_delta - t_crit * se,
    ci_upper   = mean_delta + t_crit * se
  )
