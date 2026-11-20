library(ape)  
library(betapart) 
library(betareg)
library(corrplot)
library(dplyr)
library(ecodist)
library(emmeans)
library(FD) 
library(hillR)
library(iNEXT)
library(glmmTMB)
library(ggplot2)
library(ggpubr)
library(ggrastr)
library(gridExtra)
library(lme4)
library(lmerTest)
library(mFD)
library(missForest)
library(MuMIn)
library(patchwork)
library(performance)
library(picante)
library(purrr)
#library(rethinking)
library(scales)
library(statmod)
library(tibble)
library(tidyr)
library(vegan)  
library(VennDiagram)



# GUBIC inside/outside prelim analysis ----
# Read in data sets ====
city<-read.csv('data/urban_sf_sgubik.csv')
dat5<-readRDS('data/gbif_gubic_list_per_city_final2_20241015.rds')
dat30<-readRDS('data/gbif_gubic_list_per_city_5_30km_buffer_20241003.rds')
env_diff<-readRDS('data/env_diff.rds')

dat5$ID_HDC_G0_2<-as.factor(dat5$ID_HDC_G0_2)

dat30$ID_HDC_G0_2<-as.factor(dat30$ID_HDC_G0_2)

##similarity between buffers globally, and cities globally, for the taxa groups
spat.dist<-readRDS("data/urban_center_distance.rds")

#merge buffers
dat <- merge(dat5,dat30[, c("ID_HDC_G0_2","accepted_taxon_name_binomial",
  "n_records_buffer_5_30km")],by = c("ID_HDC_G0_2", "accepted_taxon_name_binomial"),
  all = TRUE,sort = FALSE)

##replace NAs with 0s
dat[, c("n_records_urban","n_records_buffer_5km","n_records_buffer_5_30km")] <- 
  lapply(dat[, c("n_records_urban","n_records_buffer_5km","n_records_buffer_5_30km")],
  function(x) replace(x, is.na(x), 0))

B_size<-read.csv('data/city_buffer_areas.csv')

B_env <- merge(env_diff,B_size[, c("ID_HDC_G0_2","area_urban_land_km2",
  "area_buffer_0_5km_land_km2","area_buffer_5_30km_land_km2")], 
    by.x = "city_id",by.y = "ID_HDC_G0_2",all.x = TRUE,sort = FALSE)

#create Pielou's evenness for forest, grassland, crop, builtup land classes for 5 adn 30
# Pielou's evenness for four land-cover categories
landuse_evenness <- function(x) {
  if (anyNA(x) || any(x < 0) || sum(x) <= 0) {
    return(NA_real_)
  }
  
  # Normalize to proportions summing to 1
  p <- x / sum(x)
  # Exclude zeros from logarithm
  p <- p[p > 0]
  # Shannon diversity
  H <- -sum(p * log(p))
  # Pielou's evenness (four possible categories)
  H / log(4)
}

B_env <- B_env %>%rowwise() %>%
  mutate(Evenness_5km = landuse_evenness(c(buffer5km_forestland,buffer5km_cropland,
      buffer5km_builtup,buffer5km_grassland)),
    Evenness_30km = landuse_evenness(c(buffer30km_forestland,buffer30km_cropland,
        buffer30km_builtup,buffer30km_grassland))) %>%
  ungroup()

# Examine distributions
summary(B_env[, c("Evenness_5km", "Evenness_30km")])

##read in phylogeny

phy<-read.tree("data/phylo_gubic.tre")

Xunder<-function(names) {
  as.vector(sapply(names,function(x) sub("_"," ",x)))
}

phy$tip.label<-Xunder(phy$tip.label)


#traits
traits<-readRDS("data/trait_plant_all_standardized_names.rds")
traits<-column_to_rownames(traits,var = "accepted_wcvp_name_binomial")
traits<-traits[,1:6]

##need to scale traits
traits_scaled <- traits

log_traits <- c("diaspore_mass_mg","plant_height_m","leaf_area_mm2")

# Log transformation (handles zero values)
traits_scaled[log_traits] <- lapply(traits_scaled[log_traits],log1p)

# Standardize all six traits to mean 0 and SD 1
traits_scaled[] <- lapply(traits_scaled,function(x) as.numeric(scale(x)))

# Verify
sapply(traits_scaled, sd, na.rm = TRUE)

#add column for Introduced and native
dat$I_N<-"N"
dat$I_N[is.na(dat$status_glonaf)==FALSE]<-"I"

#species found in both urban and buffer (getting rid of singletons)
dat.no.0 <- dat[!is.na(dat$n_records_buffer_5km) &!is.na(dat$n_records_buffer_5_30km) &
    !is.na(dat$n_records_urban) &dat$n_records_buffer_5km > 1 &dat$n_records_buffer_5_30km > 1 &
    dat$n_records_urban > 1,]


dim(dat.no.0) #334,710 species records in all three



#percent change from urban to buffer: positive = lower in urban
#negative = higher in urban; 0 = no difference
dat.no.0$per_change <- 100 * (dat.no.0$n_records_buffer_5km -
    dat.no.0$n_records_urban) / dat.no.0$n_records_urban


## Percent change from urban to 5-30 km buffer
dat.no.0$per_change_30km <- 100 * (dat.no.0$n_records_buffer_5_30km -
    dat.no.0$n_records_urban) / dat.no.0$n_records_urban


###need data with 0s, make new long format data
dat.all.long <- dat %>%pivot_longer(cols = c(
    n_records_buffer_5km,n_records_buffer_5_30km,n_records_urban),
    names_to = "location",values_to = "records")
## Simplify location labels
dat.all.long$location <- recode(dat.all.long$location,
  n_records_buffer_5km = "buffer_5km",n_records_buffer_5_30km = "buffer_5_30km",
  n_records_urban = "urban")

##remove species not in phylogeny
tmp<-match(dat.all.long$accepted_taxon_name_binomial,phy$tip.label)
dat.long.sp.remove<-dat.all.long[!is.na(tmp),]
