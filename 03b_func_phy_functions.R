###functional and phylogenetic functions
require(dplyr)
require(ape)
require(iNEXT)


#phylo functions
mean_MPD <- function(species, dmat) {
  
  species <- intersect(unique(species), rownames(dmat))
  
  if (length(species) < 2) return(NA_real_)
  
  d <- dmat[species, species, drop = FALSE]
  mean(d[upper.tri(d)])
}

calculate_city_MPD <- function(city_dat, phy, nrep = 100) {
  
  city_species <- unique(city_dat$accepted_taxon_name_binomial[
      !is.na(city_dat$records) & city_dat$records > 0] )
  
  city_species <- intersect(city_species, phy$tip.label)
  
  if (length(city_species) < 2) return(NULL)

    # Prune tree to city
  phy_city <- keep.tip(phy, city_species)
  dmat <- cophenetic.phylo(phy_city)
  
  # Process each location × status
  results <- city_dat %>%group_by(Location, Status) %>%
    group_modify(function(df, key) {
      
      # Observed species and occurrence counts
      abund <- df %>%filter(!is.na(records),records > 0,
          !is.na(accepted_taxon_name_binomial)) %>%
        group_by(accepted_taxon_name_binomial) %>%
        summarise(records = sum(records),.groups = "drop")
      
      species <- abund$accepted_taxon_name_binomial
      counts <- abund$records
      
      # Observed MPD
      MPD_obs <- mean_MPD(species, dmat)
      
      # Only interpolate if coverage is sufficient
      coverage <- iNEXT::DataInfo(counts,datatype = "abundance")$SC
      
      if (length(coverage) == 0 ||is.na(coverage) ||coverage < 0.90) {
          return(data.frame(
          MPD_obs = MPD_obs,
          MPD_90 = NA_real_))
      }
      
      # Estimate sample size at 90% coverage
      est <- iNEXT::estimateD(counts,datatype = "abundance",q = 0,
        base = "coverage",level = 0.90,nboot = 0)
      
      n90 <- round(est$m[1])
      
      if (!is.finite(n90) || n90 < 2 || n90 > sum(counts)) {
        
        return(data.frame(MPD_obs = MPD_obs,MPD_90 = NA_real_))
      }
      
      # Subsample without replacement
      individuals <- rep(species, counts)
      
      MPD_rep <- replicate(nrep, {
        sampled_species <- unique(sample(individuals, n90, replace = FALSE) )
        mean_MPD(sampled_species, dmat)})
      
      data.frame(MPD_obs = MPD_obs,MPD_90 = mean(MPD_rep, na.rm = TRUE))
    }) %>%
    ungroup()
  
  results
}

#functional functions (uses mean_MPD from above)
calculate_city_MFD <- function(city_dat, traits, nrep = 100) {
  
  city_species <- unique(city_dat$accepted_taxon_name_binomial[
      !is.na(city_dat$records) & city_dat$records > 0])
  
  # Retain species with trait information
  city_species <- intersect(city_species, rownames(traits))
  
  # Retain complete trait records
  city_traits <- traits[city_species, , drop = FALSE]
  city_traits <- city_traits[complete.cases(city_traits),,drop = FALSE]
  
  if (nrow(city_traits) < 2) return(NULL)
  
  # Functional distance matrix for this city only
  dmat <- as.matrix(dist(city_traits, method = "euclidean"))
  
  results <- city_dat %>%
    group_by(Location, Status) %>%
    group_modify(function(df, key) {
        abund <- df %>%
        filter(!is.na(records),records > 0,!is.na(accepted_taxon_name_binomial)) %>%
        group_by(accepted_taxon_name_binomial) %>%
        summarise(records = sum(records),.groups = "drop")
      
      species <- abund$accepted_taxon_name_binomial
      counts <- abund$records
      
      # Observed MFD
      MFD_obs <- mean_MPD(species, dmat)
      
      # Observed sample coverage
      coverage <- iNEXT::DataInfo(counts,datatype = "abundance")$SC
      
      if (length(coverage) == 0 ||is.na(coverage) ||coverage < 0.90) {
        
        return(data.frame(MFD_obs = MFD_obs,MFD_90 = NA_real_))}
      
      # Sample size corresponding to 90% coverage
      est <- iNEXT::estimateD(counts,datatype = "abundance",q = 0,
        base = "coverage",level = 0.90,nboot = 0)
      
      n90 <- round(est$m[1])
      
      if (!is.finite(n90) || n90 < 2 || n90 > sum(counts)) {
          return(data.frame(
          MFD_obs = MFD_obs,
          MFD_90 = NA_real_)) }
      
      # Random subsampling
      individuals <- rep(species, counts)
      
      MFD_rep <- replicate(nrep, {
        sampled_species <- unique(sample(individuals, n90, replace = FALSE))
        mean_MPD(sampled_species, dmat) })
      
      data.frame(MFD_obs = MFD_obs,
        MFD_90 = if (all(is.na(MFD_rep))) NA_real_
        else mean(MFD_rep, na.rm = TRUE))}) %>%
    ungroup()
  
  results
}
