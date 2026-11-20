###beta diversity
source("01_data_import.R")
source("02_sample_coverage_richness.R")

### 1. within city comparisons ====

# Calculate Sørensen partition for two species assemblages
beta_pair <- function(sp1, sp2) {
  sp1 <- unique(sp1)
  sp2 <- unique(sp2)
  a <- length(intersect(sp1, sp2))
  b <- length(setdiff(sp1, sp2))
  c <- length(setdiff(sp2, sp1))
  
  # Undefined if both assemblages are empty
  if ((2*a + b + c) == 0) {
    return(c( beta_sor = NA_real_,beta_sim = NA_real_,beta_sne = NA_real_))
  }
  
  beta_sor <- (b + c) / (2*a + b + c)
  
  # Simpson turnover
  denom <- a + min(b, c)
  
  beta_sim <- if (denom == 0) if (max(b, c) > 0) 0 else NA_real_
  else min(b, c) / denom
  
  beta_sne <- beta_sor - beta_sim
  
  c(beta_sor = beta_sor, beta_sim = beta_sim,beta_sne = beta_sne)
}

##prepare lists
community_lists <- species_counts %>%
  filter(records > 0,!is.na(accepted_taxon_name_binomial)) %>%
  group_by(ID_HDC_G0_2, Location, Status) %>%
  summarise(species = list(unique(accepted_taxon_name_binomial)),.groups = "drop" )

community_all <- species_counts %>%
  filter(records > 0,!is.na(accepted_taxon_name_binomial)) %>%
  group_by(ID_HDC_G0_2, Location) %>%
  summarise(species = list(unique(accepted_taxon_name_binomial)),.groups = "drop") %>%
  mutate(Status = "All")

community_lists <- bind_rows(community_lists,community_all)

#urban buffer comparisons
urban_lists <- community_lists %>%
  filter(Location == "Urban") %>%
  select(ID_HDC_G0_2,Status,urban_species = species)

buffer_lists <- community_lists %>%
  filter(Location != "Urban") %>%
  transmute(ID_HDC_G0_2,Status,Buffer = Location,buffer_species = species)

beta_within <- buffer_lists %>%
  left_join(urban_lists,by = c("ID_HDC_G0_2", "Status")) %>%
  rowwise() %>%
  mutate(beta = list(if (is.null(urban_species)) {
        c(beta_sor = NA_real_,beta_sim = NA_real_,beta_sne = NA_real_)
      } else beta_pair(urban_species, buffer_species)
    )
  ) %>% ungroup() %>%
  unnest_wider(beta)

saveRDS(beta_within, "data/beta_within_city.rds")
#beta_within<-readRDS("data/beta_within_city.rds")

#beta diversity components
beta_within %>%
  group_by(Buffer, Status) %>%
  summarise(n = sum(!is.na(beta_sor)),Sorensen = mean(beta_sor, na.rm = TRUE),
    Turnover = mean(beta_sim, na.rm = TRUE),Nestedness = mean(beta_sne, na.rm = TRUE),
    .groups = "drop")

## Test sig. diff. for within-city comparisons ----
# Prepare paired comparisons by city
beta_wide <- beta_within %>%
  select(ID_HDC_G0_2, Buffer, Status,beta_sor, beta_sim, beta_sne) %>%
  pivot_longer(cols = starts_with("beta_"), names_to = "Component",
    values_to = "Dissimilarity") %>%
  pivot_wider(id_cols = c(ID_HDC_G0_2, Buffer, Component),
    names_from = Status,values_from = Dissimilarity)

comparisons <- list(c("All", "Native"),c("All", "Introduced"),c("Native", "Introduced"))


###########################
#just for mean and SD
msd<-function(x){
  return(data.frame(mean = apply(x,2,mean,na.rm=TRUE),
      SD = apply(x,2,sd,na.rm=TRUE)))
}

beta_summary <- by(beta_wide[, 4:6],
  list(Buffer = beta_wide$Buffer,Component = beta_wide$Component ), msd)
###########################

# paired test function
paired_test <- function(dat, group1, group2) {
  x <- dat[[group1]]
  y <- dat[[group2]]
  keep <- complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]
  p <- if (length(x) > 0 && any(x != y)) {
    wilcox.test(x, y, paired = TRUE,exact = FALSE)$p.value
  } 
  else NA_real_
  
  data.frame(group1 = group1,group2 = group2,n = length(x),
    median_diff = if (length(x) > 0)median(x - y) else NA_real_, p = p)
}

# Run tests by buffer and component
groups <- split(beta_wide,
  interaction(beta_wide$Buffer,beta_wide$Component,drop = TRUE))

beta_tests <- bind_rows(lapply(groups, function(dat) {
  results <- bind_rows(lapply(comparisons, function(comp) {
    paired_test(dat, comp[1], comp[2]) }))
  
  results$Buffer <- dat$Buffer[1]
  results$Component <- dat$Component[1]
  
  results
})) %>%
  group_by(Buffer, Component) %>%
  mutate(p_adj = p.adjust(p, method = "holm"),significance = case_when(
  is.na(p_adj) ~ "ns",p_adj < 0.001 ~ "***",p_adj < 0.01 ~ "**",p_adj < 0.05 ~ "*",
      TRUE ~ "ns" )) %>%
  ungroup()

beta_tests

## Plot within-city sorensen comparisons for Supp Mat ----
beta_plot <- beta_within %>%
  select(ID_HDC_G0_2,Buffer,Status,beta_sor,beta_sim,beta_sne) %>%
  pivot_longer(cols = c(beta_sor, beta_sim, beta_sne),names_to = "Component",
    values_to = "Dissimilarity" ) %>%
  mutate(Component = factor(Component,levels = c("beta_sor", "beta_sim", "beta_sne"),
      labels = c("A. Sørensen dissimilarity","B. Turnover","C. Nestedness-resultant")),
    Buffer = factor(Buffer,levels = c("0–5 km", "5–30 km")),
    Status = factor(Status,levels = c("All", "Native", "Introduced")))


p_beta <- ggplot(beta_plot,aes(x = Status, y = Dissimilarity, fill = Status)) +
  geom_boxplot(width = 0.65,outlier.alpha = 0.3) +
  facet_grid(Component ~ Buffer) +
  scale_fill_manual(values = c(All = "grey65",Native = "darkolivegreen4",
      Introduced = "darkorange3")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL,y = "Urban–buffer dissimilarity") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",strip.text = element_text(face = "bold"))

p_beta

##add significance labels
sig_labels <- beta_tests %>%
  filter(!is.na(p_adj),significance != "ns") %>%
  mutate( Component = factor(Component, levels = c("beta_sor", "beta_sim", "beta_sne"),
   labels = c("A. Sørensen dissimilarity", "B. Turnover", "C. Nestedness-resultant") ),
  Buffer = factor(Buffer,levels = c("0–5 km", "5–30 km")),
  x1 = case_when(group1 == "All" ~ 1,group1 == "Native" ~ 2,group1 == "Introduced" ~ 3),
  x2 = case_when(group2 == "All" ~ 1,group2 == "Native" ~ 2, group2 == "Introduced" ~ 3),
   y = case_when( group1 == "All" & group2 == "Native" ~ 1.04,
    group1 == "Native" & group2 == "Introduced" ~ 1.12,
     group1 == "All" & group2 == "Introduced" ~ 1.20))


p_beta_sig <- ggplot(beta_plot,aes(x = Status, y = Dissimilarity, fill = Status)) +
  geom_boxplot(width = 0.65,outlier.alpha = 0.3,na.rm = TRUE) +
  facet_grid(Component ~ Buffer) +
  # Horizontal bracket
  geom_segment(data = sig_labels,aes(x = x1, xend = x2, y = y, yend = y),
    inherit.aes = FALSE,linewidth = 0.4) +
  # Left bracket tip
  geom_segment(data = sig_labels,aes(x = x1, xend = x1,y = y, yend = y - 0.015),
    inherit.aes = FALSE,linewidth = 0.4) +
  # Right bracket tip
  geom_segment(data = sig_labels,aes(x = x2, xend = x2, y = y, yend = y - 0.015),
    inherit.aes = FALSE,linewidth = 0.4) +
   # Significance stars
  geom_text(data = sig_labels,aes(x = (x1 + x2) / 2,y = y + 0.025,label = significance),
    inherit.aes = FALSE, size = 4.5) +
  scale_fill_manual(values = c(All = "grey65",Native = "darkolivegreen4",
      Introduced = "darkorange3")) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  coord_cartesian( ylim = c(0, 1.30),clip = "off") +
  labs(x = NULL, y = "Urban–buffer dissimilarity") +
  theme_classic(base_size = 13) +
  theme( legend.position = "none",strip.text = element_text(face = "bold"))

p_beta_sig

ggsave("out/FS4_beta_partition.pdf",p_beta_sig,width = 10,height = 9,device = cairo_pdf)


### distance-decays ====
beta_communities <- community_lists %>%
  filter(Status %in% c("Native", "Introduced"), Location %in% c("Urban", "0–5 km", "5–30 km")) %>%
  transmute( ID_HDC_G0_2 = as.character(ID_HDC_G0_2),Status,Comparison = factor(
      Location,levels = c("Urban", "0–5 km", "5–30 km")),species )
#need to remove two cities
beta_communities <- beta_communities %>%filter(!ID_HDC_G0_2 %in% c("1226", "1639"))


##creat distance decay function
calculate_distance_decay <- function(dat, spat.dist) {
  ids <- as.character(dat$ID_HDC_G0_2)
   # Keep cities present in the distance matrix
  keep <- ids %in% rownames(spat.dist) & ids %in% colnames(spat.dist)
  
  dat <- dat[keep, , drop = FALSE]
  ids <- ids[keep]
  n <- nrow(dat)
  
  if (n < 2) return(tibble())
  pairs <- utils::combn(seq_len(n), 2)
  results <- lapply(seq_len(ncol(pairs)), function(k) {
    i <- pairs[1, k]
    j <- pairs[2, k]
    
    beta <- beta_pair(dat$species[[i]],dat$species[[j]])
    
   tibble(City1 = ids[i],City2 = ids[j],
      distance_km = as.numeric(spat.dist[ids[i], ids[j]]),
      beta_sor = unname(beta["beta_sor"]),
      beta_sim = unname(beta["beta_sim"]),
      beta_sne = unname(beta["beta_sne"]))})
  
  bind_rows(results)
}

beta_distance <- beta_communities %>%group_by(Comparison, Status) %>%
  group_modify( ~ calculate_distance_decay(.x, spat.dist)) %>%
  ungroup()

saveRDS(beta_distance, "data/beta_distance_decay.rds")

##plot distance decays

beta_distance <- beta_distance %>%
  mutate( Comparison = factor(Comparison,
      levels = c("Urban", "0–5 km", "5–30 km"),
      labels = c( "A. Urban","B. 0–5 km buffer","C. 5–30 km buffer")),
    Status = factor(Status,levels = c("Native", "Introduced"))

p_distance <- ggplot(beta_distance,aes(x = distance_km,y = beta_sor,colour = Status)) +
  #geom_point( alpha = 0.04,size = 0.3) + #makes large file
  geom_point_rast(alpha = 0.05,size = 0.3,raster.dpi = 300)+
  geom_smooth(data = subset(beta_distance, Status == "Native"),method = "lm",
    se = FALSE,colour = "#344D20",linewidth = 1.6) +
  geom_smooth(data = subset(beta_distance, Status == "Introduced"),method = "lm",
    se = FALSE,colour = "#A34700",linewidth = 1.6) +
  scale_colour_manual(values = c(Native = "darkolivegreen4",Introduced = "darkorange3" ))+
  facet_wrap(~ Comparison,ncol = 3) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(labels = scales::label_number(scale = 0.001,suffix = "k")) +
  labs(x = "Geographic distance between cities (km)",
    y = "Sørensen dissimilarity",colour = NULL) +
  theme_classic(base_size = 13) +
  theme(strip.text = element_text(face = "bold",size = 12 ),
    legend.position = "top",legend.title = element_blank())

p_distance

ggsave("out/FS3_distance_decay.pdf",p_distance,width = 9,height = 4.5,device = cairo_pdf)

###fit beta regressions

beta_reg_dat <- beta_distance %>%
  mutate(Comparison = sub("^[ABC]\\. ", "", as.character(Comparison)),
    Comparison = sub(" buffer$", "", Comparison),
    Comparison = factor(Comparison,levels = c("Urban", "0–5 km", "5–30 km") ),
    Status = factor(Status,levels = c("Native", "Introduced") ),
    log_dist = log1p(distance_km)) %>%
  filter(!is.na(beta_sor),is.finite(log_dist)) %>%
  mutate(
    # Transform [0,1] to (0,1)
    beta_adj = (beta_sor * (n() - 1) + 0.5) / n(),
    
    # Centre and scale for numerical stability
    dist_scaled = as.numeric(scale(log_dist)) )

groups <- split(beta_reg_dat,interaction(beta_reg_dat$Comparison,
    beta_reg_dat$Status,drop = TRUE))

models <- lapply(groups, function(dat) {
  glmmTMB(beta_adj ~ dist_scaled + I(dist_scaled^2),
    family = beta_family(link = "logit"),data = dat)})


coef_df <- bind_rows(
  lapply(names(models), function(nm) {
    m <- models[[nm]]
    cf <- summary(m)$coefficients$cond
    data.frame(Group = nm,Term = rownames(cf),Estimate = cf[, "Estimate"],
      SE = cf[, "Std. Error"],Lower = cf[, "Estimate"] - 1.96 * cf[, "Std. Error"],
      Upper = cf[, "Estimate"] + 1.96 * cf[, "Std. Error"] )})) %>%
  filter(Term != "(Intercept)") %>%
  mutate(Comparison = sub("\\.(Native|Introduced)$", "", Group),
    Status = sub("^.*\\.", "", Group),Term = recode(Term,
      dist_scaled = "Geographic distance",`I(dist_scaled^2)` = "Distance²"))

p_fig3x <- ggplot(coef_df,aes(x = Estimate,y = Comparison,colour = Status)) +
  geom_vline(xintercept = 0,linetype = "dashed",colour = "grey60") +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper),
    position = position_dodge(width = 0.5),height = 0.15) +
  geom_point(position = position_dodge(width = 0.5),size = 3) +
  facet_wrap(~ Term, scales = "free_x") +
  scale_colour_manual(values = c(Native = "darkolivegreen4",Introduced = "darkorange3" ) ) +
  labs(x = "Beta regression coefficient (logit scale)",y = NULL,colour = NULL) +
  theme_classic(base_size = 14) +
  theme( legend.position = "top")


p_fig3x

ggsave("out/FSX_Sor_beta_coef.pdf",p_fig3x,width = 6,height = 6,device = cairo_pdf)

###this figure examines beta regression coefficients - and native dissimiliarity
##increases more rapidly, but difficult to interpret for reader

##using fixed distances to show differences

log_mean <- mean(beta_reg_dat$log_dist)
log_sd   <- sd(beta_reg_dat$log_dist)

# Convert geographic distance to the model's predictor scale
scale_distance <- function(km) {(log1p(km) - log_mean) / log_sd}

# Predict Sørensen dissimilarity on the response scale
predict_sor <- function(model, km) {
  newdat <- data.frame(dist_scaled = scale_distance(km))
  
  as.numeric( predict(model,newdata = newdat,type = "response" ))
}

effect_df <- bind_rows(lapply(names(models), function(nm) {
    m <- models[[nm]]
    pred_500  <- predict_sor(m, 500)
    pred_1000 <- predict_sor(m, 1000)
    pred_5000 <- predict_sor(m, 5000)
    
    data.frame(Group = nm,Sorensen_1000 = pred_1000, 
      Distance_effect = pred_5000 - pred_500)})) %>%
  mutate( Comparison = sub("\\.(Native|Introduced)$", "", Group),
    Status = sub("^.*\\.", "", Group),Comparison = factor(Comparison,
      levels = c("Urban", "0–5 km", "5–30 km")),
    Status = factor(Status,levels = c("Native", "Introduced")))

effect_df


###categorical grouping of beta diversity

glob.dat3 <- beta_distance %>%
  mutate( Comp = sub("^[ABC]\\. ", "", as.character(Comparison)),
  Comp = sub(" buffer$", "", Comp),Comp = factor(Comp,levels = c("Urban", "0–5 km", "5–30 km")),
    Taxa = factor(Status,levels = c("Native", "Introduced"))) %>%
  filter(!is.na(beta_sor)) %>%
  mutate(Sor_adj = pmin(pmax(beta_sor, 1e-6), 1 - 1e-6))

modTUB <- glmmTMB(Sor_adj ~ Taxa * Comp,data = glob.dat3,family = beta_family())

summary(modTUB)

emm <- emmeans(modTUB,~ Taxa * Comp,type = "response")

emm

pairs(emm, by = "Comp", type = "response")

emm_df <- as.data.frame(emm)

p_fig3 <- ggplot(emm_df,aes(x = Comp,y = response,colour = Taxa,group = Taxa)) +
  geom_line(linewidth = 1.4) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL,ymax = asymp.UCL),width = 0.1,linewidth = 0.8) +
  scale_colour_manual(values = c(Native = "darkolivegreen4",Introduced = "darkorange3" )) +
  labs(x = "Comparison type",y = "Predicted Sørensen dissimilarity",colour = NULL) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top")

##add sig diffs
sig_df <- as.data.frame(
pairs(emm, by = "Comp", adjust = "none")) %>%
  mutate(p.adj = p.adjust(p.value, method = "holm"),
    label = case_when(p.adj < 0.001 ~ "***",p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",TRUE ~ "ns" )) %>%
  select(Comp, p.adj, label)

sig_df

sig_pos <- emm_df %>% group_by(Comp) %>%
  summarise(y = max(response) + 0.008,.groups = "drop") %>%
  left_join(sig_df, by = "Comp")

p_fig3 <- p_fig3 +
  geom_text( data = sig_pos,aes(x = Comp, y = y, label = label),inherit.aes = FALSE,
    size = 6) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))



p_fig3
ggsave("out/Fig3_Sor_comps.pdf",p_fig3,width = 4,height = 6,device = cairo_pdf)
