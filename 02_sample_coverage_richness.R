
###Estimate sample coverage
source("01_data_import.R")


# Three sampling locations
record_cols <- c("n_records_urban","n_records_buffer_5km",
                 "n_records_buffer_5_30km")

# Long-format species occurrence data
species_long <- dat %>%
  mutate(across(all_of(record_cols), ~replace_na(.x, 0))) %>%
  pivot_longer(cols = all_of(record_cols),names_to = "Location",
               values_to = "records") %>%
  mutate(Location = recode(Location,n_records_urban = "Urban",
                           n_records_buffer_5km = "0–5 km",n_records_buffer_5_30km = "5–30 km"),
         Status = recode(I_N,I = "Introduced",N = "Native"))

# Aggregate duplicate city × species × location entries, if any
species_counts <- species_long %>%
  group_by(ID_HDC_G0_2,Location,Status,
           accepted_taxon_name_binomial ) %>%
  summarise(records = sum(records, na.rm = TRUE),.groups = "drop")

##no zeros
species_counts_N0 <-species_counts[species_counts$records>0,]

#Test differences between native and introduced for each location
#urban
ks.test(species_counts_N0$records[species_counts_N0$Status=="Introduced" & 
                                    species_counts_N0$Location=="Urban"],
        species_counts_N0$records[species_counts_N0$Status=="Native" & 
                                    species_counts_N0$Location=="Urban"])
#0-5 km
ks.test(species_counts_N0$records[species_counts_N0$Status=="Introduced" & 
                                    species_counts_N0$Location=="0–5 km"],
        species_counts_N0$records[species_counts_N0$Status=="Native" & 
                                    species_counts_N0$Location=="0–5 km"])
#5-30 km
ks.test(species_counts_N0$records[species_counts_N0$Status=="Introduced" & 
                                    species_counts_N0$Location=="5–30 km"],
        species_counts_N0$records[species_counts_N0$Status=="Native" & 
                                    species_counts_N0$Location=="5–30 km"])

##anovas for no. obs. across spatial zones for natives and introduced
m.N<-lm(records~Location, data = species_counts_N0[species_counts_N0$Status=="Native",])
anova(m.N)

m.I<-lm(records~Location, data = species_counts_N0[species_counts_N0$Status=="Introduced",])
anova(m.I)

#add area
# Convert area data to long format
area_long <- B_size %>% select(ID_HDC_G0_2,area_urban_land_km2,
    area_buffer_0_5km_land_km2,area_buffer_5_30km_land_km2) %>%
  pivot_longer(cols = starts_with("area_"),
    names_to = "Location",values_to = "Area_km2" ) %>%
  mutate(ID_HDC_G0_2 = as.character(ID_HDC_G0_2),
    Location = recode(Location,"area_urban_land_km2" = "Urban",
      "area_buffer_0_5km_land_km2" = "0–5 km","area_buffer_5_30km_land_km2" = "5–30 km" ))

# Join area to species data
species_counts_N0 <- species_counts_N0 %>%
  mutate(ID_HDC_G0_2 = as.character(ID_HDC_G0_2)) %>%
  left_join(area_long,by = c("ID_HDC_G0_2", "Location"))

city_records <- species_counts_N0 %>%
  group_by(ID_HDC_G0_2, Location, Status) %>%
  summarise(total_records = sum(records, na.rm = TRUE),
    richness = n_distinct(accepted_taxon_name_binomial),
    Area_km2 = first(Area_km2),
    .groups = "drop")

p_RA<-ggplot(city_records,aes(x=Area_km2,y=total_records))+
  geom_point(alpha = 0.25,size = 1.5,colour = "grey20" ) +
  geom_smooth(method = "lm",se = TRUE,colour = "steelblue4",fill = "steelblue2",
              alpha = 0.25,linewidth = 1.2) +
  labs(x = "Log(area)",y = "Total records") +
  scale_x_log10()+
  scale_y_log10()+
  theme_classic(base_size = 14) +
  theme(strip.text = element_text(size = 13, face = "bold"),
        panel.spacing = unit(1.2, "lines"))

p_RA


p_recs<-ggplot(city_records,aes(y=total_records,x=Location,fill = Status))+
  geom_boxplot()+
  scale_y_log10()+
  scale_fill_manual(values = c(Introduced = "darkorange3",Native = "darkolivegreen4")) +
  labs(y="log(Total recards)")+
  theme_classic(base_size = 14) +
  theme(strip.text = element_text(size = 13, face = "bold"),
        panel.spacing = unit(1.2, "lines"))

m.recs<-lm(log(total_records)~Location+Status+Location*Status,data=city_records)
summary(m.recs)
anova(m.recs)
  

##correct records/area
city_records$rec.area<-city_records$total_records/city_records$Area_km2

p_cor.rec<-ggplot(city_records,aes(y=rec.area,x=Location,fill = Status))+
  geom_boxplot()+
  scale_y_log10()+
  scale_fill_manual(values = c(Introduced = "darkorange3",Native = "darkolivegreen4")) +
  labs(y="log(Recards/Area)")+
  theme_classic(base_size = 14) +
  theme(strip.text = element_text(size = 13, face = "bold"),
        panel.spacing = unit(1.2, "lines"))

ggsave( "out/Records_by_area.pdf",p_cor.rec,width = 6, height = 5,device = cairo_pdf)

P.rec.combo<-((p_recs | p_cor.rec))+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

ggsave( "out/Records.pdf",P.rec.combo,width = 12, height = 5,device = cairo_pdf)


#Tests
locs<-unique(city_records$Location)
km5<-city_records |> filter(Location == locs[1])
summary(lm(rec.area~Status,data=km5))

km30<-city_records |> filter(Location == locs[2])
summary(lm(rec.area~Status,data=km30))

kmurb<-city_records |> filter(Location == locs[3])
summary(lm(rec.area~Status,data=kmurb))

r.int<-city_records |> filter(Status == "Introduced")
summary(lm(rec.area~Location,data=r.int))

r.nat<-city_records |> filter(Status == "Native")
summary(lm(rec.area~Location,data=r.nat))

#full model

m.cor.rec<-lm(log(rec.area)~Location+Status+Location*Status,data=city_records)
summary(m.cor.rec)
anova(m.cor.rec)

####summary stats

# Function to summarize observations and richness
summarize_species <- function(dat) {
  dat <- dat %>%
    filter(!is.na(accepted_taxon_name_binomial),!is.na(records),records > 0)
  
  tibble(observations = sum(dat$records),total_species = n_distinct(dat$accepted_taxon_name_binomial),
    
    introduced_species = n_distinct(dat$accepted_taxon_name_binomial[dat$Status == "Introduced"]),
    
    native_species = n_distinct(dat$accepted_taxon_name_binomial[dat$Status == "Native"]))
}

# Overall summary
overall <- summarize_species(species_counts) %>% mutate(Location = "All locations")

# Summaries for individual locations
by_location <- species_counts %>%
  group_by(Location) %>%
  group_modify(~ summarize_species(.x)) %>%
  ungroup()

# Combine into one table
species_summary <- bind_rows(overall,by_location) %>%
  mutate(Location = factor(Location,levels = c("All locations","Urban","0–5 km","5–30 km"))) %>%
  arrange(Location) %>%
  select(Location,observations,total_species,introduced_species,native_species)

species_summary

write.csv(species_summary,"out/full_summary.csv",row.names = FALSE)

# Plot species-sampling curves
# Prepare species occurrence data
sp_dat <- species_counts %>% filter( Location %in% c("Urban", "0–5 km", "5–30 km"),
    Status %in% c("Native", "Introduced"), records > 0,!is.na(accepted_taxon_name_binomial) ) %>%
  distinct( ID_HDC_G0_2,Location, Status, accepted_taxon_name_binomial )

# Function to calculate species accumulation curves
accum_curve <- function(dat, location, status, nperm = 500) {
  sub <- dat %>%
    filter(Location == location, Status == status) %>%
    mutate(present = 1L) %>%
    pivot_wider(id_cols = ID_HDC_G0_2,names_from = accepted_taxon_name_binomial,
                values_from = present,values_fill = 0)
  
  comm <- as.matrix(sub[, -1, drop = FALSE])
  
  set.seed(123)
  acc <- specaccum(comm,method = "random",permutations = nperm )
  
  data.frame(Cities = acc$sites,Richness = acc$richness,SD = acc$sd,Location = location,
    Status = status)
}

locations <- c("Urban", "0–5 km", "5–30 km")
statuses <- c("Native", "Introduced")

accum_dat <- bind_rows(
  lapply(statuses, function(s) {
    bind_rows(lapply(locations, function(loc) {accum_curve(sp_dat, loc, s)}))
  }))

accum_dat <- accum_dat %>%
  mutate(Location = factor(Location, levels = locations),
    Status = factor(Status, levels = statuses),
    lower = pmax(0, Richness - 1.96 * SD),
    upper = Richness + 1.96 * SD )

p_accum <- ggplot(accum_dat,aes(x = Cities, y = Richness,colour = Location, fill = Location)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),alpha = 0.15,colour = NA) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ Status, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c("Urban" = "#37658C","0–5 km" = "#219FA5",
                                 "5–30 km" = "#ACF6F0") ) +
  scale_fill_manual(values = c("Urban" = "#37658C","0–5 km" = "#219FA5",
                               "5–30 km" = "#ACF6F0")) +
  labs(x = "Number of cities",y = "Cumulative species richness",
    colour = NULL,fill = NULL) +
  theme_classic(base_size = 17) +
  theme(legend.position = "top",legend.title = element_blank(),
        axis.text.x = element_text(size = 12),plot.tag = element_text(face = "bold", size = 16) ) +
  theme(strip.text = element_text(face = "bold", size = 14),
    panel.spacing = grid::unit(1, "lines"))

p_accum

ggsave( "out/FS_species_accumulation_cities.pdf",p_accum,width = 6, height = 5,device = cairo_pdf)


# Estimate sample coverage from record frequencies
coverage_est <- function(n) {
  n <- n[n > 0]
  N <- sum(n)
  if (N < 2) return(NA_real_)
  
  f1 <- sum(n == 1)
  f2 <- sum(n == 2)
  
  A <- if (f2 > 0) (N - 1) * f1 / ((N - 1) * f1 + 2 * f2) 
  else (N - 1) * (f1 - 1) /((N - 1) * (f1 - 1) + 2)
  
  if (f1 == 0) return(1)
  
  1 - (f1 / N) * A
}

sampling_summary <- species_counts %>%
  group_by(ID_HDC_G0_2, Location, Status) %>%
  summarise(total_records = sum(records),richness = sum(records > 0),
            singletons = sum(records == 1),doubletons = sum(records == 2),
            coverage = coverage_est(records),.groups = "drop")

# Check
head(sampling_summary)
summary(sampling_summary$coverage)



##plot coverage by city
quartz()

p_covs<-ggplot( sampling_summary,aes(x = Location, y = coverage, fill = Status)) +
  geom_boxplot(position = position_dodge(0.75),outlier.alpha = 0.3) +
  annotate("text",x = 4.5, y= 0.94,label="90% coverage",colour = "grey35",hjust=1)+
  geom_hline(yintercept = 0.9, linetype = "dashed", colour = "grey35", linewidth=1.5)+
  scale_fill_manual(values = c(Introduced = "darkorange3",Native = "darkolivegreen4")) +
  labs(x = NULL,y = "Estimated sample coverage",fill = NULL) +
  theme_minimal(base_size = 14)

ggsave("out/FSX_sample_coverage_cities.pdf", p_covs,width = 7,height = 5)

##averages: total records, species richness, singletons, and doubletons for 
#urban, 5km, and 30km and split between native and introduced

# Function to calculate mean and 95% CI
mean_ci <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  
  if (n < 2) return(NA_character_)
  
  m <- mean(x)
  se <- sd(x) / sqrt(n)
  ci <- qt(0.975, df = n - 1) * se
  
  sprintf("%.1f (%.1f–%.1f)", m, m - ci, m + ci)
}

# Create summary table
sampling_table <- sampling_summary %>%
  group_by(Location, Status) %>% summarise(n_cities = n(),total_records = mean_ci(total_records),
    richness = mean_ci(richness),singletons = mean_ci(singletons),doubletons = mean_ci(doubletons),.groups = "drop" ) %>%
  mutate(Location = factor(Location,levels = c("Urban", "0–5 km", "5–30 km")),
    Status = factor(Status,levels = c("Native", "Introduced")) ) %>%
  arrange(Location, Status) %>%
  rename(`Location` = Location,`Status` = Status,`Cities (n)` = n_cities,
    `Total records` = total_records,`Species richness` = richness,
    `Singletons` = singletons,`Doubletons` = doubletons)

sampling_table
write.csv(sampling_table,"out/sampling_table.csv",row.names = FALSE)

##check number of cities with different coverages
coverage_check <- sampling_summary %>%
  group_by(Location, Status) %>%
  summarise(n = n(),pct_90 = mean(coverage >= 0.90, na.rm = TRUE) * 100,
    pct_95 = mean(coverage >= 0.95, na.rm = TRUE) * 100,
    pct_99 = mean(coverage >= 0.99, na.rm = TRUE) * 100,.groups = "drop")

coverage_check

write.csv(coverage_check,"out/coverage_check.csv",row.names = FALSE)

coverage_matched <- sampling_summary %>%
  mutate(adequate = !is.na(coverage) & coverage >= 0.90) %>%
  select(ID_HDC_G0_2, Status, Location, adequate) %>%
  pivot_wider(names_from = Location,values_from = adequate,values_fill = FALSE) %>%
  mutate(all_90 = Urban & `0–5 km` & `5–30 km`)

coverage_matched %>%
  group_by(Status) %>%
  summarise(n_cities = n(),n_all_90 = sum(all_90),pct_all_90 = round(100 * mean(all_90), 1),
    .groups = "drop")

##run rarefaction
# Eligible city × status combinations
eligible <- coverage_matched %>%filter(all_90) %>%select(ID_HDC_G0_2, Status)

# Restrict species counts to eligible communities
species_rare <- species_counts %>%semi_join(eligible,by = c("ID_HDC_G0_2", "Status"))

# Create community list
community_list <- species_rare %>%group_by(ID_HDC_G0_2, Location, Status) %>%group_split()

# Rarefaction function
rare90 <- function(records) {
  records <- records[is.finite(records) & records > 0]
  
  if (sum(records) < 2) return(NA_real_)
  
  out <- iNEXT::estimateD(x = records,datatype = "abundance",q = 0,
    base = "coverage",level = 0.90,nboot = 0)
  
  as.numeric(out$qD[1])
}

# Calculate standardized richness
system.time({
  
  rarefaction_90 <- bind_rows(
    lapply(seq_along(community_list), function(i) {
      
      if (i %% 100 == 0)
        cat("Completed", i, "of",
            length(community_list), "\n")
      
      d <- community_list[[i]]
      
      data.frame(
        ID_HDC_G0_2 = d$ID_HDC_G0_2[1],
        Location = d$Location[1],
        Status = d$Status[1],
        richness_90 = rare90(d$records)
      )
    })
  )
  
})

saveRDS(rarefaction_90, "data/rarefaction_90.rds")

#rarefaction_90<-readRDS("data/rarefaction_90.rds")

#compare richness
richness_data <- sampling_summary %>%
  left_join(rarefaction_90,by = c("ID_HDC_G0_2", "Location", "Status"))

richness_data %>%group_by(Location, Status) %>%
  summarise(n_rarefied = sum(!is.na(richness_90)),mean_observed = mean(richness, na.rm = TRUE),
    mean_rarefied = mean(richness_90, na.rm = TRUE), .groups = "drop" )

#Let's plot both obs and 90% estimated
richness_matched <- richness_data %>%
  filter(!is.na(richness_90)) %>%
  mutate( Location = factor(Location,levels = c("Urban", "0–5 km", "5–30 km")))

richness_summary <- richness_matched %>%
  pivot_longer(cols = c(richness, richness_90),names_to = "Metric",values_to = "Richness") %>%
  group_by(Location, Status, Metric) %>%
  summarise(n = n(),mean = mean(Richness),SE = sd(Richness) / sqrt(n),
    CI95.L = mean - qt(0.975, n - 1) * SE,CI95.U = mean + qt(0.975, n - 1) * SE,
    .groups = "drop") %>%
  mutate(Metric = recode(Metric,richness = "Observed",richness_90 = "90% coverage"))

obs_90 <- ggplot(richness_summary,
  aes(x = Location, y = mean, colour = Status, group = Status)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI95.L, ymax = CI95.U),width = 0.12) +
  facet_wrap(~ Metric, scales = "free_y") +
  scale_colour_manual(values = c(Introduced = "darkorange3",Native = "darkolivegreen4")) +
  labs(x = NULL,y = "Mean species richness (95% CI)",colour = NULL) +
  theme_classic(base_size = 14)

obs_90

ggsave("out/richness_observed_vs_rarefied.pdf", obs_90,width = 10,height = 5)


richness_models <- richness_data %>%mutate(City = factor(ID_HDC_G0_2),
    Location = factor(Location,levels = c("Urban", "0–5 km", "5–30 km")),
    Status = factor(Status,levels = c("Native", "Introduced")))

# Observed richness: all available cities
m_observed <- lmer(log1p(richness) ~ Location * Status + (1 | City),data = richness_models)

# 90% coverage richness: eligible communities only
m_rare90 <- lmer(log1p(richness_90) ~ Location * Status + (1 | City),
  data = richness_models,na.action = na.omit)

summary(m_observed)
anova(m_observed)
summary(m_rare90)
anova(m_rare90)

get_contrasts <- function(model) {em <- emmeans(model, ~ Location | Status)
  contrast(em,method = list("Urban vs 0–5 km" = c(1, -1, 0),
      "Urban vs 5–30 km" = c(1, 0, -1))) %>% as.data.frame()}

tests_observed <- get_contrasts(m_observed)
tests_rare90 <- get_contrasts(m_rare90)

tests_observed
tests_rare90

saveRDS(list(observed = m_observed, coverage90 = m_rare90),"data/richness_models.rds")

#plot model effects

# Estimated marginal means
emm_rare90 <- emmeans(m_rare90, ~ Location * Status)

emm_df <- as.data.frame(emm_rare90)

# Plot
p_rare90 <- ggplot(emm_df,aes(x = Location,y = emmean,colour = Status,group = Status)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL,ymax = upper.CL ),width = 0.12,linewidth = 0.7) +
  scale_x_discrete(limits = c("Urban", "0–5 km", "5–30 km")) +
  scale_colour_manual(values = c(Native = "darkolivegreen4",Introduced = "darkorange3" ) ) +
  labs(x = "Location",y = "Estimated richness at 90% coverage",colour = NULL) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top" )

p_rare90


