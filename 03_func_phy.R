source("01_data_import.R")
source("02_sample_coverage_richness.R")

# Phylogenetic distance matrix
#phy_dist <- cophenetic.phylo(phy)

# Functional distance matrix
#trait_dist <- as.matrix(dist(traits, method = "euclidean"))

##these approaches require more memory than a laptop, need to do by city instead
#created functions is separate script

source("03b_func_phy_functions.R")

#run MPD calc from functions across cities (takes about 1-2 seconds per city)
#first phylo
set.seed(123)

# Split species_counts by city
city_list <- split(species_counts,species_counts$ID_HDC_G0_2)

# Run one city at a time
set.seed(123)

mpd_results <- lapply(
  seq_along(city_list),
  function(i) { cat("MPD: city", i, "of", length(city_list), "\n")
     result <- calculate_city_MPD(city_dat = city_list[[i]],phy = phy,nrep = 100)
    
    if (!is.null(result)) result$ID_HDC_G0_2 <- names(city_list)[i]
    result
  }) %>%
  dplyr::bind_rows()

# Combine results
mpd_results <- bind_rows(mpd_results)

saveRDS(mpd_results,"data/mpd_results.rds")
#mpd_results <- readRDS("data/mpd_results.rds")

m_MPD_obs <- lmer((MPD_obs) ~ Location * Status + (1 | ID_HDC_G0_2),
                  data = mpd_results)
summary(m_MPD_obs)
anova(m_MPD_obs)

m_MPD_90 <- lmer((MPD_90) ~ Location * Status + (1 | ID_HDC_G0_2),
                  data = mpd_results)
summary(m_MPD_90)
anova(m_MPD_90)


#next functional
#set.seed(123)

#city_list <- split(species_counts,species_counts$ID_HDC_G0_2)

#mfd_results <- lapply(seq_along(city_list),function(i) {
#    cat("MFD: city", i, "of", length(city_list), "\n")
    
#    result <- calculate_city_MFD(city_dat = city_list[[i]],
#      traits = traits_scaled,nrep = 100)
    
#    if (!is.null(result)) result$ID_HDC_G0_2 <- names(city_list)[i]
    
#    result
#  })

#mfd_results <- bind_rows(mfd_results)

#saveRDS(mfd_results, "data/mfd_results.rds")

####important, trait coverage is too poor to make a fair comparison. For this analysis, 
#we will forget about MFD

#combine mpd and richness for 90% coverage
plot_data <- rarefaction_90 %>%
  mutate(ID_HDC_G0_2 = as.character(ID_HDC_G0_2)) %>%
  select(ID_HDC_G0_2, Location, Status, richness_90) %>%
  inner_join(mpd_results %>%
      mutate(ID_HDC_G0_2 = as.character(ID_HDC_G0_2)) %>%
      select(ID_HDC_G0_2, Location, Status, MPD_90),
    by = c("ID_HDC_G0_2", "Location", "Status") )

#paired differences from urban - buffer
matched_cities <- coverage_matched %>%filter(all_90) %>%transmute(
    ID_HDC_G0_2 = as.character(ID_HDC_G0_2),Status )

paired_data <- plot_data %>%semi_join(matched_cities,by = c("ID_HDC_G0_2", "Status")) %>%
  pivot_wider(id_cols = c(ID_HDC_G0_2, Status),names_from = Location,
    values_from = c(richness_90, MPD_90)) %>%
  transmute(ID_HDC_G0_2,Status,
    Richness_5 = richness_90_Urban - `richness_90_0–5 km`,
    Richness_30 = richness_90_Urban - `richness_90_5–30 km`,
    MPD_5 = MPD_90_Urban - `MPD_90_0–5 km`,
    MPD_30 = MPD_90_Urban - `MPD_90_5–30 km` )
#now 95CI
plot_summary <- paired_data %>%
  pivot_longer(cols = c(Richness_5, Richness_30, MPD_5, MPD_30),
    names_to = "Comparison",values_to = "Difference" ) %>%
  mutate(Metric = if_else(
      grepl("^Richness", Comparison),
      "A. Species richness (90% coverage)",
      "B. Phylogenetic diversity (MPD)" ),
    Buffer = if_else(grepl("_5$", Comparison), "0–5 km", "5–30 km")) %>%
  group_by(Metric, Status, Buffer) %>%
  summarise(n = n(),Mean = mean(Difference),SE = sd(Difference) / sqrt(n),
    Lower = Mean - qt(0.975, n - 1) * SE,Upper = Mean + qt(0.975, n - 1) * SE,
    .groups = "drop")

plot_summary

##plot two panels
plot_summary$Group <- paste(plot_summary$Status,plot_summary$Buffer)
plot_summary$Group <- factor(plot_summary$Group,levels = c(
    "Native 0–5 km","Native 5–30 km","Introduced 0–5 km","Introduced 5–30 km"))

delta_R <- ggplot(plot_summary[plot_summary$Metric=="A. Species richness (90% coverage)",], 
                  aes(x = Mean, y = Group, colour = Status)) +
  geom_vline(xintercept = 0,linetype = "dashed",colour = "grey50") +
  geom_errorbar(aes(xmin = Lower, xmax = Upper),orientation = "y",
    width = 0.15,linewidth = 0.8) +
  geom_point(size = 3) +
  scale_colour_manual(
    values = c( Native = "darkolivegreen4",Introduced = "darkorange3")) +
  labs( x = "Urban − buffer diversity difference",y = NULL,colour = NULL,
        title = "C. Species richness (90% coverage)") +
  theme_classic(base_size = 13) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "none",panel.spacing = grid::unit(1.5, "lines"))

delta_R

delta_P <- ggplot(plot_summary[plot_summary$Metric=="B. Phylogenetic diversity (MPD)",], 
                  aes(x = Mean, y = Group, colour = Status)) +
  geom_vline(xintercept = 0,linetype = "dashed",colour = "grey50") +
  geom_errorbar(aes(xmin = Lower, xmax = Upper),orientation = "y",
                width = 0.15,linewidth = 0.8) +
  geom_point(size = 3) +
  scale_colour_manual(
    values = c( Native = "darkolivegreen4",Introduced = "darkorange3")) +
  labs( x = "Urban − buffer diversity difference",y = NULL,colour = NULL,
        title = "D. Phylogenetic diversity (MPD)") +
  theme_classic(base_size = 13) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "none",panel.spacing = grid::unit(1.5, "lines"))

delta_P

p_A<-p_A+labs(title = "A. Number of records")
p_B<-p_B+labs(title = "B. Species richness")
# Combine all four panels
figure1 <- (p_A | p_B) /
  (delta_R | delta_P)

figure1

ggsave("out/Fig1_obs_rich_MPD.pdf",figure1,width = 10,height = 8,device = cairo_pdf)


##MPD-env analyses ----
#similar to that for richness
#add MPD 90% coverage to richness_diff from 02

# Prepare urban MPD
mpd_urban <- mpd_results %>%filter(Location == "Urban") %>%
  transmute(ID_HDC_G0_2 = as.character(ID_HDC_G0_2),Status,
    MPD_obs_Urban = MPD_obs,MPD_90_Urban = MPD_90)

# Prepare buffer MPD
mpd_buffer <- mpd_results %>%
  filter(Location != "Urban") %>%
  transmute(ID_HDC_G0_2 = as.character(ID_HDC_G0_2),Status,Buffer = Location,
    MPD_obs_buffer = MPD_obs,MPD_90_buffer = MPD_90)

# Attach MPD to richness_diff
richness_diff <- richness_diff %>%
  mutate(ID_HDC_G0_2 = as.character(ID_HDC_G0_2)) %>%
  left_join(mpd_urban,by = c("ID_HDC_G0_2", "Status")) %>%
  left_join(mpd_buffer, by = c("ID_HDC_G0_2", "Status", "Buffer")) %>%
  mutate(MPD_Urban = if_else(Metric == "richness_90",MPD_90_Urban,MPD_obs_Urban),
    buffer_MPD = if_else(Metric == "richness_90",MPD_90_buffer,MPD_obs_buffer),
    MPD_diff = MPD_Urban - buffer_MPD) %>%
  select(-MPD_obs_Urban,-MPD_90_Urban,-MPD_obs_buffer,-MPD_90_buffer)


#simple
m_env_90_mpd_s <- lmer(MPD_diff ~Buffer * Status + (1 | ID_HDC_G0_2),
                     data = filter(richness_diff, Metric == "richness_90"))

summary(m_env_90_mpd_s)
anova(m_env_90_mpd_s)

#larger model
m_env_90_mpd <- lmer(MPD_diff ~Buffer * Status + z_veg_diff+ z_elevation_diff 
                 +z_buffer_area + z_built_diff + z_evenness + (1 | ID_HDC_G0_2),
                 data = filter(richness_diff, Metric == "richness_90"))

summary(m_env_90_mpd)
anova(m_env_90_mpd)
vif_90 <- check_collinearity(m_env_90_mpd)
vif_90 ###same variables inlcuded as richness models

dat_90_mpd <- model.frame(m_env_90_mpd)

global_90_mpd <-update(m_env_90_mpd,REML=FALSE,data = dat_90_mpd,na.action=na.fail)

#Dredge ----
dd_90_mpd <- dredge(global_90_mpd,
                 fixed = c("Buffer", "Status", "Buffer:Status"),rank = "AICc")
#top models
head(dd_90_mpd)
# Models within 2 AICc units of the best model
top_models90_mpd <- get.models(dd_90_mpd, subset = delta < 2)

# Model averaging
avg_model90_mpd <- model.avg(top_models90_mpd)
summary(avg_model90_mpd)
summary(avg_model90_mpd)$coefmat.full
confint(avg_model90_mpd, full = TRUE)

# Summed Akaike weights (variable importance)
sw(avg_model90_mpd)


em_mpd <- emmeans(m_env_90_mpd, ~ Buffer * Status)

results_mpd <- as.data.frame(confint(em_mpd))

results_mpd

##plot environmental effects for mpd
# Fixed-effect estimates
effects_90_mpd <- data.frame(term = names(fixef(m_env_90_mpd)),estimate = fixef(m_env_90_mpd))

# Wald 95% confidence intervals
ci_mpd <- confint(m_env_90_mpd,parm = "beta_",method = "Wald")

effects_90_mpd$conf.low <- ci[effects_90_mpd$term, 1]
effects_90_mpd$conf.high <- ci[effects_90_mpd$term, 2]



###four panel figure
# Consistent colours
status_cols <- c("Native" = "darkolivegreen4","Introduced" = "darkorange3")

location_levels <- c("Urban", "0–5 km", "5–30 km")

theme_fig <- theme_classic(base_size = 13) +
  theme(legend.position = "top",legend.title = element_blank(),
        axis.text.x = element_text(size = 11),plot.tag = element_text(face = "bold", size = 16) )

#panel A is coefficients for richness90 model
effects_90 <- data.frame(term = names(fixef(m_env_90)),
                         estimate = unname(fixef(m_env_90)))

ci <- confint( m_env_90,parm = "beta_",method = "Wald")

effects_90$conf.low <- ci[effects_90$term, 1]
effects_90$conf.high <- ci[effects_90$term, 2]

p2_A <- ggplot(env_effects, aes(x = estimate, y = Predictor)) +
  geom_vline(xintercept = 0,linetype = "dashed",colour = "grey50" ) +
  geom_segment( aes(x = conf.low,xend = conf.high,y = Predictor,yend = Predictor),
                linewidth = 0.7) +
  geom_point(size = 3, colour = "black") +
  labs( x = "Standardized predictor coefficient (95% CI)", y = NULL) +
  theme_fig

# Panel B: Buffer × Status interaction for richness90

em_90 <- as.data.frame(
  emmeans(m_env_90, ~ Buffer * Status)) %>%
  mutate(Buffer = factor( Buffer,levels = c("0–5 km", "5–30 km")),
         Status = factor(Status, levels = names(status_cols)))

p2_B <- ggplot(em_90,aes(x = Buffer, y = emmean, colour = Status, group = Status)) +
  geom_hline(yintercept = 0,linetype = "dashed",colour = "grey50") +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar( aes(ymin = lower.CL, ymax = upper.CL),width = 0.08, linewidth = 0.7) +
  scale_colour_manual(values = status_cols) +
  labs(x = "Buffer distance",y = "Urban - buffer log-richness contrast") +
  theme_fig


#panel C: coefficients for mpd model
effects_90_mpd <- data.frame(term = names(fixef(m_env_90_mpd)),
                         estimate = unname(fixef(m_env_90_mpd)))

ci_mpd <- confint( m_env_90_mpd,parm = "beta_",method = "Wald")

effects_90_mpd$conf.low <- ci_mpd[effects_90_mpd$term, 1]
effects_90_mpd$conf.high <- ci_mpd[effects_90_mpd$term, 2]

env_effects_mpd <- effects_90_mpd %>% filter(term %in% 
      c("z_built_diff","z_elevation_diff","z_veg_diff","z_evenness",
        "z_buffer_area","Buffer5–30 km","StatusNative",
        "Buffer5–30 km:StatusNative")) %>%
  mutate(Predictor = recode(term, z_built_diff = "Built area (diff)", 
      z_veg_diff = "Vegetation cover (diff)",z_buffer_area = "Buffer area",
      z_elevation_diff = "Elevation difference",z_evenness = "Pielou's evenness",
      `Buffer5–30 km` = "Buffer_5v30", StatusNative = "StatusIvN",
      `Buffer5–30 km:StatusNative` = "BufferXstatus"),
  Predictor = factor(Predictor,levels = rev(c("Built area (diff)", "Vegetation cover (diff)",
   "Pielou's evenness", "Elevation difference", "Buffer area","Buffer_5v30",
   "StatusIvN","BufferXstatus"))))


p2_C <- ggplot(env_effects_mpd, aes(x = estimate, y = Predictor)) +
  geom_vline(xintercept = 0,linetype = "dashed",colour = "grey50" ) +
  geom_segment( aes(x = conf.low,xend = conf.high,y = Predictor,yend = Predictor),
                linewidth = 0.7) +
  geom_point(size = 3, colour = "black") +
  labs( x = "Standardized predictor coefficient (95% CI)", y = NULL) +
  theme_fig

# Panel D: Buffer × Status interaction for mpd

em_90_mpd <- as.data.frame(
  emmeans(m_env_90_mpd, ~ Buffer * Status)) %>%
  mutate(Buffer = factor( Buffer,levels = c("0–5 km", "5–30 km")),
         Status = factor(Status, levels = names(status_cols)))

p2_D <- ggplot(em_90_mpd,aes(x = Buffer, y = emmean, colour = Status, group = Status)) +
  geom_hline(yintercept = 0,linetype = "dashed",colour = "grey50") +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar( aes(ymin = lower.CL, ymax = upper.CL),width = 0.08, linewidth = 0.7) +
  scale_colour_manual(values = status_cols) +
  labs(x = "Buffer distance",y = "Urban - buffer MPD contrast") +
  theme_fig

# Combine panels A–D

p_combined4_mpd <- ((p2_A | p2_B) /(p2_C | p2_D)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

p_combined4_mpd

ggsave("out/Fig2_rich_mpd_models.pdf",plot = p_combined4_mpd,width = 12,height = 9,
       device = cairo_pdf)

###plot all independent drivers of richness diff_90
richness_diff90<-richness_diff[richness_diff$Metric=="richness_90",]

p_RD_forD<-ggplot(richness_diff90,aes(x = forest_diff, y = richness_diff))+
  geom_point(alpha = 0.25,size = 1.5,colour = "grey40" ) +
  geom_smooth(method = "lm",se = TRUE,colour = "steelblue4",fill = "steelblue2",
              alpha = 0.25,linewidth = 1.2) +
  facet_wrap(~ Status, ncol = 2, scales = "free") +
  labs(x = "Urban-buffer forest diff",y = "Urban-buffer richness difference 90%") +
  theme_fig

p_RD_agD<-ggplot(richness_diff90,aes(x = ag_diff, y = richness_diff))+
  geom_point(alpha = 0.25,size = 1.5,colour = "grey40" ) +
  geom_smooth(method = "lm",se = TRUE,colour = "steelblue4",fill = "steelblue2",
              alpha = 0.25,linewidth = 1.2) +
  facet_wrap(~ Status, ncol = 2, scales = "free") +
  labs(x = "Urban-buffer agriculture diff",y = "Urban-buffer richness difference 90%") +
  theme_fig

p_RD_grD<-ggplot(richness_diff90,aes(x = grass_diff, y = richness_diff))+
  geom_point(alpha = 0.25,size = 1.5,colour = "grey40" ) +
  geom_smooth(method = "lm",se = TRUE,colour = "steelblue4",fill = "steelblue2",
              alpha = 0.25,linewidth = 1.2) +
  facet_wrap(~ Status, ncol = 2, scales = "free") +
  labs(x = "Urban-buffer grassland diff",y = "Urban-buffer richness difference 90%") +
  theme_fig

p_RD_vgD<-ggplot(richness_diff90,aes(x = veg_diff, y = richness_diff))+
  geom_point(alpha = 0.25,size = 1.5,colour = "grey40" ) +
  geom_smooth(method = "lm",se = TRUE,colour = "steelblue4",fill = "steelblue2",
              alpha = 0.25,linewidth = 1.2) +
  facet_wrap(~ Status, ncol = 2, scales = "free") +
  labs(x = "Urban-buffer vegetation diff",y = "Urban-buffer richness difference 90%") +
  theme_fig

p_RD_blD<-ggplot(richness_diff90,aes(x = built_diff, y = richness_diff))+
  geom_point(alpha = 0.25,size = 1.5,colour = "grey40" ) +
  geom_smooth(method = "lm",se = TRUE,colour = "steelblue4",fill = "steelblue2",
              alpha = 0.25,linewidth = 1.2) +
  facet_wrap(~ Status, ncol = 2, scales = "free") +
  labs(x = "Urban-buffer built diff",y = "Urban-buffer richness difference 90%") +
  theme_fig

p_RD_elD<-ggplot(richness_diff90,aes(x = elevation_diff, y = richness_diff))+
  geom_point(alpha = 0.25,size = 1.5,colour = "grey40" ) +
  geom_smooth(method = "lm",se = TRUE,colour = "steelblue4",fill = "steelblue2",
              alpha = 0.25,linewidth = 1.2) +
  facet_wrap(~ Status, ncol = 2, scales = "free") +
  labs(x = "Urban-buffer elevation diff",y = "Urban-buffer richness difference 90%") +
  theme_fig

p_RD_urA<-ggplot(richness_diff90,aes(x = urban_area, y = richness_diff))+
  geom_point(alpha = 0.25,size = 1.5,colour = "grey40" ) +
  geom_smooth(method = "lm",se = TRUE,colour = "steelblue4",fill = "steelblue2",
              alpha = 0.25,linewidth = 1.2) +
  facet_wrap(~ Status, ncol = 2, scales = "free") +
  labs(x = "Urban area",y = "Urban-buffer richness difference 90%") +
  theme_fig

p_RD_bfA<-ggplot(richness_diff90,aes(x = buffer_area, y = richness_diff))+
  geom_point(alpha = 0.25,size = 1.5,colour = "grey40" ) +
  geom_smooth(method = "lm",se = TRUE,colour = "steelblue4",fill = "steelblue2",
              alpha = 0.25,linewidth = 1.2) +
  facet_wrap(~ Status, ncol = 2, scales = "free") +
  labs(x = "Buffer area",y = "Urban-buffer richness difference 90%") +
  theme_fig

p_env_combo <- ((p_RD_forD) /(p_RD_agD)/(p_RD_grD)/(p_RD_vgD) |
                  (p_RD_blD) /(p_RD_elD)/(p_RD_urA)/(p_RD_bfA)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

ggsave("out/FigX_rich_diff_veg_diff.pdf",plot = p_env_combo,width = 10,height = 15,
       device = cairo_pdf)

#p_env_combo <- ((p_RD_blD) /(p_RD_elD)/(p_RD_urA)/(p_RD_bfA)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

#ggsave("out/FigX_rich_diff_env_diff.pdf",plot = p_env_combo,width = 5,height = 12,
       device = cairo_pdf)

##run models

m_RD_forD<-lm(richness_diff~forest_diff*Status,data=richness_diff90)
summary(m_RD_forD)
anova(m_RD_forD)

m_RD_agD<-lm(richness_diff~ag_diff*Status,data=richness_diff90)
summary(m_RD_agD)
anova(m_RD_agD)

m_RD_grD<-lm(richness_diff~grass_diff*Status,data=richness_diff90)
summary(m_RD_grD)
anova(m_RD_grD)

m_RD_vgD<-lm(richness_diff~veg_diff*Status,data=richness_diff90)
summary(m_RD_vgD)
anova(m_RD_vgD)

m_RD_blD<-lm(richness_diff~built_diff*Status,data=richness_diff90)
summary(m_RD_blD)
anova(m_RD_blD)

m_RD_elD<-lm(richness_diff~elevation_diff*Status,data=richness_diff90)
summary(m_RD_elD)
anova(m_RD_elD)

m_RD_urD<-lm(richness_diff~urban_area*Status,data=richness_diff90)
summary(m_RD_urD)
anova(m_RD_urD)

m_RD_bfD<-lm(richness_diff~buffer_area*Status,data=richness_diff90)
summary(m_RD_bfD)
anova(m_RD_bfD)

###now for MPD

p_MD_forD<-ggplot(richness_diff90,aes(x = forest_diff, y = MPD_diff))+
  geom_point(alpha = 0.25,size = 1.5,colour = "grey40" ) +
  geom_smooth(method = "lm",se = TRUE,colour = "steelblue4",fill = "steelblue2",
              alpha = 0.25,linewidth = 1.2) +
  facet_wrap(~ Status, ncol = 2, scales = "free") +
  labs(x = "Urban-buffer forest diff",y = "Urban-buffer MPD difference 90%") +
  theme_fig

p_MD_agD<-ggplot(richness_diff90,aes(x = ag_diff, y = MPD_diff))+
  geom_point(alpha = 0.25,size = 1.5,colour = "grey40" ) +
  geom_smooth(method = "lm",se = TRUE,colour = "steelblue4",fill = "steelblue2",
              alpha = 0.25,linewidth = 1.2) +
  facet_wrap(~ Status, ncol = 2, scales = "free") +
  labs(x = "Urban-buffer agriculture diff",y = "Urban-buffer MPD difference 90%") +
  theme_fig

p_MD_grD<-ggplot(richness_diff90,aes(x = grass_diff, y = MPD_diff))+
  geom_point(alpha = 0.25,size = 1.5,colour = "grey40" ) +
  geom_smooth(method = "lm",se = TRUE,colour = "steelblue4",fill = "steelblue2",
              alpha = 0.25,linewidth = 1.2) +
  facet_wrap(~ Status, ncol = 2, scales = "free") +
  labs(x = "Urban-buffer grassland diff",y = "Urban-buffer MPD difference 90%") +
  theme_fig

p_MD_vgD<-ggplot(richness_diff90,aes(x = veg_diff, y = MPD_diff))+
  geom_point(alpha = 0.25,size = 1.5,colour = "grey40" ) +
  geom_smooth(method = "lm",se = TRUE,colour = "steelblue4",fill = "steelblue2",
              alpha = 0.25,linewidth = 1.2) +
  facet_wrap(~ Status, ncol = 2, scales = "free") +
  labs(x = "Urban-buffer vegetation diff",y = "Urban-buffer MPD difference 90%") +
  theme_fig

p_MD_blD<-ggplot(richness_diff90,aes(x = built_diff, y = MPD_diff))+
  geom_point(alpha = 0.25,size = 1.5,colour = "grey40" ) +
  geom_smooth(method = "lm",se = TRUE,colour = "steelblue4",fill = "steelblue2",
              alpha = 0.25,linewidth = 1.2) +
  facet_wrap(~ Status, ncol = 2, scales = "free") +
  labs(x = "Urban-buffer built diff",y = "Urban-buffer MPD difference 90%") +
  theme_fig

p_MD_elD<-ggplot(richness_diff90,aes(x = elevation_diff, y = MPD_diff))+
  geom_point(alpha = 0.25,size = 1.5,colour = "grey40" ) +
  geom_smooth(method = "lm",se = TRUE,colour = "steelblue4",fill = "steelblue2",
              alpha = 0.25,linewidth = 1.2) +
  facet_wrap(~ Status, ncol = 2, scales = "free") +
  labs(x = "Urban-buffer elevation diff",y = "Urban-buffer MPD difference 90%") +
  theme_fig

p_MD_urA<-ggplot(richness_diff90,aes(x = urban_area, y = MPD_diff))+
  geom_point(alpha = 0.25,size = 1.5,colour = "grey40" ) +
  geom_smooth(method = "lm",se = TRUE,colour = "steelblue4",fill = "steelblue2",
              alpha = 0.25,linewidth = 1.2) +
  facet_wrap(~ Status, ncol = 2, scales = "free") +
  labs(x = "Urban area",y = "Urban-buffer MPD difference 90%") +
  theme_fig

p_MD_bfA<-ggplot(richness_diff90,aes(x = buffer_area, y = MPD_diff))+
  geom_point(alpha = 0.25,size = 1.5,colour = "grey40" ) +
  geom_smooth(method = "lm",se = TRUE,colour = "steelblue4",fill = "steelblue2",
              alpha = 0.25,linewidth = 1.2) +
  facet_wrap(~ Status, ncol = 2, scales = "free") +
  labs(x = "Buffer area",y = "Urban-buffer MPD difference 90%") +
  theme_fig

p_env_comboM <- ((p_MD_forD) /(p_MD_agD)/(p_MD_grD)/(p_MD_vgD) |
                  (p_MD_blD) /(p_MD_elD)/(p_MD_urA)/(p_MD_bfA)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

ggsave("out/FigX_MPD_diff_veg_diff.pdf",plot = p_env_comboM,width = 10,height = 15,
       device = cairo_pdf)

#p_env_combo <- ((p_RD_blD) /(p_RD_elD)/(p_RD_urA)/(p_RD_bfA)) +
plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

#ggsave("out/FigX_rich_diff_env_diff.pdf",plot = p_env_combo,width = 5,height = 12,
device = cairo_pdf)

##run models

m_MD_forD<-lm(MPD_diff~forest_diff*Status,data=richness_diff90)
summary(m_MD_forD)
anova(m_MD_forD)

m_MD_agD<-lm(MPD_diff~ag_diff*Status,data=richness_diff90)
summary(m_MD_agD)
anova(m_MD_agD)

m_MD_grD<-lm(MPD_diff~grass_diff*Status,data=richness_diff90)
summary(m_MD_grD)
anova(m_MD_grD)

m_MD_vgD<-lm(MPD_diff~veg_diff*Status,data=richness_diff90)
summary(m_MD_vgD)
anova(m_MD_vgD)

m_MD_blD<-lm(MPD_diff~built_diff*Status,data=richness_diff90)
summary(m_MD_blD)
anova(m_MD_blD)

m_MD_elD<-lm(MPD_diff~elevation_diff*Status,data=richness_diff90)
summary(m_MD_elD)
anova(m_MD_elD)

m_MD_urD<-lm(MPD_diff~urban_area*Status,data=richness_diff90)
summary(m_MD_urD)
anova(m_MD_urD)

m_MD_bfD<-lm(MPD_diff~buffer_area*Status,data=richness_diff90)
summary(m_MD_bfD)
anova(m_MD_bfD)
