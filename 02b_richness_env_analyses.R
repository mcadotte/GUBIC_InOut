##analyzing environmental predictors
source("01_data_import.R")
source("02_sample_coverage_richness.R")


#set up data for models ----
##environmental predictors
env_5 <- B_env %>%transmute(ID_HDC_G0_2 = city_id,Buffer = "0–5 km",
                            forest = buffer5km_forestland,forest_diff = forestland_diff5km,
                            agriculture = buffer5km_cropland, ag_diff = cropland_diff5km,
                            grassland = buffer5km_grassland, grass_diff = grassland_diff5km,
                            vegetation = buffer5km_veg, veg_diff = veg_diff5km,
                            builtup = buffer5km_builtup, built_diff = builtup_diff5km,
                            elevation_diff = ele_diff5km,urban_area = area_urban_land_km2,
                            buffer_area = area_buffer_0_5km_land_km2, evenness = Evenness_5km)

env_30 <- B_env %>%transmute(ID_HDC_G0_2 = city_id,Buffer = "5–30 km",
                             forest = buffer30km_forestland,forest_diff = forestland_diff30km,
                             agriculture = buffer30km_cropland, ag_diff = cropland_diff30km,
                             grassland = buffer30km_grassland, grass_diff = grassland_diff30km,
                             vegetation = buffer30km_veg, veg_diff = veg_diff30km,
                             builtup = buffer30km_builtup, built_diff = builtup_diff30km,
                             elevation_diff = ele_diff30km,urban_area = area_urban_land_km2,
                             buffer_area = area_buffer_5_30km_land_km2, evenness = Evenness_30km)

env_buffers <- bind_rows(env_5, env_30)

##plot correlations
env_num <- env_buffers[,3:ncol(env_buffers)]

# Calculate Pearson correlations
cor_mat <- cor(env_num,use = "pairwise.complete.obs",method = "pearson")

# Plot correlation matrix
corrplot(cor_mat,method = "color",type = "upper",order = "hclust",
  addCoef.col = "black",number.cex = 0.7,
  tl.col = "black",tl.cex = 0.8,tl.srt = 45,diag = FALSE,
  col = colorRampPalette(c("#B2182B", "white", "#2166AC"))(200))

### plot diff histograms

p_for<-ggplot(env_buffers,aes(forest_diff))+
  geom_histogram(binwidth = 5,fill="steelblue",colour = "white")+
  theme_classic(base_size = 14) +
  annotate("segment",x = -30, y = 150, xend = -70, yend = 150, 
    color = "grey25",linewidth = 2,
    arrow = arrow(length = unit(0.3, "cm"), type = "closed")) +
  annotate("text",x = -30, y = 120, hjust=1, label = "Greater in buffer",
           color = "grey25", size = 4)+
  labs(x="Forest difference")+
  theme(strip.text = element_text(size = 13, face = "bold"),
        panel.spacing = unit(1.2, "lines"))

p_ag<-ggplot(env_buffers,aes(ag_diff))+
  geom_histogram(binwidth = 5,fill="steelblue",colour = "white")+
  theme_classic(base_size = 14) +
  labs(x="Agriculture difference")+
  theme(strip.text = element_text(size = 13, face = "bold"),
        panel.spacing = unit(1.2, "lines"))

p_grass<-ggplot(env_buffers,aes(grass_diff))+
  geom_histogram(binwidth = 3,fill="steelblue",colour = "white")+
  theme_classic(base_size = 14) +
  labs(x="Grassland difference")+
  theme(strip.text = element_text(size = 13, face = "bold"),
        panel.spacing = unit(1.2, "lines"))

p_veg<-ggplot(env_buffers,aes(veg_diff))+
  geom_histogram(binwidth = 5,fill="steelblue",colour = "white")+
  theme_classic(base_size = 14) +
  labs(x="Vegetation difference")+
  theme(strip.text = element_text(size = 13, face = "bold"),
        panel.spacing = unit(1.2, "lines"))

p_el<-ggplot(env_buffers,aes(elevation_diff))+
  geom_histogram(binwidth = 60,fill="steelblue",colour = "white")+
  theme_classic(base_size = 14) +
  labs(x="Elevation difference")+
  theme(strip.text = element_text(size = 13, face = "bold"),
        panel.spacing = unit(1.2, "lines"))

p_bil<-ggplot(env_buffers,aes(built_diff))+
  geom_histogram(binwidth = 5,fill="steelblue",colour = "white")+
  theme_classic(base_size = 14) +
  labs(x="Built-up difference")+
  theme(strip.text = element_text(size = 13, face = "bold"),
        panel.spacing = unit(1.2, "lines"))


p_env.all <- ((p_for | p_grass | p_ag) /(p_veg | p_el | p_bil)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

ggsave("out/FSX_env_diff_hist.pdf",plot = p_env.all,width = 10,
       height = 6,device = cairo_pdf)


##plot relationships between forest grassland agriculture builtup and vegetation

# Prepare long table with variables needed
plot_dat <- bind_rows(env_buffers %>%
    transmute(x = forest, y = agriculture,Panel = "A. Agriculture vs Forest"),
  
  env_buffers %>%
    transmute(x = forest, y = grassland,Panel = "B. Grassland vs Forest"),
  
  env_buffers %>%
    transmute(x = vegetation, y = forest,Panel = "C. Forest vs Vegetation"),
  
  env_buffers %>%
    transmute(x = vegetation, y = builtup,Panel = "D. Built-up vs Vegetation")) %>%
  mutate(Panel = factor(Panel,levels = c("A. Agriculture vs Forest",
     "B. Grassland vs Forest","C. Forest vs Vegetation","D. Built-up vs Vegetation")))


# Four-panel figure
p_env <- ggplot(plot_dat, aes(x = x, y = y)) +
  geom_point(alpha = 0.25,size = 1.5,colour = "grey40" ) +
  geom_smooth(method = "lm",se = TRUE,colour = "steelblue4",fill = "steelblue2",
    alpha = 0.25,linewidth = 1.2) +
  facet_wrap(~ Panel, ncol = 2, scales = "free") +
  labs(x = NULL,y = NULL) +
  theme_classic(base_size = 14) +
  theme(strip.text = element_text(size = 13, face = "bold"),
    panel.spacing = unit(1.2, "lines"))

p_env
ggsave("out/FSX_env_four_panel.pdf",plot = p_env,width = 6,height = 6,device = cairo_pdf)

##examine richness differences

richness_diff <- richness_data %>%
  select(ID_HDC_G0_2, Status, Location,richness, richness_90) %>%
  pivot_wider( names_from = Location,values_from = c(richness, richness_90)) %>%
  pivot_longer(cols = c(`richness_0–5 km`,`richness_5–30 km`,`richness_90_0–5 km`,
                        `richness_90_5–30 km`),names_to = c("Metric", "Buffer"),
               names_pattern = "(richness(?:_90)?)_(.*)",values_to = "buffer_richness" ) %>%
  mutate(urban_richness = if_else(Metric == "richness",richness_Urban,richness_90_Urban),
         richness_diff = log1p(urban_richness) -log1p(buffer_richness)) %>%
  left_join(env_buffers,by = c("ID_HDC_G0_2", "Buffer")  )

# Standardize predictors for comparable coefficients
richness_diff <- richness_diff %>%mutate( across( c(forest, forest_diff,agriculture, 
                                                    ag_diff, grassland, grass_diff, vegetation, veg_diff, builtup, built_diff,
                                                    elevation_diff, buffer_area,evenness), ~ as.numeric(scale(.x)),.names = "z_{.col}"))

#model selection ----
#run models for dredge
#high vif so retaining 'diff' variables, update: forst and grass high cor with veg. keep veg
m_env_observed <- lmer(richness_diff ~Buffer * Status + z_veg_diff+z_elevation_diff 
                       +z_buffer_area +z_built_diff + z_evenness + (1 | ID_HDC_G0_2),
                       data = filter(richness_diff, Metric == "richness"))

summary(m_env_observed)
anova(m_env_observed)
vif_obs <- check_collinearity(m_env_observed)
vif_obs

# Extract observations actually used in the original model
dat_obs <- model.frame(m_env_observed)

global_obs <-update(m_env_observed,REML=FALSE,data = dat_obs,na.action=na.fail)

#models for the 90% coverage
m_env_90 <- lmer(richness_diff ~Buffer * Status + z_veg_diff+ z_elevation_diff 
                 +z_buffer_area + z_built_diff + z_evenness + (1 | ID_HDC_G0_2),
                 data = filter(richness_diff, Metric == "richness_90"))

summary(m_env_90)
anova(m_env_90)
vif_90 <- check_collinearity(m_env_90)
vif_90 ###need to exclude variables. I will retain 'diff' rather than percent?

dat_90 <- model.frame(m_env_90)

global_90 <-update(m_env_90,REML=FALSE,data = dat_90,na.action=na.fail)

#Dredge ----
#first observed richness (supp Mat)
dd_obs <- dredge(global_obs,
                 fixed = c("Buffer", "Status", "Buffer:Status"),rank = "AICc")
#top models
head(dd_obs)
# Models within 2 AICc units of the best model
top_models <- get.models(dd_obs, subset = delta < 4)

# Model averaging
avg_model <- model.avg(top_models)
summary(avg_model)
summary(avg_model)$coefmat.full
confint(avg_model, full = TRUE)

# Summed Akaike weights (variable importance)
sw(avg_model)

##now 90% coverage
dd_90 <- dredge(global_90,
                fixed = c("Buffer", "Status", "Buffer:Status"),rank = "AICc")
#top models
head(dd_90)
# Models within 2 AICc units of the best model
top_models90 <- get.models(dd_90, subset = delta < 2)

get.models(dd_90, subset = 1)
best90<-lmer(richness_diff ~ z_built_diff + z_elevation_diff + Buffer + Status +  
               (1 | ID_HDC_G0_2) + Buffer:Status, data = dat_90)
summary(best90)

# Model averaging
avg_model90 <- model.avg(top_models90)
summary(avg_model90)
summary(avg_model90)$coefmat.full
confint(avg_model90, full = TRUE)

# Summed Akaike weights (variable importance)
sw(avg_model90)


em <- emmeans(m_env_90, ~ Buffer * Status)

results <- as.data.frame(confint(em))

results

##plot environmental effects for 90% coverage
# Fixed-effect estimates
effects_90 <- data.frame(term = names(fixef(m_env_90)),estimate = fixef(m_env_90))

# Wald 95% confidence intervals
ci <- confint(m_env_90,parm = "beta_",method = "Wald")

effects_90$conf.low <- ci[effects_90$term, 1]
effects_90$conf.high <- ci[effects_90$term, 2]

# Select environmental predictors
env_effects <- effects_90 %>% filter(term %in% c("z_built_diff","z_elevation_diff",
           "z_veg_diff","z_evenness","z_buffer_area", "Buffer5–30 km","StatusNative",
                   "Buffer5–30 km:StatusNative")) %>%
  mutate(Predictor = recode(term, z_built_diff = "Built area (diff)", 
                            z_veg_diff = "Vegetation cover (diff)",
                            z_buffer_area = "Buffer area",
                            z_elevation_diff = "Elevation difference",
                            z_evenness = "Pielou's evenness",
                            `Buffer5–30 km` = "Buffer_5v30",
                             StatusNative = "StatusIvN",
                            `Buffer5–30 km:StatusNative` ="BufferXstatus"),
  Predictor = factor(Predictor,levels = rev(c("Built area (diff)", "Vegetation cover (diff)",
           "Pielou's evenness", "Elevation difference", "Buffer area","Buffer_5v30",
           "StatusIvN","BufferXstatus"))))

###four panel figure
# Consistent colours
status_cols <- c("Native" = "darkolivegreen4","Introduced" = "darkorange3")

location_levels <- c("Urban", "0–5 km", "5–30 km")

theme_fig <- theme_classic(base_size = 13) +
  theme(legend.position = "top",legend.title = element_blank(),
        axis.text.x = element_text(size = 11),plot.tag = element_text(face = "bold", size = 16) )

#Panel A: GBIF occurrence records

records_summary <- sampling_summary %>%
  group_by(Location, Status) %>%summarise(
    n = sum(!is.na(total_records)),
    mean = mean(total_records, na.rm = TRUE),
    SE = sd(total_records, na.rm = TRUE) / sqrt(n),
    CI95.L = mean - qt(0.975, n - 1) * SE,
    CI95.U = mean + qt(0.975, n - 1) * SE,
    .groups = "drop") %>%
  mutate(Location = factor(Location, levels = location_levels),
         Status = factor(Status, levels = names(status_cols)))

p_A <- ggplot(records_summary,
              aes(x = Location, y = mean, colour = Status, group = Status)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI95.L, ymax = CI95.U),width = 0.10,linewidth = 0.7) +
  scale_colour_manual(values = status_cols) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Location",y = "Mean number of occurrence records" ) +
  theme_fig


# Panel B: Richness standardized to 90% sample coverage

richness_90_summary <- richness_data %>%
  filter(!is.na(richness_90)) %>%
  group_by(Location, Status) %>%
  summarise(n = n(), mean = mean(richness_90), SE = sd(richness_90) / sqrt(n),
            CI95.L = mean - qt(0.975, n - 1) * SE, CI95.U = mean + qt(0.975, n - 1) * SE,
            .groups = "drop" ) %>%
  mutate(Location = factor(Location, levels = location_levels),
         Status = factor(Status, levels = names(status_cols)))

p_B <- ggplot(richness_90_summary,
              aes(x = Location, y = mean, colour = Status, group = Status)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI95.L, ymax = CI95.U), width = 0.10,linewidth = 0.7) +
  scale_colour_manual(values = status_cols) +
  labs(x = "Location",y = "Mean species richness (90% coverage)") +
  theme_fig


#Panel C: Environmental predictor effect
effects_90 <- data.frame(term = names(fixef(m_env_90)),
                         estimate = unname(fixef(m_env_90)))

ci <- confint( m_env_90,parm = "beta_",method = "Wald")

effects_90$conf.low <- ci[effects_90$term, 1]
effects_90$conf.high <- ci[effects_90$term, 2]

p_C <- ggplot(env_effects, aes(x = estimate, y = Predictor)) +
  geom_vline(xintercept = 0,linetype = "dashed",colour = "grey50" ) +
  geom_segment( aes(x = conf.low,xend = conf.high,y = Predictor,yend = Predictor),
                linewidth = 0.7) +
  geom_point(size = 3, colour = "black") +
  labs( x = "Standardized predictor coefficient (95% CI)", y = NULL) +
  theme_fig



# Panel D: Buffer × Status interaction

em_90 <- as.data.frame(
  emmeans(m_env_90, ~ Buffer * Status)) %>%
  mutate(Buffer = factor( Buffer,levels = c("0–5 km", "5–30 km")),
         Status = factor(Status, levels = names(status_cols)))

p_D <- ggplot(em_90,aes(x = Buffer, y = emmean, colour = Status, group = Status)) +
  geom_hline(yintercept = 0,linetype = "dashed",colour = "grey50") +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar( aes(ymin = lower.CL, ymax = upper.CL),width = 0.08, linewidth = 0.7) +
  scale_colour_manual(values = status_cols) +
  labs(x = "Buffer distance",y = "Urban - buffer log-richness contrast") +
  theme_fig

# Combine panels A–D

p_combined4 <- ((p_A | p_B) /(p_C | p_D)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

p_combined4

ggsave("out/richness_four_panel.pdf",plot = p_combined4,width = 12,height = 9,
       device = cairo_pdf)

##plot environmental effects for observed (Supp Mat)
# Fixed-effect estimates
effects_obs <- data.frame(term = names(fixef(m_env_observed)),estimate = fixef(m_env_observed))

# Wald 95% confidence intervals
ci <- confint(m_env_observed,parm = "beta_",method = "Wald")

effects_obs$conf.low <- ci[effects_obs$term, 1]
effects_obs$conf.high <- ci[effects_obs$term, 2]

# Select environmental predictors
env_effects <- effects_obs %>% filter(term %in% c("z_built_diff","z_elevation_diff",
                                                  "z_veg_diff","z_evenness","z_buffer_area")) %>%
  mutate(Predictor = recode(term, z_built_diff = "Built area (diff)", 
                            z_veg_diff = "Vegetation cover (diff)",z_buffer_area = "Buffer area",
                            z_elevation_diff = "Elevation difference",z_evenness = "Pielou's evenness"),
         Predictor = factor(Predictor,levels = rev(c("Built area (diff)", "Vegetation cover (diff)",
                                                     "Pielou's evenness", "Elevation difference", "Buffer area"))))

###four panel figure for observed
# Consistent colours
status_cols <- c("Native" = "darkolivegreen4","Introduced" = "darkorange3")

location_levels <- c("Urban", "0–5 km", "5–30 km")

theme_fig <- theme_classic(base_size = 13) +
  theme(legend.position = "top",legend.title = element_blank(),
        axis.text.x = element_text(size = 11),plot.tag = element_text(face = "bold", size = 16) )

#Panel A: GBIF occurrence records

records_summary <- sampling_summary %>%
  group_by(Location, Status) %>%summarise(
    n = sum(!is.na(total_records)),
    mean = mean(total_records, na.rm = TRUE),
    SE = sd(total_records, na.rm = TRUE) / sqrt(n),
    CI95.L = mean - qt(0.975, n - 1) * SE,
    CI95.U = mean + qt(0.975, n - 1) * SE,
    .groups = "drop") %>%
  mutate(Location = factor(Location, levels = location_levels),
         Status = factor(Status, levels = names(status_cols)))

p_Ao <- ggplot(records_summary,
               aes(x = Location, y = mean, colour = Status, group = Status)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI95.L, ymax = CI95.U),width = 0.10,linewidth = 0.7) +
  scale_colour_manual(values = status_cols) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Location",y = "Mean number of occurrence records" ) +
  theme_fig


# Panel B: Richness observed

richness_obs_summary <- richness_data %>%
  filter(!is.na(richness)) %>%
  group_by(Location, Status) %>%
  summarise(n = n(), mean = mean(richness), SE = sd(richness) / sqrt(n),
            CI95.L = mean - qt(0.975, n - 1) * SE, CI95.U = mean + qt(0.975, n - 1) * SE,
            .groups = "drop" ) %>%
  mutate(Location = factor(Location, levels = location_levels),
         Status = factor(Status, levels = names(status_cols)))

p_Bo <- ggplot(richness_obs_summary,
               aes(x = Location, y = mean, colour = Status, group = Status)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI95.L, ymax = CI95.U), width = 0.10,linewidth = 0.7) +
  scale_colour_manual(values = status_cols) +
  labs(x = "Location",y = "Mean species richness") +
  theme_fig


#Panel C: Environmental predictor effect
effects_obs <- data.frame(term = names(fixef(m_env_observed)),
                          estimate = unname(fixef(m_env_observed)))

ci <- confint( m_env_observed,parm = "beta_",method = "Wald")

effects_obs$conf.low <- ci[effects_obs$term, 1]
effects_obs$conf.high <- ci[effects_obs$term, 2]

p_Co <- ggplot(env_effects, aes(x = estimate, y = Predictor)) +
  geom_vline(xintercept = 0,linetype = "dashed",colour = "grey50" ) +
  geom_segment( aes(x = conf.low,xend = conf.high,y = Predictor,yend = Predictor),
                linewidth = 0.7) +
  geom_point(size = 3, colour = "black") +
  labs( x = "Standardized predictor coefficient (95% CI)", y = NULL) +
  theme_fig



# Panel D: Buffer × Status interaction

em_obs <- as.data.frame(
  emmeans(m_env_observed, ~ Buffer * Status)) %>%
  mutate(Buffer = factor( Buffer,levels = c("0–5 km", "5–30 km")),
         Status = factor(Status, levels = names(status_cols)))

p_Do <- ggplot(em_obs,aes(x = Buffer, y = emmean, colour = Status, group = Status)) +
  geom_hline(yintercept = 0,linetype = "dashed",colour = "grey50") +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar( aes(ymin = lower.CL, ymax = upper.CL),width = 0.08, linewidth = 0.7) +
  scale_colour_manual(values = status_cols) +
  labs(x = "Buffer distance",y = "Urban - buffer log-richness contrast") +
  theme_fig

# Combine panels A–D

p_combined4.obs <- ((p_Ao | p_Bo) /(p_Co | p_Do)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

p_combined4.obs

#ggsave("out/OBS_richness_four_panel_SUPP_MAT.pdf",plot = p_combined4.obs,width = 12,height = 9,
#       device = cairo_pdf)

##will change order of figures in 03

#plot species ratios for Fig. S1
ratio_dat <- species_counts %>%
  filter(Location %in% c("Urban", "0–5 km", "5–30 km"),Status %in% c("Native", "Introduced"),
         records > 0,!is.na(accepted_taxon_name_binomial)) %>%
  group_by(ID_HDC_G0_2, Location, Status) %>%
  summarise(Richness = n_distinct(accepted_taxon_name_binomial),.groups = "drop") %>%
  pivot_wider(names_from = Status,values_from = Richness,values_fill = 0) %>%
  filter(Native > 0) %>%
  mutate(IN_ratio = Introduced / Native,Location = factor(Location,
                                                          levels = c("Urban", "0–5 km", "5–30 km")) )



p_FS1 <- ggplot(ratio_dat,aes(x = Location, y = IN_ratio, fill = Location)) +
  geom_boxplot(width = 0.6,alpha = 0.8,outlier.size = 1,outlier.alpha = 0.5) +
  scale_fill_manual(values = c("Urban" = "#37658C","0–5 km" = "#219FA5",
                               "5–30 km" = "#ACF6F0")) +
  labs(x = "Location",y = "Introduced : Native richness ratio") +

  theme_classic(base_size = 14) +
  theme(legend.position = "none")

#p_FS1

#ggsave("out/FS1_IN_ratio.pdf",p_FS1,width = 6,height = 4.5,device = cairo_pdf)
#test if sig different
m_IN<-lmer(IN_ratio~Location+(1|ID_HDC_G0_2),data=ratio_dat)
summary(m_IN)
anova(m_IN)

L_F<-"F[2*','*1102] == 310.06"
L_P<-"italic(P) < 0.001~~'***'"
p_FS1<-p_FS1 + annotate("text",x=2.4,y=2.8,label=L_F,parse=TRUE,hjust=0)+
  annotate("text",x=2.4,y=2.58,label=L_P,parse=TRUE,hjust=0)
p_FS1

ggsave("out/FS1_IN_ratio.pdf",p_FS1,width = 6,height = 4.5,device = cairo_pdf)
