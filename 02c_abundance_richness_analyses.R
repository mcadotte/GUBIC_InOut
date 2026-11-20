#examine abundance richness curves

source("01_data_import.R")
source("02_sample_coverage_richness.R")

# 1. Prepare records and location labels
dat <- species_counts %>%
  filter(Location %in% c("Urban", "0–5 km", "5–30 km"),
    Status %in% c("Introduced", "Native")) %>%
  mutate(ID_HDC_G0_2 = as.character(ID_HDC_G0_2),location = recode(Location,
      "Urban" = "urban","0–5 km" = "5km","5–30 km" = "30km" ),
    status = recode(Status,"Introduced" = "intro","Native" = "native" ) )

# 2. Richness and abundance for each city/location/status
by_location <- dat %>% group_by(ID_HDC_G0_2, location, status) %>%
  summarise(richness = n_distinct(
    accepted_taxon_name_binomial[records > 0 &!is.na(records) &
          !is.na(accepted_taxon_name_binomial)]),
    abundance = sum(records, na.rm = TRUE),.groups = "drop") %>%
  pivot_wider(names_from = c(status, location),
    values_from = c(richness, abundance),
    names_glue = "{.value}_{status}_{location}",values_fill = 0)

# 3. Totals across all locations
by_total <- dat %>%group_by(ID_HDC_G0_2, status) %>%
  summarise(richness = n_distinct(
      accepted_taxon_name_binomial[records > 0 &!is.na(records) &
          !is.na(accepted_taxon_name_binomial) ]),
    abundance = sum(records, na.rm = TRUE),.groups = "drop") %>%
  pivot_wider(names_from = status,values_from = c(richness, abundance),
    names_glue = "{.value}_{status}_total",values_fill = 0)

# 4. Combined totals (native + introduced)
all_total <- dat %>% group_by(ID_HDC_G0_2) %>%
  summarise(richness_all_total = n_distinct(
      accepted_taxon_name_binomial[records > 0 &!is.na(records) &
          !is.na(accepted_taxon_name_binomial) ] ),
    abundance_all_total = sum(records, na.rm = TRUE),.groups = "drop" )

# 5. Join and arrange columns
city_diversity <- by_location %>%
  full_join(by_total, by = "ID_HDC_G0_2") %>%
  full_join(all_total, by = "ID_HDC_G0_2") %>%
  select(ID_HDC_G0_2,richness_intro_urban,abundance_intro_urban,
    richness_native_urban,abundance_native_urban,
    
    richness_intro_5km,abundance_intro_5km,
    richness_native_5km,abundance_native_5km,
    
    richness_intro_30km,abundance_intro_30km,
    richness_native_30km,abundance_native_30km,
    
    richness_intro_total,abundance_intro_total,
    richness_native_total, abundance_native_total,
    
    richness_all_total,abundance_all_total)

# Examine result
head(city_diversity)
dim(city_diversity)

theme_fig <- theme_classic(base_size = 13) +
  theme(legend.position = "top",legend.title = element_blank(),
        axis.text.x = element_text(size = 11),plot.tag = element_text(face = "bold", size = 16) )


p_RA_I_U<-ggplot(city_diversity,aes(x=abundance_intro_urban,y=richness_intro_urban))+
  geom_point(size = 3,alpha=0.6,shape = 21, color = "grey20", fill = "steelblue4")+
  scale_x_log10()+
  scale_y_log10()+
  labs(x="Log(Abundance)",y="Log(Species richness)",title = "A) Introduced:Urban")+
  theme_fig
  
p_RA_N_U<- ggplot(city_diversity,aes(x=abundance_native_urban,y=richness_native_urban))+
  geom_point(size = 3,alpha=0.6, shape = 21, color = "grey20", fill = "steelblue4")+
  scale_x_log10()+
  scale_y_log10()+
  labs(x="Log(Abundance)",y="Log(Species richness)",title = "B) Native:Urban")+
  theme_fig 

p_RA_I_5<-ggplot(city_diversity,aes(x=abundance_intro_5km,y=richness_intro_5km))+
  geom_point(size = 3,alpha=0.6, shape = 21, color = "grey20", fill = "steelblue4")+
  scale_x_log10()+
  scale_y_log10()+
  labs(x="Log(Abundance)",y="Log(Species richness)",title = "C) Introduced:0-5km")+
  theme_fig  
  
p_RA_N_5<-ggplot(city_diversity,aes(x=abundance_native_5km,y=richness_native_5km))+
  geom_point(size = 3,alpha=0.6, shape = 21, color = "grey20", fill = "steelblue4")+
  scale_x_log10()+
  scale_y_log10()+
  labs(x="Log(Abundance)",y="Log(Species richness)",title = "D) Native:0-5km")+
  theme_fig      
  
p_RA_I_30<- ggplot(city_diversity,aes(x=abundance_intro_30km,y=richness_intro_30km))+
  geom_point(size = 3,alpha=0.6, shape = 21, color = "grey20", fill = "steelblue4")+
  scale_x_log10()+
  scale_y_log10()+
  labs(x="Log(Abundance)",y="Log(Species richness)",title = "E) Introduced:5-30km")+
  theme_fig 
  
p_RA_N_30<-ggplot(city_diversity,aes(x=abundance_native_30km,y=richness_native_30km))+
  geom_point(size = 3,alpha=0.6, shape = 21, color = "grey20", fill = "steelblue4")+
  scale_x_log10()+
  scale_y_log10()+
  labs(x="Log(Abundance)",y="Log(Species richness)",title = "F) Native:5-30km")+
  theme_fig   

p_RA_I_Tot<-  ggplot(city_diversity,aes(x=abundance_intro_total,y=richness_intro_total))+
  geom_point(size = 3,alpha=0.6, shape = 21, color = "grey20", fill = "steelblue4")+
  scale_x_log10()+
  scale_y_log10()+
  labs(x="Log(Abundance)",y="Log(Species richness)",title = "GE) Introduced:total")+
  theme_fig 
  
p_RA_N_Tot<- ggplot(city_diversity,aes(x=abundance_native_total,y=richness_native_total))+
  geom_point(size = 3,alpha=0.6, shape = 21, color = "grey20", fill = "steelblue4")+
  scale_x_log10()+
  scale_y_log10()+
  labs(x="Log(Abundance)",y="Log(Species richness)",title = "H) Native:total")+
  theme_fig  

p_RA_ALL_Tot<- ggplot(city_diversity,aes(x=abundance_all_total,y=richness_all_total))+
  geom_point(size = 3,alpha=0.6, shape = 21, color = "grey20", fill = "steelblue4")+
  scale_x_log10()+
  scale_y_log10()+
  labs(x="Log(Abundance)",y="Log(Species richness)",title = "I) All:total")+
  theme_fig   

p_comb9 <- ((p_RA_I_U | p_RA_N_U | p_RA_I_5) /
              (p_RA_N_5 | p_RA_I_30 | p_RA_N_30) /
              (p_RA_I_Tot | p_RA_N_Tot | p_RA_ALL_Tot))

ggsave("out/FSx_abund-richness.pdf",p_comb9,width = 10,height = 10,device = cairo_pdf)
