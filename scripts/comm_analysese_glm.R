########################## LOADING LIBRARIES ##################################

### packages
if (!require("tidyverse")) install.packages("tidyverse"); require("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("ggpubr")) install.packages("ggpubr"); library("ggpubr")
if (!require("ape")) install.packages("ape"); library("ape")
if (!require("vegan")) install.packages("vegan"); library("vegan")
if (!require("spdep")) install.packages("spdep"); library("spdep")
source("scripts/function_plot_glm.R")

############################### LOADING DATA ###################################

### environmental data
comm_data = read.csv("0_data/comm_data.csv", sep=",", h=T)
head(comm_data)

### species data
ind_data = read.csv("0_data/ind_data.csv", sep=",", h=T, stringsAsFactors = T)
head(ind_data)

############################ PROCESSING DATA ##################################

### processing individual data
ind_data  =  ind_data %>% 
  filter(!is.na(sp)) %>% 
  unite("plot_id", site, plot, sep="_", remove = FALSE)

### processing community data
comm_data  =  comm_data %>% 
  unite("plot_id", site, plot, sep="_", remove = FALSE)

### processing variables
comm_data = comm_data %>% 
  mutate(
    burned = burned * 10, # response fire
    fire_frequency = burned * 10 # explanatory fire
  ) %>% 
  mutate(
    fire_frequency = scale(fire_frequency),
    seasonal_precipitation = scale(seasonal_precipitation),
    soil_PC1 = scale(soil_PC1),
    soil_PC2 = scale(soil_PC2)
  )

### renaming to avoid conflict predictors
p_comm_data = comm_data %>% 
  rename(
    p_burned = burned,
    p_richness = richness,
    p_fisher = fisher,
    p_floristic_PCo1 = floristic_PCo1,
    p_floristic_PCo2 = floristic_PCo2,
    p_fire_frequency = fire_frequency,
    p_seasonal_precipitation = seasonal_precipitation,
    p_soil_PC1 = soil_PC1,
    p_soil_PC2 = soil_PC2
  )

### summarize by site
s_comm_data = p_comm_data %>% 
  group_by(site) %>%  
  summarise(
    burned = round(median(p_burned, na.rm = T)),
    burned_min = min(p_burned, na.rm = T),
    burned_max = max(p_burned, na.rm = T),
    richness = round(mean(p_richness, na.rm = T)),
    richness_min = round(mean(p_richness, na.rm = T)) - sd(p_richness, na.rm = T),
    richness_max = round(mean(p_richness, na.rm = T)) + sd(p_richness, na.rm = T),
    fisher = mean(p_fisher, na.rm = T),
    fisher_min = min(p_fisher, na.rm = T),
    fisher_max = max(p_fisher, na.rm = T),
    floristic_PCo1 = mean(p_floristic_PCo1, na.rm = T),
    floristic_PCo1_min = min(p_floristic_PCo1, na.rm = T),
    floristic_PCo1_max = max(p_floristic_PCo1, na.rm = T),
    floristic_PCo2 = mean(p_floristic_PCo2, na.rm = T),
    floristic_PCo2_min = min(p_floristic_PCo2, na.rm = T),
    floristic_PCo2_max = max(p_floristic_PCo2, na.rm = T),
    fire_frequency = round(median(p_fire_frequency, na.rm = T)),
    seasonal_precipitation =  mean(p_seasonal_precipitation, na.rm = T),
    soil_PC1 =  mean(p_soil_PC1, na.rm = T),
    soil_PC2 =  mean(p_soil_PC2, na.rm = T),
    longitude = mean(longitude, na.rm = T),
    latitude = mean(latitude, na.rm = T)
  )

########################## OVERALL PARAMETERS ##########################

### all explanatory vars
all_explanatory = c("fire_frequency", 
                    "seasonal_precipitation",
                    "soil_PC1",
                    "soil_PC2"
)

### all x-axis labels
all_xlabels = c("Fire frequency",
                "Seasonal precipitation",
                "Soil PC1",
                "Soil PC2"
                
)

############################## FLORISTIC ORDINATION ############################

### sp presence matrix
plot_mtx = ind_data %>% 
  group_by(plot_id, sp) %>% 
  reframe(n = n()) %>% 
  mutate(n = case_when(n >1 ~ 1, T ~ n)) %>% 
  pivot_wider(
    names_from =  sp,
    names_expand = T,
    values_from = n,
    values_fill = 0
  )

### dissimilarity
distance = vegdist(plot_mtx[,-1], method = "bray")
exdistance = stepacross(dis = distance, path = "extended")

### pcoa
ord = pcoa(exdistance)
ord$values
sum(ord$values$Relative_eig)

################################# AUTORERESSION ###############################

### spatial neighbors list
coords = cbind(s_comm_data$longitude, s_comm_data$latitude)
neighbors = knearneigh(coords, k = 2)  # k-nearest neighbors
listw = nb2listw(knn2nb(neighbors), style = "W")

################################# FIRE REGIME #################################

### calculate sac
s_comm_data$sac = lag.listw(listw, s_comm_data$burned)

### glm model
glm_fire1 = glm(
  data = s_comm_data,
  burned ~ seasonal_precipitation + soil_PC1 + soil_PC2, 
  family = poisson("log")
)

# identity AIC 27.524
# log AIC 27.08

summary(glm_fire1)
plot(glm_fire1)
shapiro.test( residuals(glm_fire1) )

### plot model? 
show_model = c(FALSE,FALSE,FALSE, FALSE)
## graphical parameter
tiff("1_plots/fire_plots.tiff", 
     units="cm", width=14, height=14, res=600)
par(mfrow = c(2,2))
par(mar = c(4.5, 4, 1, 1))
### plots 
for(i in 1:length(all_explanatory) ){
  ## variables names
  x = all_explanatory[i]
  y = "burned"
  ## get variables
  pred = as.numeric(s_comm_data[[x]])
  resp = as.numeric(s_comm_data[[y]])
  ## simple model for each predictor
  model = glm(resp ~ pred, poisson("log"))
  ## plot model
  plot_glm(
    data = s_comm_data,
    x = x,
    y = y, 
    model = model,
    show_limits = F, 
    show_model = show_model[i],
    x_label = all_xlabels[i], 
    y_label = "Fire frequency"
  )
}
dev.off()

### export fire map
tiff("1_plots/fire_map.tiff", 
     units="cm", width=6.5, height=6.5, res=600)
ggplot(data = comm_data) +
  geom_point(aes(x = longitude, 
                 y = latitude,
                 size = fire_frequency
  ),
  alpha = 0.20,
  position = position_jitter(height= 0.25, width = 0.25)
  ) +
  scale_size(
    breaks = c(1,3,6),
    labels = c(1,3,6),
    guide = "legend"
  )+
  xlab("Longitude") +
  ylab("Latitude") +
  theme(panel.background=element_rect(fill="white"),
        panel.grid=element_line(colour=NULL),
        panel.border=element_rect(fill=NA,colour="black"),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(size= 10, angle = 0),
        axis.text.y = element_text(size= 10, angle = 0),
        legend.title= element_blank(),
        legend.key.size =  unit(0.05, 'cm'),
        legend.background = element_rect(colour = "black", linetype='solid', fill = NA),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.direction="vertical",
        legend.position = c(0.8,0.70)
  )
dev.off()

############################## SPECIES RICHNESS ###############################

##### species richeness
hist(log(s_comm_data$richness) )
### compute SAC
s_comm_data$sac = scale(lag.listw(listw, s_comm_data$richness))

### species richness model
glm_rich1 = glm(
  data = s_comm_data,
  richness ~ fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2, 
  family = poisson(link = "sqrt")
)

### model summary
summary(glm_rich1)
plot(glm_rich1)
shapiro.test(resid(glm_rich1))

### plot model? 
show_model = c(TRUE,TRUE, FALSE, FALSE)
## graphical parameter
tiff("1_plots/richness_plots.tiff", 
     units="cm", width=14, height=14, res=600)
par(mfrow = c(2,2))
par(mar = c(4.5, 4, 1, 1))
### plots 
for(i in 1:length(all_explanatory) ){
  ## variables names
  x = all_explanatory[i]
  y = "richness"
  ## get variables
  pred = as.numeric(s_comm_data[[x]])
  resp = as.numeric(s_comm_data[[y]])
  ## simple model for each predictor
  model = glm(resp ~ pred, poisson("sqrt"))
  ## plot model
  plot_glm(
    data = s_comm_data,
    x = x,
    y = y, 
    model = model,
    show_model = show_model[i],
    x_label = all_xlabels[i], 
    y_label = "Species richness"
  )
}
dev.off()

#################################### FISHER ALPHA ##############################

### compute SAC
s_comm_data$sac <- scale(lag.listw(listw, s_comm_data$fisher))

### species diversity model
glm_fish1 = glm(
  data = s_comm_data,
  fisher ~ fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2, 
  family = inverse.gaussian(link = "identity")
)

### model summary
summary(glm_fish1)
plot(glm_fish1)
shapiro.test(resid(glm_fish1))

### plot model? 
show_model = c(T, T, F, F)
## graphical parameter
tiff("1_plots/fisher_plots.tiff", 
     units="cm", width=14, height=14, res=600)
par(mfrow = c(2,2))
par(mar = c(4.5, 4, 1, 1))
### plots 
for(i in 1:length(all_explanatory) ){
  ## variables names
  x = all_explanatory[i]
  y = "fisher"
  ## get variables
  pred = as.numeric(s_comm_data[[x]])
  resp = as.numeric(s_comm_data[[y]])
  ## simple model for each predictor
  model = glm(resp ~ pred, inverse.gaussian("identity"))
  ## plot model
  plot_glm(
    data = s_comm_data,
    x = x,
    y = y, 
    model = model,
    show_model = show_model[i],
    x_label = all_xlabels[i], 
    y_label = "Fisher's alpha"
  )
}
dev.off()

############################## SPECIES COMPOSITION #############################

###### FIRST AXIS
hist(s_comm_data$floristic_PCo1)

### compute SAC
s_comm_data$sac <- scale(lag.listw(listw, s_comm_data$floristic_PCo1) )

### species compostion model 1
glm_comp1 = glm(
  data = s_comm_data,
  floristic_PCo1 ~ fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2, 
  family = gaussian(link = "identity")
)

## link function
# identity AIC: -10.17.. significant
# inverse AIC: -17.45

summary(glm_comp1)
plot(glm_comp1)
shapiro.test(resid(glm_comp1) )

### plot model? 
show_model = c(T, F, F, F)
## graphical parameter
tiff("1_plots/pco1_plots.tiff", 
     units="cm", width=14, height=14, res=600)
par(mfrow = c(2,2))
par(mar = c(4.5, 4, 1, 1))
### plots 
for(i in 1:length(all_explanatory) ){
  ## variables names
  x = all_explanatory[i]
  y = "floristic_PCo1"
  ## get variables
  pred = as.numeric(s_comm_data[[x]])
  resp = as.numeric(s_comm_data[[y]])
  ## simple model for each predictor
  model = glm(resp ~ pred, gaussian("identity"))
  ## plot model
  plot_glm(
    data = s_comm_data,
    x = x,
    y = y, 
    model = model,
    show_model = show_model[i],
    x_label = all_xlabels[i], 
    y_label = "Floristic PCo1"
  )
}
dev.off()

###### SECOND AXIS
hist(comm_data$floristic_PCo2)

### compute SAC
s_comm_data$sac = scale(lag.listw(listw, s_comm_data$floristic_PCo2))

### species compostion model 2
glm_comp2 = glm(
  data = s_comm_data ,
  floristic_PCo2 ~ fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2, 
  family = gaussian(link = "identity")
)

### summary model
summary(glm_comp2)
297plot(glm_comp2)
shapiro.test(resid(glm_comp2) )

### plot model? 
show_model = c(F, F, F, F)
## graphical parameter
tiff("1_plots/pco2_plots.tiff", 
     units="cm", width=14, height=14, res=600)
par(mfrow = c(2,2))
par(mar = c(4.5, 4, 1, 1))
### plots 
for(i in 1:length(all_explanatory) ){
  ## variables names
  x = all_explanatory[i]
  y = "floristic_PCo2"
  ## get variables
  pred = as.numeric(s_comm_data[[x]])
  resp = as.numeric(s_comm_data[[y]])
  ## simple model for each predictor
  model = glm(resp ~ pred, gaussian(link = "identity"))
  ## plot model
  plot_glm(
    data = s_comm_data,
    x = x,
    y = y, 
    model = model,
    show_model = show_model[i],
    x_label = all_xlabels[i], 
    y_label = "Floristic PCo2"
  )
}
dev.off()
