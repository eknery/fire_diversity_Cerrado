########################## LOADING LIBRARIES ##################################

### maintained packages
if (!require("tidyverse")) install.packages("tidyverse"); require("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("ggpubr")) install.packages("ggpubr"); library("ggpubr")
if (!require("ape")) install.packages("ape"); library("ape")
if (!require("vegan")) install.packages("vegan"); library("vegan")
if (!require("spdep")) install.packages("spdep"); library("spdep")

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

### transforming fire
comm_data = comm_data %>% 
  mutate(fire_frequency = burned * 10)

### scaled predictors
s_comm_data = comm_data %>% 
  mutate(
    s_fire_frequency = scale(fire_frequency),
    seasonal_precipitation = scale(seasonal_precipitation),
    soil_PC1 = scale(soil_PC1),
    soil_PC2 = scale(soil_PC2)
  )

hist(comm_data$fire_frequency)

########################## OVERALL PARAMETERS ##########################

### all explanatory vars
all_explanatory = c("s_fire_frequency", 
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

################################# ORDINATION ###################################

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

#################################### MANTEL TEST ################################

### all species names
spp_names = as.character(unique(ind_data$sp))

### species presence per plot
plot_spp = ind_data %>% 
  mutate(plot_id = as.character(plot_id)) %>% 
  group_by(plot_id) %>% 
  reframe(sp = as.character(unique(sp)), presence = 1) %>% 
  pivot_wider(names_from = sp,
              names_expand = T,
              values_from = presence,
              values_fill = as.integer(0)
  )

### diversity per plot
plot_div = comm_data %>% 
  unite("plot_id", site, plot, sep="_", remove = FALSE) %>% 
  mutate(plot_id = as.character(plot_id)) %>% 
  dplyr::select(any_of(c("plot_id", 
                         "richness", 
                         "fisher",
                         "floristic_PCo1",
                         "floristic_PCo2"
  ) 
  ) 
  )

### environment per plot
plot_env = comm_data %>% 
  unite("plot_id", site, plot, sep="_", remove = FALSE) %>% 
  mutate(plot_id = as.character(plot_id)) %>% 
  dplyr::select(any_of(c("plot_id", 
                         "fire_frequency", 
                         "seasonal_precipitation",
                         "soil_PC1",
                         "soil_PC2"
                          ) 
                      ) 
  )

### coordinates per plot
plot_coords = comm_data %>% 
  unite("plot_id", site, plot, sep="_", remove = FALSE) %>% 
  select("plot_id", "longitude", "latitude")

### transforming in matrix
plot_spp_mtx = as.matrix(plot_spp[,-1])
plot_env_mtx = as.matrix(plot_env[,-1])
plot_rich_mtx = as.matrix(plot_div[,"richness"])
plot_fish_mtx = as.matrix(plot_div[,"fisher"])
plot_flo1_mtx = as.matrix(plot_div[,"floristic_PCo1"])
plot_flo2_mtx = as.matrix(plot_div[,"floristic_PCo2"])
plot_fire_mtx = as.matrix(plot_env[,"fire_frequency"])
plot_precip_mtx = as.matrix(plot_env[,"seasonal_precipitation"])
plot_soil1_mtx = as.matrix(plot_env[,"soil_PC1"])
plot_soil2_mtx = as.matrix(plot_env[,"soil_PC2"])
plot_coords_mtx = as.matrix(plot_coords[,-1])

### naming rows
rownames(plot_spp_mtx) = plot_spp$plot_id
rownames(plot_rich_mtx) = plot_div$plot_id
rownames(plot_fish_mtx) = plot_div$plot_id
rownames(plot_flo1_mtx) = plot_div$plot_id
rownames(plot_flo2_mtx) = plot_div$plot_id
rownames(plot_env_mtx) = plot_env$plot_id
rownames(plot_fire_mtx) =  plot_env$plot_id
rownames(plot_precip_mtx) =  plot_env$plot_id
rownames(plot_soil1_mtx) =  plot_env$plot_id
rownames(plot_soil2_mtx) =  plot_env$plot_id
rownames(plot_coords_mtx) = plot_coords$plot_id

### dissimilarity matrices
distance = vegdist(plot_mtx[,-1], method = "bray")
spp_dist = stepacross(dis = distance, path = "extended")
rich_dist = dist(plot_rich_mtx, method = "euclidean")
fish_dist = dist(plot_fish_mtx, method = "euclidean")
flo1_dist = dist(plot_flo1_mtx, method = "euclidean")
flo2_dist = dist(plot_flo2_mtx, method = "euclidean")
env_dist = dist(plot_env_mtx, method = "euclidean")
fire_dist = dist(plot_fire_mtx, method = "euclidean")
precip_dist = dist(plot_precip_mtx, method = "euclidean")
soil1_dist = dist(plot_soil1_mtx, method = "euclidean")
soil2_dist = dist(plot_soil2_mtx, method = "euclidean")
geo_dist = dist(plot_coords_mtx, method = "euclidean")

### MANTEL
mantel.test(m1 = as.matrix(geo_dist), 
            m2 = as.matrix(flo2_dist), 
            graph = FALSE, 
            nperm = 999, 
            alternative = "greater"
            )

### PARTIAL MANTEL
mantel_test = vegan::mantel.partial(
  xdis = fire_dist,
  ydis = spp_dist,
  zdis = geo_dist,
  method="pearson", 
  permutations=999
)

### spatial neighbors list
coords <- cbind(s_comm_data$longitude, s_comm_data$latitude)
neighbors <- knearneigh(coords, k = 9)  # k-nearest neighbors
listw <- nb2listw(knn2nb(neighbors), style = "W")

################################# FIRE REGIME #################################

### calculate sac
s_comm_data$sac <- lag.listw(listw, s_comm_data$fire_frequency)

### glm model
glm_fire1 = glm(
  data = s_comm_data,
  #fire_frequency ~ seasonal_precipitation + soil_PC1 + soil_PC2,
  fire_frequency ~ seasonal_precipitation + soil_PC1 + soil_PC2 + sac, 
  family = poisson("log")
)

summary(glm_fire1)
plot(glm_fire1)
shapiro.test( residuals(glm_fire1) )

### plot model? 
show_model = c(FALSE,FALSE,TRUE, FALSE)
## graphical parameter
tiff("1_plots/fire_plots.tiff", 
     units="cm", width=14, height=14, res=600)
par(mfrow = c(2,2))
par(mar = c(4.5, 4, 1, 1))
### plots 
for(i in 1:length(all_explanatory) ){
  ## variables names
  x = all_explanatory[i]
  y = "fire_frequency"
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
s_comm_data$sac <- lag.listw(listw, s_comm_data$richness)
### species richness model
glm_rich1 = glm(
  data = s_comm_data,
  #richness ~ s_fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2, 
  richness ~ s_fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2 + sac, 
  family = poisson(link = "identity")
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
  model = glm(resp ~ pred, poisson("identity"))
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
s_comm_data$sac <- lag.listw(listw, s_comm_data$fisher)

### species diversity model
glm_fish1 = glm(
  data = s_comm_data,
  #fisher ~ s_fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2, 
  fisher ~ s_fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2 + sac, 
  family = Gamma(link = "identity")
)

### model summary
summary(glm_fish1)
plot(glm_fish1)
shapiro.test(resid(glm_fish1))

### plot model? 
show_model = c(TRUE, TRUE, FALSE, FALSE)
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
  model = glm(resp ~ pred, Gamma("identity"))
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
s_comm_data$sac <- lag.listw(listw, s_comm_data$floristic_PCo1)

### species compostion model 1
glm_comp1 = glm(
  data = s_comm_data,
  # floristic_PCo1 ~ s_fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2, 
  floristic_PCo1 ~ s_fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2 + sac, 
  family = gaussian(link = "identity")
)

summary(glm_comp1)
plot(glm_comp1)
shapiro.test(resid(glm_comp1) )

### plot model? 
show_model = c(F, F, T, F)
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
s_comm_data$sac <- lag.listw(listw, s_comm_data$floristic_PCo2)

### species compostion model 2
glm_comp2 = glm(
  data = s_comm_data ,
  # floristic_PCo2 ~ s_fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2, 
  floristic_PCo2 ~ s_fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2 + sac, 
  family = gaussian(link = "identity")
)

### summary model
summary(glm_comp2)
plot(glm_comp2)
shapiro.test(resid(glm_comp2) )

### plot model? 
show_model = c(F, F, T, F)
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
