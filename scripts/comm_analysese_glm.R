########################## LOADING LIBRARIES ##################################

### maintained packages
if (!require("tidyverse")) install.packages("tidyverse"); require("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("ggpubr")) install.packages("ggpubr"); library("ggpubr")
if (!require("ape")) install.packages("ape"); library("ape")
if (!require("vegan")) install.packages("vegan"); library("vegan")

### my functions
source("scripts/function_model_plot.R")

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

### transforming in matrix
plot_spp_mtx = as.matrix(plot_spp[,-1])
rownames(plot_spp_mtx) = plot_spp$plot_id

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

### transforming in matrix
plot_env_mtx = as.matrix(plot_env[,-1])
rownames(plot_env_mtx) = plot_env$plot_id

### each environ var per plot
plot_fire_mtx = as.matrix(plot_env[,"fire_frequency"])
plot_precip_mtx = as.matrix(plot_env[,"seasonal_precipitation"])
plot_soil1_mtx = as.matrix(plot_env[,"soil_PC1"])
plot_soil2_mtx = as.matrix(plot_env[,"soil_PC2"])

rownames(plot_fire_mtx) =  plot_env$plot_id
rownames(plot_precip_mtx) =  plot_env$plot_id
rownames(plot_soil1_mtx) =  plot_env$plot_id
rownames(plot_soil2_mtx) =  plot_env$plot_id

### coordinates per plot
plot_coords = comm_data %>% 
  unite("plot_id", site, plot, sep="_", remove = FALSE) %>% 
  select("plot_id", "longitude", "latitude")

### transforming in matrix
plot_coords_mtx = as.matrix(plot_coords[,-1])
rownames(plot_coords_mtx) = plot_coords$plot_id

### dissimilarity matrices
env_dist = dist(plot_env_mtx, method = "euclidean")
fire_dist = dist(plot_fire_mtx, method = "euclidean")
precip_dist = dist(plot_precip_mtx, method = "euclidean")
soil1_dist = dist(plot_soil1_mtx, method = "euclidean")
soil2_dist = dist(plot_soil2_mtx, method = "euclidean")
geo_dist = dist(plot_coords_mtx, method = "euclidean")

### MANTEL
mantel.test(m1 = as.matrix(geo_dist), 
            m2 = as.matrix(fire_dist), 
            graph = FALSE, 
            nperm = 999, 
            alternative = "greater"
            )


### PARTIAL MANTEL
mantel_test = vegan::mantel.partial(
  xdis = env_dist,
  ydis = spp_dist,
  zdis = geo_dist,
  method="pearson", 
  permutations=999
)

plot(env_dist, spp_dist)

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

################################# FIRE REGIME #################################

### glm model
glm_fire1 = glm(
  data = s_comm_data,
  fire_frequency ~ seasonal_precipitation + soil_PC1 + soil_PC2, 
  family = poisson()
)

summary(glm_fire1)
plot(glm_fire1)
shapiro.test( residuals(glm_fire1) )

###### fire plots
fire_plots = list()
fire_relationships = c("linear","linear", "linear", "linear")

### plots 
for(i in 2:length(all_explanatory) ){
  
  fire_plots[[i]] = model_plot(data = s_comm_data,
                               x = all_explanatory[i],
                               y = "fire_frequency", 
                               model = glm_fire1,
                               relationship = fire_relationships[i],
                               x_label = all_xlabels[i], 
                               y_label = "Fire frequency")
}

### export plots
tiff("1_plots/fire_plots.tiff", units="cm", width=14, height=14, res=600)
ggarrange(fire_plots[[2]], fire_plots[[2]],fire_plots[[3]], fire_plots[[4]],
          labels = c("", "A", "B", "C"),
          ncol = 2, nrow = 2)
dev.off()

############################## SPECIES RICHNESS ###############################

##### species richeness
hist(log(s_comm_data$richness) )

### species richness model
glm_rich1 = glm(
  data = s_comm_data,
  richness ~ s_fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2, 
  family = poisson(link = "identity")
)

### model summary
summary(glm_rich1)
plot(glm_rich1)
shapiro.test(resid(glm_rich1))

### plot model? 
show_plot = c(TRUE,TRUE, FALSE, FALSE)

## graphical parameter
par(mfrow = c(2,2))
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
    show = show_plot[i],
    x_label = all_xlabels[i], 
    y_label = "Species richness"
  )
}

tiff("plots/richness_plots.tiff", units="cm", width=14, height=14, res=600)
dev.off()

############################### SPECIES ABUNDANCE #############################

glm_abund1 = glm(
  data = s_comm_data,
  density ~ s_fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2, 
  family = poisson(link = "identity")
)

summary(glm_abund1)
plot(glm_abund1)
shapiro.test(resid(glm_abund1))

#################################### FISHER ALPHA ##############################

### species diversity model
glm_fish1 = glm(
  data = s_comm_data,
  fisher ~ s_fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2, 
  family = Gamma(link = "log")
)

### model summary
summary(glm_fish1)
plot(glm_fish1)
shapiro.test(resid(glm_fish1))

###### richness plots
fish_plots = list()
fish_relationships = c("none","none", "none", "none")

### plots 
for(i in 1:length(all_explanatory) ){
  
  fish_plots[[i]] = model_plot(data = s_comm_data %>% mutate(fisher = log(fisher)),
                               x = all_explanatory[i],
                               y = "fisher", 
                               model = gls_fish1,
                               relationship = fish_relationships[i],
                               x_label = all_xlabels[i], 
                               y_label = "ln(Fisher's alpha)")
}

tiff("plots/fisher_plots.tiff", units="cm", width=14, height=14, res=600)
ggarrange(fish_plots[[1]], fish_plots[[2]], fish_plots[[3]], fish_plots[[4]],
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
dev.off()

############################## SPECIES COMPOSITION #############################

hist(s_comm_data$floristic_PCo1)

### species compostion model 1
glm_comp1 = glm(
  data = s_comm_data,
  floristic_PCo1 ~ s_fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2, 
  family = gaussian(link = "inverse")
)

summary(glm_comp1)
plot(glm_comp1)
shapiro.test(resid(glm_comp1) )

###### composition2 plots
comp1_plots = list()
comp1_relationships = c("none","none", "none", "none")

### plots 
for(i in 1:length(all_explanatory) ){
  
  comp1_plots[[i]] = model_plot(data = s_comm_data,
                                x = all_explanatory[i],
                                y = "floristic_PCo1", 
                                model = gls_comp1,
                                relationship = comp1_relationships[i],
                                x_label = all_xlabels[i], 
                                y_label = "Floristic PCo1")
}

tiff("plots/comp1_plots.tiff", units="cm", width=14, height=14, res=600)
ggarrange(comp1_plots[[1]], comp1_plots[[2]], comp1_plots[[3]], comp1_plots[[4]],
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
dev.off()


###### SECOND AXIS
hist(comm_data$floristic_PCo2)

### species compostion model 2
glm_comp2 = glm(
  data = s_comm_data ,
  floristic_PCo2 ~ s_fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2, 
  family = gaussian(link = "identity")
)

### summary model
summary(glm_comp2)
plot(glm_comp2)
check_resid(model = gls_comp2)

###### composition2 plots
comp2_plots = list()
comp2_relationships = c("linear","none", "none", "linear")

### plots 
for(i in 1:length(all_explanatory) ){
  
  comp2_plots[[i]] = model_plot(data = s_comm_data,
                                x = all_explanatory[i],
                                y = "floristic_PCo2", 
                                model = gls_comp2,
                                relationship = comp2_relationships[i],
                                x_label = all_xlabels[i], 
                                y_label = "Floristic PCo2")
}

tiff("plots/comp2_plots.tiff", units="cm", width=14, height=14, res=600)
ggarrange(comp2_plots[[1]], comp2_plots[[2]], comp2_plots[[3]], comp2_plots[[4]],
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
dev.off()



