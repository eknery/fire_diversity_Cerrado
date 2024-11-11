########################## LOADING LIBRARIES ##################################

### maintained packages
if (!require("tidyverse")) install.packages("tidyverse"); require("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("ggpubr")) install.packages("ggpubr"); library("ggpubr")
if (!require("ape")) install.packages("ape"); library("ape")
if (!require("nlme")) install.packages("nlme"); library("nlme")
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
  mutate(fire_frequency = burned*10)

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

### check residual
check_resid = function(model){
  N = model$dims$N
  resid_n = resid(model)[1:N]
  shapiro.test(resid_n)
}
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
##
ord$values

sum(ord$values$Relative_eig)

################################# FIRE REGIME #################################

### gls model
glm_fire1 = glm(
  data = s_comm_data,
  fire_frequency ~ seasonal_precipitation + soil_PC1 + soil_PC2, 
  family = poisson(link = "identity")
)

summary(glm_fire1)
plot(gls_fire1)
shapiro.test( residuals(glm_fire1) )

###### fire plots
fire_plots = list()
fire_relationships = c("none","linear", "linear", "none")

### plots 
for(i in 2:length(all_explanatory) ){
  
  fire_plots[[i]] = model_plot(data = s_comm_data,
                               x = all_explanatory[i],
                               y = "fire_frquency", 
                               model = glm_fire1,
                               relationship = fire_relationships[i],
                               x_label = all_xlabels[i], 
                               y_label = "ln(Fire frequency)")
}

