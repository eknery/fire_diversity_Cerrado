########################## LOADING LIBRARIES ##################################
if (!require("tidyverse")) install.packages("tidyverse"); require("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("nlme")) install.packages("nlme"); library("nlme")
if (!require("philentropy")) install.packages("philentropy"); library("philentropy")
if (!require("ape")) install.packages("ape"); library("ape")

############################### LOADING DATA ###################################

### environmental data
comm_data = read.csv("0_data/comm_data.csv", sep=",", h=T)
head(comm_data)

### species data
spp_data = read.csv("0_data/spp_data.csv", sep=",", h=T, stringsAsFactors = T)
head(spp_data)

############################ PROCESSING DATA ##################################

### removing non-used info from species data
spp_data  =  spp_data %>% 
  filter(!is.na(sp))

### transforming fire
comm_data = comm_data %>% 
  mutate(trans_fire = fire_frequency^(1/3) )

################################# FIRE REGIME #################################

### check residual
check_resid = function(model){
  N = model$dims$N
  resid_n = resid(model)[1:N]
  shapiro.test(resid_n)
}

##### GLS

### my model
gls_fire1 = gls(data = comm_data,
    trans_fire ~ seasonal_precipitation + soil_PC1 + soil_PC2, 
    method = "REML",
    correlation = corGaus(form = ~ longitude + latitude, nugget = TRUE))

summary(gls_fire1)
plot(gls_fire1)
check_resid(model = gls_fire1)

### alternative model
gls_fire2 = gls(data = comm_data,
               trans_fire ~ site, 
               method = "REML",
               correlation = corGaus(form = ~ longitude + latitude, nugget = TRUE))

summary(gls_fire2)
plot(gls_fire2)
check_resid(model = gls_fire2)

############################## COMMUNITY METRICS ###############################

##### richeness

hist(comm_data$richeness)

### simple model
gls_rich1 = gls(data = comm_data,
                richness ~ fire_frequency, 
                method = "REML",
                correlation = corGaus(form = ~ longitude + latitude, nugget = TRUE))

summary(gls_rich1)
plot(gls_rich1)
check_resid(model = gls_rich1)

### full model
gls_rich2 = gls(data = comm_data,
    richness ~ fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2, 
    method = "REML",
    correlation = corGaus(form = ~ longitude + latitude, nugget = TRUE))

summary(gls_rich2)
plot(gls_rich2)
check_resid(model = gls_rich1)


##### species diversity

hist(comm_data$fisher)

### simple model
gls_div1 = gls(data = comm_data,
                fisher ~ fire_frequency, 
                method = "REML",
                correlation = corGaus(form = ~ longitude + latitude, nugget = TRUE))

summary(gls_div1)
plot(gls_div1)
check_resid(model = gls_div1)

### full model
gls_div2 = gls(data = comm_data,
                fisher ~ fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2, 
                method = "REML",
                correlation = corGaus(form = ~ longitude + latitude, nugget = TRUE))

summary(gls_div2)
plot(gls_div2)
check_resid(model = gls_div2)



########################### SPECIES COMPOSITION ############################

### all species names
spp_names = as.character(unique(spp_data$sp))

### species presence per plot
plot_spp = spp_data %>% 
  unite("plot_id", site, plot, sep="_", remove = FALSE) %>% 
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

### fire per plot
plot_fire = comm_data %>% 
  unite("plot_id", site, plot, sep="_", remove = FALSE) %>% 
  mutate(plot_id = as.character(plot_id)) %>% 
  dplyr::select(any_of(c("plot_id", "fire_frequency") ) )

### transforming in matrix
plot_fire_mtx = as.matrix(plot_fire[,-1])
rownames(plot_fire_mtx) = plot_fire$plot_id

### fire per plot
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


### dissimilarity matrices
spp_dist = distance(plot_spp_mtx, "sorensen")
fire_dist = distance(plot_fire_mtx, "euclidean")
env_dist = distance(plot_env_mtx, "euclidean")

ape::mantel.test(spp_dist, env_dist, 
                 nperm = 999,
                 alternative = "greater",
                 graph = T)




