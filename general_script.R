########################## LOADING LIBRARIES ##################################

### old package
install.packages("0_data/mvpart_1.6-2.tar.gz", repos = NULL, type = "source")
require("mvpart")

### maintained packages
if (!require("tidyverse")) install.packages("tidyverse"); require("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("nlme")) install.packages("nlme"); library("nlme")
if (!require("vegan")) install.packages("vegan"); library("vegan")

############################### LOADING DATA ###################################

### environmental data
comm_data = read.csv("0_data/comm_data.csv", sep=",", h=T)
head(comm_data)

### species data
ind_data = read.csv("0_data/ind_data.csv", sep=",", h=T, stringsAsFactors = T)
head(ind_data)

############################ PROCESSING DATA ##################################

### removing non-used info from species data
ind_data  =  ind_data %>% 
  filter(!is.na(sp))

### transforming fire
comm_data = comm_data %>% 
  mutate(trans_fire = fire_frequency^(1/3) )

### scaled predictors
scale_comm = comm_data %>% 
  mutate(fire_frequency =scale(fire_frequency),
         seasonal_precipitation= scale(seasonal_precipitation),
         soil_PC1 = scale(soil_PC1),
         soil_PC2 = scale(soil_PC2)
  )

################################# FIRE REGIME #################################

### check residual
check_resid = function(model){
  N = model$dims$N
  resid_n = resid(model)[1:N]
  shapiro.test(resid_n)
}

##### GLS

### my model
gls_fire1 = gls(
    data = scale_comm,
    trans_fire ~ seasonal_precipitation + soil_PC1 + soil_PC2, 
    method = "REML",
    correlation = corLin(form = ~ longitude + latitude))

summary(gls_fire1)
plot(gls_fire1)
check_resid(model = gls_fire1)


############################## COMMUNITY METRICS ###############################

##### species richeness

hist(comm_data$richeness)

### species richness model
gls_rich1 = gls(
    data = scale_comm,
    richness ~ fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2, 
    method = "REML",
    correlation = corLin(form = ~ longitude + latitude))

summary(gls_rich1)
plot(gls_rich1)
check_resid(model = gls_rich1)

### species abundance model
gls_abund1 = gls(data = comm_data,
    density ~ fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2, 
    method = "REML",
    correlation = corLin(form = ~ longitude + latitude))

summary(gls_abund1)
plot(gls_abund1)
check_resid(model = gls_abund1)

##### species diversity

hist(comm_data$fisher)

### species diversity model
gls_div1 = gls(data = scale_comm,
                fisher ~ fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2, 
                method = "REML",
                correlation = corLin(form = ~ longitude + latitude))

summary(gls_div1)
plot(gls_div1)
check_resid(model = gls_div1)

##### species composition

hist( comm_data$floristic_PCo2 )

### species compostion model
gls_comp1 = gls(data = scale_comm,
               floristic_PCo2 ~ fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2, 
               method = "REML",
               correlation = corLin(form = ~ longitude + latitude))

summary(gls_comp1)
plot(gls_comp1)
check_resid(model = gls_comp1)


ggplot(comm_data)+
geom_point(aes(x= floristic_PCo1, y=floristic_PCo2, size= fire_frequency))

########################### SPECIES COMPOSITION ############################

### all species names
spp_names = as.character(unique(ind_data$sp))

### species presence per plot
plot_spp = ind_data %>% 
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
plot_env = comm_data %>% 
  unite("plot_id", site, plot, sep="_", remove = FALSE) %>% 
  mutate(plot_id = as.character(plot_id)) %>% 
  dplyr::select(any_of(c("plot_id", 
                         "fire_frequency", 
                         "seasonal_precipitation",
                         "soil_PC1"
                         ) 
                       ) 
                )

### transforming in matrix
plot_env_mtx = as.matrix(plot_env[,-1])
rownames(plot_env_mtx) = plot_env$plot_id

### coordinates per plot
plot_coords = comm_data %>% 
  unite("plot_id", site, plot, sep="_", remove = FALSE) %>% 
  select("plot_id", "longitude", "latitude")

### transforming in matrix
plot_coords_mtx = as.matrix(plot_coords[,-1])
rownames(plot_coords_mtx) = plot_coords$plot_id

### dissimilarity matrices
spp_dist = xdiss(plot_spp_mtx, method = "bray")
env_dist = dist(plot_env_mtx, method = "euclidean")
geo_dist = dist(plot_coords_mtx, method = "euclidean")

### MANTEL
mantel_test = vegan::mantel.partial(
                      xdis = env_dist,
                      ydis = spp_dist,
                      zdis = geo_dist,
                      method="pearson", 
                      permutations=9999
                      )

plot(env_dist, spp_dist)

### dist by dist plot

