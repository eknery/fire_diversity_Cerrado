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
  mutate(fire_frequency = case_when(
    fire_frequency == 0 ~ 0.1,
    TRUE                ~ fire_frequency
    )
  ) %>% 
  mutate(trans_fire = log(fire_frequency) )

### scaled predictors
s_comm_data = comm_data %>% 
  mutate(
         fire_frequency = scale(fire_frequency),
         seasonal_precipitation = scale(seasonal_precipitation),
         soil_PC1 = scale(soil_PC1),
         soil_PC2 = scale(soil_PC2)
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
gls_fire1 = gls(
    data = s_comm_data,
    trans_fire ~ seasonal_precipitation + soil_PC1 + soil_PC2, 
    method = "ML",
    correlation = corLin(form = ~ longitude + latitude))

summary(gls_fire1)
plot(gls_fire1)
check_resid(model = gls_fire1)

###### fire plots
fire_plots = list()
fire_relationships = c("none","linear", "linear", "none")

### plots 
for(i in 2:length(all_explanatory) ){
  
  fire_plots[[i]] = model_plot(data = s_comm_data,
             x = all_explanatory[i],
             y = "trans_fire", 
             model = gls_fire1,
             relationship = fire_relationships[i],
             x_label = all_xlabels[i], 
             y_label = "ln(Fire frequency)")
}

### export plots
tiff("plots/fire_plots.tiff", units="cm", width=14, height=14, res=600)
  ggarrange(fire_plots[[2]], fire_plots[[2]],fire_plots[[3]], fire_plots[[4]],
            labels = c("", "A", "B", "C"),
            ncol = 2, nrow = 2)
dev.off()

### export fire map
tiff("plots/fire_map.tiff", units="cm", width=7, height=7, res=600)
  ggplot(data = comm_data) +
    geom_point(aes(x = longitude, 
                   y = latitude,
                   size = fire_frequency*enhance,
                   ),
               alpha = 0.20
    ) +
    scale_size(
      breaks = c(0.2,1,5, 25)*3,
      labels = c(0.2,1,5, 25),
      guide = "legend"
    )+
    xlab("Longitude") +
    ylab("Latitude") +
    theme(panel.background=element_rect(fill="white"),
          panel.grid=element_line(colour=NULL),
          panel.border=element_rect(fill=NA,colour="black"),
          axis.title.x = element_text(size=12, face="bold"),
          axis.title.y = element_text(size=12, face="bold"),
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
hist(log(comm_data$richness) )

### species richness model
gls_rich1 = gls(
    data = s_comm_data,
    log(richness) ~ fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2, 
    method = "REML",
    correlation = corLin(form = ~ longitude + latitude))

### model summary
summary(gls_rich1)
plot(gls_rich1)
check_resid(model = gls_rich1)

###### richness plots
rich_plots = list()
rich_relationships = c("linear","linear", "none", "none")

### plots 
for(i in 1:length(all_explanatory) ){
  
  rich_plots[[i]] = model_plot(data = s_comm_data %>% mutate(richness = log(richness)),
                               x = all_explanatory[i],
                               y = "richness", 
                               model = gls_rich1,
                               relationship = rich_relationships[i],
                               x_label = all_xlabels[i], 
                               y_label = "ln(N species per plot)")
}

tiff("plots/richness_plots.tiff", units="cm", width=14, height=14, res=600)
ggarrange(rich_plots[[1]], rich_plots[[2]], rich_plots[[3]], rich_plots[[4]],
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
dev.off()


############################### SPECIES ABUNDANCE #############################

gls_abund1 = gls(data = s_comm_data,
    density ~ fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2, 
    method = "REML",
    correlation = corLin(form = ~ longitude + latitude))

summary(gls_abund1)
plot(gls_abund1)
check_resid(model = gls_abund1)

#################################### FISHER ALPHA ##############################

hist(comm_data$fisher)

### species diversity model
gls_fish1 = gls(data = s_comm_data,
                log(fisher) ~ fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2, 
                method = "REML",
                correlation = corLin(form = ~ longitude + latitude))

### model summary
summary(gls_fish1)
plot(gls_fish1)
check_resid(model = gls_fish1)

###### richness plots
fish_plots = list()
fish_relationships = c("linear","linear", "none", "none")

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

###### FISRT AXIS
hist( comm_data$floristic_PCo1 )

### species compostion model 1
gls_comp1 = gls(data = s_comm_data,
               (floristic_PCo1) ~ fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2, 
               method = "REML",
               correlation = corLin(form = ~ longitude + latitude))

summary(gls_comp1)
plot(gls_comp1)
check_resid(model = gls_comp1)

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
gls_comp2 = gls(data = s_comm_data ,
                floristic_PCo2 ~ fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2, 
                method = "REML",
                correlation = corLin(form = ~ longitude + latitude))

### summary model
summary(gls_comp2)
plot(gls_comp2)
check_resid(model = gls_comp2)

###### composition2 plots
comp2_plots = list()
comp2_relationships = c("linear","none", "none", "none")

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

############################### OTHER MODELS ##################################

ggplot(data = subset(s_comm_data, fire_frequency < 4 )  ) +
  geom_point(aes(
    x = fire_frequency,
    y = floristic_PCo2
      )
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

