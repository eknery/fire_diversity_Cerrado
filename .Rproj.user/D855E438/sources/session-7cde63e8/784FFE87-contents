
### maintained packages
if (!require("tidyverse")) install.packages("tidyverse"); require("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("ggpubr")) install.packages("ggpubr"); library("ggpubr")
if (!require("nlme")) install.packages("nlme"); library("nlme")
if (!require("vegan")) install.packages("vegan"); library("vegan")

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

### fire frequency
comm_data = comm_data %>% 
 mutate(fire_frequency = burned * 10)

############################## PLOTS and SPP ###########################

### overall spp abundance
ind_data %>% 
  group_by(sp) %>% 
  reframe( n = n()) %>% 
  arrange(desc(n))

####### FIRE
hist(comm_data$fire_frequency)
### high fire plots
high_fire_plot = comm_data %>% 
  filter(fire_frequency >= 5) %>% 
  select(plot_id) 

high_fire_sp = ind_data %>% 
  filter(plot_id %in% high_fire_plot$plot_id) %>% 
  group_by(sp) %>% 
  reframe(n = n()) %>% 
  arrange(desc(n))

### low fire plots
low_fire_plot = comm_data %>% 
  filter(fire_frequency <= 1) %>% 
  select(plot_id) 

low_fire_sp = ind_data %>% 
  filter(plot_id %in% low_fire_plot$plot_id) %>% 
  group_by(sp) %>% 
  reframe(n = n()) %>% 
  arrange(desc(n))

### comparing low and high fire
which(low_fire_sp$sp == high_fire_sp$sp[[1]])

##### SOIL NUTRIENTS
hist(comm_data$soil_PC1)
### high nutri plots
high_nutri_plot = comm_data %>% 
  filter(soil_PC1 >= 1) %>% 
  select(plot_id) 

high_nutri_sp = ind_data %>% 
  filter(plot_id %in% high_nutri_plot$plot_id) %>% 
  group_by(sp) %>% 
  reframe(n = n()) %>% 
  arrange(desc(n))

### high nutri plots
low_nutri_plot = comm_data %>% 
  filter(soil_PC1 <= -3) %>% 
  select(plot_id) 

low_nutri_sp = ind_data %>% 
  filter(plot_id %in% low_nutri_plot$plot_id) %>% 
  group_by(sp) %>% 
  reframe(n = n()) %>% 
  arrange(desc(n))

################################# other plots #################################

ggplot(data = comm_data) +
  geom_point(aes(x = fire_frequency, 
                 y = density
                ),
                alpha = 0.20
  ) +
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
