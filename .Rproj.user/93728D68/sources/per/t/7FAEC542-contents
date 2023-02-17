### input data
## env data
env_data = read.csv("0_data/env_data.csv", sep=",", h=T)

## spp data
spp_data = read.csv("0_data/spp_data.csv", sep=",", h=T, stringsAsFactors = T)
spp_data = spp_data[1:15]
head(spp_data)


##  create site_plot variable
# for env data
site_plot = paste(env_data$site, env_data$plot, sep="_")
env_data = data.frame(site_plot, env_data)
# for spp_data
site_plot = paste(spp_data$site, spp_data$plot, sep="_")
spp_data = data.frame(site_plot, spp_data)

unique(site_plot)

## abundance per plot spp df
abundace_plot = data.frame(table(spp_data$site_plot, spp_data$sp))
# remove unusal plot naming
abundace_plot = abundace_plot[abundace_plot$Var1 %in% env_data$site_plot,]

env_data$site_plot[!env_data$site_plot %in% abundace_plot$Var1]

unique(abundace_plot$Var1)

# add fire values 
fire = c()
for (i in 1:nrow(abundace_plot)){
  plot_index = which(abundace_plot$Var1[i] == env_data$site_plot)
  fire = c(fire, env_data$fire_frequency[plot_index] )
}

data.frame(abundace_plot, fire)
str(abundace_plot)

# iti vs ca