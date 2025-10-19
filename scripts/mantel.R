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
