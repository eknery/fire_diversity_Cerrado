############################## SPECIES RICHNESS ###############################

##### species richeness
hist(log(comm_data$richness) )

### species richness model
glm_rich1 = glm(
  data = s_comm_data,
  richness ~ s_fire_frequency, 
  family = poisson(link = "identity")
)

glm_rich2 = glm(
  data = s_comm_data,
  richness ~ s_fire_frequency + seasonal_precipitation, 
  family = poisson(link = "identity")
)

glm_rich3 = glm(
  data = s_comm_data,
  richness ~ s_fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2, 
  family = poisson(link = "identity")
)

glm_rich4 = glm(
  data = s_comm_data,
  richness ~ (s_fire_frequency + seasonal_precipitation + soil_PC1 + soil_PC2)/site, 
  family = poisson(link = "identity")
)

aic_vals = c(
  model1 =glm_rich1$aic,
  model2 =glm_rich2$aic,
  model3 =glm_rich3$aic,
  model4 =glm_rich4$aic
)
param_vals = c(
  length(coef(glm_rich1)),
  length(coef(glm_rich2)),
  length(coef(glm_rich3)),
  length(coef(glm_rich4))
)
anova(
  glm_rich1,
  glm_rich2,
  glm_rich3,
  glm_rich4, 
  test = "Chisq")

### model summary
summary(glm_rich1)
plot(glm_rich1)
shapiro.test(resid(glm_rich1))
