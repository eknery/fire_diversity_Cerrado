plot_glm = function(data, x, y, model, show_limits = F, show_model = F, x_label, y_label){
  ### predictor and response variables
  pred = as.numeric(data[[x]])
  resp = as.numeric(data[[y]])
  ### limits for bars
  y_min = paste0(y, "_min")
  minimal = data[[y_min]]
  y_max = paste0(y, "_max")
  maximum = data[[y_max]]
  ### get extreme values from predicor
  xlim = range(pred)
  xrange = xlim[2] - xlim[1]
  xweight = seq(xlim[1],xlim[2], xrange/100)
  newdata = as.data.frame(xweight)
  colnames(newdata) = "pred"
  ### get predicted values
  yweight = predict(
    object = model, 
    newdata = newdata, 
    type="response"
  )
  ### predicted data
  pdata = cbind(xweight, yweight)
  ### plot
  plot(
    x = data[[x]],
    y = data[[y]],
    pch = 21,
    col = "black",
    bg="gray",
    cex = 1.5,
    xlab = x_label,
    ylab = y_label
    )
  if(show_limits){
    segments(
      data[[x]], minimal, 
      data[[x]], maximum, 
      col = "gray50",
      lwd = 2
    )
  }
  if(show_model){
    lines(xweight, yweight, lty=1, lwd= 2.5)
  }
    
}


