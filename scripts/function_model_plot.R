
model_plot = function(data, 
                      x, 
                      y,
                      model, 
                      relationship = "none", 
                      x_label = NULL, 
                      y_label = NULL){

intercept = coef(model)[["(Intercept)"]]
slope = coef(model)[[x]]

if(relationship == "linear"){
  model_func = function(x){intercept + slope*x }
}
if(relationship == "exponential"){
  model_func = function(x){exp(intercept) * exp(slope)^x }
}

if( relationship == "none" ) {
  
  myplot = ggplot(data = data) +
    
    geom_point(aes(x = .data[[x]], 
                   y = .data[[y]]
    ),
    size = 3,
    alpha = 0.5
    )+
    
    xlab(x_label) +
    ylab(y_label) +
    
    theme(panel.background=element_rect(fill="white"),
          panel.grid=element_line(colour=NULL),
          panel.border=element_rect(fill=NA,colour="black"),
          axis.title.x = element_text(size=12, face="bold"),
          axis.title.y = element_text(size=12, face="bold"),
          axis.text.x = element_text(size= 10, angle = 0),
          axis.text.y = element_text(size= 10, angle = 0),
          legend.position = "none"
    )
  
} else {

  myplot = ggplot(data = data) +

    geom_point(aes(x = .data[[x]],
                   y = .data[[y]]
              ),
              size = 3,
              alpha = 0.5
    )+


    geom_line(aes(x = .data[[x]],
                  y = model_func(.data[[x]])
              ),
              linewidth = 1,
              linetype = "solid"
    ) +

    xlab(x_label) +
    ylab(y_label) +

    theme(panel.background=element_rect(fill="white"),
          panel.grid=element_line(colour=NULL),
          panel.border=element_rect(fill=NA,colour="black"),
          axis.title.x = element_text(size=12, face="bold"),
          axis.title.y = element_text(size=12, face="bold"),
          axis.text.x = element_text(size= 10, angle = 0),
          axis.text.y = element_text(size= 10, angle = 0),
          legend.position = "none"
    )

}

### return my plot
return(myplot)
  
}
