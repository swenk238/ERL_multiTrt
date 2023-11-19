## plotting DR with dual treatment


wd <- "~/GitHub/ERL_multiTrt"
setwd(wd)

load(file="ERL_dualtrt_DR.RData")


# plotting ----------------------------------------------------------------
## axis titles
axx <- list(title = "Specified MML days")
axy <- list(title = "Specified GGHE values")
axz <- list(title = "Median dose response values")

# Average dose response for all countries ---------------------------------------------------------
plot_ly(x = seqT1, y = seqT2, z = DR_mean, scene='scene1') %>% 
  add_surface(showscale=F, zmin=-0.2, zmax=0.2, 
              contours=list(z=list(show=TRUE,
                                   usecolormap=TRUE,
                                   highlightcolor="#ff0000",
                                   project=list(z=TRUE)))) %>% 
  layout(title = list(text = "Dose response curve for dual treatment", y = 0.9),
         scene = list(xaxis=axx,
                      yaxis=axy,
                      zaxis=list(title = "Average dose response values",
                                 range = c(-0.15,0.15))))

DR_list = list(seqT1, seqT2, cbind(DR_mean, DR_mean))
names(DR_list) = c('seqT1', 'seqT2', 'DR_mean')

plot(x=seqT1, y=DR_mean[,1])
plot(x=seqT2, y=DR_mean[,2])


