library(tidyverse)

modelPredictions <- readRDS("data/modelPredictions.rds")

#total nationwide abundance is just:
sum(modelPredictions$occ_mean * modelPredictions$density_mean)#yep, about right
#1168194

### occupancy ####

resampleOcc <- function(modelPredictions,focalGrid){
  
  nonfocalDensity <- modelPredictions$occ_mean[-focalGrid]*modelPredictions$density_mean[-focalGrid]
  focalDensity <- modelPredictions$density_mean[focalGrid]
    
  totalNu <- as.numeric()
  for(i in 1:1000){
    newOccSample <- rnorm(1,modelPredictions$occ_mean[focalGrid],modelPredictions$occ_sd[focalGrid])
    totalNu[i] <- sum(nonfocalDensity,newOccSample*focalDensity)
  }
  
  data.frame(grid=focalGrid,meanPop=mean(totalNu),medianPop=median(totalNu),sdPop=sd(totalNu))
}

#apply function to each grid
sensitivityOcc <- 1:nrow(modelPredictions) %>%
                     map_dfr(~resampleOcc(modelPredictions,focalGrid=.x))
sensitivityOcc$x <- varDF$x
sensitivityOcc$y <- varDF$y

qplot(x,y, data=sensitivityOcc, colour = meanPop)+ scale_colour_viridis_c()
qplot(x,y, data=sensitivityOcc, colour = sdPop)+ scale_colour_viridis_c()
qplot(x,y, data=modelPredictions, colour = occ_sd)+ scale_colour_viridis_c()


mygrid[] <- NA
mygrid[varDF$grid] <- sensitivityOcc$sdPop
crs(mygrid) <- equalM
t1 <- tm_shape(mygrid)+ tm_raster(title="Occupancy SD",palette="Spectral",
                                  style="cont", breaks=c(0,50,75,100,125,150))


modelPredictions$sdPop <- sensitivityOcc$sdPop
qplot(occ_mean, occ_sd, data=modelPredictions,colour=log(sdPop))+ 
  scale_colour_viridis_c()
qplot(occ_mean, occ_sd, data=subset(modelPredictions,log(sdPop)>3.9),colour=log(sdPop))+ 
  scale_colour_viridis_c()

### abundance model ####

resampleDensity <- function(modelPredictions,focalGrid){
  
  nonfocalDensity <- modelPredictions$occ_mean[-focalGrid]*modelPredictions$density_mean[-focalGrid]
  focalOcc<- modelPredictions$occ_mean[focalGrid]
  
  totalNu <- as.numeric()
  for(i in 1:1000){
    newDensitySample <- rnorm(1,modelPredictions$density_mean[focalGrid],modelPredictions$density_sd[focalGrid])
    totalNu[i] <- sum(nonfocalDensity,newDensitySample*focalOcc)
  }
  
  data.frame(grid=focalGrid,meanPop=mean(totalNu),medianPop=median(totalNu),sdPop=sd(totalNu))
}

#apply function to each grid
sensitivityDensity <- 1:nrow(modelPredictions) %>%
  map_dfr(~resampleDensity(modelPredictions,focalGrid=.x))
sensitivityDensity$x <- varDF$x
sensitivityDensity$y <- varDF$y

qplot(x,y, data=sensitivityDensity, colour = meanPop)+ scale_colour_viridis_c()
qplot(x,y, data=sensitivityDensity, colour = sdPop)+ scale_colour_viridis_c()
qplot(x,y, data=modelPredictions, colour = density_sd)+ scale_colour_viridis_c()

mygrid[] <- NA
mygrid[varDF$grid] <- sensitivityDensity$sdPop
crs(mygrid) <- equalM
t2 <- tm_shape(mygrid)+ tm_raster(title="Density SD  ",palette="Spectral", 
                                  style="cont", breaks=c(0,50,75,100,125,150))

### combine ####

tmaptools::palette_explorer()
tmap_arrange(t1,t2,nrow=1)

median(sensitivityDensity$sdPop)
median(sensitivityOcc$sdPop)

temp <- tmap_arrange(t1,t2,nrow=1)
temp

tmap_save(temp, "plots/Fig_6.png",width = 6, height = 4)

### end #######################################

