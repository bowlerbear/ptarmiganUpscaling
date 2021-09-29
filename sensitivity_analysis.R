library(tidyverse)

### occupancy and abundance

#get grid info
varDF <- readRDS("data/varDF_allEnvironData_5km_idiv.rds")

#get model predictions
modelPredictions <- readRDS("data/modelPredictions.rds")
nrow(modelPredictions)==nrow(varDF)

#add x and y
modelPredictions$x <- varDF$x
modelPredictions$y <- varDF$y

#total nationwide abundance is just:
sum(modelPredictions$occ_mean * modelPredictions$density_mean)#yep, about right

#first check how much uncertainty in each prediction affects the overall size:

#occupancy model
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

#using tmap - run start of modelSummaries.R
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

#abundance model
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

#using tmap - run start of modelSummaries.R
mygrid[] <- NA
mygrid[varDF$grid] <- sensitivityDensity$sdPop
crs(mygrid) <- equalM
t2 <- tm_shape(mygrid)+ tm_raster(title="Density SD  ",palette="Spectral", 
                                  style="cont", breaks=c(0,50,75,100,125,150))

#check relationships with model predictions
modelPredictions$densitySens <- sensitivityDensity$sdPop
modelPredictions$occSens <- sensitivityOcc$sdPop
pairs(modelPredictions)

qplot(occ_mean,occSens, data=modelPredictions)
qplot(density_mean,densitySens, data=modelPredictions)
qplot(occ_sd,occSens, data=modelPredictions)
qplot(density_sd,densitySens, data=modelPredictions)
qplot(occSens,densitySens, data=modelPredictions)

#analyse the uncertainty as a function of land-use variables and density and occupancy estimates??
hist(modelPredictions$densitySens)
summary(lm(densitySens ~ density_mean, data = modelPredictions))#higher density, lower the sensitivity
hist(modelPredictions$occSens)
summary(lm(occSens ~ occ_mean, data = modelPredictions))#higher occupancy, lower the sensitivity

#combine the two
sensitivityDensity$both <- sensitivityDensity$sdPop * sensitivityOcc$sdPop
qplot(x,y, data=sensitivityDensity, colour = both)+ scale_colour_viridis_c()

#or by simulation
resampleBoth<- function(modelPredictions,focalGrid){
  
  nonfocalDensity <- modelPredictions$occ_mean[-focalGrid]*modelPredictions$density_mean[-focalGrid]
  totalNu <- as.numeric()
  
  for(i in 1:1000){
    newDensitySample <- rnorm(1,modelPredictions$density_mean[focalGrid],modelPredictions$density_sd[focalGrid])
    newOccSample <- rnorm(1,modelPredictions$occ_mean[focalGrid],modelPredictions$occ_sd[focalGrid])
    totalNu[i] <- sum(nonfocalDensity,newDensitySample*newOccSample)
  }
  
  data.frame(grid=focalGrid,meanPop=mean(totalNu),medianPop=median(totalNu),sdPop=sd(totalNu))
}

sensitivityBoth <- 1:nrow(modelPredictions) %>%
  map_dfr(~resampleBoth(modelPredictions,focalGrid=.x))
sensitivityBoth$x <- varDF$x
sensitivityBoth$y <- varDF$y
qplot(x,y, data=sensitivityBoth, colour = sdPop)+ scale_colour_viridis_c()

#using tmap - run start of modelSummaries.R
mygrid[] <- NA
mygrid[varDF$grid] <- sensitivityBoth$sdPop
crs(mygrid) <- equalM
t3 <- tm_shape(mygrid)+ tm_raster(title="Both SD     ",palette="Spectral", 
                                  style="cont", breaks=c(0,50,75,100,125,150))

#combine all
tmaptools::palette_explorer()
tmap_arrange(t1,t2,nrow=1)

median(sensitivityDensity$sdPop)
median(sensitivityOcc$sdPop)

### simulate all grids #######################################

#just abundance
totalNu <- as.numeric()

for(i in 1:1000){
  
  newDensitySample <- rnorm(nrow(modelPredictions),modelPredictions$density_mean,modelPredictions$density_sd)
  totalNu[i] <- sum(newDensitySample*modelPredictions$occ_mean)
  
}
data.frame(grid=focalGrid,meanPop=mean(totalNu),medianPop=median(totalNu),sdPop=sd(totalNu))


#just occupancy
totalNu <- as.numeric()

for(i in 1:1000){
  
  newOccSample <- rnorm(nrow(modelPredictions),modelPredictions$occ_mean,modelPredictions$occ_sd)
  totalNu[i] <- sum(modelPredictions$density_mean*newOccSample)
  
}
data.frame(grid=focalGrid,meanPop=mean(totalNu),medianPop=median(totalNu),sdPop=sd(totalNu))


#both density and occupancy
totalNu <- as.numeric()

for(i in 1:1000){
  
newDensitySample <- rnorm(nrow(modelPredictions),modelPredictions$density_mean,modelPredictions$density_sd)
newOccSample <- rnorm(nrow(modelPredictions),modelPredictions$occ_mean,modelPredictions$occ_sd)
totalNu[i] <- sum(newDensitySample*newOccSample)

}
data.frame(grid=focalGrid,meanPop=mean(totalNu),medianPop=median(totalNu),sdPop=sd(totalNu))


### occupancy only #############################################
 
tempDF$fitsL <- predict(glm1,type="link",newdata=tempDF)
tempDF$fitsL.se <- predict(glm1,type="link",newdata=tempDF,se.fit=TRUE)$se.fit

#sensitivity analysis
allPreds <- tempDF[,c("grid","x","y","nuVisits",
                      "fits","fits.se","fitsL","fitsL.se")]

#total predicted number of grid cells
sum(allPreds$fits)
sum(allPreds$fits.se)

#for each grid
#all current occupancies for all other grids
#take new samples for focus grid
#get new predicted total sample size

myGridSurveys <- function(allPreds,mygrid){
  
  currentFits <- allPreds$fits[allPreds$grid!=mygrid]
  
  totalNu <- as.numeric()
  require(boot)
  for(i in 1:100){
    newSample <- invlogit(rnorm(1,allPreds$fitsL[allPreds$grid==mygrid],
                                allPreds$fitsL.se[allPreds$grid==mygrid]))
    
    #new predicted total number of grids
    totalNu[i] <- sum(currentFits,newSample)
  }
  
  data.frame(grid=mygrid,mean(totalNu),sd(totalNu))
}

#apply function to each grid
sensitivityOutput <- ldply(allPreds$grid,function(x){
  myGridSurveys(allPreds,mygrid=x)
})

allPreds <- merge(allPreds,sensitivityOutput,by="grid")
ggplot(allPreds)+
  geom_point(aes(x,y,colour=sd.totalNu.),shape=15,size=rel(1))+
  scale_colour_gradient2(low="grey",mid="pink",high="purple",
                         midpoint=0.004)

qplot(fits,sd.totalNu.,data=allPreds)
qplot(sd.totalNu.,nuVisits,data=allPreds)

summary(lm(sqrt(sd.totalNu.) ~ scale(fits) + scale(nuVisits),data=allPreds))

#simulate change in binomial variance
#with different values of n

#variation only partly caused by number of visits
#pull out the predicted effect of another visit?
summary(lm(sqrt(sd.totalNu.) ~ scale(fits)*nuVisits,data=allPreds))

allPreds$fitsB <- ifelse(allPreds$fits < median(allPreds$fits),"L","H")
summary(lm(sqrt(sd.totalNu.) ~ fitsB + fitsB:nuVisits,data=allPreds))
#effect of visit more important in higher density areas


### binomial #########################################################

n <- 1:100

p1 <- 0.1
p2 <- 0.5
p3 <- 0.9

varBin <- function(n,p){
  n*p*(1-p)
}

varBin(n=n,p=p1)

df1 <- data.frame(n,varBin(n,p1))
df1$P <- p1
names(df1)<-c("n","var","p")

df2 <- data.frame(n,varBin(n,p2))
df2$P <- p2
names(df2)<-c("n","var","p")

df3 <- data.frame(n,varBin(n,p3))
df3$P <- p3
names(df3)<-c("n","var","p")

allDF <- rbind(df1,df2,df3)
qplot(n,var,data=allDF,facets=~factor(p))

# ###influence measures#######################################################################
# 
# #https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/influence.measures
# 
# #Bar Plot of Cook's distance to detect observations that strongly influence fitted values of the mode
# #cooks.distance
# 
# tempDF$cooks <- cooks.distance(glm1)
# ggplot(tempDF)+
#   geom_point(aes(x,y,colour=cooks),shape=15,size=rel(1))+
#   scale_colour_gradient(low="steelblue",high="red")
# 
# #DFBETA measures the difference in each parameter estimate with and without the influential point
# #dfbeta - for the overall mean
# temp <- dfbeta(glm1)
# 
# tempDF$dfbeta <- temp[,1]
# ggplot(tempDF)+
#   geom_point(aes(x,y,colour=abs(dfbeta)),shape=15,size=rel(1))+
#   scale_colour_gradient2(low="grey",mid="pink",high="purple")