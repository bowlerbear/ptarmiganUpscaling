### occupancy models ############################################

#see Broms et al. 2006
#we want 3000 posterior samples for each parameter

### goodness of fit #############################################

#Bayesian p-values 
#Pearson's residuals - Tobler et al. 2015
#plot them against the residuals too
#aggregate detections per site
#se also MacKenzie and Bailey (2004).
#compare them across sites

### model selection ############################################

#If the goal is prediction, then one could avoid model
#selection by including all the covariates in the model and
#using regularization to restrict the parameter space
#through strong priors

#In a similar vein, we used weakly informative
#priors as a regularization on the model. The weakly
#informative priors helped with collinearity between
#covariates and separation related to indicator variables.

#WAIC??
#indicator variable
#selection
#Rushing

### out-of-sample ##############################################

#k-fold cross validation - use 5...

#Withholding one of the k∗ subsets, one fits the models to
#the rest of the data, and calculates a scoring rule to measure
#how close the models’ predictions are to the hold-out
#data. This process is repeated for each subset of data,
#then the scoring function is summed or averaged over
#the k∗ folds

#After the final model is chosen, it is refit using the entire
#data set so that inference is based on all the data

#https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13107
#https://cran.r-project.org/web/packages/blockCV/vignettes/BlockCV_for_SDM.html

library(blockCV)
#spatially or environmentally separated folds

# loading raster library
library(raster)
library(sf)
library(plyr)
library(ggthemes)
library(ggplot2)

#create grid
equalM<-"+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

library(raster)
library(maptools)
library(rgeos)
data(wrld_simpl)
Norway <- subset(wrld_simpl,NAME=="Norway")
NorwayOrig <- Norway
Norway <- gBuffer(Norway,width=1)
Norway <- spTransform(Norway,crs(equalM))
plot(Norway)

#create grid
newres = 5000#5 km grid
mygrid <- raster(extent(projectExtent(Norway,equalM)))
res(mygrid) <- newres 
mygrid[] <- 1:ncell(mygrid)
plot(mygrid)
gridTemp <- mygrid
plot(Norway,add=T)
myGridDF <- as.data.frame(mygrid,xy=T)

#### for occu model #########

siteInfo <- readRDS("data/siteInfo_ArtsDaten.rds")

# import raster data
#make each variable column a grid
myVars <- names(siteInfo)[c(14:26,32:34,40:42)]
myRasters <- list()
for(i in 1:length(myVars)){
  temp <- mygrid
  temp[] <- NA
  temp[siteInfo$grid] <- siteInfo[,myVars[i]]
  myRasters[i] <- temp
}
myRasters <- stack(myRasters)
awt<- raster::brick(myRasters)
projection(awt) <- equalM
plot(awt)

# make a SpatialPointsDataFrame object from data.frame
PA <- siteInfo[,c("x","y","species")]
PA$species[is.na(PA$species)] <- 0
pa_data <- st_as_sf(PA, coords = c("x", "y"), crs = st_crs(awt,asText=TRUE))
# see the first few rows
pa_data
st_crs(pa_data)==st_crs(awt)

# plot species data on the map
plot(awt[[1]]) # plot raster data
plot(pa_data[which(pa_data$species==1), ], pch = 16, col="red", add=TRUE) 
plot(pa_data[which(pa_data$species==0), ], pch = 16, col="blue", add=TRUE)
plot(pa_data)

# spatial blocking

sb <- spatialBlock(speciesData = pa_data,
                   species = "species",
                   rasterLayer = awt,
                   theRange = 75000, # 150 km works 
                   k = 5,
                   selection = "random")

sb <- spatialBlock(speciesData = pa_data, 
                    species = "species",
                    rasterLayer = awt,
                    rows = 25,
                    cols = 1,
                    k = 5,
                    selection = "systematic")


# buffering with presence-absence data
#sb <- buffering(speciesData= pa_data,
#                 species= "species",
#                 theRange= 70000,
#                 spDataType = "PA",
#                 progress = TRUE)

# # environmental clustering
# 
# eb <- envBlock(rasterLayer = awt,
#                speciesData = pa_data,
#                species = "species",
#                k = 5,
#                standardization = "standard", # rescale variables between 0 and 1
#                rasterBlock = TRUE,
#                numLimit = 50)

sac <- spatialAutoRange(rasterLayer = awt,
                        sampleNumber = 5000,
                        doParallel = TRUE,
                        showPlots = TRUE)

#interactive

# explore generated folds
foldExplorer(blocks = sb, 
             rasterLayer = awt, 
             speciesData = pa_data)

# explore the block size
rangeExplorer(rasterLayer = awt) # the only mandatory input

# add species data to add them on the map
rangeExplorer(rasterLayer = awt,
              speciesData = pa_data,
              species = "species",
              rangeTable = NULL,
              minRange = 30000, # limit the search domain
              maxRange = 100000)

### retrieve each fold in the data
mydata <- raster::extract(awt, pa_data,df=TRUE)
nrow(mydata) #sites with NA species data are excluded
nrow(siteInfo)
folds <- sb$foldID
length(folds)

#plot the folds
siteInfo$folds <- factor(sb$foldID)
g1 <- qplot(x, y, data=siteInfo, colour=folds) + theme_void()+
  theme(legend.position = "top")
g1

#random
saveRDS(siteInfo[,c("grid","siteIndex","folds")],
        file="data/folds_occModel.rds")

#systematic bands
saveRDS(siteInfo[,c("grid","siteIndex","folds")],
        file="data/folds_occModel_bands.rds")

#how many are in each fold

#### for distance model ####

library(blockCV)

siteInfo <- readRDS("data/siteInfo_ArtsDaten.rds")

# import raster data for whole country
myVars <- names(siteInfo)[c(14:26,32:34,40:42)]
myRasters <- list()
for(i in 1:length(myVars)){
  temp <- mygrid
  temp[] <- NA
  temp[siteInfo$grid] <- siteInfo[,myVars[i]]
  myRasters[i] <- temp
}
myRasters <- stack(myRasters)
awt<- raster::brick(myRasters)
projection(awt) <- equalM
plot(awt)

# get transect data
siteInfo <- readRDS("data/lines_to_grids.rds")

# import presence-absence species data
pa_data <- siteInfo
pa_data$species <- 1
coordinates(pa_data) <- c("x","y")
proj4string(pa_data) <- CRS(equalM)

# plot species data on the map
plot(awt[[1]]) # plot raster data
plot(pa_data,add=TRUE)

# spatial blocking
sb <- spatialBlock(speciesData = pa_data,
                   species = "species",
                   rasterLayer = awt,
                   theRange = 100000,
                   #theRange = 150000, # size of the blocks
                   k = 5,
                   selection = "random")


sb <- spatialBlock(speciesData = pa_data, 
                   species = "species",
                   rasterLayer = awt,
                   rows = 25,
                   cols = 1,
                   k = 5,
                   selection = "systematic")

sac <- spatialAutoRange(rasterLayer = awt,
                        sampleNumber = 5000,
                        doParallel = TRUE,
                        showPlots = TRUE)

### retrieve each fold in the data
mydata <- raster::extract(awt, pa_data,df=TRUE)
nrow(mydata) #sites with NA species data are excluded
nrow(siteInfo)
folds <- sb$foldID
length(folds)

#plot the folds
siteInfo$folds <- sb$foldID
ggplot2::qplot(x, y, data=siteInfo, colour=factor(folds))

#random
saveRDS(siteInfo[,c("grid","siteIndex","LinjeID","folds")],
        file="data/folds_distanceModel.rds")

#bands
saveRDS(siteInfo[,c("LinjeID","grid","folds")],
        file="data/folds_distanceModel_bands.rds")

### AUC ##################################################

#choice for binary data models is AUC, area under
#the receiver operator characteristic curve

#AUC, Brier, Logarithmic, and 0–1 scoring rules
#https://www.r-bloggers.com/2013/07/occupancy-model-fit-auc/

#Now that we have our posteriors for $\psi$ at each site in the final
#year, we can fit a single-year model to the final year’s data to
#estimate $Z$. 

#To set up the data for AUC calculations, produce site by
#iteration arrays for $\psi$ and $Z$

#Now generate the posterior for AUC and store data on the true and false
#positive rates to produce ROC curves

#https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/11-1936.1
#https://rviews.rstudio.com/2019/03/01/some-r-packages-for-roc-curves/

#AUC is equal to the probability that a true positive is scored greater than a true negative
#In a ROC plot the true positive rate (sensitivity) is plotted against the false positive rate (1.0-specificity) as the threshold varies from 0 to 1.

library(ROCR)

#using point estimate
siteInfo_obsSites <- subset(siteInfo,!is.na(species))
auc.tmp <- performance(pred,"auc"); auc <- as.numeric(auc.tmp@y.values)

library(data.table)
library(mltools)
preds <- siteInfo_obsSites$preds
actuals <- siteInfo_obsSites$species
  
auc_roc(preds, actuals)  
df <- auc_roc(preds, actuals, returnDT=TRUE)
qplot(CumulativeFPR,CumulativeTPR,data=df)

pred <- prediction(preds, actuals)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE)

#using full posterior

#Rather than use average values to determine a single point estimate, we used the full posterior distribution (3000 draws) and the R package ROCR

#make data frame into a matrix - site by iterations matrix
#calculate auc for each iteration - need both z and psi for this bit
library(ggmcmc)
ggd2 <- ggs(out2)
nstore <- niter / nthin
psi.pred <- array(dim=c(nsite, nstore*nchain))
Z.est <- array(dim=c(nsite, nstore*nchain))
for (j in 1:nsite){
  indx <- which(ggd$Parameter == paste("psi[", j, ",", nyear, "]", sep=""))
  psi.pred[j, ] <- ggd$value[indx]
  indx <- which(ggd2$Parameter == paste("z[", j, "]", sep=""))
  Z.est[j,] <- ggd2$value[indx]
}

require(ROCR)
AUC1 <- rep(NA, nstore * nchain)
fpr <- array(dim=c(nstore * nchain, nsite + 1))
tpr <- array(dim=c(nstore * nchain, nsite + 1))
for (i in 1:(nstore * nchain)){
  psi.vals <- psi.pred[, i]
  Z.vals <- Z.est[, i]
  pred <- prediction(psi.vals, factor(Z.vals, levels=c("0", "1")))
  perf <- performance(pred, "auc")
  AUC1[i] <- perf@y.values[[1]]
  perf <- performance(pred, "tpr","fpr")
  fpr[i, ] <- perf@x.values[[1]]
  tpr[i, ] <- perf@y.values[[1]]
}
require(reshape2)
fprm <- melt(fpr, varnames=c("iter", "site"))
tprm <- melt(tpr, varnames = c("iter", "site"))
ROC1 <- data.frame(fpr = fprm$value,
                   tpr = tprm$value,
                   iter = rep(fprm$iter, 2))

#actually compares z and psi...

#https://rdrr.io/github/jdyen/occupancy/man/validate.html
library(occupancy)

occupancy#object just Jags objects
calculate_metrics
r2_calc
validate
https://rdrr.io/github/jdyen/occupancy/src/R/validate.R
#fitted is p_obs???


