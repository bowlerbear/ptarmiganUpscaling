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
siteInfo <- readRDS("data/siteInfo_ArtsDaten.rds")

# import raster data
#make each variable column a grid
myVars <- names(siteInfo)[c(12:21,25:31)]
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

# import presence-absence species data
speciesSummary <- ddply(siteInfo,.(grid,siteIndex),summarise,PA=max(y))
siteInfo$species <- speciesSummary$PA[match(siteInfo$grid,speciesSummary$grid)]
siteInfo <- merge(siteInfo,myGridDF,by.x="grid",by.y="layer")

# make a SpatialPointsDataFrame object from data.frame
PA <- siteInfo[!is.na(siteInfo$species),c("x","y.y","species")]
pa_data <- st_as_sf(PA, coords = c("x", "y.y"), crs = st_crs(awt,asText=TRUE))
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
                   theRange = 150000, # size of the blocks
                   k = 5,
                   selection = "random")

sb2 <- spatialBlock(speciesData = pa_data, 
                    species = "species",
                    rasterLayer = awt,
                    rows = 10,
                    cols = 10,
                    k = 5,
                    selection = "systematic")

# environmental clustering

eb <- envBlock(rasterLayer = awt,
               speciesData = pa_data,
               species = "species",
               k = 5,
               standardization = "standard", # rescale variables between 0 and 1
               rasterBlock = TRUE,
               numLimit = 50)

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
              species = "Species",
              rangeTable = NULL,
              minRange = 30000, # limit the search domain
              maxRange = 100000)

### retrieve each fold in the data
mydata <- raster::extract(awt, pa_data,df=TRUE)
nrow(mydata) #sites with NA species data are excluded
nrow(siteInfo)
folds <- eb$foldID
length(folds)

#plot the folds
siteInfo$folds <- sb$foldID
qplot(x, y.y, data=siteInfo, colour=factor(folds))
saveRDS(siteInfo[,c("grid","siteIndex","folds")],
        file="data/folds_occModel.rds")

#same order as in the siteInfo filer, as long as data with only presence/absence
#records are used
siteInfo$fold <- 6
siteInfo$fold[!is.na(siteInfo$species)] <- folds
table(siteInfo$fold,siteInfo$adm)

# this way only works with foldID
# for k = 1 to 5
k = 1
trainData = mydata[which(folds != k), ]
testData = mydata[which(folds == k), ]  
  
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