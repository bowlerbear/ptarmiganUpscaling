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

# import raster data
#make each variable column a grid
myVars <- names(siteInfo)[c(12:21,25:27)]
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
PA <- siteInfo[!is.na(siteInfo$species),c("x","y","species")]
# make a SpatialPointsDataFrame object from data.frame
pa_data <- st_as_sf(PA, coords = c("x", "y"), crs = st_crs(awt,asText=TRUE))
# see the first few rows
pa_data
st_crs(pa_data)==st_crs(awt)
st_crs(pa_data)
st_crs(awt)
st_crs(pa_data) = 9001

# plot species data on the map
plot(awt[[1]]) # plot raster data
plot(pa_data[which(pa_data$species==1), ], pch = 16, col="red", add=TRUE) 
plot(pa_data[which(pa_data$species==0), ], pch = 16, col="blue", add=TRUE)
plot(pa_data)

# spatial blocking by specified range with random assignment
sb <- spatialBlock(speciesData = pa_data,
                   species = "species",
                   rasterLayer = awt,
                   theRange = 200000, # size of the blocks
                   k = 5,
                   selection = "random")

# spatial blocking by rows and columns with checkerboard assignment
sb2 <- spatialBlock(speciesData = pa_data, 
                    species = "species",
                    rasterLayer = awt,
                    rows = 5,
                    cols = 6,
                    k = 5,
                    selection = "systematic")

# spatial blocking by rows with systematic assignment
sb3 <- spatialBlock(speciesData = pa_data,
                    species = "species",
                    rasterLayer = awt,
                    rows = 6,
                    selection = "checkerboard",
                    biomod2Format = TRUE)

# adding points on spatialBlock plot
library(ggplot2)

sb$plots + geom_sf(data = pa_data, alpha = 0.5)

# environmental clustering
eb <- envBlock(rasterLayer = awt,
               speciesData = pa_data,
               species = "species",
               k = 5,
               standardization = "standard", # rescale variables between 0 and 1
               rasterBlock = FALSE,
               numLimit = 50)

sac <- spatialAutoRange(rasterLayer = awt,
                        sampleNumber = 5000,
                        doParallel = TRUE,
                        showPlots = TRUE)

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


#check that each fold contains data from all the adms

## apply

# extract the raster values for the species points as a dataframe
mydata <- raster::extract(awt, pb_data)
mydata <- as.data.frame(mydata)
# create a vector of 1 (for presence) and 0 (for background samples)
pb <- pb_data$Species

# extract the folds in spatialBlock object created 
# in the previous section (with presence-background data)
# the foldID only works for spatialBlock and envBlock folds
folds <- sb2$foldID

# create an empty vector to store the AUC of each fold
AUCs <- vector(mode = "numeric")
for(k in seq_len(5)){
  # extracting the training and testing indices
  # this way only works with foldID
  trainSet <- which(folds != k) # training set indices
  testSet <- which(folds == k) # testing set indices

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