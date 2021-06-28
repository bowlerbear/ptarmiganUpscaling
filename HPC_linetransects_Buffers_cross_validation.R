#script to analysis line transect data on the HPC
library(tidyverse)
library(sp)
library(rgeos)
library(raster)
library(maptools)

#specify top level folder
#myfolder <- "Data" #on local PC
myfolder <- "/data/idiv_ess/ptarmiganUpscaling" #HPC

### ptarmigan data ###############################################################

#read in data frame
allData <- readRDS(paste(myfolder,"allData.rds",sep="/"))

#subset to years of interest
allData <- subset(allData,Year>2006 & Year<2018)

#remove hyphens for help with subsetting
allData$Fylkesnavn <- gsub("-"," ",allData$Fylkesnavn)
allData$Fylkesnavn[which(allData$Rapporteringsniva=="Indre Troms")] <- "Troms"

#mistake with 1405 - transect length
allData$LengdeTaksert[which(allData$LinjeID==1405&allData$LengdeTaksert==1100)] <- 11000

### aggregate data to the lines ######################################

#Get statistics per year and line
tlDF <- allData %>%
  dplyr::group_by(LinjeID,Year) %>%
  dplyr::summarise(nuGroups=length(totalIndiv[!is.na(totalIndiv) & totalIndiv!=0]),
                   totalsInfo=sum(totalIndiv,na.rm=T),
                   groupSize=mean(totalIndiv[!is.na(totalIndiv) & totalIndiv!=0]),
                   length = mean(LengdeTaksert,na.rm=T))
sum(tlDF$totalsInfo,na.rm=T)
#67946

#insert NA when there is no transect 
tlDF$length[is.na(tlDF$length)] <- 0
tlDF$nuGroups[tlDF$length==0] <- NA
tlDF$totalsInfo[tlDF$length==0] <- NA
tlDF$groupSize[tlDF$length==0] <- NA
summary(tlDF)
sum(tlDF$length==0)#1302

### get environ data #################################################

bufferData <- readRDS(paste(myfolder,"varDF_allEnvironData_buffers_idiv.rds",sep="/"))
tlDF <- subset(tlDF, LinjeID %in% bufferData$LinjeID)

### make siteInfo ######################################

tlDF$siteIndex <- as.numeric(as.factor(tlDF$LinjeID))
siteInfo <- unique(tlDF[,c("LinjeID","siteIndex")])
siteInfo <- arrange(siteInfo,siteIndex)
siteInfo$adm <- bufferData$adm[match(siteInfo$LinjeID,bufferData$LinjeID)]
siteInfo$admN <- as.numeric(as.factor(siteInfo$adm))

### folds ########################################################

folds <- readRDS(paste(myfolder,"folds_distanceModel_bands.rds",sep="/"))
siteInfo$fold <- folds$fold[match(siteInfo$LinjeID,folds$LinjeID)]

#select fold of this task
fold.id = as.integer(Sys.getenv("SGE_TASK_ID", "1"))

#split into test and train
siteInfo_test <- subset(siteInfo,fold == fold.id)
siteInfo_train<- subset(siteInfo,fold != fold.id)
tlDF_train <- subset(tlDF,siteIndex %in% siteInfo_train$siteIndex)
tlDF_test <- subset(tlDF,siteIndex %in% siteInfo_test$siteIndex)

### make arrays ###################################################

#cast into arrays
groupInfo <- reshape2::acast(tlDF_train,siteIndex~Year,value.var="nuGroups")
totalsInfo <- reshape2::acast(tlDF_train,siteIndex~Year,value.var="totalsInfo")
groupSizes <- reshape2::acast(tlDF_train,siteIndex~Year,value.var="groupSize")
transectLengths <- reshape2::acast(tlDF_train,siteIndex~Year,value.var="length")

#put NAs where transect length is zero
groupInfo[transectLengths==0] <- NA
totalsInfo[transectLengths==0] <- NA
groupSizes[groupSizes==0] <- NA
sum(as.numeric(totalsInfo),na.rm=T)

#check alignment with other datasets
all(row.names(groupInfo)==siteInfo_train$siteIndex)

### detection data ################################################

allDetections <- subset(allData, !is.na(totalIndiv) & totalIndiv!=0 
                        & LinjeID %in% siteInfo_train$LinjeID)
allDetections$yearIndex <- as.numeric(factor(allDetections$Year))
allDetections$siteIndex <- siteInfo$siteIndex[match(allDetections$LinjeID,
                                                    siteInfo$LinjeID)]

#use raw data as predictor in model on sigma
groupSizes[is.na(groupSizes)] <- median(groupSizes,na.rm=T)

#same for testing dataset
groupSizes_test <- reshape2::acast(tlDF_test,siteIndex~Year,value.var="groupSize")
groupSizes_test[is.na(groupSizes_test)] <- median(groupSizes_test,na.rm=T)
transectLengths_test <- reshape2::acast(tlDF_test,siteIndex~Year,value.var="length")
totalsInfo_test <- reshape2::acast(tlDF_test,siteIndex~Year,value.var="totalsInfo")

### make bugs objects ###########################################

bugs.data <- list(#For the state model
  nsite_train = length(unique(siteInfo_train$siteIndex)),
  nsite_test = length(unique(siteInfo_test$siteIndex)),
  nyear = length(unique(allData$Year)),
  nadm = length(unique(siteInfo$admN)),
  site = siteInfo$siteIndex,
  adm = siteInfo$admN,
  year = (1:length(unique(allData$Year))),
  NuIndivs = totalsInfo,
  TransectLength_train = transectLengths,
  TransectLength_test = transectLengths_test,
  #For the distance model
  W = 200,
  ndetections_train = nrow(allDetections),
  y = allDetections$LinjeAvstand,
  ln_GroupSize = log(allDetections$totalIndiv),
  groupSizes_train = groupSizes,
  groupSizes_test = groupSizes_test,
  zeros.dist = rep(0,nrow(allDetections)))

names(bugs.data)

### choose model ##############################################

modelTaskID <- read.delim(paste(myfolder,"modelTaskID_distanceModel_CV.txt",sep="/"),as.is=T)

#get task id
task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

#get model for this task
mymodel <- modelTaskID$Model[which(modelTaskID$TaskID==task.id)]

### scale vars ################################################

bufferData[,-1] <- plyr::numcolwise(scale)(bufferData[,-1])

#train dataset
bufferData_train <- subset(bufferData,LinjeID %in% siteInfo_train$LinjeID)
all(siteInfo_train$LinjeID==bufferData_train$LinjeID)

#test dataset
bufferData_test <- subset(bufferData,LinjeID %in% siteInfo_test$LinjeID)
all(siteInfo_test$LinjeID==bufferData_test$LinjeID)


### standard model  ###########################################

#variables selected based on first simple analyses

if(mymodel == "linetransectModel_variables_CV.txt"){
  
#add new variables to the bugs data
bugs.data$occDM_train <- model.matrix(~ bufferData_train$y +
                                    bufferData_train$bio6 +
                                    I(bufferData_train$bio6^2) +
                                    bufferData_train$distCoast +
                                    I(bufferData_train$distCoast^2) +
                                    bufferData_train$bio5 +
                                    I(bufferData_train$bio5^2) +
                                    bufferData_train$tree_line +
                                    I(bufferData_train$tree_line^2) +
                                    bufferData_train$OSF +
                                    bufferData_train$SnowBeds)[,-1]
  
#predictions to full grid
bugs.data$occDM_test <- model.matrix(~ bufferData_test$y +
                                     bufferData_test$bio6 +
                                     I(bufferData_test$bio6^2) +
                                     bufferData_test$distCoast +
                                     I(bufferData_test$distCoast^2) +
                                     bufferData_test$bio5 +
                                     I(bufferData_test$bio5^2) +
                                     bufferData_test$tree_line +
                                     I(bufferData_test$tree_line^2) +
                                     bufferData_test$OSF +
                                     bufferData_test$SnowBeds)[,-1]
  
### indicator model selection ################################
#or lasso
#all linear and select quadratic

} else { 
  
  bugs.data$occDM_train <- model.matrix(~ bufferData_train$bio6 +
                                    bufferData_train$bio5 +
                                    bufferData_train$y +
                                    bufferData_train$distCoast +
                                    bufferData_train$tree_line +
                                    bufferData_train$MountainBirchForest +
                                    bufferData_train$Bog +
                                    bufferData_train$ODF +
                                    bufferData_train$Meadows +
                                    bufferData_train$OSF +
                                    bufferData_train$Mire +
                                    bufferData_train$SnowBeds +
                                    I(bufferData_train$bio6^2) +
                                    I(bufferData_train$bio5^2) +
                                    I(bufferData_train$y^2) +
                                    I(bufferData_train$distCoast^2) +
                                    I(bufferData_train$tree_line^2) +
                                    I(bufferData_train$MountainBirchForest^2) +
                                    I(bufferData_train$Bog^2) +
                                    I(bufferData_train$ODF^2) +
                                    I(bufferData_train$Meadows^2) +
                                    I(bufferData_train$OSF^2) +
                                    I(bufferData_train$Mire^2) +
                                    I(bufferData_train$SnowBeds^2))[,-1]
  
  # #predictions to full grid
  bugs.data$occDM_test <- model.matrix(~ bufferData_test$bio6 +
                                     bufferData_test$bio5 +
                                     bufferData_test$y +
                                     bufferData_test$distCoast +
                                     bufferData_test$tree_line +
                                     bufferData_test$MountainBirchForest +
                                     bufferData_test$Bog +
                                     bufferData_test$ODF +
                                     bufferData_test$Meadows +
                                     bufferData_test$OSF +
                                     bufferData_test$Mire +
                                     bufferData_test$SnowBeds +
                                     I(bufferData_test$bio6^2) +
                                     I(bufferData_test$bio5^2) +
                                     I(bufferData_test$y^2) +
                                     I(bufferData_test$distCoast^2) +
                                     I(bufferData_test$tree_line^2) +
                                     I(bufferData_test$MountainBirchForest^2) +
                                     I(bufferData_test$Bog^2) +
                                     I(bufferData_test$ODF^2) +
                                     I(bufferData_test$Meadows^2) +
                                     I(bufferData_test$OSF^2) +
                                     I(bufferData_test$Mire^2) +
                                     I(bufferData_test$SnowBeds^2))[,-1]
} 

bugs.data$n.covs <- ncol(bugs.data$occDM_train)

### fit model #################################################

library(rjags)
library(jagsUI)

params <- c("int.d","line.d.sd","year.d.sd","beta")

modelfile <- paste(myfolder,mymodel,sep="/")

n.iterations <- 40000
n.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")) 

out1 <- jags(bugs.data, 
             inits=NULL, 
             params, 
             modelfile, 
             n.thin=10,
             n.chains=n.cores, 
             n.burnin=n.iterations/2,
             n.iter=n.iterations,
             parallel=T)

#saveRDS(out1$summary,file="outSummary_linetransectModel_CV.rds")

#update to extract cross validation comparisons
out2 <- update(out1,
               parameters.to.save = c("mean.expNuIndivs_train","mean.expNuIndivs_test",
                                      "mean.Density_train","mean.Density_test"),
               n.iter = 1000)

saveRDS(out2,file=paste0("out_update_linetransectModel_CV_",fold.id,".rds"))

#and save test and training dataset
saveRDS(totalsInfo,file=paste0("totalsInfo_train_CV_",fold.id,".rds"))
saveRDS(totalsInfo_test,file=paste0("totalsInfo_test_CV_",fold.id,".rds"))

### end ########################################################################

