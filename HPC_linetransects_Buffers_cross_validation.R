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

folds <- readRDS(paste(myfolder,"folds_distanceModel.rds",sep="/"))
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

### get environ data #########################################

#scale it 
bufferData[,-1] <- plyr::numcolwise(scale)(bufferData[,-1])

#train dataset
bufferData_train <- subset(bufferData,LinjeID %in% siteInfo_train$LinjeID)
all(siteInfo_train$LinjeID==bufferData_train$LinjeID)

#add new variables to the bugs data
bugs.data$occDM_train <- model.matrix(~ bufferData_train$tree_line +
                                  bufferData_train$bio1 +
                                  bufferData_train$bio1^2 +
                                  bufferData_train$bio6 +
                                  bufferData_train$Meadows +   
                                  bufferData_train$elevation +
                                  bufferData_train$Bog+
                                  bufferData_train$Forest)[,-1]

bugs.data$n.covs <- ncol(bugs.data$occDM)

#test dataset
bufferData_test <- subset(bufferData,LinjeID %in% siteInfo_test$LinjeID)
all(siteInfo_test$LinjeID==bufferData_test$LinjeID)

#add new variables to the bugs data
bugs.data$occDM_test <- model.matrix(~ bufferData_test$tree_line +
                                        bufferData_test$bio1 +
                                        bufferData_test$bio1^2 +
                                        bufferData_test$bio6 +
                                        bufferData_test$Meadows +
                                        bufferData_test$elevation +
                                        bufferData_test$Bog +
                                        bufferData_test$Forest)[,-1]

### fit model #################################################

library(rjags)
library(jagsUI)

params <- c("int.d","line.d.sd","year.d.sd","beta")

modelfile <- paste(myfolder,"linetransectModel_variables_CV.txt",sep="/")

#add adm to model eventually

n.iterations <- 40000

out1 <- jags(bugs.data, 
             inits=NULL, 
             params, 
             modelfile, 
             n.thin=10,
             n.chains=3, 
             n.burnin=n.iterations/2,
             n.iter=n.iterations,
             parallel=T)

saveRDS(out1$summary,file="outSummary_linetransectModel_CV.rds")


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

