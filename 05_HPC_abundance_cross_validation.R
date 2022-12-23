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

#subset to years of interest - 2008 onwards
allData <- subset(allData,Year>2007 & Year<2018)

#remove hyphens for help with subsetting
allData$Fylkesnavn <- gsub("-"," ",allData$Fylkesnavn)
allData$Fylkesnavn[which(allData$Rapporteringsniva=="Indre Troms")] <- "Troms"

#mistake with 1405 - transect length
allData$LengdeTaksert[which(allData$LinjeID==1405&allData$LengdeTaksert==1100)] <- 11000

### remove outliers (see section below) #############################

#LinjeID 1925 has twice as high counts as all others
allData <- subset(allData, LinjeID!=1925)

#remove LinjeID 131?? -only 503 m long - smallest transect

#drop lines visited in less than 5 years - see below
allData <- subset(allData, !LinjeID %in% 
                    c(935,874,876,882,884,936,2317,2328,2338,878,886,1250,1569,2331,2339))

### aggregate data to the lines ######################################

#Get statistics per year and line
tlDF <- allData %>%
  dplyr::group_by(LinjeID,Year) %>%
  dplyr::summarise(nuGroups = length(totalIndiv[!is.na(totalIndiv)]),
                   totalsInfo = sum(totalIndiv,na.rm=T),
                   groupSize = mean(totalIndiv,na.rm=T),
                   length = mean(LengdeTaksert,na.rm=T))
sum(tlDF$totalsInfo,na.rm=T)

#insert NA when there is no transect but evidence of a survey
tlDF$length[is.na(tlDF$length)] <- 0
tlDF$nuGroups[tlDF$length==0 ] <- NA
tlDF$totalsInfo[tlDF$length==0] <- NA
tlDF$groupSize[tlDF$length==0] <- NA
summary(tlDF)
sum(tlDF$length==0)

### get environ data #################################################

bufferData <- readRDS(paste(myfolder,
                            "varDF_allEnvironData_buffers_idiv.rds",sep="/"))
bufferData <- subset(bufferData, !LinjeID %in% 
                       c(935,874,876,882,884,936,2317,2328,2338,878,886,1250,1569,2331,2339,1925))

tlDF <- subset(tlDF, LinjeID %in% bufferData$LinjeID)

siteInfo_ArtsDaten <- readRDS(paste(myfolder,
                                    "siteInfo_ArtsDaten.rds",sep="/"))

### make siteInfo ######################################

tlDF$siteIndex <- as.numeric(as.factor(tlDF$LinjeID))
siteInfo <- unique(tlDF[,c("LinjeID","siteIndex")])
siteInfo <- arrange(siteInfo,siteIndex)
siteInfo$adm <- bufferData$adm[match(siteInfo$LinjeID,bufferData$LinjeID)]
siteInfo$admN <- as.numeric(as.factor(siteInfo$adm))

### choose model ##############################################

modelTaskID <- read.delim(paste(myfolder,"modelTaskID_distanceModel_CV.txt",sep="/"),as.is=T)

#get task id
task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

#get model for this task
mymodel <- modelTaskID$Model[which(modelTaskID$TaskID==task.id)]
mymodel

### folds ########################################################

folds <- readRDS(paste(myfolder,"folds_distanceModel_bands.rds",sep="/"))
siteInfo$fold <- folds$fold[match(siteInfo$LinjeID,folds$LinjeID)]

#select fold of this task
fold.id = modelTaskID$Fold[which(modelTaskID$TaskID==task.id)]

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

#check alignment with other datasets
all(row.names(groupInfo)==siteInfo_train$siteIndex)

### detection data ################################################

allDetections <- subset(allData, !is.na(totalIndiv) & totalIndiv!=0 
                        & LinjeID %in% siteInfo_train$LinjeID)
allDetections$yearIndex <- as.numeric(factor(allDetections$Year))
allDetections$siteIndex <- siteInfo$siteIndex[match(allDetections$LinjeID,
                                                    siteInfo$LinjeID)]

#same for testing dataset
groupSizes_test <- reshape2::acast(tlDF_test,siteIndex~Year,value.var="groupSize")
transectLengths_test <- reshape2::acast(tlDF_test,siteIndex~Year,value.var="length")
totalsInfo_test <- reshape2::acast(tlDF_test,siteIndex~Year,value.var="totalsInfo")

### group sizes ###################################################

#predict possible ESW for all transects - impute for mean value when it is missing
meanGS = apply(groupSizes,1,function(x)median(x[!is.na(x)]))
for(i in 1:nrow(groupSizes)){
  for(j in 1:ncol(groupSizes)){
    groupSizes[i,j] <- ifelse(is.na(groupSizes[i,j]),
                              meanGS[i],
                              groupSizes[i,j])
  }
}

meanGS_test = apply(groupSizes_test,1,function(x)median(x[!is.na(x)]))
for(i in 1:nrow(groupSizes_test)){
  for(j in 1:ncol(groupSizes_test)){
    groupSizes_test[i,j] <- ifelse(is.na(groupSizes_test[i,j]),
                              meanGS_test[i],
                              groupSizes_test[i,j])
  }
}

### make bugs objects ###########################################

bugs.data <- list(#For the state model
  nsite_train = length(unique(siteInfo_train$siteIndex)),
  nsite_test = length(unique(siteInfo_test$siteIndex)),
  nyear = length(unique(allData$Year)),
  nadm = length(unique(siteInfo$admN)),
  site = siteInfo$siteIndex,
  adm = siteInfo$admN,
  NuIndivs = totalsInfo,
  TransectLength_train = transectLengths,
  TransectLength_test = transectLengths_test,
  #For the distance model
  W = 200,
  ndetections_train = nrow(allDetections),
  y = allDetections$LinjeAvstand,
  ln_GroupSize = log(allDetections$totalIndiv+1),
  groupSizes_train = groupSizes,
  groupSizes_test = groupSizes_test,
  zeros.dist = rep(0,nrow(allDetections)))

names(bugs.data)

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
bugs.data$occDM_train <- model.matrix(~ bufferData_train$bio6 +
                                    bufferData_train$bio5 +
                                    bufferData_train$tree_line +
                                    I(bufferData_train$tree_line^2))[,-1]
  
#predictions to full grid
bugs.data$occDM_test <- model.matrix(~ bufferData_test$bio6 +
                                     bufferData_test$bio5 +
                                     bufferData_test$tree_line +
                                     I(bufferData_test$tree_line^2))[,-1]
  
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

params <- c("int.d","line.d.sd","year.d.sd","beta",
            "mid.expNuIndivs_train","mid.expNuIndivs_test",
            "mean.expNuIndivs_train","mean.expNuIndivs_test")

modelfile <- paste(myfolder,mymodel,sep="/")

n.iterations <- 50000
n.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")) 

Sys.time()

out1 <- jags(bugs.data, 
             inits=NULL, 
             params, 
             modelfile, 
             n.thin=50,
             n.chains=n.cores, 
             n.burnin=n.iterations/2,
             n.iter=n.iterations,
             parallel=T)

saveRDS(out1$summary,file="outSummary_linetransectModel_CV.rds")

Sys.time()
print("Done main model")

### summary annual ###################################

temp <- data.frame(out1$summary)
temp$Param <- row.names(temp)

#train
temp_train <- subset(temp, grepl("mid.expNuIndivs_train", temp$Param))$mean
data_train <- bugs.data$NuIndivs[,6]
mean(abs(data_train[!is.na(data_train)] - temp_train[!is.na(data_train)]))
cor.test(data_train[!is.na(data_train)],temp_train[!is.na(data_train)])

#test
temp_test <- subset(temp, grepl("mid.expNuIndivs_test", temp$Param))$mean
data_test <- totalsInfo_test[,6]
mean(abs(data_test[!is.na(data_test)] - temp_test[!is.na(data_test)]))
cor.test(data_test[!is.na(data_test)],temp_test[!is.na(data_test)])

print("Simple annual stats done now")

### summary mean #############################

temp_train <- subset(temp, grepl("mean.expNuIndivs_train", temp$Param))$mean
data_train <- as.numeric(apply(bugs.data$NuIndivs,1,mean,na.rm=T))
mean(abs(data_train[!is.na(data_train)] - temp_train[!is.na(data_train)]))
cor.test(data_train[!is.na(data_train)],temp_train[!is.na(data_train)])

#test
temp_test <- subset(temp, grepl("mean.expNuIndivs_test", temp$Param))$mean
data_test <- as.numeric(apply(totalsInfo_test,1,mean,na.rm=T))
mean(abs(data_test[!is.na(data_test)] - temp_test[!is.na(data_test)]))
cor.test(data_test[!is.na(data_test)],temp_test[!is.na(data_test)])

print("Simple mean stats done now")

### annual samples ###################################

library(ggmcmc)
ggd <- ggs(out1$samples)

#train
out1_dataset <- subset(ggd,grepl("mid.expNuIndivs_train",ggd$Parameter))
out1_dataset$index <- as.numeric(interaction(out1_dataset$Iteration,out1_dataset$Chain))

#get actual NuIndiv
totalsInfo_mid <- as.numeric(totalsInfo[,6])
useData <- !is.na(totalsInfo_mid)

#get difference between this value and the simulated values
mad_dataset <- as.numeric()
rmse_dataset <- as.numeric()
n.index <- max(out1_dataset$index)

for(i in 1:n.index){
  
  mad_dataset[i] <- mean(abs(totalsInfo_mid[useData] - 
                               out1_dataset$value[out1_dataset$index==i][useData]))
  
  rmse_dataset[i] <- sqrt(mean((totalsInfo_mid[useData] - 
                                  out1_dataset$value[out1_dataset$index==i][useData])^2))
  
}

summary(mad_dataset)
summary(rmse_dataset)

saveRDS(summary(mad_dataset),file=paste0("MAD_train_linetransectModel_CV_",task.id,".rds"))
saveRDS(summary(rmse_dataset),file=paste0("RMSE_train_linetransectModel_CV_",task.id,".rds"))

#test
out1_dataset <- subset(ggd,grepl("mid.expNuIndivs_test",ggd$Parameter))
out1_dataset$index <- as.numeric(interaction(out1_dataset$Iteration,out1_dataset$Chain))

#get actual NuIndiv
totalsInfo_mid <- as.numeric(totalsInfo_test[,6])
useData <- !is.na(totalsInfo_mid)

#get difference between this value and the simulated values
mad_dataset <- as.numeric()
rmse_dataset <- as.numeric()
n.index <- max(out1_dataset$index)

for(i in 1:n.index){
  mad_dataset[i] <- mean(abs(totalsInfo_mid[useData] - 
                               out1_dataset$value[out1_dataset$index==i][useData]))
  
  rmse_dataset[i] <- sqrt(mean((totalsInfo_mid[useData] - 
                                  out1_dataset$value[out1_dataset$index==i][useData])^2))
  
}

summary(mad_dataset)
summary(rmse_dataset)

saveRDS(summary(mad_dataset),file=paste0("MAD_test_linetransectModel_CV_",task.id,".rds"))
saveRDS(summary(rmse_dataset),file=paste0("RMSE_test_linetransectModel_CV_",task.id,".rds"))

### mean samples ###################################

library(ggmcmc)
ggd <- ggs(out1$samples)

#train
out1_dataset <- subset(ggd,grepl("mean.expNuIndivs_train",ggd$Parameter))
out1_dataset$index <- as.numeric(interaction(out1_dataset$Iteration,out1_dataset$Chain))

#get actual NuIndiv
totalsInfo_mid <- as.numeric(apply(totalsInfo,1,mean,na.rm=T))
useData <- !is.na(totalsInfo_mid)

#get difference between this value and the simulated values
mad_dataset <- as.numeric()
rmse_dataset <- as.numeric()
n.index <- max(out1_dataset$index)

for(i in 1:n.index){
  
  mad_dataset[i] <- mean(abs(totalsInfo_mid[useData] - 
                               out1_dataset$value[out1_dataset$index==i][useData]))
  
  rmse_dataset[i] <- sqrt(mean((totalsInfo_mid[useData] - 
                                  out1_dataset$value[out1_dataset$index==i][useData])^2))
  
}

summary(mad_dataset)
summary(rmse_dataset)

saveRDS(summary(mad_dataset),file=paste0("mMAD_train_linetransectModel_CV_",task.id,".rds"))
saveRDS(summary(rmse_dataset),file=paste0("mRMSE_train_linetransectModel_CV_",task.id,".rds"))

#test
out1_dataset <- subset(ggd,grepl("mean.expNuIndivs_test",ggd$Parameter))
out1_dataset$index <- as.numeric(interaction(out1_dataset$Iteration,out1_dataset$Chain))

#get actual NuIndiv
totalsInfo_mid <- as.numeric(apply(totalsInfo_test,1,mean,na.rm=T))
useData <- !is.na(totalsInfo_mid)

#get difference between this value and the simulated values
mad_dataset <- as.numeric()
rmse_dataset <- as.numeric()
n.index <- max(out1_dataset$index)

for(i in 1:n.index){
  mad_dataset[i] <- mean(abs(totalsInfo_mid[useData] - 
                               out1_dataset$value[out1_dataset$index==i][useData]))
  
  rmse_dataset[i] <- sqrt(mean((totalsInfo_mid[useData] - 
                                  out1_dataset$value[out1_dataset$index==i][useData])^2))
  
}

summary(mad_dataset)
summary(rmse_dataset)

saveRDS(summary(mad_dataset),file=paste0("mMAD_test_linetransectModel_CV_",task.id,".rds"))
saveRDS(summary(rmse_dataset),file=paste0("mRMSE_test_linetransectModel_CV_",task.id,".rds"))


print("End")

### end ########################################################################

