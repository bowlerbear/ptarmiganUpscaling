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

#subset to years of interest - 2008 onwards since many not visited in 2007
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


### plot data #######################################################

#get lines as a spatial object
# library(sf)
# library(tmap)
# 
# Lines_spatial <- readRDS("data/Lines_spatial.rds")
# Lines_spatial <- subset(Lines_spatial, LinjeID %in% allData$LinjeID)
# Lines_spatial <-st_as_sf(Lines_spatial)
# Lines_spatial <- st_transform(Lines_spatial, st_crs(NorwayOrigProj))
# 
# occ_tmap <- tm_shape(NorwayOrigProj) + 
#    tm_borders() +
#  tm_shape(Lines_spatial)+
#    tm_lines(col="skyblue4",lwd=2)
# occ_tmap
# 
# saveRDS(occ_tmap,"plots/transects.RDS")


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

#### outlier check ##################################################

#row/siteIndex 423 is an outlier/LinejeID 1925
# summaryData <- tlDF %>%
#     group_by(LinjeID) %>%
#     summarise(med = median(totalsInfo,na.rm=T),
#               medDensity = median(totalsInfo/length,na.rm=T),
#               transectlength=mean(length,na.rm=T),
#               nuObsYears = length(unique(Year[!totalsInfo==0 & !is.na(totalsInfo)])),
#               nuYears = length(unique(Year[!is.na(totalsInfo)])),
#               propObsYears = nuObsYears/nuYears) %>%
#     arrange(desc(med))
# 
# qplot(summaryData$transectlength,summaryData$medDensity)
# qplot(summaryData$nuYears,summaryData$medDensity)
# 
# subset(summaryData,propObsYears<0.3)
# 
# tlDF %>%
#   group_by(LinjeID) %>%
#   summarise(nuYears = sum(!is.na(totalsInfo))) %>%
#   arrange(nuYears) %>%
#   filter(nuYears <5)
# #15 line

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

### make arrays ###################################################

#cast into arrays
groupInfo <- reshape2::acast(tlDF,siteIndex~Year,value.var="nuGroups")
totalsInfo <- reshape2::acast(tlDF,siteIndex~Year,value.var="totalsInfo")
groupSizes <- reshape2::acast(tlDF,siteIndex~Year,value.var="groupSize")
transectLengths <- reshape2::acast(tlDF,siteIndex~Year,value.var="length")

#put NAs where transect length is zero
groupInfo[transectLengths==0] <- NA
totalsInfo[transectLengths==0] <- NA
groupSizes[groupSizes==0] <- NA
sum(as.numeric(totalsInfo),na.rm=T)

#where there is a NA for transect length - put the mean for the line
#just for imputation purposes
meanTL = apply(transectLengths,1,median)
for(i in 1:nrow(transectLengths)){
  for(j in 1:ncol(transectLengths)){
    transectLengths[i,j] <- ifelse(transectLengths[i,j]==0,
                                   meanTL[i],
                                   transectLengths[i,j])
  }
}

#check alignment with other datasets
all(row.names(groupInfo)==siteInfo$siteIndex)

### site abundances ##############################################

# siteSummary <- tlDF %>%
#                 filter(!is.na(totalsInfo)) %>%
#                 group_by(siteIndex) %>%
#                 summarise(nuZeros = sum(totalsInfo==0),meanCount = mean(totalsInfo))
# table(siteSummary$nuZeros)  
# summary(siteSummary$meanCount) #all above zero 


### detection data ################################################

allDetections <- subset(allData, !is.na(totalIndiv) & totalIndiv!=0 
                        & LinjeID %in% bufferData$LinjeID)
allDetections$yearIndex <- as.numeric(factor(allDetections$Year))
allDetections$siteIndex <- siteInfo$siteIndex[match(allDetections$LinjeID,
                                                    siteInfo$LinjeID)]
allDetections$admN <- siteInfo$admN[match(allDetections$LinjeID,
                                                    siteInfo$LinjeID)]

#add on admN index
siteInfo_ArtsDaten$admNgrouped <- siteInfo$admN[match(siteInfo_ArtsDaten$admGrouped,
                                                      siteInfo$adm)]
  
#use raw data as predictor in model on sigma
groupSizes[is.na(groupSizes)] <- median(groupSizes,na.rm=T)

### line-transect index ########################################

#remember: some lines are given the same siteIndex (when they overlap in the same grid)

#get mapping
siteIndex_linetransects <- readRDS(paste(myfolder,"siteIndex_linetransects.rds",sep="/"))
siteIndex_linetransects <- siteIndex_linetransects %>% ungroup() %>% filter(LinjeID %in% siteInfo$LinjeID)
siteIndex_linetransects$siteIndex_All <- as.numeric(as.factor(siteIndex_linetransects$siteIndex_All))
summary(siteIndex_linetransects$siteIndex_All)
#302 grids are sampled

#map to siteInfo
siteInfo$siteIndex_All <- siteIndex_linetransects$siteIndex_All[match(siteInfo$LinjeID,
                                                                      siteIndex_linetransects$LinjeID)]
summary(siteInfo$siteIndex_All)

#map to siteInfo_ArtsDaten (and adding indices to unsampled grids)
siteInfo_ArtsDaten$siteIndex_All <- siteIndex_linetransects$siteIndex_All[match(siteInfo_ArtsDaten$grid,siteIndex_linetransects$grid)]
summary(siteInfo_ArtsDaten$siteIndex_All)
siteInfo_ArtsDaten <- plyr::arrange(siteInfo_ArtsDaten,siteIndex_All)

#fill in unsampled ones with indices
nuMissing <- sum(is.na(siteInfo_ArtsDaten$siteIndex_All))
maxIndex <- max(siteIndex_linetransects$siteIndex_All)
siteInfo_ArtsDaten$siteIndex_All[is.na(siteInfo_ArtsDaten$siteIndex_All)] <- (maxIndex+1):(maxIndex + nuMissing)
summary(siteInfo_ArtsDaten$siteIndex_All)

### make bugs objects ###########################################

bugs.data <- list(#For the state model
  nsite = length(unique(siteInfo$siteIndex)),
  nsiteAll = length(unique(siteInfo_ArtsDaten$siteIndex_All)),
  nyear = length(unique(allData$Year)),
  nadm = length(unique(siteInfo$admN)),
  site = siteInfo$siteIndex,
  siteAll = siteInfo$siteIndex_All,
  adm = siteInfo$admN,
  pred.adm = siteInfo_ArtsDaten$admNgrouped,
  year = (1:length(unique(allData$Year))),
  NuIndivs = totalsInfo,
  TransectLength = transectLengths,
  #For the distance model
  W = 200,
  ndetections = nrow(allDetections),
  y = allDetections$LinjeAvstand,
  ln_GroupSize = log(allDetections$totalIndiv+1),
  GroupSizes =  groupSizes,
  detectionYear = allDetections$yearIndex,
  detectionSite = allDetections$siteIndex,
  detectionAdm = allDetections$admN,
  zeros.dist = rep(0,nrow(allDetections)))

names(bugs.data)

bugs.data$GroupSizes[is.na(bugs.data$GroupSizes)] <- 0

### get environ data #########################################

all(bufferData$LinjeID==siteInfo$LinjeID)

myVars <- c("bio1", "bio5", "bio6","MountainBirchForest", "Bog","ODF",
            "Meadows","OSF","Mire","SnowBeds",
            "tree_line_position","tree_line",
            "y","distCoast")

# scale them
bufferData <- bufferData[,c("LinjeID",myVars)]
bufferMeans <- as.numeric(apply(bufferData,2,mean))
bufferSD <- as.numeric(apply(bufferData,2,sd))
for(i in 2:ncol(bufferData)){
  bufferData[,i] <- (bufferData[,i] - bufferMeans[i])/bufferSD[i]
}

#also for the siteInfo_ArtsDaten with same scaling
siteInfo_ArtsDaten <- siteInfo_ArtsDaten[,c("grid",myVars)]
for(i in 2:ncol(siteInfo_ArtsDaten)){
  siteInfo_ArtsDaten[,i] <- (siteInfo_ArtsDaten[,i] - bufferMeans[i])/bufferSD[i]
}

### choose model ##############################################

modelTaskID <- read.delim(paste(myfolder,"modelTaskID_distanceModel.txt",sep="/"),as.is=T)

#get task id
task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

#get model for this task
mymodel <- modelTaskID$Model[which(modelTaskID$TaskID==task.id)]

### standard model  ###########################################

#variables selected based on first simple analyses

if(mymodel == "linetransectModel_variables.txt"){
  
#add new variables to the bugs data
bugs.data$occDM <- model.matrix(~ bufferData$bio6 +
                                  bufferData$bio5 +
                                  bufferData$tree_line +
                                  I(bufferData$tree_line^2))[,-1]

#predictions to full grid
bugs.data$predDM <- model.matrix(~ siteInfo_ArtsDaten$bio6 +
                                   siteInfo_ArtsDaten$bio5 +
                                   siteInfo_ArtsDaten$tree_line +
                                   I(siteInfo_ArtsDaten$tree_line^2))[,-1]

### indicator model selection ################################

#all linear and select quadratic
} else { 
  
bugs.data$occDM <- model.matrix(~ bufferData$bio6 +
                                  bufferData$bio5 +
                                  bufferData$y +
                                  bufferData$distCoast +
                                  bufferData$tree_line +
                                  bufferData$MountainBirchForest +
                                  bufferData$Bog +
                                  bufferData$ODF +
                                  bufferData$Meadows +
                                  bufferData$OSF +
                                  bufferData$Mire +
                                  bufferData$SnowBeds +
                                  I(bufferData$bio6^2) +
                                  I(bufferData$bio5^2) +
                                  I(bufferData$y^2) +
                                  I(bufferData$distCoast^2) +
                                  I(bufferData$tree_line^2) +
                                  I(bufferData$MountainBirchForest^2) +
                                  I(bufferData$Bog^2) +
                                  I(bufferData$ODF^2) +
                                  I(bufferData$Meadows^2) +
                                  I(bufferData$OSF^2) +
                                  I(bufferData$Mire^2) +
                                  I(bufferData$SnowBeds^2))[,-1]

# #predictions to full grid
bugs.data$predDM <- model.matrix(~ siteInfo_ArtsDaten$bio6 +
                                  siteInfo_ArtsDaten$bio5 +
                                  siteInfo_ArtsDaten$y +
                                  siteInfo_ArtsDaten$distCoast +
                                  siteInfo_ArtsDaten$tree_line +
                                  siteInfo_ArtsDaten$MountainBirchForest +
                                  siteInfo_ArtsDaten$Bog +
                                  siteInfo_ArtsDaten$ODF +
                                  siteInfo_ArtsDaten$Meadows +
                                  siteInfo_ArtsDaten$OSF +
                                  siteInfo_ArtsDaten$Mire +
                                  siteInfo_ArtsDaten$SnowBeds +
                                  I(siteInfo_ArtsDaten$bio6^2) +
                                  I(siteInfo_ArtsDaten$bio5^2) +
                                  I(siteInfo_ArtsDaten$y^2) +
                                  I(siteInfo_ArtsDaten$distCoast^2) +
                                  I(siteInfo_ArtsDaten$tree_line^2) +
                                  I(siteInfo_ArtsDaten$MountainBirchForest^2) +
                                  I(siteInfo_ArtsDaten$Bog^2) +
                                  I(siteInfo_ArtsDaten$ODF^2) +
                                  I(siteInfo_ArtsDaten$Meadows^2) +
                                  I(siteInfo_ArtsDaten$OSF^2) +
                                  I(siteInfo_ArtsDaten$Mire^2) +
                                  I(siteInfo_ArtsDaten$SnowBeds^2))[,-1]

} 


# myvars <- c("y",'bio6',"bio6_2","distCoast","distCoast_2",
#             "bio5","bio5_2","tree_line","tree_line_2","OSF","SnowBeds")
# 
# myvars <- c("bio6","bio5","y","distCoast","tree_line","MountainBirchForest",
#             "Bog","ODF","Meadows","OSF","Mire","SnowBeds",
#             "bio6_2","bio5_2","y_2","distCoast_2","tree_line_2","MountainBirchForest_2",
#             "Bog_2","ODF_2","Meadows_2","OSF_2","Mire_2","SnowBeds_2")

bugs.data$n.covs <- ncol(bugs.data$occDM)
bugs.data$n.preds <- dim(bugs.data$predDM)[1]


#saveRDS(bugs.data, file="data/bugs.data_linetransects.rds")

### fit model #################################################

library(rjags)
library(jagsUI)

params <- c("int.d","beta","g",
            "b.group.size","meanESW","random.d.line",
            "meanDensity","Density.p",
            "bpv",
            "NuIndivs.j","NuIndivs.new.j","exp.j",
            "mid.expNuIndivs")

#choose model - already done above now
#modelfile <- paste(myfolder,"linetransectModel_variables.txt",sep="/")
#modelfile <- paste(myfolder,"linetransectModel_variables_LASSO.txt",sep="/")
#modelfile <- paste(myfolder,"linetransectModel_variables_ModelSelection.txt",sep="/")
modelfile <- paste(myfolder,mymodel,sep="/")

#n.cores = as.integer(Sys.getenv("NSLOTS", "1")) 
n.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")) 
#n.cores = 3

n.iterations = 8000

out1 <- jags(bugs.data, 
             inits=NULL, 
             params, 
             modelfile, 
             n.thin=10,
             n.chains=n.cores, 
             n.burnin=n.iterations/2,
             n.iter=n.iterations,
             parallel=T)

saveRDS(out1$summary,file=paste0("outSummary_linetransectModel_variables_",task.id,".rds"))

### fit #######################################################

library(ggmcmc)

ggd <- ggs(out1$samples)
out1_dataset <- subset(ggd,grepl("mid.expNuIndivs",ggd$Parameter))
out1_dataset$index <- as.numeric(interaction(out1_dataset$Iteration,out1_dataset$Chain))

#get actual NuIndiv
totalsInfo_mid <- bugs.data$NuIndivs[,6]

#get difference between this value and the simulated values
mad_dataset <- as.numeric()
rmse_dataset <- as.numeric()
n.index <- max(out1_dataset$index)

for(i in 1:n.index){
  mad_dataset[i] <- mean(abs(totalsInfo_mid[!is.na(totalsInfo_mid)] - 
                               out1_dataset$value[out1_dataset$index==i][!is.na(totalsInfo_mid)]))
  
  rmse_dataset[i] <- sqrt(mean((totalsInfo_mid[!is.na(totalsInfo_mid)] - 
                                  out1_dataset$value[out1_dataset$index==i][!is.na(totalsInfo_mid)])^2))
  
}

summary(mad_dataset)
summary(rmse_dataset)

saveRDS(summary(mad_dataset),file=paste0("MAD_linetransectModel_variables_",task.id,".rds"))
saveRDS(summary(rmse_dataset),file=paste0("RMSE_linetransectModel_variables_",task.id,".rds"))

### get site and year predictions ############################

out2 <- update(out1, parameters.to.save = c("Density.pt"),n.iter=2000)
ggd <- ggs(out2$samples)
saveRDS(ggd,file=paste0("Density.pt_linetransectModel_variables_",task.id,".rds"))
        
### end #######################################################
