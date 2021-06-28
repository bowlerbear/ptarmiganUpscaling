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
  dplyr::summarise(nuGroups=length(totalIndiv[!is.na(totalIndiv)]),
                   totalsInfo=sum(totalIndiv,na.rm=T),
                   groupSize=mean(totalIndiv,na.rm=T),
                   length = mean(LengdeTaksert,na.rm=T))
sum(tlDF$totalsInfo,na.rm=T)
#67946

#insert NA when there is no transect but evidence of a survey
tlDF$length[is.na(tlDF$length)] <- 0
tlDF$nuGroups[tlDF$length==0 ] <- NA
tlDF$totalsInfo[tlDF$length==0] <- NA
tlDF$groupSize[tlDF$length==0] <- NA
summary(tlDF)
sum(tlDF$length==0)#1302

### get environ data #################################################

bufferData <- readRDS(paste(myfolder,
                            "varDF_allEnvironData_buffers_idiv.rds",sep="/"))
#why are some missing???
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
sum(as.numeric(totalsInfo),na.rm=T)
#67946

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

siteInfo_ArtsDaten$admNgrouped <- siteInfo$admN[match(siteInfo_ArtsDaten$admGrouped,
                                                      siteInfo$adm)]
  
### make bugs objects ###########################################

bugs.data <- list(#For the state model
  nsite = length(unique(siteInfo$siteIndex)),
  nyear = length(unique(allData$Year)),
  nadm = length(unique(siteInfo$admN)),
  site = siteInfo$siteIndex,
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
bugs.data$occDM <- model.matrix(~ bufferData$y +
                                  bufferData$bio6 +
                                  I(bufferData$bio6^2) +
                                  bufferData$distCoast +
                                  I(bufferData$distCoast^2) +
                                  bufferData$bio5 +
                                  I(bufferData$bio5^2) +
                                  bufferData$tree_line +
                                  I(bufferData$tree_line^2) +
                                  bufferData$OSF +
                                  bufferData$SnowBeds)[,-1]

#predictions to full grid
bugs.data$predDM <- model.matrix(~ siteInfo_ArtsDaten$y +
                                   siteInfo_ArtsDaten$bio6 +
                                   I(siteInfo_ArtsDaten$bio6^2) +
                                   siteInfo_ArtsDaten$distCoast +
                                   I(siteInfo_ArtsDaten$distCoast^2) +
                                   siteInfo_ArtsDaten$bio5 +
                                   I(siteInfo_ArtsDaten$bio5^2) +
                                   siteInfo_ArtsDaten$tree_line +
                                   I(siteInfo_ArtsDaten$tree_line^2) +
                                   siteInfo_ArtsDaten$OSF +
                                   siteInfo_ArtsDaten$SnowBeds)[,-1]

### indicator model selection ################################

#all linear and select quadratic
} else if (mymodel == "linetransectModel_variables_ModelSelection.txt"){ 
  
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

### lasso model ##############################################

} else if (mymodel == "linetransectModel_variables_LASSO.txt"){

#all linear and select quadratic
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

myvars <- c("bio6","bio5","y","distCoast","tree_line","MountainBirchForest",
            "Bog","ODF","Meadows","OSF","Mire","SnowBeds",
            "bio6_2","bio5_2","y_2","distCoast_2","tree_line_2","MountainBirchForest_2",
            "Bog_2","ODF_2","Meadows_2","OSF_2","Mire_2","SnowBeds_2")
bugs.data$n.covs <- ncol(bugs.data$occDM)
bugs.data$n.preds <- dim(bugs.data$predDM)[1]

myvars <- c("y",'bio6',"bio6_2","distCoast","distCoast_2",
            "bio5","bio5_2","tree_line","tree_line_2","OSF","SnowBeds")

### fit model #################################################

library(rjags)
library(jagsUI)

params <- c("int.d","beta","meanDensity","meanExpNu","fit","fit.new",
            "Density","g","expNuIndivs","exp","NuIndivs.new")

#choose model - already done above now
#modelfile <- paste(myfolder,"linetransectModel_variables.txt",sep="/")
#modelfile <- paste(myfolder,"linetransectModel_variables_LASSO.txt",sep="/")
#modelfile <- paste(myfolder,"linetransectModel_variables_ModelSelection.txt",sep="/")
modelfile <- paste(myfolder,mymodel,sep="/")

#n.cores = as.integer(Sys.getenv("NSLOTS", "1")) 
n.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")) 
#n.cores = 3

out1 <- jags(bugs.data, 
             inits=NULL, 
             params, 
             modelfile, 
             n.thin=10,
             n.chains=n.cores, 
             n.burnin=10000,
             n.iter=20000,
             parallel=T)

saveRDS(out1,file=paste0("out_linetransectModel_variables_",task.id,".rds"))
saveRDS(out1$summary,file=paste0("outSummary_linetransectModel_variables_",task.id,".rds"))

### end #######################################################

