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

bufferData <- readRDS(paste(myfolder,"varDF_allEnvironData_buffers_idiv.rds",sep="/"))
#why are some missing???
tlDF <- subset(tlDF, LinjeID %in% bufferData$LinjeID)

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

### detection data ################################################

allDetections <- subset(allData, !is.na(totalIndiv) & totalIndiv!=0)
allDetections$yearIndex <- as.numeric(factor(allDetections$Year))
allDetections$siteIndex <- siteInfo$siteIndex[match(allDetections$LinjeID,
                                                    siteInfo$LinjeID)]

### make bugs objects ###########################################

bugs.data <- list(#For the state model
  nsite = length(unique(siteInfo$siteIndex)),
  nyear = length(unique(allData$Year)),
  nadm = length(unique(siteInfo$admN)),
  site = siteInfo$siteIndex,
  adm = siteInfo$admN,
  year = (1:length(unique(allData$Year))),
  NuIndivs = totalsInfo,
  TransectLength = transectLengths,
  #For the distance model
  W = 200,
  ndetections = nrow(allDetections),
  y = allDetections$LinjeAvstand,
  ln_GroupSize = log(allDetections$totalIndiv+1),
  GroupSize = (allDetections$totalIndiv-1),#so it can start at 0
  detectionYear = allDetections$yearIndex,
  detectionSite = allDetections$siteIndex,
  zeros.dist = rep(0,nrow(allDetections)))

names(bugs.data)

### get environ data #########################################

#check everything aligns
all(bufferData$LinjeID==siteInfo$LinjeID)

#add new variables to the bugs data
bugs.data$occDM <- model.matrix(~ bufferData$tree_line +
                                  bufferData$bio1 +
                                  bufferData$bio5 +
                                  bufferData$elevation +
                                  bufferData$Bog +
                                  bufferData$Forest)[,-1]

bugs.data$n.covs <- ncol(bugs.data$occDM)

siteInfo_ArtsDaten <- readRDS(paste(myfolder,
                                    "siteInfo_ArtsDaten.rds",sep="/"))

environData <- siteInfo_ArtsDaten[,c(12:33)]
bugs.data$predDM <- model.matrix(~ environData$tree_line +
                                   environData$bio1 +
                                   environData$bio5 +
                                   environData$elevation + 
                                   environData$Bog +
                                   environData$Forest)[,-1]

bugs.data$npreds <- nrow(environData)

### fit model #################################################

library(rjags)
library(jagsUI)

#params <- c("int.d","line.d.sd","year.d.sd",
#            "beta","bpv","totalPop","Density","Dens_lt")

params <- c("Dens_lt","Density","meanExpNu","ExpNu_5")

modelfile <- paste(myfolder,"linetransectModel_variables.txt",sep="/")

#add adm to model eventually

out1 <- jags(bugs.data, 
             inits=NULL, 
             params, 
             modelfile, 
             n.thin=10,
             n.chains=3, 
             n.burnin=20000,
             n.iter=40000,
             parallel=T)

saveRDS(out1$summary,file="outSummary_linetransectModel_variables.rds")

### end #######################################################

