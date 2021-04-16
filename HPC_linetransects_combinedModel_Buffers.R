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

#make the number smaller

siteInfo_ArtsDaten$distCoast <- siteInfo_ArtsDaten$distCoast/10000
siteInfo_ArtsDaten$elevation <- siteInfo_ArtsDaten$elevation/100  
siteInfo_ArtsDaten$bio1 <- siteInfo_ArtsDaten$bio1/10
siteInfo_ArtsDaten$bio5 <- siteInfo_ArtsDaten$bio5/100  
siteInfo_ArtsDaten$bio6 <- siteInfo_ArtsDaten$bio6/100    
siteInfo_ArtsDaten$tree_line <- siteInfo_ArtsDaten$tree_line/100
siteInfo_ArtsDaten$y <- siteInfo_ArtsDaten$y/1000000


bufferData$distCoast <- bufferData$distCoast/10000
bufferData$elevation <- bufferData$elevation/100  
bufferData$bio1 <- bufferData$bio1/10
bufferData$bio5 <- bufferData$bio5/100  
bufferData$bio6 <- bufferData$bio6/100    
bufferData$tree_line <- bufferData$tree_line/100
bufferData$y <- bufferData$y/1000000

#add new variables to the bugs data
all(bufferData$LinjeID==siteInfo$LinjeID)
bugs.data$occDM <- model.matrix(~ bufferData$y +
                                  bufferData$y^2 +#insig
                                  bufferData$bio6 +
                                  bufferData$distCoast +
                                  bufferData$bio5 +
                                  bufferData$tree_line +
                                  bufferData$elevation +#insig
                                  bufferData$Meadows +#insig
                                  bufferData$Bog +#insig
                                  bufferData$Forest)[,-1]
bugs.data$n.covs <- ncol(bugs.data$occDM)


#predictions to full grid
bugs.data$predDM <- model.matrix(~ siteInfo_ArtsDaten$y +
                                  siteInfo_ArtsDaten$y^2 +#insig
                                  siteInfo_ArtsDaten$bio6 +
                                  siteInfo_ArtsDaten$distCoast +
                                  siteInfo_ArtsDaten$bio5 +
                                  siteInfo_ArtsDaten$tree_line +
                                  siteInfo_ArtsDaten$elevation +#insig
                                  siteInfo_ArtsDaten$Meadows +#insig
                                  siteInfo_ArtsDaten$Bog +#insig
                                  siteInfo_ArtsDaten$Forest)[,-1]
bugs.data$n.preds <- dim(bugs.data$predDM)[1]
 

### fit model #################################################

library(rjags)
library(jagsUI)

#params <- c("int.d","line.d.sd","year.d.sd",
#            "beta","bpv","totalPop","Density","Dens_lt")

params <- c("int.d","beta","meanDensity","meanExpNu","bpv","expNuIndivs","Density")

modelfile <- paste(myfolder,"linetransectModel_variables.txt",sep="/")

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

#v 1
# bugs.data$occDM <- model.matrix(~ bufferData$tree_line +
#                                   bufferData$bio1 +
#                                   bufferData$bio5 +
#                                   bufferData$elevation +
#                                   bufferData$Bog +
#                                   bufferData$Forest)[,-1]



#v2
# bugs.data$occDM <- model.matrix(~ bufferData$tree_line +
#                                   bufferData$bio1 +
#                                   bufferData$bio1^2 +
#                                   bufferData$bio6 +
#                                   bufferData$Meadows +
#                                   bufferData$OSF +
#                                   bufferData$ODF +
#                                   bufferData$elevation +
#                                   bufferData$Bog +
#                                   bufferData$Forest)[,-1]

#v3
#as v2 except with random adm effect

#v4
#as V3 except with ESW model

#5
#as V3 except with random adm effect on the ESW model and state model

#6 
# bugs.data$occDM <- model.matrix(~ bufferData$y +
#                                   bufferData$y^2 +
#                                   bufferData$bio6 +
#                                   bufferData$distCoast +
#                                   bufferData$bio5 +
#                                   bufferData$tree_line +
#                                   bufferData$elevation +
#                                   bufferData$Meadows +
#                                   bufferData$Bog +
#                                   bufferData$Forest)[,-1]
