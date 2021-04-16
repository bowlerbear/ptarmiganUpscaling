#script to analysis line transect data on the HPC
library(tidyverse)
library(sp)
library(rgeos)
library(raster)
library(maptools)

source('C:/Users/db40fysa/Dropbox/ptarmigan Upscaling/generalFunctions.R', encoding = 'UTF-8')

#specify top level folder
myfolder <- "Data" #on local PC
#myfolder <- "/data/idiv_ess/ptarmiganUpscaling" #HPC

### get norway ##############################################

equalM<-"+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

library(raster)
library(maptools)
library(rgeos)
data(wrld_simpl)
Norway <- subset(wrld_simpl,NAME=="Norway")
NorwayOrig <- Norway
Norway <- gBuffer(Norway,width=1)
Norway <- spTransform(Norway,crs(equalM))
NorwayOrig <- spTransform(NorwayOrig,crs(equalM))
plot(Norway)

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

### get environ data ##############################################

#get mean across all years
siteMeans <- tlDF %>%
  group_by(LinjeID) %>%
  summarise(meanNu = mean(totalsInfo,na.rm=T),
            meanTL = mean(length[length!=0])) %>%
  filter(!is.na(meanNu)) %>%
  filter(!is.infinite(meanNu))

bufferData <- readRDS("data/varDF_allEnvironData_buffers_idiv.rds")

siteMeans <- merge(siteMeans,bufferData,by="LinjeID")
siteMeans <- subset(siteMeans,!is.na(tree_line))

siteMeans$adm <- mergeCounties(siteMeans$adm,further=TRUE)
table(siteMeans$adm)

### get spatial polygons ###########################################

Polys_spatial <- readRDS("data/Polys_spatial.rds")
Polys_spatial <- merge(Polys_spatial,siteMeans,by="LinjeID")
Polys_spatial <- subset(Polys_spatial,!is.na(meanNu))
Polys_spatial$meanDensity <- Polys_spatial$meanNu/Polys_spatial$meanTL
head(Polys_spatial)
summary(Polys_spatial)

library(tmap)
tm_shape(NorwayOrig) +
  tm_borders()+
tm_shape(Polys_spatial)+
  tm_fill(col="adm")

tm_shape(NorwayOrig) +
  tm_fill("grey")+
  tm_shape(Polys_spatial)+
  tm_fill(col="meanDensity",style="pretty",n=7)+
  tm_legend(legend.position=c("left","top"))

### glm model ######################################################

hist(siteMeans$meanNu)
summary(siteMeans$meanNu)
hist(log(siteMeans$meanNu+1))

glm1 <- lm(log(meanNu+1) ~ scale(bio5) + scale(bio6) + 
            Forest + Bog + ODF + OSF + Mire + SnowBeds + Human + 
            scale(y) +
            scale(distCoast) +
            scale(tree_line) + 
            scale(elevation),
            offset=log(meanTL),
            data=siteMeans)
summary(glm1)

library(MuMIn)
options(na.action = "na.fail")
dd <- dredge(glm1)
subset(dd, delta < 2)
#bog, forest, mire, ODF, OSF,bio5, bio6, distCoast,treeline, y

### brt #########################################################################

#Boosted regression tree
library(dismo)
library(gbm)

siteMeans$meanNu <- log(siteMeans$meanNu/siteMeans$meanTL+1)

brt1 <- gbm.step(data=siteMeans, 
                 gbm.x = c(4:16,18:20,26,27), 
                 gbm.y = 2,
                 family = 'gaussian')

summary(brt1)
#                                     var    rel.inf
# y                                     y 26.2817454#nonlinear
# bio6                               bio6 17.2973079
# distCoast                     distCoast 12.2791788
# bio5                               bio5 11.5308334
# tree_line                     tree_line  8.5324660
# elevation                     elevation  4.1847467
# Bog                                 Bog  3.7565941
# Meadows                         Meadows  3.0434379
# Forest                           Forest  2.4127264

#plot main effects
gbm.plot(brt1, n.plots=12, write.title = TRUE)

#interactions?
find.int <- gbm.interactions(brt1)
find.int$interactions#none!
find.int$rank.list

### end ################################################################
