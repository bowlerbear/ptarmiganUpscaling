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

### spatial line transects #######################################################

spatialLines <- readRDS(paste(myfolder,"spatialLineTransects.rds",sep="/"))

spatialLines <- subset(spatialLines,LinjeID%in%unique(allData$LinjeID))

LongLat = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

for (i in seq(nrow(spatialLines))) {
  if (i == 1) {
    spTemp = readWKT(spatialLines$STAsText[i], spatialLines$LinjeID[i], p4s=LongLat)
    
  }
  else {
    spTemp = rbind(
      spTemp, readWKT(spatialLines$STAsText[i], spatialLines$LinjeID[i], p4s=LongLat)
      
    )
  }
}

#Make Lines spatial data frame
mydata <- spatialLines[,c("LinjeID","Fylkesnavn","Region","Rapporteringsniva","OmradeID", "OmradeNavn")]
rownames(mydata) <- paste(mydata$LinjeID)
Lines_spatial <- SpatialLinesDataFrame(spTemp, mydata, match.ID=T)
#plot(Lines_spatial)

#### quality check ###########################################################

#get all observations
allDataObs <- subset(allData,!is.na(totalIndiv))
allDataObs <- subset(allDataObs,totalIndiv!=0)

#Remove obs without coords
allDataObs <- subset(allDataObs,!is.na(Latitude))
allDataObs <- subset(allDataObs,Latitude!=0)
summary(allDataObs$Latitude)
summary(allDataObs$Longitude)

#Remove obs whose coords differ from predicted coords by 500m
allDataObs$diff<-abs(allDataObs$LinjeAvstand-allDataObs$N_dist)
summary(allDataObs$diff)
sum(allDataObs$diff>500)#lose 818 observations..check with Erlend
allDataObs<-subset(allDataObs,diff<=500)

#make into a spatial object
coordinates(allDataObs)<-c("Longitude","Latitude")
#plot(allDataObs)

### common grid ###########################################################

#using a m grid
equalM<-"+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

#get Norway
data(wrld_simpl)
Norway <- subset(wrld_simpl,NAME=="Norway")
NorwayOrig <- Norway
Norway <- gBuffer(Norway,width=1)
Norway <- spTransform(Norway,crs(equalM))
NorwayOrigProj <- spTransform(NorwayOrig,crs(equalM))

#create grid
newres = 5000#5 km grid
mygrid<-raster(extent(projectExtent(Norway,equalM)))
res(mygrid) <- newres
mygrid[] <- 1:ncell(mygrid)
plot(mygrid)
gridTemp <- mygrid
myGridDF <- as.data.frame(mygrid,xy=T)

### define grids of interest ###############################################

focusGrids <- readRDS(paste(myfolder,"focusGrids.rds",sep="/"))
varDF <- readRDS(paste(myfolder,"varDF_allEnvironData_5km_idiv.rds",sep="/"))
focusGrids <- focusGrids[focusGrids %in% varDF$grid]

### map bird obs to grid ###################################################

#convert observations to this CRS as well
proj4string(allDataObs) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
allDataObs <- sp::spTransform(allDataObs, crs(equalM))
allDataObs$grid <- raster::extract(mygrid,allDataObs)
allDataObs <- subset(allDataObs, grid %in% focusGrids)

#plot
#plot(mygrid)
#plot(Norway,add=T)
#plot(allDataObs,add=T)

### map transects to grid ###############################################

#polygonise grid
projection(mygrid) <- CRS(equalM)
rsp <- rasterToPolygons(mygrid)

#convert lines into a m grid
Lines_spatial <-spTransform(Lines_spatial,CRS("+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

#get overlap
rp <- raster::intersect(Lines_spatial,rsp)
rp$length <- gLength(rp, byid=TRUE) 

head(rp)
#LinjeID           Fylkesnavn   Region Rapporteringsniva OmradeID OmradeNavn layer   length
#1     122 Oppland              Statskog         Kongsvoll      285    Gåvålia 56471 2158.865
#2     122 Oppland              Statskog         Kongsvoll      285    Gåvålia 56708 1521.135
#3      93 Oppland              Statskog         Kongsvoll      286  Kongsvoll 56707 3741.265
#4     103 Oppland              Statskog         Kongsvoll      286  Kongsvoll 56470 1808.584
#5     103 Oppland              Statskog         Kongsvoll      286  Kongsvoll 56471 1485.609
#6    1389 Oppland              Statskog         Kongsvoll      286  Kongsvoll 56470  928.865

#checking this works
# rp2<-subset(rp,LinjeID==122)
# rsp2<-subset(rsp,layer%in%c("56471","56708"))
# plot(rsp2)
# plot(rp2,add=T)#seem to overlap 2 cells indeed

#get information on each year in which each grid was sampled
allSurveys <- unique(allData[!is.na(allData$TakseringID),c("Year","LinjeID")])
rp <- merge(allSurveys,rp@data,by=c("LinjeID"),all.x=T)
newgrid <- expand.grid(Year=sort(unique(rp$Year)), layer=unique(rp$layer))
rp <- merge(newgrid,rp,by=c("Year","layer"),all.x=T)

#get total length per grid cell per year
transectLengths <- rp %>% 
                  dplyr::group_by(layer,Year) %>%
                  dplyr::summarise(length=sum(length,na.rm=T))

#cast transect lengths into an array
tlDF <- transectLengths
names(tlDF)[1] <- "grid"
transectLengths <- reshape2::acast(tlDF, grid~Year,value.var="length")
#transectLengths[1:10,1:10]

### aggregate data to the grid ######################################

#Get statistics per year and grid
obsDF <- allDataObs@data %>%
          dplyr::group_by(grid,Year) %>%
          dplyr::summarise(nuGroups=length(totalIndiv),totalsInfo=sum(totalIndiv),groupSize=mean(totalIndiv))
#all observations are above zero

sum(obsDF$totalsInfo,na.rm=T)
#[1] 62044

#add info to the transects data frame
tlDF <- merge(tlDF,obsDF,by=c("grid","Year"),all.x=T)

#insert zeros when there is missing data but evidence of a survey
tlDF$nuGroups[tlDF$length>0 & is.na(tlDF$nuGroups)] <- 0
tlDF$totalsInfo[tlDF$length>0 & is.na(tlDF$totalsInfo)] <- 0
tlDF$groupSize[tlDF$length>0 & is.na(tlDF$groupSize)] <- 0

### artdendaten bank site info ######################################

siteInfo_ArtsDaten <- readRDS(paste(myfolder,
                                    "siteInfo_ArtsDaten.rds",sep="/"))

### extend to whole grid ############################################

#add all the grids we wish to impute for
newgrid <- expand.grid(Year = unique(tlDF$Year),
                     grid = siteInfo_ArtsDaten$grid)
tlDF <- merge(newgrid,tlDF,by=c("Year","grid"),all.x=T)

#add zero transect length info
tlDF$length[is.na(tlDF$length)] <- 0

#add site index data
tlDF$siteIndex <- siteInfo_ArtsDaten$siteIndex[match(tlDF$grid,siteInfo_ArtsDaten$grid)]
tlDF <- arrange(tlDF,siteIndex,Year)

#plot where we have data
#dataPresent <- subset(myGridDF,layer%in%tlDF$grid[!is.na(tlDF$totalsInfo)])
#qplot(x,y,data=dataPresent)

### make siteInfo ##################################################

siteInfo <- unique(tlDF[,c("grid","siteIndex")])

#make indices the same as in the Artsdatenbank data
siteInfo$admN <- siteInfo_ArtsDaten$admN[match(siteInfo$grid,siteInfo_ArtsDaten$grid)]
siteInfo$admN2 <- siteInfo_ArtsDaten$admN2[match(siteInfo$grid,siteInfo_ArtsDaten$grid)]

#give a unique name
siteInfo_LineTransects<-siteInfo
saveRDS(siteInfo, "data/siteInfo_LineTransects.rds")

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
#61993

#check alignment with other datasets
all(row.names(groupInfo)==siteInfo_ArtsDaten$siteIndex)

### glm model ######################################################

#get mean across all years
siteMeans <- tlDF %>%
              group_by(grid,siteIndex) %>%
              summarise(meanNu = mean(totalsInfo,na.rm=T)) %>%
              filter(!is.na(meanNu)) %>%
              filter(!is.infinite(meanNu))

siteMeans <- merge(siteMeans,varDF,by="grid")
  
hist(siteMeans$meanNu)
summary(siteMeans$meanNu)

glm1 <- glm(log(meanNu+1) ~ scale(bio1) + scale(bio5) + scale(bio6) + 
            Open + PrefOpen + PrefClosed +
            scale(tree_line_position) + 
            scale(I(tree_line_position^2))+
            scale(elevation)+
            scale(elevation^2),data=siteMeans)
summary(glm1)

library(MuMIn)
options(na.action = "na.fail")
dd <- dredge(glm1)
subset(dd, delta < 2)


#look at predictions
siteInfo_ArtsDaten$preds <- exp(predict(glm1,newdata=siteInfo_ArtsDaten,type="response"))
mygrid[] <- NA
mygrid[siteInfo_ArtsDaten$grid] <- siteInfo_ArtsDaten$preds
plot(mygrid)

### brt #########################################################################

#Boosted regression tree
library(dismo)
library(gbm)

brt1 <- gbm.step(data=siteMeans, 
                 gbm.x = c(4:13,17:23), 
                 gbm.y = 3,
                 family = 'gaussian')

summary(brt1)
# var    rel.inf
# Top                               Top 31.3139235
# tree_line                   tree_line 20.0760114
# bio5                             bio5  7.6325236
# Open                             Open  7.6268486
# elevation                   elevation  7.5250045
# tree_line_position tree_line_position  6.4396742
# PrefOpen                     PrefOpen  6.1405622
# bio6                             bio6  4.7662491
# Forest                         Forest  2.6996189
# alpine_habitat2       alpine_habitat2  2.4619149
# bio1                             bio1  1.4815438
# PrefClosed                 PrefClosed  1.1766827
# alpine_habitat1       alpine_habitat1  0.3399723
# Bottom                         Bottom  0.3194701
# Agriculture               Agriculture  0.0000000
# alpine_habitat3       alpine_habitat3  0.0000000
# alpine_habitat4       alpine_habitat4  0.0000000

#plot main effects
gbm.plot(brt1, n.plots=11, write.title = TRUE)
gbm.plot.fits(brt1)
#non-linear plot for tree line position and bio1

#interactions?
find.int <- gbm.interactions(brt1)
find.int$interactions#none!
find.int$rank.list
### end ################################################################
