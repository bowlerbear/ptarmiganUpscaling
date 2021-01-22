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

#check grid 57656
rp2<-subset(rp,LinjeID == 1506)
rsp2<-subset(rsp,layer %in% unique(rp2$layer))
plot(rsp2)
plot(rp2,add=T)
myObs <- subset(allDataObs,LinjeID==1506)
plot(myObs,add=T)
#they are in the next grid??

#maybe create 5 km grid polygons around each transect
#each transect is its own site then
#cluster transects that are nearby

### checking transect and obs overlap ###############

myLines <- unique(rp$LinjeID)

line <- myLines[1]

for(i in 1:length(myLines)){
  
  line <- myLines[i]
  
#draw transect
png(filename=paste0("obs/LinjeID",line,".png"))
rp2<-subset(rp,LinjeID == line)
rsp2<-subset(rsp,layer %in% unique(rp2$layer))
plot(rsp2)
plot(rp2,add=T)#seem to overlap 2 cells indeed

#get observations for this line
myObs <- subset(allDataObs,LinjeID==line)
plot(myObs,add=T)
dev.off()

}

#grid 56707 has 16 transects crossing it??
rp2<-subset(rp,layer == 56707)
rsp2<-subset(rsp,layer == 56707)
plot(rsp2)
plot(rp2,add=T)#seem to overlap 2 cells indeed

#grid 32112 has a super high density of individuals because of some survey area
myObs <- subset(allDataObs, grid == 32112)
rp2<-subset(rp,layer == 32112)
rsp2<-subset(rsp,layer == 32112)
plot(rsp2)
plot(rp2,add=T)
plot(myObs,add=T)

#61680 has a high density buth also survey Prop>0.005
myObs <- subset(allDataObs, grid == 61680)
rp2<-subset(rp,layer == 61680)
rsp2<-subset(rsp,layer == 61680)
plot(rsp2)
plot(rp2,add=T)
plot(myObs,add=T)

#57656 has a high density ???? only see one group here
myObs <- subset(allDataObs, grid == 57656)
rp2<-subset(rp,layer == 57656)
rsp2<-subset(rsp,layer == 57656)
plot(rsp2)
plot(rp2,add=T)
plot(myObs,add=T)

#threshold 54823
myObs <- subset(allDataObs, grid == 51514)
rp2<-subset(rp,layer == 51514)
rsp2<-subset(rsp,layer == 51514)
plot(rsp2)
plot(rp2,add=T)
plot(myObs,add=T)

### fix issues #######################################

#issues with 96, 112, 131, 417, 484, 689, 714, 715, 724
#757, 783, 805, 808, 952, 953, 1458, 1496, 1502, 
# 1504, 1506, 1530, 1560, 1571, 1581, 1591, 1598, 1609
# 1640, 1667, 1690, 1772, 1776, 1825, 1883, 1918,
#1944, 2176, 2187, 2189, 2224, 2330, 2336, 2764, 2803, 2833,
#2838, 2843, 2850, 2869, 

## make a rule that if transect overlap is small and 
#no individuals ever seen assume grid not surveys

#fix the transect and move the dataobs
#loop through each rp
rp <- subset(rp, layer %in% focusGrids)
transectGrids <- unique(rp@data[,c("LinjeID","layer")])

#get all observations where there is an observation in a grid,
# not covered by a transect grid for that line
x <- subset(allDataObs,LinjeID==2263)

allDataObs <- plyr::ddply(allDataObs@data,"LinjeID",function(x){
  
  surveysGrids <- unique(x$grid)
  
  #check whether they all match with true grids for that transects
  checkMatch <- !all(surveysGrids %in% transectGrids$layer[transectGrids$LinjeID %in% x$LinjeID])

  #get number of true grids
  grids <- transectGrids$layer[transectGrids$LinjeID %in% x$LinjeID]
    
  if(checkMatch){
  #if there is only one true grid, assign all obs at that
    if(length(grids)==1){
      
      x$grid  <- grids
      
    }else if(length(grids)>1){ #if there is more than one possible grid, assign obs to grid with more obs
      
      Summary <- plyr::ddply(x,"grid",summarise,Sum = sum(totalIndiv))
      Summary <- subset(Summary, grid %in% grids)
      Summary <- arrange(Summary,desc(Sum))
      maxGrid <- Summary$grid[1]
      
      #which observations outside of the grid to these
      outsideGrid <- surveysGrids[!surveysGrids %in% grids] 
      x$grid[x$grid %in% outsideGrid] <- maxGrid
    }
    
    return(x)
  
  }
  
  else{
    
    return(x)
  }
    
})

#some grid have multiple line transect surveys
x <- subset(allDataObs,LinjeID==2263)
unique(x$grid)#should be two

### merging TL ###########################################

#get information on each year in which each grid was sampled
allSurveys <- unique(allData[!is.na(allData$TakseringID),c("Year","LinjeID")])
rp <- merge(allSurveys,rp@data,by=c("LinjeID"),all.x=T)

#expand to all years (NAs when there was no transect)
newgrid <- expand.grid(Year=sort(unique(rp$Year)), layer=unique(rp$layer))
rp <- merge(newgrid,rp,by=c("Year","layer"),all.x=T)
rp$length[is.na(rp$length)] <- 0

#get total length per grid cell per year
transectLengths <- rp %>% 
                  dplyr::group_by(layer,Year) %>%
                  dplyr::summarise(length=sum(length),nuTransects = length(unique(LinjeID,na.rm=T)))

#rename
tlDF <- transectLengths
names(tlDF)[1] <- "grid"

#which grid has the most transects
transectLengths <- arrange(transectLengths,desc(nuTransects))
head(transectLengths)
#56707 has the longest length covering it
#guess area surveyd for this grid
#length of transects is 24.762 km
#width is 106 *2 m
#area is 5.2495444 km 2

#area is whole grid is 25 km

#proportion of area surveyed is 5.2495444/25
#about 0.2

### aggregate data to the grid ######################################

#Get statistics per year and grid
obsDF <- allDataObs %>%
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

### check data ######################################################

tlDF$surveyArea <- tlDF$length/1000 * 2 * (106/1000)
tlDF$surveyProp <- tlDF$surveyArea/25
summary(tlDF$surveyProp)
# on average only 0.2 % of a grid is surveyed

qplot(surveyArea,totalsInfo,data=tlDF)
qplot(surveyArea,totalsInfo/surveyArea,data=tlDF)
qplot(surveyProp,totalsInfo/surveyArea,data=tlDF)
qplot(surveyProp,totalsInfo/surveyArea,data=subset(tlDF,surveyProp>0.05))
#smaller survey areas tend to see more individuals

#few extreme values to check
#tlDF$Density <- tlDF$totalsInfo/tlDF$surveyArea
#tlDF <- arrange(tlDF,desc(Density))
#grid Year   length nuTransects nuGroups totalsInfo groupSize  surveyArea   surveyProp  Density
#1 32112 2014  12.5664           1        1          2  2.000000 0.002664078 0.0001065631 750.7289
#2 51514 2008 138.5417           1        1         11 11.000000 0.029370835 0.0011748334 374.5212
#3 57656 2016 521.1351           1        4         37  9.250000 0.110480643 0.0044192257 334.9003

#ignore observations in a grid with only small area covered? place in other grid for these transects
#less than 5%
nrow(subset(tlDF,surveyProp<0.05))#3483
#look at threshold example
nrow(subset(tlDF,surveyProp<0.005))#1166

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
              summarise(meanNu = mean(totalsInfo,na.rm=T),
                        meanTL = mean(length),
                        meanSP = mean(surveyProp)) %>%
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
            scale(elevation^2),
            data=siteMeans)
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


#corrected by transect length
summary(siteMeans$meanNu/siteMeans$meanSP)
glm1 <- glm((meanNu+1)/meanSP ~ scale(bio1) + scale(bio5) + scale(bio6) + 
              Open + PrefOpen + PrefClosed +
              scale(tree_line_position) + 
              scale(elevation),
              data=siteMeans)
summary(glm1)

#look at predictions
siteInfo_ArtsDaten$preds <- exp(predict(glm1,newdata=siteInfo_ArtsDaten,type="response"))
mygrid[] <- NA
mygrid[siteInfo_ArtsDaten$grid] <- siteInfo_ArtsDaten$preds
plot(mygrid)

### gams ########################################################################

library(mgcv)

glm1 <- glm(log(meanNu+1) ~ scale(bio1) + scale(bio5) + scale(bio6) + 
              Open + PrefOpen + PrefClosed +
              scale(tree_line_position) + 
              scale(I(tree_line_position^2))+
              scale(elevation)+
              scale(elevation^2),data=siteMeans)
summary(glm1)

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
