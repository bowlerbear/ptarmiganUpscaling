########################################################################
setwd("/data/dbowler/ptarmiganUpscaling")

#load in general functions

source('generalFunctions.R')

library(sp)
library(raster)
library(maptools)
library(ggplot2)
library(rgeos)
library(plyr)

#get norway##############################################################

#using a m grid
equalM<-"+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
data(wrld_simpl)
Norway <- subset(wrld_simpl,NAME=="Norway")
NorwayOrig <- Norway
Norway <- gBuffer(Norway,width=1)
Norway <- spTransform(Norway,crs(equalM))
NorwayOrigProj <- spTransform(NorwayOrig,crs(equalM))

#########################################################################

#create grid
newres=5000#5 km grid
mygrid<-raster(extent(projectExtent(Norway,equalM)))
res(mygrid) <- newres
mygrid[] <- 1:ncell(mygrid)
plot(mygrid)
gridTemp <- mygrid
myGridDF <- as.data.frame(mygrid,xy=T)

#########################################################################

#within this grid, define cells of interest:

# plot(mygrid)
# plot(NorwayOrigProj,add=T)
# 
# #get cells at least 50% covering Norway
# myGridMask <- data.frame(extract(mygrid,NorwayOrigProj,weights=T,normalizeWeights=F)[[1]])
# myGridMask <- myGridMask$value[myGridMask$weight>0.5]
# myGrid2 <- mygrid
# myGrid2[!getValues(mygrid) %in% myGridMask]<-NA
# plot(myGrid2)
# plot(NorwayOrigProj,add=T)
# #looks good!
# 
# #get these grid numbers
# focusGrids <- getValues(myGrid2)[!is.na(getValues(myGrid2))]
# focusGridsDF <- as.data.frame(myGrid2,xy=T)
focusGrids <- readRDS("focusGrids.rds")

#########################################################################

#get willow ptarmigan data - occurrence data

#get GBIF data
#tdir <- "C:/Users/db40fysa/Dropbox/Alpine/GBIF"
#lirype <- read.delim(paste(tdir,"willow_ptarmigan_GBIF.txt",sep="/"),as.is=T,row.names=NULL)
tdir <- "GBIF"
lirype <- read.delim(paste(tdir,"willow_ptarmigan_GBIF.txt",sep="/"),as.is=T,row.names=NULL,fileEncoding="latin1",quote = "")
names(lirype) <- gsub("X.","",names(lirype))
names(lirype) <- gsub("\\.","",names(lirype))

#clean
lirype$decimalLatitude<-as.numeric(lirype$decimalLatitude)
lirype$decimalLongitude<-as.numeric(lirype$decimalLongitude)
lirype<-subset(lirype,!is.na(lirype$decimalLatitude)|!is.na(lirype$decimalLongitude))

#which column has the dataset information
temp <- subset(lirype,samplingProtocol=="\"Distance sampling based on line transects\"")
unique(temp$datasetName)
unique(temp$ownerInstitutionCode)
table(temp$ownerInstitutionCode)
#see also datasetName
#Statskog and FeFo are both in the line transect surveys
#might need to remove these....???

lirype<-lirype[,c("decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters",
                  "name","year","month","day","recordedBy")]

#get artsdatenbanken data
#tdir<-"C:/Users/db40fysa/Dropbox/Alpine/SpeciesMapServices/Ptarmigan"
tdir<-"SpeciesMapServices/Ptarmigan"
lirype2 <- read.delim(paste(tdir,"DataExportFraArtskart.txt",sep="/"),skipNul=T,dec=",")
lirype2$Species<-apply(lirype2,1,function(x)paste(x["Genus"],x["Species"],sep=" "))

#check dataset names
unique(lirype2$InstitutionName)

lirype2<-lirype2[,c("Latitude","Longitude","CoordinatePrecision",
                    "Species","YearCollected","MonthCollected","DayCollected","Collector")]
names(lirype2)<-names(lirype)

#combine them
lirype<-rbind(lirype,lirype2)
lirypeD<-lirype

#remove duplicates
lirype$year <- as.numeric(lirype$year)
lirype$date<-paste(lirype$year,lirype$month,lirype$day,sep="-")
lirype$id<-paste(lirype$date,round(lirype$decimalLatitude,digits=4),round(lirype$decimalLongitude,digits=4),sep="_")
lirype<-subset(lirype,!duplicated(id))

#remove dodgy records
summary(lirype$coordinateUncertaintyInMeters)
lirype$coordinateUncertaintyInMeters <- as.numeric(lirype$coordinateUncertaintyInMeters)
lirype <- subset(lirype, coordinateUncertaintyInMeters<5000|is.na(coordinateUncertaintyInMeters))

#remove those without 2 decimal places
lirype$latDP<-getDecimalPlaces(lirype$decimalLatitude)
lirype$lonDP<-getDecimalPlaces(lirype$decimalLongitude)
lirype<-subset(lirype,latDP>2 & lonDP>2)

#examine recorders - choose those with a name?? ut some are just numbers....
#mean(is.na(lirype$recordedBy))#0.3178803
#lirype <- subset(lirype,!is.na(recordedBy))
#out <- ddply(lirype,.(recordedBy),summarise,nuObs=length(name))

#subsetting
lirype<-subset(lirype,year>2006&year<2018)
#exclude Nov to March  - rock and willow both have white coats#
lirype<-subset(lirype,month >4 & month <10)
table(lirype$year)

#make spatial and convert CRS
coordinates(lirype)<-c("decimalLongitude","decimalLatitude")
proj4string(lirype) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#overlay to the grid
plot(NorwayOrig)
plot(lirype,add=T,pch=16,cex=0.1,col= alpha("red",0.1))
#all points look pretty good!!

#pull out the grid cells
lirype <- sp::spTransform(lirype, crs(equalM))
plot(mygrid)
plot(lirype,add=T,pch=16,cex=0.5,col= alpha("red",0.1))
lirype$grid <- raster::extract(mygrid,lirype)

#subset to focal grids
lirype@data <- subset(lirype@data,grid %in% focusGrids)

#get number of points per grid cell
gridSummary<-ddply(lirype@data,.(grid),summarise,nuRecs=length(name))
mygrid[]<-0
mygrid[gridSummary$grid] <- gridSummary$nuRecs
length(gridSummary$grid)==length(gridSummary$nuRecs)
plot(mygrid)

#how many records per site are there
hist(gridSummary$nuRecs)
gridSummary$RepeatedVisits<-sapply(gridSummary$nuRecs,function(x)ifelse(x>1,1,0))
table(gridSummary$RepeatedVisits)


#grid cell per year
gridSummary<-ddply(lirype@data,.(grid,year),summarise,nuRecs=length(name))
gridSummary$RepeatedVisits<-sapply(gridSummary$nuRecs,function(x)ifelse(x>1,1,0))
table(gridSummary$RepeatedVisits)


lirype <- lirype@data
nrow(lirype)
#12135

##########################################################################

#get data for all birds - as an effort layer

#for first time:

# #get GBIF data
#tdir<-"C:/Users/diana.bowler/OneDrive - NINA/Alpine/GBIF"
tdir<-"GBIF"
speciesFiles<-list.files(tdir)[grepl("all_birds",list.files(tdir))]
library(plyr)
allbirds<-ldply(speciesFiles,function(x){
  read.delim(paste(tdir,x,sep="/"),as.is=T,fileEncoding="latin1",quote = "")
})

#clean names
names(allbirds) <- gsub("X.","",names(allbirds))
names(allbirds) <- gsub("\\.","",names(allbirds))

#filter
allbirds<-allbirds[,c("decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters",
                      "name","year","month","day","recordedBy")]
allbirds$decimalLatitude<-as.numeric(allbirds$decimalLatitude)
allbirds$decimalLongitude<-as.numeric(allbirds$decimalLongitude)
allbirds$year<-as.numeric(allbirds$year)
allbirds$month<-as.numeric(allbirds$month)
allbirds$day<-as.numeric(allbirds$day)
allbirds <- subset(allbirds,!is.na(year))
allbirds <- subset(allbirds,!is.na(month))
allbirds<-subset(allbirds,!is.na(allbirds$decimalLatitude))
allbirds<-subset(allbirds,!is.na(allbirds$decimalLongitude))

#get artsdatenbanken data
#tdir<-"C:/Users/diana.bowler/OneDrive - NINA/Alpine/SpeciesMapServices/Birds"
tdir<-"SpeciesMapServices/Birds"
allbirds2 <- read.delim(paste(tdir,"DataExportFraArtskart.txt",sep="/"),skipNul=T,dec=",",as.is=T)
allbirds2$Species<-apply(allbirds2,1,function(x)paste(x["Genus"],x["Species"],sep=" "))
allbirds2<-allbirds2[,c("Latitude","Longitude","CoordinatePrecision",
                        "Species","YearCollected","MonthCollected","DayCollected","Collector")]
names(allbirds2)<-names(allbirds)
allbirds2 <- subset(allbirds2,!is.na(day))

#combine them
allbirds<-rbind(allbirds,allbirds2)

#also with ptarmigan data
allbirds<-rbind(allbirds,lirypeD)
allbirds$year <- as.numeric(allbirds$year)

#tidy species names
allbirds$name <- gsub('[\"]','',allbirds$name)
sort(unique(allbirds$name))

#remove records that arent species level
allbirds <- subset(allbirds, grepl(" ",allbirds$name))

#remove duplicates
allbirds$date<-paste(allbirds$year,allbirds$month,allbirds$day,sep="-")
allbirds$id<-paste(allbirds$name,allbirds$date,round(allbirds$decimalLatitude,digits=4),
                   round(allbirds$decimalLongitude,digits=4),sep="_")
allbirds<-subset(allbirds,!duplicated(id))

#remove dodgy records
allbirds$coordinateUncertaintyInMeters <- as.numeric(allbirds$coordinateUncertaintyInMeters)
summary(allbirds$coordinateUncertaintyInMeters)
allbirds <- subset(allbirds, coordinateUncertaintyInMeters<5000|is.na(coordinateUncertaintyInMeters))
allbirds<-subset(allbirds,!is.na(allbirds$decimalLatitude)|!is.na(allbirds$decimalLongitude))

#remove those without more than 2 decimal places
allbirds$latDP<-getDecimalPlaces(allbirds$decimalLatitude)
allbirds$lonDP<-getDecimalPlaces(allbirds$decimalLongitude)
allbirds<-subset(allbirds,latDP>2 & lonDP>2)

#extract survyed grids within time period of interest
allbirds<-subset(allbirds,month >4 & month <10)
allbirds<-subset(allbirds,year>2006&year<2018)
table(allbirds$year)

#rename
names(allbirds)[which(names(allbirds)=="name")]<-"Species"
nrow(allbirds)#7068029

saveRDS(allbirds,file="allbirds_uniqueRecords.rds")
 
#####################################################################

#subsequent times:
load("C:/Users/db40fysa/Dropbox/ptarmigan Upscaling/data/allbirds_uniqueRecords.RData")

#####################################################################################

#make spatial and get grid info
coordinates(allbirds)<-c("decimalLongitude","decimalLatitude")
proj4string(allbirds)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#plotting original
#plot(NorwayOrig)
#plot(allbirds[sample(1:nrow(allbirds),591174),],add=T,pch=16,cex=0.5,col= alpha("red",0.1))#10%

#plotting
allbirds <- spTransform(allbirds, CRS(equalM))
mygrid[]<-1:ncell(mygrid)
plot(mygrid)
plot(Norway,add=T)
plot(allbirds[1:1000,],add=T)
#all points look pretty good!! 

#exclude points beyond the mask
#out <- over(allbirds,Norway)
#allbirds <- allbirds[!is.na(out),]

#overlay to the grid
allbirds$grid <- extract(mygrid,allbirds)

#subset to focus grids
allbirds <- subset(allbirds, grid %in% focusGrids)

#reorganise
allbirds <-subset(allbirds,!is.na(Species))
allbirds <- subset(allbirds,!is.na(grid))
all_data <- data.frame(allbirds@data) 
all_data <- all_data[,c("year","month","day","grid","Species","recordedBy")]
table(all_data$year,all_data$month)

########################################################################

#downsample the dataset
nrow(all_data)#4276862

dataSummary <- ddply(all_data,.(grid,year,month),summarise,nuObs=length(unique(day)))
summary(dataSummary$nuObs)

all_data
#cap at 25?
all_dataS <- ddply(all_data,.(grid,year,month),function(x){
  
  #get number of observations
  uniq = unique(x$day)
  nuObs = length(uniq)
  
  if(nuObs > 4){
    mySample = sample(uniq,4)
    subset(x, day %in% mySample)
  }else{
    x
  }})
  
nrow(all_dataS)#1883262

#reduces the dataset by over a half

#########################################################################

#are all grid cells in focus grids included here?

length(focusGrids)#11846
length(unique(all_dataS$grid))#9209
obsGrids <- sort(unique(all_dataS$grid))
#data from 78% of grids

#where are we missing data?
missingGrids <- focusGrids[!focusGrids %in% obsGrids]
mygrid[] <- 0
mygrid[missingGrids] <- 1
plot(mygrid)

#where we have data?
mygrid[] <- 0
mygrid[obsGrids] <- 1
plot(mygrid)

##################################################################################################

#get list length info

#add list length of all birds on each sampling visit
listlengthDF <- ddply(all_dataS,.(year,month,day,grid),summarise,
                    L = length(unique(Species)), #number of species
                    L2 = length(Species))
                    
nrow(listlengthDF)#322363

#get total number of species ever seen per grid
# gridRichness <- ddply(all_dataS,.(grid),summarise,nuSpecies=length(unique(Species)))
# mygrid[] <- 0
# mygrid[gridRichness$grid] <- gridRichness$nuSpecies
# plot(mygrid)
# 
# #add species richness to the listlength df
# listlengthDF$nuSpecies <- gridRichness$nuSpecies[match(listlengthDF$grid,gridRichness$grid)]
# summary(listlengthDF$nuSpecies)

##################################################################################

#add on ptarmigan obs

listlengthDF$visit <- paste(listlengthDF$year,listlengthDF$month,
                          listlengthDF$day,listlengthDF$grid,sep="-")
lirype$visit <- paste(lirype$year,lirype$month,
                    lirype$day,lirype$grid,sep="-")

listlengthDF$y <- sapply(listlengthDF$visit,function(x)ifelse(x%in%lirype$visit,1,0))

table(listlengthDF$y)
#0      1 
#314459   7904

####################################################################################################

#add NAs for all years and grids

newgrid <- expand.grid(year=unique(listlengthDF$year),
                     grid=focusGrids)

listlengthDF <- merge(listlengthDF,newgrid,by=c("year","grid"),all.y=T)
summary(listlengthDF)
#table(listlengthDF$grid,listlengthDF$year)

######################################################################################################

#adding non-detections??
#listlengthDF$L[is.na(listlengthDF$L)] <- 0
#listlengthDF$y[listlengthDF$L==0] <- 0

#########################################
#######################################################################

#order data by site and year
listlengthDF<-arrange(listlengthDF,year,grid)

#add indices
listlengthDF$siteIndex<-as.numeric(as.factor(listlengthDF$grid))
listlengthDF$yearIndex<-as.numeric(as.factor(listlengthDF$year))

saveRDS(listlengthDF,file="listlength_iDiv.rds")

######################################################################################################


