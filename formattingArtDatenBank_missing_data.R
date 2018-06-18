########################################################################

#load in general functions

source('C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/generalFunctions.R')

#get norway##############################################################

library(raster)
library(maptools)
library(rgeos)
data(wrld_simpl)
Norway<- subset(wrld_simpl,NAME=="Norway")
Norway <- gBuffer(Norway,width=1)

########################################################################

#choose modelling resolution:

#newres=1

newres=0.05

#########################################################################

#get willow ptarmigan data - occurrence data

#get GBIF data
tdir<-"C:/Users/diana.bowler/OneDrive - NINA/Alpine/GBIF"
lirype <- read.delim(paste(tdir,"willow_ptarmigan_GBIF.txt",sep="/"),as.is=T,row.names=NULL)
lirype$decimalLatitude<-as.numeric(lirype$decimalLatitude)
lirype$decimalLongitude<-as.numeric(lirype$decimalLongitude)
lirype<-subset(lirype,!is.na(lirype$decimalLatitude)|!is.na(lirype$decimalLongitude))
lirype<-lirype[,c("decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters",
                  "name","year","month","day","recordedBy")]

#get artsdatenbanken data
tdir<-"C:/Users/diana.bowler/OneDrive - NINA/Alpine/SpeciesMapServices/Ptarmigan"
lirype2 <- read.delim(paste(tdir,"DataExportFraArtskart.txt",sep="/"),skipNul=T,dec=",")
lirype2$Species<-apply(lirype2,1,function(x)paste(x["Genus"],x["Species"],sep=" "))
lirype2<-lirype2[,c("Latitude","Longitude","CoordinatePrecision",
                    "Species","YearCollected","MonthCollected","DayCollected","Collector")]
names(lirype2)<-names(lirype)

#combine them
lirype<-rbind(lirype,lirype2)
lirypeD<-lirype

#remove duplicates
lirype$date<-paste(lirype$year,lirype$month,lirype$day,sep="-")
lirype$id<-paste(lirype$date,round(lirype$decimalLatitude,digits=4),round(lirype$decimalLongitude,digits=4),sep="_")
lirype<-subset(lirype,!duplicated(id))

#remove dodgy records
summary(lirype$coordinateUncertaintyInMeters)
lirype <- subset(lirype, !coordinateUncertaintyInMeters>10000|is.na(coordinateUncertaintyInMeters))

#remove those without 2 decimal places
lirype$latDP<-getDecimalPlaces(lirype$decimalLatitude)
lirype$lonDP<-getDecimalPlaces(lirype$decimalLongitude)
lirype<-subset(lirype,latDP>2 & lonDP>2)

#subsetting
lirype<-subset(lirype,year>2006&year<2018)
#exclude Nov to March  - rock and willow both have white coats#
lirype<-subset(lirype,month >3 & month <11)
table(lirype$year)

#make spatial
coordinates(lirype)<-c("decimalLongitude","decimalLatitude")
plot(Norway)
plot(lirype,add=T)
#all points look pretty good!! few in the sea

#overlay to the grid
library(raster)
mygrid<-raster(extent(Norway))
res(mygrid)<-newres
mygrid[]<-1:ncell(mygrid)
#plot(mygrid)
gridTemp<-mygrid
#get grid number for each point
lirype$grid<-extract(mygrid,lirype)

#get number of points per grid cell
library(plyr)
gridSummary<-ddply(lirype@data,.(grid),summarise,nuRecs=length(name))
mygrid[]<-0
mygrid[gridSummary$grid]<-gridSummary$nuRecs
length(gridSummary$grid)==length(gridSummary$nuRecs)
#plot(mygrid)

#how many records per site are there
hist(gridSummary$nuRecs)
gridSummary$RepeatedVisits<-sapply(gridSummary$nuRecs,function(x)ifelse(x>1,1,0))
#table(gridSummary$RepeatedVisits)
#0    1 
#1654 1662

#grid cell per year
gridSummary<-ddply(lirype@data,.(grid,year),summarise,nuRecs=length(name))
gridSummary$RepeatedVisits<-sapply(gridSummary$nuRecs,function(x)ifelse(x>1,1,0))
#table(gridSummary$RepeatedVisits)
#0    1 
#4726 2435

##########################################################################

#get data for all birds - as an effort layer

#for first time:

# #get GBIF data
# tdir<-"C:/Users/diana.bowler/OneDrive - NINA/Alpine/GBIF"
# speciesFiles<-list.files(tdir)[grepl("all_birds",list.files(tdir))]
# library(plyr)
# allbirds<-ldply(speciesFiles,function(x){
#   read.delim(paste(tdir,x,sep="/"),as.is=T,row.names=NULL)
# })
# 
# allbirds<-allbirds[,c("decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters",
#                       "name","year","month","day","recordedBy")]
# allbirds$decimalLatitude<-as.numeric(allbirds$decimalLatitude)
# allbirds$decimalLongitude<-as.numeric(allbirds$decimalLongitude)
# allbirds$year<-as.numeric(allbirds$year)
# allbirds$month<-as.numeric(allbirds$month)
# allbirds<-subset(allbirds,!is.na(allbirds$decimalLatitude))
# allbirds<-subset(allbirds,!is.na(allbirds$decimalLongitude))
# 
# #get artsdatenbanken data
# tdir<-"C:/Users/diana.bowler/OneDrive - NINA/Alpine/SpeciesMapServices/Birds"
# allbirds2 <- read.delim(paste(tdir,"DataExportFraArtskart.txt",sep="/"),skipNul=T,dec=",",as.is=T)
# allbirds2$Species<-apply(allbirds2,1,function(x)paste(x["Genus"],x["Species"],sep=" "))
# allbirds2<-allbirds2[,c("Latitude","Longitude","CoordinatePrecision",
#                         "Species","YearCollected","MonthCollected","DayCollected","Collector")]
# names(allbirds2)<-names(allbirds)
# 
# #combine them
# allbirds<-rbind(allbirds,allbirds2)
# 
# #also with ptarmigan data
# allbirds<-rbind(allbirds,lirypeD)
# 
# #remove duplicates
# allbirds$date<-paste(allbirds$year,allbirds$month,allbirds$day,sep="-")
# allbirds$id<-paste(allbirds$name,allbirds$date,round(allbirds$decimalLatitude,digits=4),
#                    round(allbirds$decimalLongitude,digits=4),sep="_")
# allbirds<-subset(allbirds,!duplicated(id))
# 
# #remove dodgy records
# summary(allbirds$coordinateUncertaintyInMeters)
# allbirds <- subset(allbirds, !coordinateUncertaintyInMeters>10000|is.na(coordinateUncertaintyInMeters))
# allbirds<-subset(allbirds,!is.na(allbirds$decimalLatitude)|!is.na(allbirds$decimalLongitude))
# 
# #remove those without 2 decimal places
# allbirds$latDP<-getDecimalPlaces(allbirds$decimalLatitude)
# allbirds$lonDP<-getDecimalPlaces(allbirds$decimalLongitude)
# allbirds<-subset(allbirds,latDP>2 & lonDP>2)
# 
# #extract survyed grids within time period of interest
# allbirds<-subset(allbirds,month >3 & month <11)
# allbirds<-subset(allbirds,year>2006&year<2018)
# table(allbirds$year)
# save(allbirds,file="allbirds_uniqueRecords.RData")
# 
# #sort out observer ids
# allRecorders<-names(table(allbirds$recordedBy))
# #if comma separate the names
# allRecorders<-sort(unique(trim(unlist(lapply(allRecorders,function(x)strsplit(x,",")[[1]])))))
# #work out how many times each has contributed
# effort<-as.numeric()
# for(i in 1:length(allRecorders)){
#   temp<-allbirds[grepl(allRecorders[i],allbirds$recordedBy),]
#   effort[i]<-length(unique(interaction(temp$month,temp$year)))
# }
# bestRecorders<-data.frame(allRecorders,effort)
# save(bestRecorders,file="bestRecorders.RData")

#####################################################################

#subsequent times:

load("allbirds_uniqueRecords.RData")
names(allbirds)[which(names(allbirds)=="name")]<-"Species"

#####################################################################################

#if removing poor recorders? 
# load("bestRecorders.RData")
# bestRecorders<-bestRecorders[order(bestRecorders$effort,decreasing=TRUE),]
# head(bestRecorders)
# #remove those without spaces in the name (i.e., first name and surname)
# bestRecorders<-bestRecorders[grepl("\\ ",bestRecorders$allRecorders),]
# summary(bestRecorders$effort)
# bestRecorders<-subset(bestRecorders,effort>52)#upper quartile
# allbirds$Best<-sapply(allbirds$recordedBy,function(x){
#   all<-strsplit(x,",")[[1]]
#   all<-trim(all)
#   any(is.element(all,bestRecorders$allRecorders))
# })
# table(allbirds$Best)
# allbirds<-subset(allbirds,Best=="TRUE")

#####################################################################################

#make spatial
coordinates(allbirds)<-c("decimalLongitude","decimalLatitude")
proj4string(allbirds)<-crs(Norway)
#plot(Norway)
#plot(allbirds,add=T)
#all points look pretty good!! few in the sea...

#exclude points beyond the mask
out<-over(allbirds,Norway)
allbirds<-allbirds[!is.na(out),]

#overlay to the grid
library(raster)
mygrid<-raster(extent(Norway))
res(mygrid)<-newres
mygrid[]<-1:ncell(mygrid)
plot(mygrid)
#get grid number for each point
allbirds$grid<-extract(mygrid,allbirds)

#identify all grid cells within the buffer region
library(rgeos)
mygridMask<-mask(mygrid,Norway)
plot(mygridMask)
mygridMaskDF<-as.data.frame(mygridMask,xy=T)
mygridMaskDF<-subset(mygridMaskDF,!is.na(layer))

#recoganise
allbirds<-subset(allbirds,!is.na(Species))
allbirds <- subset(allbirds,!is.na(grid))
all_data<-data.frame(allbirds@data,allbirds@coords) 
all_data <- all_data[,c("year","month","day","grid","Species")]
table(all_data$year,all_data$month)

#########################################################################

#get info on administrative names for the grid

library(rgdal)
NorwayADM<-readOGR(dsn="C:/Users/diana.bowler/OneDrive - NINA/Alpine/NOR_adm",layer="NOR_adm2")
crs(NorwayADM)
mygridDF<-as.data.frame(mygrid,xy=T)
mygridPoints<-mygridDF
coordinates(mygridPoints)<-c("x","y")
proj4string(mygridPoints)<-crs(NorwayADM)
myAdm<-over(mygridPoints,NorwayADM)
myAdm$grid<-mygridPoints@data$layer

myAdm$VARNAME_1<-as.character(myAdm$NAME_1)
myAdm$VARNAME_1[is.na(myAdm$VARNAME_1)]<-"outside"
myAdm$VARNAME_2<-as.character(myAdm$NAME_2)
myAdm$VARNAME_2[is.na(myAdm$VARNAME_2)]<-"outside"
table(myAdm$VARNAME_1)
table(myAdm$VARNAME_2)
mygrid[myAdm$grid]<-as.numeric(as.factor(myAdm$VARNAME_1))
plot(mygrid)#looks good!
mygrid[myAdm$grid]<-as.numeric(as.factor(myAdm$VARNAME_2))
plot(mygrid)#looks good!

#########################################################################

#get list length info

#add list length of all birds on each sampling visit
listlengthDF<-ddply(all_data,.(year,month,day,grid),summarise,
                    L = length(unique(Species[!is.na(Species)])), #number of species
                    L2 = length(Species[!is.na(Species)]), #number of records
                    L3 = L2/L) #records per species

summary(listlengthDF$L)
#cap list length
table(listlengthDF$L)
listlengthDF$L[listlengthDF$L>quantile(listlengthDF$L,0.975,na.rm=T)]<-quantile(listlengthDF$L,0.975,na.rm=T)

#other measures
hist(listlengthDF$L3)
listlengthDF$L2[listlengthDF$L2>quantile(listlengthDF$L2,0.975,na.rm=T)]<-quantile(listlengthDF$L2,0.975,na.rm=T)
listlengthDF$L3[listlengthDF$L3>quantile(listlengthDF$L3,0.975,na.rm=T)]<-quantile(listlengthDF$L3,0.975,na.rm=T)

#get total number of species ever seen per grid
gridRichness <- ddply(all_data,.(grid),summarise,nuSpecies=length(unique(Species[!is.na(Species)])),nuRecs=length(Species[!is.na(Species)]))
mygrid[]<-0
gridRichness$nuRecs[gridRichness$nuRecs>quantile(gridRichness$nuRecs,0.975)]<-quantile(gridRichness$nuRecs,0.975)
mygrid[gridRichness$grid]<-gridRichness$nuRecs/gridRichness$nuSpecies
plot(mygrid)
mygrid[gridRichness$grid]<-gridRichness$nuRecs
plot(mygrid)

#add species richnes to the listlength df
listlengthDF$nuSpecies <- gridRichness$nuSpecies[match(listlengthDF$grid,gridRichness$grid)]
summary(listlengthDF$nuSpecies)

##################################################################################

#sort out occupancy matrix of ptarmigan

listlengthDF$visit<-paste(listlengthDF$year,listlengthDF$month,
                          listlengthDF$day,listlengthDF$grid,sep="-")
lirype$visit<-paste(lirype$year,lirype$month,
                    lirype$day,lirype$grid,sep="-")
listlengthDF$y<-sapply(listlengthDF$visit,function(x)ifelse(x%in%lirype@data$visit,1,0))

table(listlengthDF$y)
#0      1 
#784145  10770 
#removing poor observers
#0      1 
#604912   6063

#add NAs for all years and grids
newgrid<-expand.grid(year=unique(listlengthDF$year),
                     grid=sort(unique(mygridMaskDF$layer)))
listlengthDF<-merge(listlengthDF,newgrid,by=c("year","grid"),all.y=T)
summary(listlengthDF)
table(listlengthDF$grid,listlengthDF$year)

#########################################
#adding non-detections??
#listlengthDF$L[is.na(listlengthDF$L)]<-0
#listlengthDF$y[listlengthDF$L==0]<-0
#########################################
#######################################################################

#order data by site and year
listlengthDF<-arrange(listlengthDF,year,grid)

#add indices
listlengthDF$siteIndex<-as.numeric(as.factor(listlengthDF$grid))
listlengthDF$yearIndex<-as.numeric(as.factor(listlengthDF$year))

#add adm
listlengthDF$adm<-myAdm$VARNAME_1[match(listlengthDF$grid,myAdm$grid)]
table(listlengthDF$adm)
listlengthDF$admN<-as.numeric(factor(listlengthDF$adm))

listlengthDF$adm2<-myAdm$VARNAME_2[match(listlengthDF$grid,myAdm$grid)]
table(listlengthDF$adm2)
listlengthDF$admN2<-as.numeric(factor(listlengthDF$adm2))

######################################################################################################

#combine environmental and pop data:

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling")
load("varDF_allEnvironData_5km.RData")
varDF<-subset(varDF,!duplicated(grid))
listlengthDF<-merge(listlengthDF,varDF,by="grid",sort=F)

#check where we have data
mygrid[]<-0
mygrid[listlengthDF$grid]<-listlengthDF$admN
plot(mygrid)

#if testing, reduce the number of sites - get number of times each grid was visited
#gridSummary <- ddply(listlengthDF,.(grid),summarize,nuVisits=length(y[!is.na(y)]))
#listlengthDF <- subset(listlengthDF,grid%in%gridSummary$grid[gridSummary$nuVisits>0])


##################################################################################################

#can we limit by the tree line?

summary(listlengthDF$tree_line_position)
summary(listlengthDF$tree_line_position[listlengthDF$y==1&!is.na(listlengthDF$y)])
#seen at most at -535m

#nrow(subset(listlengthDF,tree_line_position<(-536)))#3475
#nrow(subset(listlengthDF,tree_line_position<(-536) & !is.na(y)))#175...

#exclude very high sites
listlengthDF$High <- ifelse(listlengthDF$tree_line_position<(-540),1,0)
listlengthDF <- subset(listlengthDF, High==0)

##################################################################################################

#add indices
listlengthDF$site<-paste(listlengthDF$adm2,listlengthDF$grid,sep="-")
listlengthDF$siteIndex<-as.numeric(factor(listlengthDF$grid))
listlengthDF$yearIndex<-as.numeric(factor(listlengthDF$year))

#order data by site and year
listlengthDF<-arrange(listlengthDF,yearIndex,siteIndex)

#extract site data
listlengthDF_SiteCovs<-subset(listlengthDF,!duplicated(grid))

#adm data
listlengthDF$admN<-as.numeric(factor(listlengthDF$adm))
listlengthDF$admN2<-as.numeric(factor(listlengthDF$adm2))
siteInfo<-unique(listlengthDF[,c("siteIndex","admN","grid","admN2")])

############################################################################

#also get number of unique records for each year/site as another measure of effort

siteyearInfo<-ddply(all_data,.(year,grid),summarise,
              L = length(unique(Species[!is.na(Species)])), #number of species
              L2 = length(Species[!is.na(Species)]), #number of records
              L3 = L2/L,
              nuDays = length(unique(interaction(month,day))))#number of sampling days

summary(siteyearInfo$nuDays)
#up to 214!!! almost every days

averagesiteInfo<-ddply(siteyearInfo,.(grid),summarise,nuDays=mean(nuDays))
mygrid[]<-0
mygrid[averagesiteInfo$grid]<-log(averagesiteInfo$nuDays)
plot(mygrid,main="log nu sampling days")

############################################################################