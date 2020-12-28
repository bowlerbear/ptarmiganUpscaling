########################################################################

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

plot(mygrid)
plot(NorwayOrigProj,add=T)

#get cells at least 50% covering Norway
myGridMask <- data.frame(extract(mygrid,NorwayOrigProj,weights=T,normalizeWeights=F)[[1]])
myGridMask <- myGridMask$value[myGridMask$weight>0.5]
myGrid2 <- mygrid
myGrid2[!getValues(mygrid) %in% myGridMask]<-NA
plot(myGrid2)
plot(NorwayOrigProj,add=T)
#looks good!

#get these grid numbers
focusGrids <- getValues(myGrid2)[!is.na(getValues(myGrid2))]
focusGridsDF <- as.data.frame(myGrid2,xy=T)

#########################################################################

#get willow ptarmigan data - occurrence data

#get GBIF data
tdir <- "C:/Users/db40fysa/Dropbox/Alpine/GBIF"
lirype <- read.delim(paste(tdir,"willow_ptarmigan_GBIF.txt",sep="/"),as.is=T,row.names=NULL)
lirype$decimalLatitude<-as.numeric(lirype$decimalLatitude)
lirype$decimalLongitude<-as.numeric(lirype$decimalLongitude)
lirype<-subset(lirype,!is.na(lirype$decimalLatitude)|!is.na(lirype$decimalLongitude))
lirype<-lirype[,c("decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters",
                  "name","year","month","day","recordedBy")]

#get artsdatenbanken data
tdir<-"C:/Users/db40fysa/Dropbox/Alpine/SpeciesMapServices/Ptarmigan"
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
lirype <- subset(lirype, !coordinateUncertaintyInMeters>5000|is.na(coordinateUncertaintyInMeters))

#remove those without 2 decimal places
lirype$latDP<-getDecimalPlaces(lirype$decimalLatitude)
lirype$lonDP<-getDecimalPlaces(lirype$decimalLongitude)
lirype<-subset(lirype,latDP>2 & lonDP>2)

#subsetting
lirype<-subset(lirype,year>2006&year<2018)
#exclude Nov to March  - rock and willow both have white coats#
lirype<-subset(lirype,month >3 & month <11)
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
lirype$grid <- extract(mygrid,lirype)

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
#0    1 
#974 1175

#grid cell per year
gridSummary<-ddply(lirype@data,.(grid,year),summarise,nuRecs=length(name))
gridSummary$RepeatedVisits<-sapply(gridSummary$nuRecs,function(x)ifelse(x>1,1,0))
table(gridSummary$RepeatedVisits)
#0    1 
#3217 1815

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
# allbirds <- subset(allbirds, !coordinateUncertaintyInMeters>5000|is.na(coordinateUncertaintyInMeters))
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
load("C:/Users/db40fysa/Dropbox/ptarmigan Upscaling/data/allbirds_uniqueRecords.RData")
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

#random sample
allbirds <- allbirds[sample(1:nrow(allbirds),3000000),]

#make spatial
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
out <- over(allbirds,Norway)
allbirds <- allbirds[!is.na(out),]

#overlay to the grid
allbirds$grid <- extract(mygrid,allbirds)

#identify all grid cells within the buffer region
mygridMask <- mask(mygrid,Norway)
plot(mygridMask)
mygridMaskDF <- as.data.frame(mygridMask,xy=T)
mygridMaskDF <- subset(mygridMaskDF,!is.na(layer))

#use existing terminology 
all(focusGrids %in% mygridMaskDF$layer)#TRUE
all(mygridMaskDF$layer %in% focusGrids)#FALSE
#focusGrids <- sort(mygridMaskDF$layer)
#focusGridsDF <- mygridMaskDF

#subset to focus grids
allbirds <- subset(allbirds, grid %in% focusGrids)

#reorganise
allbirds <-subset(allbirds,!is.na(Species))
allbirds <- subset(allbirds,!is.na(grid))
all_data <- data.frame(allbirds@data) 
all_data <- all_data[,c("year","month","day","grid","Species")]
table(all_data$year,all_data$month)

#########################################################################

#are all grid cells in focus grids included here?

length(focusGrids)#11846
length(unique(all_data$grid))

#plot focusGrids
plot(myGrid2)

#plot where we have data
focusGridsDF_Data <- subset(focusGridsDF, (layer %in% unique(all_data$grid)))
dataRaster <- rasterize(focusGridsDF_Data [,1:2],field=focusGridsDF_Data [,3],mygrid)
par(mfrow=c(1,2))
plot(dataRaster)

#plot where we dont have data
focusGridsDF_nonData <- subset(focusGridsDF,!is.na(layer))
focusGridsDF_nonData <- subset(focusGridsDF_nonData, !(layer %in% unique(all_data$grid)))
focusGridsDF_nonData <- subset(focusGridsDF_nonData,layer %in% focusGrids)
nondataRaster <- rasterize(focusGridsDF_nonData[,1:2],field=focusGridsDF_nonData [,3],mygrid)
plot(nondataRaster)

#identify these cells
missingGrids <- focusGridsDF_nonData$layer
length(missingGrids)#3307
#total grids
length(focusGrids)#11846
length(missingGrids)/length(focusGrids)#0.28

par(mfrow=c(1,1))

#########################################################################

#get list length info

#add list length of all birds on each sampling visit
listlengthDF <- ddply(all_data,.(year,month,day,grid),summarise,
                    L = length(unique(Species[!is.na(Species)])), #number of species
                    L2 = length(Species[!is.na(Species)]), #number of records
                    L3 = L2/L) #records per species
summary(listlengthDF$L)

#cap list length
# table(listlengthDF$L)
# listlengthDF$L[listlengthDF$L>quantile(listlengthDF$L,0.975,na.rm=T)]<-quantile(listlengthDF$L,0.975,na.rm=T)
# 
# #other measures
# hist(listlengthDF$L3)
# listlengthDF$L2[listlengthDF$L2>quantile(listlengthDF$L2,0.975,na.rm=T)]<-quantile(listlengthDF$L2,0.975,na.rm=T)
# listlengthDF$L3[listlengthDF$L3>quantile(listlengthDF$L3,0.975,na.rm=T)]<-quantile(listlengthDF$L3,0.975,na.rm=T)

#get total number of species ever seen per grid
gridRichness <- ddply(all_data,.(grid),summarise,nuSpecies=length(unique(Species[!is.na(Species)])),nuRecs=length(Species[!is.na(Species)]))
mygrid[]<-0
mygrid[gridRichness$grid] <- gridRichness$nuSpecies
plot(mygrid)

#add species richnes to the listlength df
listlengthDF$nuSpecies <- gridRichness$nuSpecies[match(listlengthDF$grid,gridRichness$grid)]
summary(listlengthDF$nuSpecies)

##################################################################################

#sort out occupancy matrix of ptarmigan
listlengthDF$visit <- paste(listlengthDF$year,listlengthDF$month,
                          listlengthDF$day,listlengthDF$grid,sep="-")
lirype$visit <- paste(lirype$year,lirype$month,
                    lirype$day,lirype$grid,sep="-")

listlengthDF$y <- sapply(listlengthDF$visit,function(x)ifelse(x%in%lirype@data$visit,1,0))

table(listlengthDF$y)
#0      1 
#399092  6477 

#add NAs for all years and grids
newgrid <- expand.grid(year=unique(listlengthDF$year),
                     grid=focusGrids)
listlengthDF <- merge(listlengthDF,newgrid,by=c("year","grid"),all.y=T)
summary(listlengthDF)
#table(listlengthDF$grid,listlengthDF$year)

#########################################

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

######################################################################################################

#get environmental data
load("C:/Users/db40fysa/Dropbox/ptarmigan Upscaling/data/varDF_allEnvironData_5km.RData")

#check all focus grids present
missingenvironData <- focusGrids[!focusGrids %in% varDF$grid]#62 are missing - drop these sites
length(missingenvironData)#61 grids

#where are they?
mygrid[] <- 0
mygrid[missingenvironData]<-1
plot(mygrid)#at the edge

#check bio1
mygrid[] <- NA
mygrid[varDF$grid]<-varDF$bio1
plot(mygrid)

#check open
mygrid[] <- NA
mygrid[varDF$grid]<-varDF$Open
plot(mygrid)

#check adm
mygrid[] <- NA
mygrid[varDF$grid[varDF$adm=="Finnmark"]] <- 1
plot(mygrid)
mygrid[varDF$grid[varDF$adm=="outside"]] <- 1
plot(mygrid)#along the coast

#check adm2
mygrid[varDF$grid[varDF$adm2=="outside"]] <- 1
plot(mygrid)#along the coast

#tree line
mygrid[varDF$grid]<-varDF$tree_line_position
plot(mygrid)

###############################################################################

#exclude all those "ouside" with adm or adm2
varDF <- subset(varDF,adm!="outside")
varDF <- subset(varDF,adm2!="outside")

#check plot again
mygrid[] <- NA
mygrid[varDF$grid] <- varDF$bio1
plot(mygrid)

###sort adm out################################################################

#not done - in the sea

# mygrid[] <- 1:ncell(mygrid)
# 
# norwayAdm <- getData('GADM',country='NOR',level=1)
# norwayAdm <- spTransform(norwayAdm,CRS(equalM))
# plot(mygrid)
# plot(norwayAdm,add=T)
# 
# #extract adm data for each grid
# out1 <- extract(mygrid,norwayAdm,df=T)
# nrow(out1)
# unique(out1$ID)#i guess in same order as data frame
# norwayAdm@data$ID <- unique(out1$ID)
# out1$Adm <- norwayAdm@data$NAME_1[match(out1$ID,norwayAdm@data$ID)]
# 
# #add to varDF
# varDF$newAdm <- out1$Adm[match(varDF$grid,out1$layer)]#still many blanks!!
# varDF$x <- mygridMaskDF$x[match(varDF$grid,mygridMaskDF$layer)]
# varDF$y <- mygridMaskDF$y[match(varDF$grid,mygridMaskDF$layer)]
# 
# #fill in nas with mode n 3 x 3
# varDF$admIndex <- as.numeric(as.factor(varDF$newAdm))
# varDF_noNA <- subset(varDF,!is.na(admIndex))
# admRaster <- rasterize(varDF_noNA[,c("x","y")],field=varDF_noNA[,"admIndex"],mygrid)
# 
# Mode <- function(x, na.rm = FALSE) {
#   if(na.rm){
#     x = x[!is.na(x)]
#   }
#   
#   ux <- unique(x)
#   return(ux[which.max(tabulate(match(x, ux)))])
# }
# 
# fill.na <- function(x, i=85) {
#   if( is.na(x)[i]) {
#     return(Mode(x, na.rm=TRUE))
#   } else {
#     return(x)
#   }
# }
#   
# r2 <- focal(admRaster, w = matrix(1,13,13), fun = fill.na, 
#             pad = TRUE)
#   
# #align new adm and grid numbers
# r2DF <- as.data.frame(r2,xy=T)
# myGridDF$admIndex <- r2DF$layer
# varDF$newAdmIndex <- myGridDF$admIndex[match(varDF$grid,myGridDF$layer)]
# varDF$admIndex[is.na(varDF$admIndex)] <- varDF$newAdmIndex[is.na(varDF$admIndex)]
# summary(varDF$admIndex)
# #check
# admRaster <- rasterize(varDF[,c("x","y")],field=varDF[,"admIndex"],mygrid)
# plot(admRaster)
# admMap <- unique(varDF_noNA[,c("adm","admIndex")])
# varDF$adm <- admMap$adm[match(varDF$admIndex,admMap$admIndex)]
# table(varDF$adm)
# 
# ###sort adm2###################################################################
# 
# mygrid[] <- 1:ncell(mygrid)
# 
# norwayAdm <- getData('GADM',country='NOR',level=2)
# norwayAdm <- spTransform(norwayAdm,CRS(equalM))
# plot(mygrid)
# plot(norwayAdm,add=T)
# 
# #extract adm data for each grid
# out1 <- extract(mygrid,norwayAdm,df=T)
# nrow(out1)
# unique(out1$ID)#i guess in same order as data frame
# norwayAdm@data$ID <- unique(out1$ID)
# out1$Adm2 <- norwayAdm@data$NAME_2[match(out1$ID,norwayAdm@data$ID)]
# 
# #add to varDF
# varDF$newAdm2 <- out1$Adm2[match(varDF$grid,out1$layer)]#still many blanks!!
# 
# #fill in nas with mode n 3 x 3
# varDF$admIndex2 <- as.numeric(as.factor(varDF$newAdm2))
# varDF_noNA <- subset(varDF,!is.na(admIndex2))
# admRaster <- rasterize(varDF_noNA[,c("x","y")],field=varDF_noNA[,"admIndex2"],mygrid)
# 
# Mode <- function(x, na.rm = FALSE) {
#   if(na.rm){
#     x = x[!is.na(x)]
#   }
#   
#   ux <- unique(x)
#   return(ux[which.max(tabulate(match(x, ux)))])
# }
# 
# fill.na <- function(x, i=13) {
#   if(is.na(x)[i]) {
#     return(Mode(x, na.rm=TRUE))
#   } else {
#     return(x)
#   }
# }
# 
# r2 <- focal(admRaster, w = matrix(1,5,5), fun = fill.na, 
#             pad = TRUE)
# 
# r2 <- focal(r2, w = matrix(1,5,5), fun = fill.na, 
#             pad = TRUE)
# 
# #align new adm and grid numbers
# r2DF <- as.data.frame(r2,xy=T)
# myGridDF$admIndex2 <- r2DF$layer
# varDF$newAdmIndex2 <- myGridDF$admIndex2[match(varDF$grid,myGridDF$layer)]
# varDF$admIndex2[is.na(varDF$admIndex)] <- varDF$newAdmIndex2[is.na(varDF$admIndex2)]
# summary(varDF$admIndex2)
# #check
# admRaster <- rasterize(varDF[,c("x","y")],field=varDF[,"admIndex"],mygrid)
# plot(admRaster)
# admMap <- unique(varDF_noNA[,c("adm","admIndex")])
# varDF$adm <- admMap$adm[match(varDF$admIndex,admMap$admIndex)]
# table(varDF$adm)


######################################################################################

#cut data

#missingenvironData - check environ data eventually

#if testing, reduce the number of sites - get number of times each grid was visited
#gridSummary <- ddply(listlengthDF,.(grid),summarize,nuVisits=length(y[!is.na(y)]))
#listlengthDF <- subset(listlengthDF,grid%in%gridSummary$grid[gridSummary$nuVisits>0])

#merge data
listlengthDF <- merge(listlengthDF,varDF,by="grid",sort=F)

##################################################################################################

#can we limit by the tree line?

summary(listlengthDF$tree_line_position)
summary(listlengthDF$tree_line_position[listlengthDF$y==1&!is.na(listlengthDF$y)])
#seen at most at -553

#nrow(subset(listlengthDF,tree_line_position<(-536)))#3475
#nrow(subset(listlengthDF,tree_line_position<(-560) & !is.na(y)))#119...

#exclude very high sites
listlengthDF$High <- ifelse(listlengthDF$tree_line_position<(-560),1,0)
listlengthDF <- subset(listlengthDF, High==0)

##################################################################################################

#add adm
table(listlengthDF$adm)
listlengthDF$admN<-as.numeric(factor(listlengthDF$adm))
table(listlengthDF$adm2)
listlengthDF$admN2<-as.numeric(factor(listlengthDF$adm2))

#add indices
listlengthDF<-arrange(listlengthDF,admN,admN2)
listlengthDF$site<-paste(listlengthDF$admN,listlengthDF$admN2,listlengthDF$grid,sep="-")
listlengthDF$siteIndex<-as.numeric(factor(listlengthDF$site))
listlengthDF$yearIndex<-as.numeric(factor(listlengthDF$year))

#order data by site and year
listlengthDF<-arrange(listlengthDF,yearIndex,siteIndex)

#extract site data
listlengthDF_SiteCovs<-subset(listlengthDF,!duplicated(grid))
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