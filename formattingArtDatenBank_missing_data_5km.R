#get norway##############################################################

library(raster)
library(maptools)
library(rgeos)
data(wrld_simpl)
Norway<- subset(wrld_simpl,NAME=="Norway")
plot(Norway)
Norway <- gBuffer(Norway,width=1)

########################################################################

#choose modelling resolution:

#newres=1

newres=0.05

#########################################################################

#get willow ptarmigan data - occurrence data

tdir<-"C:/Users/diana.bowler/OneDrive - NINA/Alpine/SpeciesMapServices/Ptarmigan"
lirype <- read.delim(paste(tdir,"DataExportFraArtskart.txt",sep="/"),skipNul=T,dec=",")
summary(lirype$CoordinatePrecision)
lirype <- subset(lirype, !CoordinatePrecision>10000)
coordinates(lirype)<-c("Longitude","Latitude")
plot(Norway)
plot(lirype,add=T)
#all points look pretty good!! few in the sea

###############################################################
#exclude Nov to March  - rock and willow both have white coats#
###############################################################

#overlay to the grid
library(raster)
mygrid<-raster(extent(Norway))
res(mygrid)<-newres
mygrid[]<-1:ncell(mygrid)
plot(mygrid)
gridTemp<-mygrid

#get grid number for each point
lirype$grid<-extract(mygrid,lirype)

#get number of points per grid cell
library(plyr)
gridSummary<-ddply(lirype@data,.(grid),summarise,nuRecs=length(NorskNavn))
mygrid[]<-0
mygrid[gridSummary$grid]<-gridSummary$nuRecs
length(gridSummary$grid)==length(gridSummary$nuRecs)
plot(mygrid)

#plot for each year
lirype<-subset(lirype,YearCollected>2006&YearCollected<2018)
#gridSummary<-ddply(lirype@data,.(grid,YearCollected),summarise,nuRecs=length(NorskNavn))
#sort(unique(gridSummary$YearCollected))
#par(mfrow=c(3,4))
#for(i in 2007:2017){
#  mygrid[]<-0
#  mygrid[gridSummary$grid[gridSummary$YearCollected==i]]<-gridSummary$nuRecs[gridSummary$YearCollected==i]
#  plot(mygrid,main=i)
#}

#how many records per site are there, per year
hist(gridSummary$nuRecs)
gridSummary$RepeatedVisits<-sapply(gridSummary$nuRecs,function(x)ifelse(x>1,1,0))
table(gridSummary$RepeatedVisits)
#0    1 
#4754 2173 

##########################################################################

#get data for all birds - as an effort layer

tdir<-"C:/Users/diana.bowler/OneDrive - NINA/Alpine/SpeciesMapServices/Ptarmigan"
lirypeDF <- read.delim(paste(tdir,"DataExportFraArtskart.txt",sep="/"),skipNul=T,dec=",",as.is=T)
tdir<-"C:/Users/diana.bowler/OneDrive - NINA/Alpine/SpeciesMapServices/Birds"
allbirds <- read.delim(paste(tdir,"DataExportFraArtskart.txt",sep="/"),skipNul=T,dec=",",as.is=T)
allbirds <- rbind(allbirds,lirypeDF)
unique(allbirds$Class)
summary(allbirds$CoordinatePrecision)
allbirds <- subset(allbirds, !CoordinatePrecision>10000)
coordinates(allbirds)<-c("Longitude","Latitude")
plot(Norway)
plot(allbirds,add=T)
#all points look pretty good!! few in the sea...

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
mygridMaskDF<-as.data.frame(mygridMask,xy=T)
mygridMaskDF<-subset(mygridMaskDF,!is.na(layer))

#survyed grids
allbirds<-subset(allbirds,YearCollected>2006&YearCollected<2018)
allbirds@data$Species<-apply(allbirds@data,1,function(x)paste(x["Genus"],x["Species"],sep=" "))
allbirds_data<-data.frame(allbirds@data,allbirds@coords) 
allbirds_data <- allbirds_data[,c("YearCollected","MonthCollected","DayCollected","grid","Species")]
allbirds_data <- subset(allbirds_data,DayCollected!=0)

#unsurveyed grids
nobirds_data <- expand.grid(YearCollected=2017,
                           grid=mygridMaskDF$layer[!mygridMaskDF$layer%in%allbirds_data$grid],#268 grids
                           MonthCollected=NA,
                           DayCollected=NA,
                           Species=NA)

#combine them
all_data <- rbind(allbirds_data,nobirds_data)
all_data <- subset(all_data,!is.na(grid))

#########################################################################

#get info on administrative names for the grid
library(rgdal)
NorwayADM<-readOGR(dsn="C:/Users/diana.bowler/OneDrive - NINA/Alpine/NOR_adm",layer="NOR_adm1")
unique(NorwayADM$VARNAME_1)
crs(NorwayADM)
mygridDF<-as.data.frame(mygrid,xy=T)
mygridPoints<-mygridDF
coordinates(mygridPoints)<-c("x","y")
proj4string(mygridPoints)<-crs(NorwayADM)
myAdm<-over(mygridPoints,NorwayADM)
myAdm$grid<-mygridPoints@data$layer
myAdm$VARNAME_1<-as.character(myAdm$VARNAME_1)
myAdm$VARNAME_1[is.na(myAdm$VARNAME_1)]<-"outside"
table(myAdm$VARNAME_1)
mygrid[myAdm$grid]<-as.numeric(as.factor(myAdm$VARNAME_1))
plot(mygrid)#looks good!

#########################################################################

#Run occupancy model for willow ptarmigan - using the sparta base code

#get total number of species ever seen per grid
gridRichness <- ddply(all_data,.(grid),summarise,nuSpecies=length(unique(Species[!is.na(Species)])),nuRecs=length(Species[!is.na(Species)]))
mygrid[]<-0
gridRichness$nuRecs[gridRichness$nuRecs>quantile(gridRichness$nuRecs,0.975)]<-quantile(gridRichness$nuRecs,0.975)
mygrid[gridRichness$grid]<-gridRichness$nuRecs/gridRichness$nuSpecies
plot(mygrid)
mygrid[gridRichness$grid]<-gridRichness$nuRecs
plot(mygrid)

#add list length of all birds on each sampling visit
listlengthDF<-ddply(all_data,.(YearCollected,MonthCollected,DayCollected,grid),summarise,
                    L = length(unique(Species[!is.na(Species)])), #number of species
                    L2 = length(Species[!is.na(Species)]), #number of records
                    L3 = L2/L) #records per species

#sort out occupancy matrix of ptarmigan
listlengthDF$visit<-paste(listlengthDF$YearCollected,listlengthDF$MonthCollected,
                          listlengthDF$DayCollected,listlengthDF$grid,sep="-")
lirype$visit<-paste(lirype$YearCollected,lirype$MonthCollected,
                     lirype$DayCollected,lirype$grid,sep="-")
listlengthDF$y<-sapply(listlengthDF$visit,function(x)ifelse(x%in%lirype$visit,1,0))

sum(listlengthDF$y,na.rm=T)
#10927

#add NAs for all years
newgrid<-expand.grid(YearCollected=unique(listlengthDF$YearCollected),grid=unique(listlengthDF$grid))
listlengthDF<-merge(listlengthDF,newgrid,by=c("YearCollected","grid"),all=T)

#add NAs for 2017
#newgrid<-expand.grid(YearCollected=2017,grid=unique(listlengthDF$grid))
#listlengthDF<-merge(listlengthDF,newgrid,by=c("YearCollected","grid"),all=T)

#########################################
#listlengthDF$L[is.na(listlengthDF$L)]<-0
#listlengthDF$y[listlengthDF$L==0]<-0
#########################################

#hist(listlengthDF$L)
summary(listlengthDF$L)
#cap list length
table(listlengthDF$L)
listlengthDF$L[listlengthDF$L>quantile(listlengthDF$L,0.975,na.rm=T)]<-quantile(listlengthDF$L,0.975,na.rm=T)

#other measures
#table(listlengthDF$L3)
#hist(listlengthDF$L3)
#hist(log(listlengthDF$L3+1))

#add species richnes for each grid
listlengthDF$nuSpecies <- gridRichness$nuSpecies[match(listlengthDF$grid,gridRichness$grid)]
summary(listlengthDF$nuSpecies)

#order data by site and year
listlengthDF<-arrange(listlengthDF,YearCollected,grid)

#add indices
listlengthDF$siteIndex<-as.numeric(as.factor(listlengthDF$grid))
listlengthDF$yearIndex<-as.numeric(as.factor(listlengthDF$YearCollected))

#add adm
listlengthDF$adm<-myAdm$VARNAME_1[match(listlengthDF$grid,myAdm$grid)]
table(listlengthDF$adm)
listlengthDF$adm2<-as.numeric(as.factor(listlengthDF$adm))

bugs.data <- list(nsite = length(unique(listlengthDF$siteIndex)),
                  nyear = length(unique(listlengthDF$yearIndex)),
                  nvisit = nrow(listlengthDF),
                  grid = listlengthDF$grid,
                  site = listlengthDF$siteIndex,
                  year = listlengthDF$yearIndex,
                  L = listlengthDF$L - median(listlengthDF$L),
                  y = listlengthDF$y)

########################################################################

#general functions:

#get function to extract raster into for these grids:

getEnvironData<-function(myraster,mygridTemp){
  require(maptools)
  
  #crop raster to Norway extent
  rasterCRS<-crs(myraster)
  Norway<-spTransform(Norway,rasterCRS)
  myraster<-crop(myraster,extent(Norway))
  
  #convert raster into points and convert to lon lat
  myrasterDF<-as.data.frame(myraster,xy=T)
  names(myrasterDF)[3]<-"myraster"
  coordinates(myrasterDF)<-c("x","y")
  proj4string(myrasterDF)<-rasterCRS
  myrasterDF<-spTransform(myrasterDF,"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
  
  #get general grid
  mygrid<-gridTemp
  projection(mygrid)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0") 
  
  #get mean myraster values per grid cell
  mygrid[]<-1:ncell(mygrid)
  variable<-extract(mygrid,myrasterDF)
  myrasterDF<-data.frame(myrasterDF@data)
  myrasterDF$grid<-variable
  myrasterDF<-ddply(myrasterDF,.(grid),summarise,myraster=mean(myraster,na.rm=T))
  return(myrasterDF)
}

plotBUGSData<-function(myvar){
  temp<-listlengthDF
  temp<-subset(temp,siteIndex %in% bugs.data$site)
  temp$variable<-as.numeric(bugs.data[myvar][[1]])[match(temp$siteIndex,bugs.data$site)]
  temp$variablePA<-sapply(temp$variable,function(x)ifelse(is.na(x),0,1))
  temp<-subset(temp,!duplicated(siteIndex))
  mygrid<-gridTemp
  par(mfrow=c(1,2))
  mygrid[]<-0
  mygrid[temp$grid]<-temp$variable
  plot(mygrid)
  mygrid[]<-0
  mygrid[temp$grid]<-temp$variablePA
  plot(mygrid)
}

plotZ<-function(model,param="z"){
  zSummary<-data.frame(model$summary[grepl(param,row.names(model$summary)),])
  zSummary$ParamNu <- as.character(sub(".*\\[([^][]+)].*", "\\1", row.names(zSummary)))
  zSummary$Site<-sapply(zSummary$ParamNu,function(x)strsplit(x,",")[[1]][1])
  zSummary$Year<-sapply(zSummary$ParamNu,function(x)strsplit(x,",")[[1]][2])
  zSummary$grid<-listlengthDF$grid[match(zSummary$Site,listlengthDF$siteIndex)]
  zSummary_year<-subset(zSummary,Year==11)
  
  #predicted occupancy across all years
  #zSummary_year<-ddply(zSummary,.(grid),summarise,prop=mean(mean))
  #plot mean prop occupancy
  par(mfrow=c(1,1))
  mygrid[]<-0
  mygrid[zSummary_year$grid]<-zSummary_year$mean
  plot(mygrid)
  plot(Norway,add=T)
}

plotZerror<-function(model){
  zSummary<-data.frame(model$summary[grepl("z",row.names(model$summary)),])
  zSummary$ParamNu <- as.character(sub(".*\\[([^][]+)].*", "\\1", row.names(zSummary)))
  zSummary$Site<-sapply(zSummary$ParamNu,function(x)strsplit(x,",")[[1]][1])
  zSummary$Year<-sapply(zSummary$ParamNu,function(x)strsplit(x,",")[[1]][2])
  zSummary$grid<-listlengthDF$grid[match(zSummary$Site,listlengthDF$siteIndex)]
  
  #predicted occupancy across all years
  zSummary_year<-ddply(zSummary,.(grid),summarise,prop=mean(mean),prop_sd=mean(sd),prop_cov=prop_sd/prop)
  
  #plot sd
  par(mfrow=c(2,1))
  mygrid[]<-0
  mygrid[zSummary_year$grid]<-zSummary_year$prop_sd
  plot(mygrid)
  plot(Norway,add=T)
  mygrid[]<-0
  #plot cov
  mygrid[zSummary_year$grid]<-zSummary_year$prop_cov
  plot(mygrid)
  plot(Norway,add=T)
}

#######################################################################

#Code to get environmental data:

#get accessibility map
#https://www.nature.com/articles/nature25181
setwd("C:/Users/diana.bowler/OneDrive - NINA/maps/accessibility/accessibility_to_cities_2015_v1.0")
access<-raster("accessibility_to_cities_2015_v1.0.tif")
out<-getEnvironData(access,mygrid)

#check data
hist(out$myraster)
summary(out$myraster)

#Add to the bugs data
bugs.data$access = out$myraster[match(listlengthDF$grid,out$grid)]

#plotting
plot(access)
plotBUGSData("access")

rm(access)

######################################################################################################

#habitats:
setwd("R:/GeoSpatialData/LandCover/Norway_Satveg_NORUT/Original/Satveg_deling_nd_partnere_09_12_2009/tiff")
library(raster)

#project other rasters onto this
habitatRasterTop <- raster("NNred25-30-t1.tif")#30 x 30 m
habitatRasterTop <- aggregate(habitatRasterTop,fact=10,fun=modal,na.rm=T)
habitatRasterBot <- raster("sn25_geocorr.tif")
habitatRasterBot <- aggregate(habitatRasterBot,fact=10,fun=modal,na.rm=T)
extent(habitatRasterTop)
extent(habitatRasterBot)

#merge each dataset
totalExtent <- c(246285,1342485,6414184,8014594)
habitatRasterTop<-extend(habitatRasterTop,extent(totalExtent))
habitatRasterBot<-extend(habitatRasterBot,extent(totalExtent))
origin(habitatRasterBot) <-origin(habitatRasterTop)
habitatRaster <- merge(habitatRasterTop,habitatRasterBot)
plot(habitatRaster)
#yeah!!!

#Set NAs for irrelevant habitats
habitatRaster[habitatRaster>24] <- NA
habitatRaster[habitatRaster==22] <- NA
habitatRaster[habitatRaster==0] <- NA
plot(habitatRaster)

#plot each grid, extract what????
#crop raster to Norway extent
myraster<-habitatRaster
rasterCRS<-crs(myraster)

#convert raster into points and convert to lon lat
myrasterDF<-as.data.frame(myraster,xy=T)
names(myrasterDF)[3]<-"myraster"
myrasterDF<-subset(myrasterDF,!is.na(myraster))
coordinates(myrasterDF)<-c("x","y")
proj4string(myrasterDF)<-rasterCRS
myrasterDF<-spTransform(myrasterDF,"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#get general grid
mygrid<-gridTemp
projection(mygrid)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0") 

#get myraster values per grid cell
mygrid[]<-1:ncell(mygrid)
variable<-extract(mygrid,myrasterDF)
mygrid[]<-0
mygrid[variable]<-1
myrasterDF<-data.frame(myrasterDF@data)
myrasterDF$grid<-variable

library(reshape2)
myrasterDF<-melt(table(myrasterDF$grid,myrasterDF$myraster))
names(myrasterDF)<-c("grid","raster","count")

#get rid of 0s
myrasterDF<-subset(myrasterDF,raster!=0)

#simplify habitat counts
myrasterDF$Forest<-0
myrasterDF$Forest[myrasterDF$raster%in%c(1:8)]<-myrasterDF$count[myrasterDF$raster%in%c(1:8)]
myrasterDF$Open<-0
myrasterDF$Open[myrasterDF$raster%in%c(9:21)]<-myrasterDF$count[myrasterDF$raster%in%c(9:21)]
myrasterDF$Top<-0
myrasterDF$Top[myrasterDF$raster%in%c(14,17)]<-myrasterDF$count[myrasterDF$raster%in%c(14,17)]
#Heather-rich alpine ridge vegetation, Fresh heather and dwarf-shrub communities
myrasterDF$Bottom<-0
myrasterDF$Bottom[myrasterDF$raster%in%c(23:24)]<-myrasterDF$count[myrasterDF$raster%in%c(23:24)]
#Agricultural areas, Cities and built-up areas
myrasterDF$Agriculture[myrasterDF$raster%in%c(23)]<-myrasterDF$count[myrasterDF$raster%in%c(23)]
#Agricultural areas
myrasterDF<-ddply(myrasterDF,.(grid),summarise,Forest=sum(Forest,na.rm=T),Open=sum(Open,na.rm=T),
                  Bottom=sum(Bottom,na.rm=T),Top=sum(Top,na.rm=T),Agriculture=sum(Agriculture,na.rm=T))


#Simplify into factors
myrasterDF$Habitat<-apply(myrasterDF,1,function(x)ifelse(as.numeric(x["Forest"])>as.numeric(x["Open"]),
                                                           "Forest","Open"))
myrasterDF$Habitat[myrasterDF$Top>150]<-"Top"
sum(is.na(myrasterDF$Habitat))#none!

#Add to the bugs data
bugs.data$habitat = myrasterDF$Habitat[match(listlengthDF$grid,myrasterDF$grid)]
bugs.data$Forest = myrasterDF$Forest[match(listlengthDF$grid,myrasterDF$grid)]
bugs.data$Open = myrasterDF$Open[match(listlengthDF$grid,myrasterDF$grid)]
bugs.data$Top = myrasterDF$Top[match(listlengthDF$grid,myrasterDF$grid)]
bugs.data$Agriculture = myrasterDF$Agriculture[match(listlengthDF$grid,myrasterDF$grid)]

table(bugs.data$habitat)
sum(is.na(bugs.data$habitat))
#132

unique(listlengthDF$grid[is.na(bugs.data$habitat)])
#130!

# get coordinates of grid cell with NA habitat
emptygrids<-unique(listlengthDF$grid[is.na(bugs.data$habitat)])
pts <- data.frame(xyFromCell(mygrid,emptygrids))
coordinates(pts) <- c("x","y")
proj4string(pts) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
pts<-spTransform(pts,crs(habitatRaster))
plot(habitatRaster)
plot(pts,col="red",add=T)#the points are in the sea!!!

rm(habitatRaster)
rm(myraster)

######################################################################################################

#average climatic conditions

#try the EuroLST dataset
setwd("R:/GeoSpatialData/Meteorology/Europe_EuroLST_bioclim")
#-BIO1: Annual mean temperature (°C*10): eurolst_clim.bio01.zip (MD5) 72MB
#-BIO2: Mean diurnal range (Mean monthly (max - min tem)): eurolst_clim.bio02.zip (MD5) 72MB
#-BIO3: Isothermality ((bio2/bio7)*100): eurolst_clim.bio03.zip (MD5) 72MB
#-BIO4: Temperature seasonality (standard deviation * 100): eurolst_clim.bio04.zip (MD5) 160MB
#-BIO5: Maximum temperature of the warmest month (°C*10): eurolst_clim.bio05.zip (MD5) 106MB
#-BIO6: Minimum temperature of the coldest month (°C*10): eurolst_clim.bio06.zip (MD5) 104MB
#-BIO7: Temperature annual range (bio5 - bio6) (°C*10): eurolst_clim.bio07.zip (MD5) 132MB
#-BIO10: Mean temperature of the warmest quarter (°C*10): eurolst_clim.bio10.zip (MD5) 77MB
#-BIO11: Mean temperature of the coldest quarter (°C*10): eurolst_clim.bio11.zip (MD5) 78MB

######
#bio1#
######
temp_bio1<-raster("eurolst_clim.bio01.tif")#res 250 m
temp_bio1<-aggregate(temp_bio1,fact=4,fun=mean)

#extract the data
out<-getEnvironData(temp_bio1,mygrid)

#check data
hist(out$myraster)
summary(out$myraster)

#Add to the bugs data
bugs.data$bio1 = out$myraster[match(listlengthDF$grid,out$grid)]
summary(bugs.data$bio1)
plotBUGSData("bio1")

rm(temp_bio1)

##############
#maximum temp#
##############
temp_bio5<-raster("eurolst_clim.bio05.tif")
temp_bio5<-aggregate(temp_bio5,fact=4,fun=mean)
out<-getEnvironData(temp_bio5,mygrid)
bugs.data$bio5 = out$myraster[match(listlengthDF$grid,out$grid)]
summary(bugs.data$bio5)
plotBUGSData("bio5")
rm(temp_bio5)

#minimum temp
temp_bio6<-raster("eurolst_clim.bio06.tif")
temp_bio6<-aggregate(temp_bio6,fact=4,fun=mean)
out<-getEnvironData(temp_bio6,mygrid)
bugs.data$bio6 = out$myraster[match(listlengthDF$grid,out$grid)]
summary(bugs.data$bio6)
plotBUGSData("bio6")

#####################################################################################################

#Examine correlations among bugs variables

varDF<-data.frame(do.call(cbind,bugs.data[9:length(names(bugs.data))]),stringsAsFactors = FALSE)
varDF[,-which(names(varDF)=="habitat")]<-sapply(varDF[,-which(names(varDF)=="habitat")],as.numeric)
save(varDF,file="varDF_missing_5km.RData")
library(GGally)
ggpairs(varDF)
table(varDF$habitat)

#bio1 and bio6 strongly related
#bio1 and elevation strongly related
#access and bio6 are related
#open and elevation are correlated
#bio1 and open are correlated
#top and open are correlated
#use habitat, bio 1 and bio 5

#plot in space
gridDF<-as.data.frame(gridTemp,xy=T)
varDF$grid<-bugs.data$grid
varDF<-merge(varDF,gridDF,by.x="grid",by.y="layer",all.x=T)

qplot(x,y,data=varDF,color=habitat)
qplot(x,y,data=varDF,color=Forest)
qplot(x,y,data=varDF,color=Open)
qplot(x,y,data=varDF,color=Top)
qplot(x,y,data=varDF,color=Agriculture)
qplot(x,y,data=varDF,color=bio1)
qplot(x,y,data=varDF,color=bio5)
qplot(x,y,data=varDF,color=bio6)

#which grid cells do we not have information for all these
NAgrid<-unique(c(varDF$grid[is.na(varDF$habitat)],varDF$grid[is.na(varDF$bio1)]))
NAgrid
varDF<-subset(varDF,!grid%in%NAgrid)
#replot using code above
#still good!

#####################################################################################################

setwd( "C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling")
#later loading
load("varDF_missing_5km.RData")
varDF$grid<-bugs.data$grid
NAgrid<-unique(c(varDF$grid[is.na(varDF$habitat)],varDF$grid[is.na(varDF$bio1)]))
varDF<-subset(varDF,!grid%in%NAgrid)
nrow(varDF)
length(unique(varDF$grid))#28102

#####################################################################################################

#also get alpine data

load("alpineData_5km.RData")
#1     2     3     4 
#14928  8310  2313   589

#plot it
qplot(x,y,data=alpineData,color=tree_line)
qplot(x,y,data=alpineData,color=elevation)
qplot(x,y,data=alpineData,color=alpine_habitat)
#look ok

#summarise it
alpineData$tree_line_position <- alpineData$tree_line-alpineData$elevation
alpineData$alpine_habitat1 <- sapply(alpineData$alpine_habitat,function(x)ifelse(x==1,1,0))
alpineData$alpine_habitat2 <- sapply(alpineData$alpine_habitat,function(x)ifelse(x==2,1,0))
alpineData$alpine_habitat3 <- sapply(alpineData$alpine_habitat,function(x)ifelse(x==3,1,0))
alpineData$alpine_habitat4 <- sapply(alpineData$alpine_habitat,function(x)ifelse(x==4,1,0))

#aggregate to site level
alpineData <- ddply(alpineData,.(site),summarise,
                    tree_line_position = median(tree_line_position,na.rm=T),
                    alpine_habitat1 = mean(alpine_habitat1,na.rm=T), 
                    alpine_habitat2 = mean(alpine_habitat2,na.rm=T), 
                    alpine_habitat3 = mean(alpine_habitat3,na.rm=T),
                    alpine_habitat4 = mean(alpine_habitat4,na.rm=T)) 

names(alpineData)[which(names(alpineData)=="site")]<-"grid"

#any more NAs??
moreNAs <- alpineData$grid[!complete.cases(alpineData)]
NAgrid<-append(NAgrid,moreNAs)
alpineData<-subset(alpineData,complete.cases(alpineData))
moreNAs<-setdiff(varDF$grid,alpineData$grid)
NAgrid<-unique(append(NAgrid,moreNAs))

#remove from both datasets
varDF<-subset(varDF,!grid%in%NAgrid)
alpineData<-subset(alpineData,!grid%in%NAgrid)
varDF<-merge(varDF,alpineData,by="grid",all.x=T,sort=FALSE)
summary(varDF)
length(unique(varDF$grid))#28081
save(varDF,file="varDF_allEnvironData_5km.RData")

######################################################################################################

#combine environmental and pop data:

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling")
load("varDF_allEnvironData_5km.RData")
varDF<-subset(varDF,!duplicated(grid))
listlengthDF<-merge(listlengthDF,varDF,by="grid")

#add indices
listlengthDF$site<-paste(listlengthDF$adm2,listlengthDF$grid,sep="-")
listlengthDF$siteIndex<-as.numeric(factor(listlengthDF$site))
listlengthDF$yearIndex<-as.numeric(factor(listlengthDF$YearCollected))

#order data by site and year
listlengthDF<-arrange(listlengthDF,yearIndex,adm2,siteIndex)

#extract site data
listlengthDF_SiteCovs<-subset(listlengthDF,!duplicated(grid))

####################################################################################################

#fit as glm with explanatory variables

occupancyGrid<-ddply(listlengthDF,.(grid),summarise,species=max(y,na.rm=T))
varDF<-merge(varDF,occupancyGrid,by="grid",all.x=T)
varDF$species[is.infinite(varDF$species)]<-NA
varDF2<-subset(varDF,!is.na(species))

summary(glm(species~alpine_habitat2,data=varDF2))
summary(glm(species~bio1,data=varDF2))
summary(glm(species~bio5,data=varDF2))
summary(glm(species~tree_line_position + I(tree_line_position^2),data=varDF2))
summary(glm(species~Top,data=varDF2))

#all together
summary(glm(species~bio1 + Open + tree_line_position + I(tree_line_position^2),data=varDF))
#all significant...
glm1<-glm(species~bio1 + Open + tree_line_position + I(tree_line_position^2),data=varDF2)
varDF2$fits<-predict(glm1,type="response")


#need to get x and y coords below
library(ggplot2)
ggplot(varDF2)+
  geom_point(aes(x,y,colour=fits),shape=15,size=rel(1),alpha=0.5)+
  scale_colour_gradient(low="steelblue",high="red")+
  geom_point(data=subset(varDF2,species==1),aes(x,y))

###################################################################################################

#Boosted regression tree

library(dismo)
library(gbm)

str(varDF)

brt1 <- gbm.step(data=varDF2, gbm.x = 4:15, gbm.y = 16,family = "bernoulli")

summary(brt1)
                                var    rel.inf
tree_line_position tree_line_position 66.5219809
bio1                             bio1 14.6543613
Open                             Open  6.0688441
bio5                             bio5  2.5385401
Agriculture               Agriculture  2.3005598
bio6                             bio6  2.1910442
Forest                         Forest  1.4662198
alpine_habitat3       alpine_habitat3  1.4248853
alpine_habitat1       alpine_habitat1  1.2022284
Top                               Top  1.0933468
alpine_habitat4       alpine_habitat4  0.3625442
alpine_habitat2       alpine_habitat2  0.175444


#plot main effects
gbm.plot(brt1, n.plots=12, write.title = TRUE)
gbm.plot.fits(brt1)

#interactions?
find.int <- gbm.interactions(brt1)
find.int$interactions#none!
find.int$rank.list

####################################################################################################

#for BUGS

bugs.data <- list(nsite = length(unique(listlengthDF$siteIndex)),
                  nyear = length(unique(listlengthDF$yearIndex)),
                  nvisit = nrow(listlengthDF),
                  grid = listlengthDF$grid,
                  site = listlengthDF$siteIndex,
                  year = listlengthDF$yearIndex,
                  L = listlengthDF$L,
                  nuSpecies = listlengthDF$nuSpecies,
                  y = listlengthDF$y,
                  #add an adm effect
                  adm = listlengthDF$adm2,
                  n.adm = length(unique(listlengthDF$adm2)),
                  habitat = ifelse(listlengthDF_SiteCovs$habitat=="Forest",1,0),
                  bio1 = listlengthDF_SiteCovs$bio1,
                  bio1_2 = listlengthDF_SiteCovs$bio1^2,
                  bio6 = listlengthDF_SiteCovs$bio6,
                  bio5 = listlengthDF_SiteCovs$bio5,
                  bio5_2 = listlengthDF_SiteCovs$bio5^2,
                  forest = listlengthDF_SiteCovs$Forest,
                  open = listlengthDF_SiteCovs$Open,
                  top = log(listlengthDF_SiteCovs$Top+1),
                  alpine_habitat1 = listlengthDF_SiteCovs$alpine_habitat1,
                  alpine_habitat2 = log(listlengthDF_SiteCovs$alpine_habitat2+1),
                  alpine_habitat3 = log(listlengthDF_SiteCovs$alpine_habitat3+1),
                  alpine_habitat4 = log(listlengthDF_SiteCovs$alpine_habitat4+1),
                  tree_line_position = listlengthDF_SiteCovs$tree_line_position,
                  tree_line_position2 = listlengthDF_SiteCovs$tree_line_position^2)

#alpine_habitat:
#1= Open lowland, 
#2 = Low alpine zone, 
#3 = intermediate alpine zone, 
#4 = high alpine zone 

#need to specify initial values
library(reshape2)
zst <- acast(listlengthDF, siteIndex~YearCollected, value.var="y",fun=max,na.rm=T)
zst [is.infinite(zst)] <- 0
inits <- function(){list(z = zst)}

#get BUGS functions
source('C:/Users/diana.bowler/OneDrive - NINA/methods/models/bugsFunctions.R')

##########################################################################################

#run model with random effects

#specify parameters to monitor
params <- c("int","psi.fs","dtype.p","mu.lp","a","random.adm.sd")

#specify model structure
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "BUGS_sparta_random_missing.txt", n.thin=10,
             n.chains=3, n.burnin=300,n.iter=1000,parallel=T)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_OM_random_missing_5km.RData")

#########################################################################################

#run model including explanatory variables

#specify parameters to monitor
params <- c("int","psi.fs","dtype.p","mu.lp","beta","a","det.elevation","det.elevation2")

##########################################################################################

#specify model structure
bugs.data$occDM <- model.matrix(~ scale(bugs.data$bio1)+
                                  scale(bugs.data$tree_line_position) + scale(bugs.data$tree_line_position2)+
                                  scale(bugs.data$open))[,-1]

bugs.data$n.covs <- ncol(bugs.data$occDM)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "BUGS_sparta_variables_missing.txt", n.thin=10,
             n.chains=3, n.burnin=5000,n.iter=10000,parallel=T)

print(out1,2)
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_OM_variables_missing_5km.RData")
colnames(bugs.data$occDM)

##########################################################################################

#plot occupancy model
out2<-update(out1,parameters.to.save=c("muZ","z"),n.iter=100)
plotZ(out2)
plotZerror(out2)

##########################################################################################
 
#also using jagam

#add to the dataset, the coordinates of the grid
gridDF<-as.data.frame(gridTemp,xy=T)
varDF<-merge(varDF,gridDF,by.x="grid",by.y="layer",all.x=T)
 
#get info on whether there was an observation within each grid
occupancyGrid<-ddply(occupancyDF,.(grid),summarise,species=max(y))
varDF<-merge(varDF,occupancyGrid,by="grid",all.x=T)

#fit as gam
library(mgcv)
gam1 <- gam(species~ 1 + s(x, y,k=10), 
                    data=varDF, 
                    family="binomial")

#plot it
varDF$fits<-gam1$fitted.values
library(ggplot2)
ggplot(varDF)+
  geom_point(aes(x,y,colour=fits),shape=15,size=rel(4))+
  scale_colour_gradient(low="blue",high="red")+
  geom_point(data=subset(varDF,species==1),aes(x,y))


#Use the JAGAM function to get the BUGS code for the spline
#in the BUGS model
jags.ready <- jagam(species ~ 1 + s(x, y,k=10), 
                    data = varDF, 
                    family="binomial", 
                    sp.prior="log.uniform",
                    file="jagam.txt")

#get the data bits we need from jags data
bugs.data$X = jags.ready$jags.data$X
bugs.data$S1 = jags.ready$jags.data$S1
bugs.data$zero = jags.ready$jags.data$zero

#specify parameters to monitor
params <- c("dtype.p","mu.lp","rho")

#run model
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "BUGS_sparta_JAGAM.txt", n.thin=nt,
             n.chains=3, n.burnin=10000,n.iter=50000)

traceplot(out1)
print(out1,2)

#########################################################################################
