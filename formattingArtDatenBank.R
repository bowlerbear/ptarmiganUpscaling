#get norway##############################################################

library(maptools)
data(wrld_simpl)
Norway<- subset(wrld_simpl,NAME=="Norway")
plot(Norway)

########################################################################

#choose modelling resolution
newres=1
newres=0.5

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

#overlay to a grid
#make 10 km (0.1) x 10(0.1) km grid
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
gridSummary<-ddply(lirype@data,.(grid,YearCollected),summarise,nuRecs=length(NorskNavn))
sort(unique(gridSummary$YearCollected))
par(mfrow=c(3,4))
for(i in 2007:2017){
  mygrid[]<-0
  mygrid[gridSummary$grid[gridSummary$YearCollected==i]]<-gridSummary$nuRecs[gridSummary$YearCollected==i]
  plot(mygrid,main=i)
}

#how many records per site are there, per year
hist(gridSummary$nuRecs)
gridSummary$RepeatedVisits<-sapply(gridSummary$nuRecs,function(x)ifelse(x>1,1,0))
table(gridSummary$RepeatedVisits)
#0   1 
#121 705

##########################################################################

#get data for all birds - as an effort layer

tdir<-"C:/Users/diana.bowler/OneDrive - NINA/Alpine/SpeciesMapServices/Birds"
allbirds <- read.delim(paste(tdir,"DataExportFraArtskart.txt",sep="/"),skipNul=T,dec=",")
unique(allbirds$Class)
summary(allbirds$CoordinatePrecision)
allbirds <- subset(allbirds, !CoordinatePrecision>10000)
coordinates(allbirds)<-c("Longitude","Latitude")
plot(Norway)
plot(allbirds,add=T)
#all points look pretty good!! few in the sea...

#overlay to a grid
#make 10 km (0.1) x 10(0.1) km grid
library(raster)
mygrid<-raster(extent(Norway))
res(mygrid)<-newres
mygrid[]<-1:ncell(mygrid)

#get grid number for each point
allbirds$grid<-extract(mygrid,allbirds)

#get number of points per grid cell
library(plyr)
gridSummary<-ddply(allbirds@data,.(grid),summarise,nuRecs=length(NorskNavn))
mygrid[]<-0
mygrid[gridSummary$grid]<-gridSummary$nuRecs
length(gridSummary$grid)==length(gridSummary$nuRecs)
plot(mygrid)

#plot for each year
allbirds<-subset(allbirds,YearCollected>2006&YearCollected<2018)
gridSummary<-ddply(allbirds@data,.(grid,YearCollected),summarise,nuRecs=length(NorskNavn))
sort(unique(gridSummary$YearCollected))
par(mfrow=c(3,4))
for(i in 2007:2017){
  mygrid[]<-0
  mygrid[gridSummary$grid[gridSummary$YearCollected==i]]<-gridSummary$nuRecs[gridSummary$YearCollected==i]
  plot(mygrid,main=i)
}

#get list length per cell
allbirds<-subset(allbirds,YearCollected>2006&YearCollected<2018)
allbirds@data$Species<-apply(allbirds@data,1,function(x)paste(x["Genus"],x["Species"],sep=" "))
gridSummary<-ddply(allbirds@data,.(grid,YearCollected),summarise,nuRecs=length(unique(Species)))
sort(unique(gridSummary$YearCollected))
par(mfrow=c(3,4))
for(i in 2007:2017){
  mygrid[]<-0
  mygrid[gridSummary$grid[gridSummary$YearCollected==i]]<-gridSummary$nuRecs[gridSummary$YearCollected==i]
  plot(mygrid,main=i)
}

#########################################################################

#Run occupancy model for willow ptarmigan -using the sparta base code

#add list length of all birds on each sampling visit
listlengthDF<-ddply(allbirds@data,.(YearCollected,MonthCollected,DayCollected,grid),summarise,L = length(Species))
listlengthDF$visit<-paste(listlengthDF$YearCollected,listlengthDF$MonthCollected,
                          listlengthDF$DayCollected,listlengthDF$grid,sep="-")
listlengthDF$siteIndex<-as.numeric(factor(listlengthDF$grid))
listlengthDF$yearIndex<-as.numeric(factor(listlengthDF$YearCollected))
listlengthDF<-subset(listlengthDF,!is.na(grid))

#sort out occupancy matrix of ptarmigan
occupancyDF<-data.frame(visit=listlengthDF$visit,year=listlengthDF$YearCollected,grid=listlengthDF$grid)
lirype$visit<-paste(lirype$YearCollected,lirype$MonthCollected,
                    lirype$DayCollected,lirype$grid,sep="-")

occupancyDF$y<-sapply(occupancyDF$visit,function(x)ifelse(x%in%lirype$visit,1,0))
sum(occupancyDF$y)

bugs.data <- list(nsite = length(unique(listlengthDF$siteIndex)),
                  nyear = length(unique(listlengthDF$yearIndex)),
                  nvisit = nrow(listlengthDF),
                  site = listlengthDF$siteIndex,
                  year = listlengthDF$yearIndex,
                  L = listlengthDF$L - median(listlengthDF$L),
                  y = occupancyDF$y)


#specify parameters to monitor
source('C:/Users/diana.bowler/OneDrive - NINA/methods/models/bugsFunctions.R')
params <- c("psi.fs","dtype.p","mu.lp")

#need to specify initial values
library(reshape2)
zst <- acast(occupancyDF, grid~year, value.var="y",fun=max)
zst [is.infinite(zst)] <- 0
inits <- function(){list(z = zst)}

#run model
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "BUGS_sparta.txt", n.thin=nt,
             n.chains=3, n.burnin=500,n.iter=10000)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_0.5_BUGS_sparta.RData")
print(out1,2)

#plot time series of occupancy
psiSummary<-data.frame(out1$summary[grepl("psi.fs",row.names(out1$summary)),])
psiSummary$Year <- 2007:2017
ggplot(psiSummary)+
  geom_line(aes(x=Year,y=mean))+
  geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.5)+
  theme_bw()
ggsave("ts.png")

#plot occupancy model
out2<-update(out1,parameters.to.save="z",n.iter=2000)
zSummary<-data.frame(out2$summary[grepl("z",row.names(out2$summary)),])
zSummary$ParamNu <- as.character(sub(".*\\[([^][]+)].*", "\\1", row.names(zSummary)))
zSummary$Site<-sapply(zSummary$ParamNu,function(x)strsplit(x,",")[[1]][1])
zSummary$Year<-sapply(zSummary$ParamNu,function(x)strsplit(x,",")[[1]][2])
zSummary$grid<-listlengthDF$grid[match(zSummary$Site,listlengthDF$siteIndex)]

#predicted occupancy across all years
zSummary_year<-ddply(zSummary,.(grid),summarise,prop=mean(mean))
mygrid[]<-0
mygrid[zSummary_year$grid]<-zSummary_year$prop
plot(mygrid)
plot(Norway,add=T)
ggsave("occupancyMap_0.5.png")

########################################################################

#get function to extract raster into for these grids:

getEnvironData<-function(myraster,mygridTemp){
  require(maptools)
  #convert into points
  rasterCRS<-crs(myraster)
  Norway<-spTransform(Norway,rasterCRS)
  myraster<-crop(myraster,extent(Norway))
  myrasterDF<-as.data.frame(myraster,xy=T)
  names(myrasterDF)[3]<-"myraster"
  coordinates(myrasterDF)<-c("x","y")
  proj4string(myrasterDF)<-rasterCRS
  myrasterDF<-spTransform(myrasterDF,"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
  grid<-gridTemp
  projection(mygrid)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0") 
  
  #get myraster values per grid cell
  mygrid[]<-1:ncell(mygrid)
  variable<-extract(mygrid,myrasterDF)
  myrasterDF<-data.frame(myrasterDF@data)
  myrasterDF$grid<-variable
  myrasterDF<-ddply(myrasterDF,.(grid),summarise,myraster=mean(myraster,na.rm=T))
  return(myrasterDF)
}

#######################################################################
#get accessibility map
#https://www.nature.com/articles/nature25181
setwd("C:/Users/diana.bowler/OneDrive - NINA/maps/accessibility/accessibility_to_cities_2015_v1.0")
library(raster)
library(maptools)
access<-raster("accessibility_to_cities_2015_v1.0.tif")
out<-getEnvironData(access,mygrid)

#check data
hist(out$myraster)
summary(out$myraster)

#Add to the bugs data
bugs.data$access = out$myraster[match(listlengthDF$grid,out$grid)]

#Fill in blanks with mean for the moment...
bugs.data$access[is.na(bugs.data$access)]<-mean(bugs.data$access,na.rm=T)

#run model including the accessibility term

#specify parameters to monitor
source('C:/Users/diana.bowler/OneDrive - NINA/methods/models/bugsFunctions.R')
params <- c("psi.fs","dtype.p","mu.lp","beta.access")

#need to specify initial values
library(reshape2)
zst <- acast(occupancyDF, grid~year, value.var="y",fun=max)
zst [is.infinite(zst)] <- 0
inits <- function(){list(z = zst)}

#run model
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/SpeciesMapServices")
out1 <- jags(bugs.data, inits=inits, params, "BUGS_sparta_access.txt", n.thin=nt,
             n.chains=3, n.burnin=500,n.iter=2000)

print(out1,2)

#############################################################################################

#get human density map

#Europe_GEOSTAT_1km2_population
setwd("R:/GeoSpatialData/Population_demography/Europe_GEOSTAT_1km2_population/Original")
library(rgdal)
popdensity<-readOGR(dsn=getwd(),layer="Grid_ETRS89_LAEA_1K_ref_GEOSTAT_POP_2011_V2_0_1")
#plot(popdensity)

#read the data file
popData<-read.csv("GEOSTAT_grid_POP_1K_2011_V2_0_1.csv")
#crop to Norway
popData<-subset(popData,CNTR_CODE=="NO")
summary(popData$TOT_P)

#subset the polygon grid
popdensity<-subset(popdensity,GRD_ID%in%popData$GRD_ID)
plot(popdensity)
#get the coordinates for each grid
out<-sapply(popdensity@polygons, function(x) x@labpt)
grids<-data.frame(x=out[1,],y=out[2,],GRD_ID=popData$GRD_ID)

#add coords to pop data file
popData<-merge(popData,grids,by="GRD_ID")
qplot(x,y,data=popData,colour=log(TOT_P))

#need to fill in the blank grids
full.grid<-expand.grid(x=seq(from=min(popData$x),to=max(popData$x),by=1000),
                       y=seq(from=min(popData$y),to=max(popData$y),by=1000))
full.grid$values<-popData$TOT_P[match(interaction(full.grid$x,full.grid$y),
                                      interaction(popData$x,popData$y))]
full.grid$values[is.na(full.grid$values)]<-0
qplot(x,y,data=full.grid,colour=values)

#turn it into a raster
full.grid$values<-log(full.grid$values+1)
rasterPop<-rasterFromXYZ(full.grid)
projection(rasterPop)<-"+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
plot(rasterPop)

#mask to Norway
Norway<-spTransform(Norway,CRS(projection(rasterPop)))
rasterPop<-crop(rasterPop,extent(Norway))
rasterPop<-mask(rasterPop,Norway)
plot(rasterPop)

#######################################################################################################

#elevation map
tdir<-"R:/GeoSpatialData/Elevation/Europe_EU_DSM_25m"#there are higher resolution ones
elevation<-raster(paste(tdir,"eudem_dem_3035_europe.tif",sep="/"))
#its really high resolution to lower it
rasterCRS<-crs(elevation)
Norway<-spTransform(Norway,rasterCRS)
elevation<-crop(elevation,extent(Norway))
elevation<-aggregate(elevation,fact=100,fun=mean)

#extract the data
out<-getEnvironData(elevation,mygrid)

#check data
hist(out$myraster)
summary(out$myraster)

#Add to the bugs data
bugs.data$elevation = out$myraster[match(listlengthDF$grid,out$grid)]

######################################################################################################

#habitats
setwd("R:/GeoSpatialData/Habitats_biotopes/Europe_EU_DSM_25m")
Norway_Naturbase??

######################################################################################################

#climate
setwd("R:/GeoSpatialData/Habitats_biotopes/Europe_EU_DSM_25m")
Norway_Naturbase??

######################################################################################################
