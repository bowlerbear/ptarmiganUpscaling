#get norway##############################################################

library(raster)
library(maptools)
data(wrld_simpl)
Norway<- subset(wrld_simpl,NAME=="Norway")
plot(Norway)

########################################################################

#choose modelling resolution:

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

########################################################################

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
  grid<-gridTemp
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

#######################################################################

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

#############################################################################################

#get human density map???

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
elevation<-aggregate(elevation,fact=50,fun=mean)

#extract the data
out<-getEnvironData(elevation,mygrid)

#check data
hist(out$myraster)
summary(out$myraster)

#Add to the bugs data
bugs.data$elevation = out$myraster[match(listlengthDF$grid,out$grid)]
summary(bugs.data$elevation)

#plotting
plot(elevation)
plotBUGSData("elevation")

rm(elevation)
######################################################################################################

#habitats:
setwd("R:/GeoSpatialData/LandCover/Norway_Satveg_NORUT/Original/Satveg_deling_nd_partnere_09_12_2009/tiff")
library(raster)
habitatRaster<-raster("NNred25-30-t1.tif")#30 x 30 m)
habitatRaster[habitatRaster>24]<-NA
habitatRaster[habitatRaster==22]<-NA
habitatRaster2<-aggregate(habitatRaster,fact=10,fun=modal,na.rm=T)
plot(habitatRaster2)

#plot each grid, extract what????
#crop raster to Norway extent
myraster<-habitatRaster2
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
myrasterDF<-data.frame(myrasterDF@data)
myrasterDF$grid<-variable
library(reshape2)
myrasterDF2<-melt(table(myrasterDF$grid,myrasterDF$myraster))
names(myrasterDF2)<-c("grid","raster","count")

#get rid of 0s
myrasterDF2<-subset(myrasterDF2,raster!=0)

#simplify habitat counts
myrasterDF2$Forest<-0
myrasterDF2$Forest[myrasterDF2$raster%in%c(1:8)]<-myrasterDF2$count[myrasterDF2$raster%in%c(1:8)]
myrasterDF2$Open<-0
myrasterDF2$Open[myrasterDF2$raster%in%c(9:21)]<-myrasterDF2$count[myrasterDF2$raster%in%c(9:21)]
myrasterDF2$Top<-0
myrasterDF2$Top[myrasterDF2$raster%in%c(14,17)]<-myrasterDF2$count[myrasterDF2$raster%in%c(14,17)]
#Heather-rich alpine ridge vegetation, Fresh heather and dwarf-shrub communities
myrasterDF2$Bottom<-0
myrasterDF2$Bottom[myrasterDF2$raster%in%c(23:24)]<-myrasterDF2$count[myrasterDF2$raster%in%c(23:24)]
#Agricultural areas, Cities and built-up areas
myrasterDF2<-ddply(myrasterDF2,.(grid),summarise,Forest=sum(Forest),Open=sum(Open),Bottom=sum(Bottom),Top=sum(Top))


#Simplify into factors
myrasterDF2$Habitat<-apply(myrasterDF2,1,function(x)ifelse(as.numeric(x["Forest"])>as.numeric(x["Open"]),
                                                           "Forest","Open"))
myrasterDF2$Habitat[myrasterDF2$Top>150]<-"Top"
sum(is.na(myrasterDF2$Habitat))#none!

#Add to the bugs data
bugs.data$habitat = myrasterDF2$Habitat[match(listlengthDF$grid,myrasterDF2$grid)]
table(bugs.data$habitat)
sum(is.na(bugs.data$habitat))

#855

rm(habitatRaster)
rm(habitatRaster2)

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
temp_bio1<-raster("eurolst_clim.bio01.tif")
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

varDF<-do.call(cbind,bugs.data[8:13])


######################################################################################################

#basic model

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

##########################################################################################

#run model including explanatory variables

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

##########################################################################################