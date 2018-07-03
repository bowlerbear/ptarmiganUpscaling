#using a m grid
equalM<-"+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

#get norway##############################################################

library(raster)
library(maptools)
library(rgeos)
data(wrld_simpl)
Norway<- subset(wrld_simpl,NAME=="Norway")
Norway <- gBuffer(Norway,width=1)
Norway<-spTransform(Norway,crs(equalM))
plot(Norway)

########################################################################

#create grid
newres=5000#5 km grid
mygrid<-raster(extent(projectExtent(Norway,equalM)))
res(mygrid)<-newres
mygrid[]<-1:ncell(mygrid)
plot(mygrid)
gridTemp<-mygrid
plot(Norway,add=T)

#########################################################################

#identify all grid cells within the buffer region

library(rgeos)
mygridMask<-mask(mygrid,Norway)
mygridMaskDF<-as.data.frame(mygridMask,xy=T)
head(mygridMaskDF)
plot(mygridMask)

#########################################################################

#get info on administrative names for the grid
library(rgdal)

#make both spatial objects in the same crs
NorwayADM<-readOGR(dsn="C:/Users/diana.bowler/OneDrive - NINA/Alpine/NOR_adm",layer="NOR_adm1")
NorwayADM<-spTransform(NorwayADM,crs(equalM))
mygridPoints<-mygridMaskDF
coordinates(mygridPoints)<-c("x","y")
proj4string(mygridPoints)<-crs(NorwayADM)

#check they overlap
plot(mygridPoints)
plot(NorwayADM,add=T,col="red")

#extract the data
myAdm<-over(mygridPoints,NorwayADM)
myAdm$grid<-mygridPoints$layer
myAdm$VARNAME_1<-as.character(myAdm$VARNAME_1)
myAdm$VARNAME_1[is.na(myAdm$VARNAME_1)]<-"outside"
myAdm<-subset(myAdm,!is.na(grid))
table(myAdm$VARNAME_1)

#check the results
mygrid[]<-0
mygrid[myAdm$grid]<-as.numeric(as.factor(myAdm$VARNAME_1))
plot(mygrid)#looks good!

########################################################################

#general functions:

#get function to extract raster into for these grids:

getEnvironData<-function(myraster,mygridTemp){
  require(maptools)
  require(plyr)
  
  #crop raster to Norway extent
  rasterCRS<-crs(myraster)
  NorwayB<-spTransform(Norway,rasterCRS)
  myraster<-crop(myraster,extent(NorwayB))
  
  #convert raster into points and convert
  myrasterDF<-as.data.frame(myraster,xy=T)
  names(myrasterDF)[3]<-"myraster"
  coordinates(myrasterDF)<-c("x","y")
  proj4string(myrasterDF)<-rasterCRS
  myrasterDF<-spTransform(myrasterDF,crs(equalM))
  
  #get general grid
  mygrid<-gridTemp
  projection(mygrid)<- equalM
  
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

#Code to get environmental data:

#get accessibility map
#https://www.nature.com/articles/nature25181
setwd("C:/Users/diana.bowler/OneDrive - NINA/maps/accessibility/accessibility_to_cities_2015_v1.0")
access<-raster("accessibility_to_cities_2015_v1.0.tif")
out<-getEnvironData(access,mygrid)

#check data
hist(out$myraster)
summary(out$myraster)
outAccess<-out

rm(access)
rm(out)

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
habitatRaster[habitatRaster==0] <- NA
plot(habitatRaster)

#plot each grid, extract what????
#crop raster to Norway extent
myraster<-habitatRaster
rasterCRS<-crs(myraster)

#convert raster into points and convert 
myrasterDF<-as.data.frame(myraster,xy=T)
names(myrasterDF)[3]<-"myraster"
myrasterDF<-subset(myrasterDF,!is.na(myraster))
coordinates(myrasterDF)<-c("x","y")
proj4string(myrasterDF)<-rasterCRS
myrasterDF<-spTransform(myrasterDF,crs(equalM))
                        
#get general grid
mygrid<-gridTemp
projection(mygrid)<-CRS(equalM) 

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
myrasterDF$Open[myrasterDF$raster%in%c(9,10,12:21)]<-myrasterDF$count[myrasterDF$raster%in%c(9,10,12:21)]
myrasterDF$Top<-0
myrasterDF$Top[myrasterDF$raster%in%c(14,17)]<-myrasterDF$count[myrasterDF$raster%in%c(14,17)]
#Heather-rich alpine ridge vegetation, Fresh heather and dwarf-shrub communities
myrasterDF$PrefOpen[myrasterDF$raster%in%c(10,17,18)]<-myrasterDF$count[myrasterDF$raster%in%c(10,17,18)]
myrasterDF$PrefClosed[myrasterDF$raster%in%c(6,7)]<-myrasterDF$count[myrasterDF$raster%in%c(6,7)]
myrasterDF$Bottom<-0
myrasterDF$Bottom[myrasterDF$raster%in%c(23:24)]<-myrasterDF$count[myrasterDF$raster%in%c(23:24)]
#Agricultural areas, Cities and built-up areas
myrasterDF$Agriculture[myrasterDF$raster%in%c(23)]<-myrasterDF$count[myrasterDF$raster%in%c(23)]
#Agricultural areas

#aggregate
myrasterDF<-ddply(myrasterDF,.(grid),summarise,Forest=sum(Forest,na.rm=T),
                  Open=sum(Open,na.rm=T),
                  Bottom=sum(Bottom,na.rm=T),
                  Top=sum(Top,na.rm=T),
                  PrefOpen=sum(PrefOpen,na.rm=T),
                  PrefClosed=sum(PrefClosed,na.rm=T),
                  Agriculture=sum(Agriculture,na.rm=T))


#Simplify into factors
myrasterDF$Habitat<-apply(myrasterDF,1,
                          function(x)
                            ifelse(as.numeric(x["Forest"])>as.numeric(x["Open"]),
                                                         "Forest","Open"))
myrasterDF$Habitat[myrasterDF$Top>150]<-"Top"

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

out_Bio1<-out
rm(temp_bio1)
rm(out)

##############
#maximum temp#
##############
temp_bio5<-raster("eurolst_clim.bio05.tif")
temp_bio5<-aggregate(temp_bio5,fact=4,fun=mean)
out<-getEnvironData(temp_bio5,mygrid)

out_Bio5<-out
rm(temp_bio5)
rm(out)

#minimum temp
temp_bio6<-raster("eurolst_clim.bio06.tif")
temp_bio6<-aggregate(temp_bio6,fact=4,fun=mean)
out<-getEnvironData(temp_bio6,mygrid)

out_Bio6<-out
rm(temp_bio6)
rm(out)

###################################################################
#combine all

#temp data
names(out_Bio1)[2]<-"bio1"
names(out_Bio5)[2]<-"bio5"
names(out_Bio6)[2]<-"bio6"
out_Bio<-cbind(out_Bio1,bio5=out_Bio5[,2],bio6=out_Bio6[,2])

#combine others
varDF<-merge(out_Bio,myrasterDF,by="grid")
varDF<-merge(varDF,myAdm,by="grid",all.x=T)
varDF<-subset(varDF,is.finite(bio1))

save(varDF,file="varDF_missing_5km.RData")

#####################################################################################################

#Examine correlations among bugs variables
library(GGally)
ggpairs(varDF[,2:11])
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
NAgrid<-unique(c(varDF$grid[is.na(varDF$Habitat)],varDF$grid[is.na(varDF$bio1)]))
varDF<-subset(varDF,!grid%in%NAgrid)
nrow(varDF)
length(unique(varDF$grid))#16064

#check coverage
mygrid[]<-0
mygrid[varDF$grid]<-varDF$bio1
plot(mygrid)
#looks good
mygrid[]<-0
mygrid[varDF$grid]<-varDF$PrefOpen
plot(mygrid)
#ok
mygrid[]<-0
mygrid[varDF$grid]<-varDF$bio6
plot(mygrid)
#good!

#####################################################################################################

#also get alpine data

load("alpineData_5kmGrid.RData")

#plot it
library(ggplot2)
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
library(plyr)
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
length(unique(varDF$grid))#15636
mygrid[]<-0
mygrid[varDF$grid]<-varDF$bio6
plot(mygrid)

save(varDF,file="varDF_allEnvironData_5km.RData")
########################################################################################