library(tidyverse)
library(sp)
library(rgeos)
library(raster)
library(maptools)

myfolder <- "Data"
equalM<-"+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

source('generalFunctions.R')

### get norway ##############################################

data(wrld_simpl)
Norway <- subset(wrld_simpl,NAME=="Norway")
NorwayOrig <- Norway
Norway <- gBuffer(Norway,width=1)
Norway <- spTransform(Norway,crs(equalM))
plot(Norway)

#create grid
newres = 5000#5 km grid
mygrid <- raster(extent(projectExtent(Norway,equalM)))
res(mygrid) <- newres 
mygrid[] <- 1:ncell(mygrid)
plot(mygrid)
gridTemp <- mygrid
plot(Norway,add=T)

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
saveRDS(Lines_spatial,file="data/Lines_spatial.rds")

### create buffers ##############################################################

#convert to utm 
spTemp <- spTransform(Lines_spatial,crs(equalM))

#get centroids
lineCentres <- gCentroid(spTemp,byid=TRUE)@coords

#put buffers around each one and make into a polygons
makePolygon <- function(x){
  myCentre <- data.frame(t(x))
  coordinates(myCentre) <- c("x","y")
  myPoly <- gBuffer(myCentre,width=2821)#sqrt(25/pi)
  return(myPoly)
}
p <- apply(lineCentres,1,makePolygon)

#make into spatial polygons object
library(purrr)
allPolys <- list(p, makeUniqueIDs = T) %>% 
  flatten() %>% 
  do.call(rbind, .)
plot(allPolys)

#make into a spatial object
Polys_spatial <- SpatialPolygonsDataFrame(Sr=allPolys, 
                                          data=mydata,
                                          match.ID=FALSE)

proj4string(Polys_spatial) <- equalM
saveRDS(Polys_spatial,file="data/Polys_spatial.rds")

### plot ####################################################

myLines <- Lines_spatial$LinjeID

Lines_spatial <- spTransform(Lines_spatial,crs(equalM))
plot(subset(Polys_spatial,LinjeID==122))
plot(subset(Lines_spatial,LinjeID==122),add=T)

for(i in 1:length(myLines)){
  
  line <- myLines[i]

  png(filename=paste0("plots/obs/circles/LinjeID",line,".png"))
  plot(subset(Polys_spatial,LinjeID==line))
  plot(subset(Lines_spatial,LinjeID==line),add=T)
  dev.off()

}

### get main grid of each polygon ##########################

centreDF <- data.frame(LinjeID = spTemp$LinjeID,
                       x = lineCentres[,1],
                       y = lineCentres[,2])

coordinates(centreDF) <- c("x","y")
proj4string(centreDF) <- CRS(equalM)

plot(mygrid)
plot(Norway,add=T)
plot(centreDF,add=T)

centreDF$grid <- extract(mygrid,centreDF)
centreDF@data$x <- centreDF@coords[,1]
centreDF@data$y <- centreDF@coords[,2]
saveRDS(centreDF@data,"data/lines_to_grids.rds")

# #### all grids for each line #################################
 
#why do some lines not overlap with the focal grids??
plot(NorwayOrig)
plot(Lines_spatial,add=T,col="red")

#for those missing (see 'mapping_lines_to_grids.R')
plot(NorwayOrig)
plot(subset(Lines_spatial,LinjeID %in% missing),add=T,col="red")
#these are in the north - in islands and in Sweden (not in focal grids)

#put buffer around each line
Lines_spatial <- spTransform(Lines_spatial,crs(equalM))

#Lines_spatial_buffer <- gBuffer(Lines_spatial,width = 5000, byid=TRUE)
Lines_spatial_buffer <- gBuffer(Lines_spatial,width = 15000, byid=TRUE)

Polys_grids <- extract(mygrid,Lines_spatial_buffer,df=T)
head(Polys_grids)

#add to data frame
Lines_spatial_buffer$ID <- 1:nrow(Lines_spatial_buffer)
Polys_grids$LinjeID <- Lines_spatial_buffer$LinjeID[match(Polys_grids$ID,Lines_spatial_buffer$ID)]
names(Polys_grids)[2]<- "grid"

#saveRDS(Polys_grids,file = "data/lineBuffers_to_grids_5km.rds")
saveRDS(Polys_grids,file = "data/lineBuffers_to_grids_15km.rds")

### get environ data for all circles ########################

Polys_spatial <- readRDS("data/Polys_spatial.rds")

### distance to coast ######################################

#get coastline data -very smooth
sldf_coast <- rnaturalearth::ne_coastline()
plot(sldf_coast)
coast <- sf::st_as_sf(sldf_coast)

#get better coastline data? too slow
# coastNorway <- getData('GADM', country='Norway', level=0)
# coastSweden <- getData('GADM', country='Sweden', level=0)
# coastFinland <- getData('GADM', country='Finland', level=0)
# coastRussia <- getData('GADM', country='Russia', level=0)
# coast <- sf::st_union(sf::st_as_sf(coastNorway), sf::st_as_sf(coastSweden))
# coast <- sf::st_union(coast, sf::st_as_sf(coastFinland))
# coast <- sf::st_union(coast, sf::st_as_sf(coastRussia))

#get points
lineCentres_xy <- lineCentres
lineCentres <- sf::st_as_sf(data.frame(lineCentres), 
                            coords = c("x", "y"), crs = equalM)
lineCentres <- sf::st_transform(lineCentres,sf::st_crs(coast))
ggplot()+
  geom_sf(data=coast)+
  geom_sf(data=lineCentres,color="red")

#get distance between them
lineCentres <- sf::st_coordinates(lineCentres)
dist <- geosphere::dist2Line(p = as.matrix(lineCentres), 
                             line = sldf_coast)

#pack into data frame
lineCentres_xy <- as.data.frame(lineCentres_xy)
lineCentres_xy$dist <- dist[,1]
lineCentres_xy$LinjeID <- Lines_spatial$LinjeID

### admin ##################################################

#get info on administrative names for the buffers
library(rgdal)
library(plyr)

#make both spatial objects in the same crs
NorwayADM <- readOGR(dsn="C:/Users/db40fysa/Dropbox/Alpine/NOR_adm",layer="NOR_adm2")
plot(NorwayADM)
NorwayADM$NAME_1 <- iconv(NorwayADM$NAME_1, "UTF-8","latin1")
table(NorwayADM$NAME_1)
NorwayADM <- spTransform(NorwayADM,crs(Polys_spatial))
NorwayADM$NAME_1 <- mergeCounties(as.character(NorwayADM$NAME_1),further=TRUE)
table(NorwayADM$NAME_1)

#overlay with the Polygons
myAdm <- over(Polys_spatial,NorwayADM)
myAdm$LinjeID <- Polys_spatial@data$LinjeID

#group??
table(myAdm$NAME_1)
myAdm <- myAdm[,c("LinjeID","NAME_1")]
names(myAdm)[2] <- "adm"

### habitat ################################################

setwd("C:/Users/db40fysa/Dropbox/Alpine/Habitat/Satveg_deling_nd_partnere_09_12_2009/Satveg_deling_nd_partnere_09_12_2009/tiff")
library(raster)

#project other rasters onto this
habitatRasterTop <- raster("NNred25-30-t1.tif")#30 x 30 m
habitatRasterTop <- aggregate(habitatRasterTop,fact=5,fun=modal,na.rm=T)
habitatRasterBot <- raster("sn25_geocorr.tif")
habitatRasterBot <- aggregate(habitatRasterBot,fact=5,fun=modal,na.rm=T)
extent(habitatRasterTop)
extent(habitatRasterBot)

#merge each dataset
totalExtent <- c(246285,1342485,6414184,8014594)
habitatRasterTop <- extend(habitatRasterTop,extent(totalExtent))
habitatRasterBot <- extend(habitatRasterBot,extent(totalExtent))
origin(habitatRasterBot) <- origin(habitatRasterTop)
habitatRaster <- merge(habitatRasterTop,habitatRasterBot)
plot(habitatRaster)
#yeah!!!

#Set NAs for irrelevant habitats
habitatRaster[habitatRaster==22] <- NA
habitatRaster[habitatRaster>24] <- NA
habitatRaster[habitatRaster==0] <- NA
plot(habitatRaster)

#plot each grid, extract what????
#crop raster to Norway extent
myraster <- habitatRaster
rasterCRS <- crs(myraster)

#get raster values for each polygon
myrasterDF <- raster::extract(myraster,Polys_spatial,weights=TRUE, 
                               normalizeWeights=FALSE, df=TRUE)

#simplify habitat counts
myrasterDF$MountainBirchForest <- ifelse(myrasterDF$layer%in%c(6,7,8),1,0)
myrasterDF$Bog <- ifelse(myrasterDF$layer%in%c(9,10),1,0)
myrasterDF$Forest <- ifelse(myrasterDF$layer%in%c(1:5),1,0)
myrasterDF$ODF <- ifelse(myrasterDF$layer%in%c(16,17),1,0)
myrasterDF$Meadows <- ifelse(myrasterDF$layer%in%c(18),1,0)
myrasterDF$OSF <- ifelse(myrasterDF$layer%in%c(12:15),1,0)
myrasterDF$Mire <- ifelse(myrasterDF$layer%in%c(11),1,0)
myrasterDF$SnowBeds <- ifelse(myrasterDF$layer%in%c(19,20),1,0)
myrasterDF$Open <- ifelse(myrasterDF$layer%in%c(11:20),1,0)
myrasterDF$Human <- ifelse(myrasterDF$layer%in%c(23,24),1,0)

#aggregate
myrasterDF<-ddply(myrasterDF,.(ID),summarise,
                  MountainBirchForest = sum(weight[MountainBirchForest==1])/sum(weight),
                  Bog = sum(weight[Bog==1])/sum(weight),
                  Forest = sum(weight[Forest=1])/sum(weight),
                  ODF = sum(weight[ODF==1])/sum(weight),
                  Meadows = sum(weight[Meadows==1])/sum(weight),
                  OSF = sum(weight[OSF==1])/sum(weight),
                  Mire = sum(weight[Mire==1])/sum(weight),
                  SnowBeds = sum(weight[SnowBeds==1])/sum(weight),
                  Open = sum(weight[Open==1])/sum(weight),
                  Human = sum(weight[Human==1])/sum(weight))

#add line ID
myrasterDF$LinjeID <- Polys_spatial$LinjeID

### climate ################################################################################

#average climatic conditions

#try the EuroLST dataset
setwd("C:/Users/db40fysa/Dropbox/Alpine/Bioclim_LST")
#-BIO1: Annual mean temperature (°C*10): eurolst_clim.bio01.zip (MD5) 72MB
#-BIO2: Mean diurnal range (Mean monthly (max - min tem)): eurolst_clim.bio02.zip (MD5) 72MB
#-BIO3: Isothermality ((bio2/bio7)*100): eurolst_clim.bio03.zip (MD5) 72MB
#-BIO4: Temperature seasonality (standard deviation * 100): eurolst_clim.bio04.zip (MD5) 160MB
#-BIO5: Maximum temperature of the warmest month (°C*10): eurolst_clim.bio05.zip (MD5) 106MB
#-BIO6: Minimum temperature of the coldest month (°C*10): eurolst_clim.bio06.zip (MD5) 104MB
#-BIO7: Temperature annual range (bio5 - bio6) (°C*10): eurolst_clim.bio07.zip (MD5) 132MB
#-BIO10: Mean temperature of the warmest quarter (°C*10): eurolst_clim.bio10.zip (MD5) 77MB
#-BIO11: Mean temperature of the coldest quarter (°C*10): eurolst_clim.bio11.zip (MD5) 78MB

#bio1#
temp_bio1 <- raster("eurolst_clim.bio01/eurolst_clim.bio01.tif")#res 250 m
temp_bio1 <- aggregate(temp_bio1,fact=4,fun=mean,na.rm=T)

#extract the data
out <- getBufferData(temp_bio1,Polys_spatial)

out_Bio1 <- out
rm(temp_bio1)
rm(out)

#maximum temp#
temp_bio5 <- raster("eurolst_clim.bio05/eurolst_clim.bio05.tif")
temp_bio5 <- aggregate(temp_bio5,fact=4,fun=mean,na.rm=T)
out <- getBufferData(temp_bio5,Polys_spatial)

out_Bio5 <- out
rm(temp_bio5)
rm(out)

#minimum temp
temp_bio6 <- raster("eurolst_clim.bio06/eurolst_clim.bio06.tif")
temp_bio6 <- aggregate(temp_bio6,fact=4,fun=mean,na.rm=T)
out <- getBufferData(temp_bio6,Polys_spatial)

out_Bio6 <- out
rm(temp_bio6)
rm(out)

### merge ################################################################

#temp data
names(out_Bio1)[2]<-"bio1"
names(out_Bio5)[2]<-"bio5"
names(out_Bio6)[2]<-"bio6"
out_Bio<-cbind(out_Bio1,bio5=out_Bio5[,2],bio6=out_Bio6[,2])

#combine others
varDF <- merge(out_Bio,myrasterDF,by="LinjeID",all=T)
varDF <- merge(varDF,myAdm,by="LinjeID",all=T)

#and with distance to coast data
names(lineCentres_xy)[3] <- "distCoast"
varDF <- merge(varDF,lineCentres_xy,by="LinjeID",all.x=T)

### alpine data #############################################################################

#also get alpine data

load("data/alpineData_5kmGrid.RData")

#tree line position
alpineData$tree_line_position <- alpineData$tree_line-alpineData$elevation
alpineData$alpine_habitat1 <- sapply(alpineData$alpine_habitat,function(x)ifelse(x==1,1,0))
alpineData$alpine_habitat2 <- sapply(alpineData$alpine_habitat,function(x)ifelse(x==2,1,0))
alpineData$alpine_habitat3 <- sapply(alpineData$alpine_habitat,function(x)ifelse(x==3,1,0))
alpineData$alpine_habitat4 <- sapply(alpineData$alpine_habitat,function(x)ifelse(x==4,1,0))


#make spatial
coordinates(alpineData) <- c("x","y")
proj4string(alpineData) <- CRS("+proj=utm +no_defs +zone=33 +a=6378137 +rf=298.257222101 +towgs84=0,0,0,0,0,0,0 +to_meter=1")
alpineData <- spTransform(alpineData,crs(equalM))

#overlay with the spatial polgons
alpineData$LinjeID <- over(alpineData, Polys_spatial)$LinjeID
alpineData <- subset(alpineData,!is.na(LinjeID))

#aggregate to site level
library(plyr)
alpineData <- ddply(alpineData@data,.(LinjeID),summarise,
                    tree_line_position = median(tree_line_position,na.rm=T),
                    tree_line = median(tree_line,na.rm=T),
                    elevation = median(elevation,na.rm=T),
                    alpine_habitat1 = mean(alpine_habitat1,na.rm=T), 
                    alpine_habitat2 = mean(alpine_habitat2,na.rm=T), 
                    alpine_habitat3 = mean(alpine_habitat3,na.rm=T),
                    alpine_habitat4 = mean(alpine_habitat4,na.rm=T)) 

### merge all ##############################################################################

varDF <- merge(varDF,alpineData,by="LinjeID")
varDF <- varDF[,-which(names(varDF)=="ID")]

### correlations ############################################################################

#Examine correlations among bugs variables
library(GGally)
ggpairs(varDF[,c(2:14,23:25)])

#correlations >0.7
#bio1 vs bio6
#x vs y
#choose bio6 for line transects

### save #############################################################################

saveRDS(varDF,file="data/varDF_allEnvironData_buffers_idiv.rds")

### end #####################################################