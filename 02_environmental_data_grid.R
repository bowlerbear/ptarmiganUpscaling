
#using a m grid
equalM<-"+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

#get helper file
source('generalFunctions.R')

#get focal grids
focusGrids <- readRDS("data/focusGrids.rds")

### get norway ##############################################################

library(raster)
library(maptools)
library(rgeos)
data(wrld_simpl)
Norway <- subset(wrld_simpl,NAME=="Norway")
NorwayOrig <- Norway
Norway <- gBuffer(Norway,width=1)
Norway <- spTransform(Norway,crs(equalM))
plot(Norway)

### grid ###################################################################

#create grid
newres = 5000#5 km grid
mygrid <- raster(extent(projectExtent(Norway,equalM)))
res(mygrid) <- newres 
mygrid[] <- 1:ncell(mygrid)
plot(mygrid)
gridTemp <- mygrid
plot(Norway,add=T)

### adm #####################################################################

#get info on administrative names for the grid
library(rgdal)
library(plyr)

#make both spatial objects in the same crs
NorwayADM <- readOGR(dsn="C:/Users/db40fysa/Dropbox/Alpine/NOR_adm",layer="NOR_adm2")

#get multiple points per grid
NorwayADM <- spTransform(NorwayADM,crs(equalM))
mygridPoints <- mygrid
mygridPoints <- disaggregate(mygridPoints,fact=10)
mygridPoints <- as.data.frame(mygridPoints,xy=T)
coordinates(mygridPoints) <- c("x","y")
proj4string(mygridPoints) <- crs(NorwayADM)

#check they overlap
plot(mygridPoints)
plot(NorwayADM,add=T,col="red")

#extract the data
myAdm <- over(mygridPoints,NorwayADM)
myAdm$grid <- mygridPoints$layer

#get mode of adm and adm2 per grid
myAdm <- ddply(myAdm,.(grid),summarise,
               adm = Mode(NAME_1),
               adm2 = Mode(NAME_2))

myAdm$adm[is.na(myAdm$adm)] <- "outside"
myAdm$adm2[is.na(myAdm$adm2)] <- "outside"
myAdm <- subset(myAdm,!is.na(grid))

saveRDS(myAdm,file="C:/Users/db40fysa/Dropbox/ptarmigan Upscaling/data/grid_Adm.rds")

#check the results
mygrid[] <- 0
mygrid[myAdm$grid] <- as.numeric(as.factor(myAdm$adm))
plot(mygrid)#looks good!

#how many outside do we have?
table(getValues(mygrid))

### accessibility ###########################################################

#Code to get environmental data:
library(raster)

#get accessibility map
#https://www.nature.com/articles/nature25181
#setwd("C:/Users/diana.bowler/OneDrive - NINA/maps/accessibility/accessibility_to_cities_2015_v1.0")
access <- raster("C:/Users/db40fysa/Nextcloud/sMon/Gis-DataData/Accessibility/2015_accessibility_to_cities_v1.0/2015_accessibility_to_cities_v1.0.tif")
out <- getEnvironData(access,mygrid)

#check data
hist(out$myraster)
summary(out$myraster)
outAccess <- out
names(outAccess)[2] <- "Accessibility"
outAccess$Accessibility[is.nan(outAccess$Accessibility)] <- NA
summary(outAccess$Accessibility)
outAccess <- subset(outAccess,!is.na(Accessibility))
outAccess <- subset(outAccess,!is.na(grid))

saveRDS(outAccess,file="data/grid_Access.rds")

#check the results
mygrid[] <- 0
mygrid[outAccess$grid] <- outAccess$Accessibility
plot(mygrid)#looks good!

rm(access)
rm(out)

### human pop density #######################################################

humanpop <- raster("C:/Users/db40fysa/Dropbox/Alpine/PopDensity/SEDAC_pop_data/gpw-v4-population-density-adjusted-to-2015-unwpp-country-totals-rev10_2000_30_sec_tif/gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev10_2000_30_sec.tif")
out <- getEnvironData(humanpop,mygrid)

#check data
hist(out$myraster)
summary(out$myraster)
outHP<- out
names(outHP)[2] <- "HumanPop"
outHP$HumanPop[is.nan(outHP$HumanPop)] <- NA
summary(outHP$HumanPop)
outHP <- subset(outHP,!is.na(HumanPop))
outHP <- subset(outHP,!is.na(grid))

saveRDS(outHP,file="data/grid_HumanDensity.rds")

#check the results
mygrid[] <- 0
mygrid[outHP$grid] <- outHP$HumanPop
plot(mygrid)#looks good!

# ## habitats ################################################################

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

#convert raster into points and convert 
myrasterDF <- as.data.frame(myraster,xy=T)
names(myrasterDF)[3] <- "myraster"
myrasterDF<-subset(myrasterDF,!is.na(myraster))
coordinates(myrasterDF) <-c("x","y")
proj4string(myrasterDF) <- rasterCRS
myrasterDF <- spTransform(myrasterDF,crs(equalM))
                        
#get general grid
mygrid<-gridTemp
projection(mygrid) <- CRS(equalM) 

#get myraster values per grid cell
mygrid[] <- 1:ncell(mygrid)
variable <- raster::extract(mygrid,myrasterDF)
mygrid[] <- 0
mygrid[variable] <- 1
myrasterDF <- data.frame(myrasterDF@data)
myrasterDF$grid <- variable

library(reshape2)
myrasterDF <- melt(table(myrasterDF$grid,myrasterDF$myraster))
names(myrasterDF)<-c("grid","raster","count")

#add on total count per grid
library(plyr)
gridTotals <- ddply(myrasterDF,.(grid),summarise,total=sum(count))
myrasterDF$total <- gridTotals$total[match(myrasterDF$grid,gridTotals$grid)]

#simplify habitat counts

#predict positive effect
myrasterDF$MountainBirchForest<-0
myrasterDF$MountainBirchForest[myrasterDF$raster%in%c(6,7,8)]<-myrasterDF$count[myrasterDF$raster%in%c(6,7,8)]

#predict positive effect
myrasterDF$Bog<-0
myrasterDF$Bog[myrasterDF$raster%in%c(9,10)]<-myrasterDF$count[myrasterDF$raster%in%c(9,10)]

#predict negative effect
myrasterDF$Forest<-0
myrasterDF$Forest[myrasterDF$raster%in%c(1:5)]<-myrasterDF$count[myrasterDF$raster%in%c(1:5)]

#predict positive effect 
myrasterDF$ODF<-0
myrasterDF$ODF[myrasterDF$raster%in%c(16,17)]<-myrasterDF$count[myrasterDF$raster%in%c(16,17)]

#predict positive effect 
myrasterDF$Meadows<-0
myrasterDF$Meadows[myrasterDF$raster%in%c(18)]<-myrasterDF$count[myrasterDF$raster%in%c(18)]

#predict negative effect
myrasterDF$OSF<-0
myrasterDF$OSF[myrasterDF$raster%in%c(12,13,14,15)]<-myrasterDF$count[myrasterDF$raster%in%c(12,13,14,15)]

#?
myrasterDF$Mire<-0
myrasterDF$Mire[myrasterDF$raster%in%c(11)]<-myrasterDF$count[myrasterDF$raster%in%c(11)]

#expect negative effect
myrasterDF$SnowBeds<-0
myrasterDF$SnowBeds[myrasterDF$raster%in%c(19,20)]<-myrasterDF$count[myrasterDF$raster%in%c(19,20)]

#open habitat
myrasterDF$Open<-0
myrasterDF$Open[myrasterDF$raster%in%c(11:20)]<-myrasterDF$count[myrasterDF$raster%in%c(11:20)]

#human habitat - agriculture and urban
myrasterDF$Human<-0
myrasterDF$Human[myrasterDF$raster%in%c(23,24)]<-myrasterDF$count[myrasterDF$raster%in%c(23,24)]


#aggregate
myrasterDF<-ddply(myrasterDF,.(grid),summarise,
                  MountainBirchForest = sum(MountainBirchForest)/unique(total),
                  Bog = sum(Bog)/unique(total),
                  Mire = sum(Mire)/unique(total),
                  Open = sum(Open)/unique(total),
                  Human = sum(Human)/unique(total),
                  Forest = sum(Forest)/unique(total),
                  ODF = sum(ODF)/unique(total),
                  Meadows = sum(Meadows)/unique(total),
                  OSF = sum(OSF)/unique(total),
                  SnowBeds = sum(SnowBeds)/unique(total),
                  total = unique(total))


saveRDS(myrasterDF,"data/grid_Habitats.rds")

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
out <- getEnvironData(temp_bio1,mygrid)

out_Bio1 <- out
rm(temp_bio1)
rm(out)

#maximum temp#
temp_bio5 <- raster("eurolst_clim.bio05/eurolst_clim.bio05.tif")
temp_bio5 <- aggregate(temp_bio5,fact=4,fun=mean,na.rm=T)
out <- getEnvironData(temp_bio5,mygrid)

out_Bio5 <- out
rm(temp_bio5)
rm(out)

#minimum temp
temp_bio6 <- raster("eurolst_clim.bio06/eurolst_clim.bio06.tif")
temp_bio6 <- aggregate(temp_bio6,fact=4,fun=mean,na.rm=T)
out <- getEnvironData(temp_bio6,mygrid)

out_Bio6 <- out
rm(temp_bio6)
rm(out)

### merge ################################################################

#temp data
names(out_Bio1)[2]<-"bio1"
names(out_Bio5)[2]<-"bio5"
names(out_Bio6)[2]<-"bio6"
out_Bio<-cbind(out_Bio1,bio5=out_Bio5[,2],bio6=out_Bio6[,2])
saveRDS(out_Bio,"data/grid_Climate.rds")

#combine others
varDF <- merge(out_Bio,myrasterDF,by="grid",all=T)
varDF <- merge(varDF,myAdm,by="grid",all=T)
varDF <- merge(varDF,outHP,all=T)
varDF <- merge(varDF,outAccess,all=T)
varDF <- subset(varDF,!is.na(grid))

saveRDS(varDF,file="data/varDF_missing_5km_idiv.rds")

### correlations ############################################################################

#Examine correlations among bugs variables
library(GGally)
ggpairs(varDF[,2:11])

#bio1 and bio6 strongly related (0.808)
#top and open are strongly related (0.726)
#pref open and top (0.901)
#Agriculture and bottom are strongly related (0.907)

# #plot in space
# gridDF<-as.data.frame(gridTemp,xy=T)
# varDF$grid<-bugs.data$grid
# varDF<-merge(varDF,gridDF,by.x="grid",by.y="layer",all.x=T)
# 
# qplot(x,y,data=varDF,color=habitat)
# qplot(x,y,data=varDF,color=Forest)
# qplot(x,y,data=varDF,color=Open)
# qplot(x,y,data=varDF,color=Top)
# qplot(x,y,data=varDF,color=Agriculture)
# qplot(x,y,data=varDF,color=bio1)
# qplot(x,y,data=varDF,color=bio5)
# qplot(x,y,data=varDF,color=bio6)

### coverage ################################################################################

#subset to focal grids
varDF <- subset(varDF,grid %in% focusGrids)
# missing habitat data for 57 grids
apply(varDF,2,function(x)sum(is.na(x)))

#check plots
mygrid[] <- NA
mygrid[varDF$grid] <- varDF$bio1
plot(mygrid)
#looks good
mygrid[] <- NA
mygrid[varDF$grid] <- varDF$Forest
plot(mygrid)
#ok
mygrid[] <- NA
mygrid[varDF$grid]<-varDF$bio6
plot(mygrid)
#good!

#where are these 57 grids with missing habitat data
mygrid[]<-0
mygrid[varDF$grid[is.na(varDF$Forest)]]<-1
plot(mygrid)
plot(Norway,add=T)
#these are on the edge of Norway

### alpine data #############################################################################

#also get alpine data
load("data/alpineData_5kmGrid.RData")

alpineData <- subset(alpineData,site %in% focusGrids)
 
#plot it
library(ggplot2)
qplot(x,y,data=alpineData,color=elevation)
qplot(x,y,data=alpineData,color=alpine_habitat)
#look ok

#summarise it
summary(alpineData$tree_line)
summary(alpineData$elevation)
apply(alpineData,2,function(x)sum(is.na(x)))

#tree line position
alpineData$tree_line_position <- alpineData$tree_line-alpineData$elevation
alpineData$alpine_habitat1 <- sapply(alpineData$alpine_habitat,function(x)ifelse(x==1,1,0))
alpineData$alpine_habitat2 <- sapply(alpineData$alpine_habitat,function(x)ifelse(x==2,1,0))
alpineData$alpine_habitat3 <- sapply(alpineData$alpine_habitat,function(x)ifelse(x==3,1,0))
alpineData$alpine_habitat4 <- sapply(alpineData$alpine_habitat,function(x)ifelse(x==4,1,0))

#aggregate to site level
library(plyr)
alpineData <- ddply(alpineData,.(site),summarise,
                    tree_line_position = median(tree_line_position,na.rm=T),
                    tree_line = median(tree_line,na.rm=T),
                    elevation = median(elevation,na.rm=T),
                    alpine_habitat1 = mean(alpine_habitat1,na.rm=T), 
                    alpine_habitat2 = mean(alpine_habitat2,na.rm=T), 
                    alpine_habitat3 = mean(alpine_habitat3,na.rm=T),
                    alpine_habitat4 = mean(alpine_habitat4,na.rm=T)) 

names(alpineData)[which(names(alpineData)=="site")]<-"grid"

#plot data 
mygrid[] <- NA
mygrid[alpineData$grid] <- alpineData$tree_line
plot(mygrid)

mygrid[] <- NA
mygrid[alpineData$grid] <- alpineData$elevation
plot(mygrid)

#plot missing data 
mygrid[] <- NA
mygrid[alpineData$grid[is.na(alpineData$tree_line)]] <- 1 
plot(mygrid)

mygrid[] <- NA
mygrid[alpineData$grid[is.na(alpineData$elevation)]] <- 1 
plot(mygrid)

#any more NAs??
apply(alpineData,2,function(x)sum(is.na(x)))
#56 are missing except elevation

#are these in addition of different to the missing data in varDF
missingVarDF <- varDF$grid[is.na(varDF$Forest)]
missingAlpine <- alpineData$grid[is.na(alpineData$tree_line)]
length(unique(missingVarDF,missingAlpine))#57 in total!!

#merge all
varDF <- merge(varDF,alpineData,by="grid",all=T)
varDF <- subset(varDF,!is.na(Forest))
varDF <- subset(varDF,!is.na(tree_line))

### check missing data again #########################################################

apply(varDF,2,function(x)sum(is.na(x)))

varDF$adm <- iconv(varDF$adm,"UTF-8","latin1")
varDF$adm2 <- iconv(varDF$adm2,"UTF-8","latin1")

#5 grids are outside - can we figure out there admin???
table(varDF$adm)
table(varDF$adm2)
subset(varDF,adm=="outside")$grid
#22415 25249 25253 25721 74673

#for each one, see what is the adm before and after the grid
subset(varDF,grid %in% c(22412,22413,22414,22415,22416,22417,22418))#none
#see about cells below
mygrid[] <- 1:ncell(mygrid)
rowColFromCell(mygrid,22415)
getValues(mygrid,row=95)[137]
getValues(mygrid,row=94)[137]
getValues(mygrid,row=96)[137]
subset(varDF,grid %in% c(22178,22652))
varDF$adm[varDF$grid==22415] <- "Troms"
varDF$adm2[varDF$grid==22415] <- "Bardu"

subset(varDF,grid %in% c(25247,25248,25250,25251))
varDF$adm[varDF$grid==25249] <- "Nordland"
varDF$adm2[varDF$grid==25249] <- "Tysfjord"

subset(varDF,grid %in% c(25251,25252,25253,25254,25255))
varDF$adm[varDF$grid==25253] <- "Nordland"
varDF$adm2[varDF$grid==25253] <- "Narvik"

subset(varDF,grid %in% c(25719,25720,25721,25722,25723))
varDF$adm[varDF$grid==25721] <- "Nordland"
varDF$adm2[varDF$grid==25721] <- "Tysfjord"

subset(varDF,grid %in% c(74671,74672,74673,74674,74675))
varDF$adm[varDF$grid==74673] <- "Rogaland"
varDF$adm2[varDF$grid==74673] <- "Klepp"

#also:
#11587, 12063, 25485
subset(varDF,grid %in% c(11585,11586,11588,11589))
varDF$adm[varDF$grid==11587] <- "Finnmark"
varDF$adm2[varDF$grid==11587] <- "Sør-Varanger"

subset(varDF,grid %in% c(12061,12062,12064,12065))
varDF$adm[varDF$grid==12063] <- "Finnmark"
varDF$adm2[varDF$grid==12063] <- "Sør-Varanger"

subset(varDF,grid %in% c(25483,25484,25486,25487))
varDF$adm[varDF$grid==25485] <- "Nordland"
varDF$adm2[varDF$grid==25485] <- "Tysfjord"


table(varDF$adm)

### save #############################################################################

saveRDS(varDF,file="data/varDF_allEnvironData_5km_idiv.rds")

### group admins ####################################################################

source('generalFunctions.R', encoding = 'UTF-8'f)

varDF$admGrouped <- mergeCounties(varDF$adm,further=TRUE)
table(varDF$adm)
table(varDF$admGrouped)
sum(is.na(varDF$admGrouped))

saveRDS(varDF,file="data/varDF_allEnvironData_5km_idiv.rds")

### latitude and distance to the coast ###############################################

#coordinates
mygridDF <- as.data.frame(mygrid,xy=T)
names(mygridDF)[3] <- "grid"

varDF <- merge(varDF, mygridDF, by="grid",all.x=T)

#distance to coast
sldf_coast <- rnaturalearth::ne_coastline()
plot(sldf_coast)
coast <- sf::st_as_sf(sldf_coast)

#get points
mygridPoints <- sf::st_as_sf(varDF[,c("grid","x","y")], 
                            coords = c("x", "y"), crs = proj4string(NorwayADM))
mygridPoints <- sf::st_transform(mygridPoints,sf::st_crs(coast))
ggplot()+
  geom_sf(data=coast)+
  geom_sf(data=mygridPoints,color="red")

#get distance between them
mygridPoints1 <- sf::st_coordinates(mygridPoints)
dist <- geosphere::dist2Line(p = as.matrix(mygridPoints1), 
                             line = sldf_coast)

#pack into data frame
dist <- as.data.frame(dist)
dist$grid <- mygridPoints$grid
names(dist)[1] <- "distCoast"
  
#merge with varDF
varDF <- merge(varDF, dist[,c(1:3,5)], by="grid", all.x=T)

saveRDS(varDF,file="data/varDF_allEnvironData_5km_idiv.rds")

### correlations ############################################

library(GGally)
ggpairs(varDF[,c(2:13,20:22,28:30)])
cor(varDF[,c(2:13,20:22,28:30)])

#cors >0.7
#bio1 and bio6
#bio1 and tree line position
#forest and open
#osf and open
#open and tree line position
#y and tree line (0.76)
#y and x 

### end #####################################################
