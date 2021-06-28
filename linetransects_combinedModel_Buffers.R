#script to analysis line transect data on the HPC
library(tidyverse)
library(sp)
library(rgeos)
library(raster)
library(maptools)
library(tmap)

source('generalFunctions.R', encoding = 'UTF-8')

#specify top level folder
myfolder <- "Data" #on local PC
#myfolder <- "/data/idiv_ess/ptarmiganUpscaling" #HPC

### get norway ##############################################

equalM<-"+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

library(raster)
library(maptools)
library(rgeos)
data(wrld_simpl)
Norway <- subset(wrld_simpl,NAME=="Norway")
NorwayOrig <- Norway
Norway <- gBuffer(Norway,width=1)
Norway <- spTransform(Norway,crs(equalM))
NorwayOrig <- spTransform(NorwayOrig,crs(equalM))
plot(Norway)

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

### aggregate data to the lines ######################################

#Get statistics per year and line
tlDF <- allData %>%
          dplyr::group_by(LinjeID,Year) %>%
          dplyr::summarise(nuGroups=length(totalIndiv[!is.na(totalIndiv)]),
                           totalsInfo=sum(totalIndiv,na.rm=T),
                           groupSize=mean(totalIndiv,na.rm=T),
                           length = mean(LengdeTaksert,na.rm=T))
sum(tlDF$totalsInfo,na.rm=T)
#67946

#insert NA when there is no transect but evidence of a survey
tlDF$length[is.na(tlDF$length)] <- 0
tlDF$nuGroups[tlDF$length==0 ] <- NA
tlDF$totalsInfo[tlDF$length==0] <- NA
tlDF$groupSize[tlDF$length==0] <- NA
summary(tlDF)
sum(tlDF$length==0)#1302

### get environ data ##############################################

#get mean across all years
siteMeans <- tlDF %>%
  mutate(density = totalsInfo/(length/1000 * 106/1000 * 2)) %>%
  group_by(LinjeID) %>%
  summarise(meanNu = mean(totalsInfo,na.rm=T),
            meanTL = mean(length[length!=0]),
            meanDensity = mean(density,na.rm=T)) %>%
  filter(!is.na(meanNu)) %>%
  filter(!is.infinite(meanNu))
summary(siteMeans)


bufferData <- readRDS("data/varDF_allEnvironData_buffers_idiv.rds")
siteMeans <- merge(siteMeans,bufferData,by="LinjeID")
siteMeans <- subset(siteMeans,!is.na(tree_line))

siteMeans$adm <- mergeCounties(siteMeans$adm,further=TRUE)
table(siteMeans$adm)

### get spatial polygons ###########################################

Polys_spatial <- readRDS("data/Polys_spatial.rds")
Polys_spatial <- merge(Polys_spatial,siteMeans,by="LinjeID")
Polys_spatial <- subset(Polys_spatial,!is.na(meanDensity))
Polys_spatial$meanDensity[Polys_spatial$meanDensity > 25] <- 25
head(Polys_spatial)
summary(Polys_spatial)

tm_shape(NorwayOrig) +
  tm_borders()+
tm_shape(Polys_spatial)+
  tm_fill(col="adm")

tm_shape(NorwayOrig) +
  tm_fill("grey")+
tm_shape(Polys_spatial)+
  tm_fill(col="meanDensity",style="pretty",n=7)+
  tm_legend(legend.position=c("left","top"))

### glm model ######################################################

hist(siteMeans$meanDensity)
summary(siteMeans$meanDensity)
hist(log(siteMeans$meanDensity+1))

glm1 <- lm(log(meanDensity+1) ~ scale(bio5) + scale(bio6) + 
             Forest + Bog + ODF + OSF + Mire + SnowBeds + Human + 
             scale(y) +
             scale(distCoast) +
             scale(tree_line) + 
             scale(elevation),
           data=siteMeans)
summary(glm1)#21%

#bio5 has negative effect
#bio6 has positive effect
#distCoast has negative effect

library(MuMIn)
options(na.action = "na.fail")
dd <- dredge(glm1)
subset(dd, delta < 2)
#bog, forest, mire, ODF, OSF,bio5, bio6, distCoast,treeline, y, snowbeds

#check variance inflation

glm1 <- lm(log(meanDensity+0.01) ~ 
             scale(bio5) + 
             scale(bio6) + 
             Forest + 
             Bog + 
             SnowBeds +
             Meadows + 
             scale(y) +
             scale(distCoast) +
             scale(tree_line) + 
             scale(elevation),
           data=siteMeans)
summary(glm1)#only 20%
car::vif(glm1)
#tree line, elevation and y

pairs(siteMeans[,c("y","distCoast","tree_line","elevation","tree_line_position")])
#tree line and elevation are highly correlated
#y and elevation are also correlated

ggplot(siteMeans,aes(y=log(meanDensity+1),
                     x=y))+
  geom_point()+stat_smooth()
#quadratic should be fine

ggplot(siteMeans,aes(y=log(meanDensity+1),
                     x=distCoast))+
  geom_point()+stat_smooth()
#nothing obvious

ggplot(siteMeans,aes(y=log(meanDensity+1),
                     x=tree_line))+
  geom_point()+stat_smooth()
#increase, but humped

ggplot(siteMeans,aes(y=log(meanDensity+1),
                     x=tree_line_position))+
  geom_point()+stat_smooth()
#nothing obvious

ggplot(siteMeans,aes(y=log(meanDensity+1),
                     x=elevation))+
  geom_point()+stat_smooth()
#general increase - humped

#look at partial plots

glm1 <- lm(log(meanDensity+0.01) ~ 
             scale(bio5) + 
             scale(bio6) + 
             Forest + 
             Bog + 
             SnowBeds +
             Meadows + 
             scale(y) +
             scale(distCoast) +
             scale(tree_line),
           data=siteMeans)
car::avPlots(glm1)

### brt #########################################################################

#Boosted regression tree
library(dismo)
library(gbm)
siteMeans$log.Density <- log(siteMeans$meanDensity+0.01)

brt1 <- gbm.step(data=siteMeans, 
                 gbm.x = c(5:17,20,26:28), 
                 gbm.y = 29,
                 family = 'gaussian')

summary(brt1)
#                                     var    rel.inf
# bio6                               bio6 22.3458169
# bio5                               bio5 14.2431986
# y                                     y 11.7659117
# x                                     x  9.9893607
# distCoast                     distCoast  9.9046832#non linear effect maybe
# tree_line                     tree_line  6.8581819
# Open                               Open  3.8421042
# bio1                               bio1  2.8542700
# MountainBirchForest MountainBirchForest  2.7410876
# Mire                               Mire  2.6665185
# OSF                                 OSF  2.6246154
# Forest                           Forest  2.4941139
# Meadows                         Meadows  2.1145808
# SnowBeds                       SnowBeds  1.9913050
# Bog                                 Bog  1.8380248

#plot main effects
gbm.plot(brt1, n.plots=12, write.title = TRUE)

#interactions?
find.int <- gbm.interactions(brt1)
find.int$interactions#none!
find.int$rank.list

#check similar glm
glm1 <- lm(log(meanDensity+0.01) ~ bio6 + y + bio5 + distCoast + tree_line+
             Open+bio1 + Mire + Forest + MountainBirchForest + Meadows + OSF +
             SnowBeds + Bog,data=siteMeans)
summary(glm1)#21%...

#### gam model #########################################################

library(mgcv)

gam1 <- gam(log(meanDensity+0.01) ~ s(bio6,k=3) + s(y,k=3) + s(bio5,k=3) + s(distCoast,k=3) + s(tree_line,k=3),data=siteMeans)
summary(gam1)#25%
plot(gam1)

#### quadratic model ##################################################

glm1 <- lm(log(meanDensity+0.01) ~ bio6 + I(bio6^2)+
                                  distCoast + I(distCoast^2)+
                                  bio5 + I(bio5^2)+
                                  tree_line + I(tree_line^2)+
                                  SnowBeds + 
                                  y +
                                  OSF,data=siteMeans)
summary(glm1)#26%

### end ################################################################

