#script to analysis line transect data on the HPC
library(tidyverse)
library(sp)
library(rgeos)
library(raster)
library(maptools)

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

### aggregate data to the grid ######################################

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
  group_by(LinjeID) %>%
  summarise(meanNu = mean(totalsInfo,na.rm=T),
            meanTL = mean(length[length!=0])) %>%
  filter(!is.na(meanNu)) %>%
  filter(!is.infinite(meanNu))

bufferData <- readRDS("data/varDF_allEnvironData_buffers_idiv.rds")

siteMeans <- merge(siteMeans,bufferData,by="LinjeID")
siteMeans <- subset(siteMeans,!is.na(tree_line))

### get spatial polygons ###########################################

Polys_spatial <- readRDS("data/Polys_spatial.rds")
Polys_spatial <- merge(Polys_spatial,siteMeans,by="LinjeID")
Polys_spatial <- subset(Polys_spatial,!is.na(meanNu))
Polys_spatial$meanDensity <- Polys_spatial$meanNu/Polys_spatial$meanTL
head(Polys_spatial)
summary(Polys_spatial)

library(tmap)
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

hist(siteMeans$meanNu)
summary(siteMeans$meanNu)

glm1 <- glm(log(meanNu+1) ~ scale(bio1) + scale(bio5) + scale(bio6) + 
            Open + PrefOpen + PrefClosed +
            scale(tree_line) + 
            scale(elevation)+
            scale(elevation^2),
            data=siteMeans)
summary(glm1)
  #                           Estimate Std. Error t value Pr(>|t|)    
  # (Intercept)                1.9298457  0.1438699  13.414  < 2e-16 ***
  # scale(bio1)                0.1892296  0.0614568   3.079 0.002171 ** 
  # scale(bio5)               -0.0766258  0.0363032  -2.111 0.035207 *  
  # scale(bio6)                0.0357897  0.0567560   0.631 0.528548    
  # Open                       0.0017925  0.0008570   2.092 0.036883 *  
  # PrefOpen                   0.0006795  0.0006042   1.125 0.261183    
  # PrefClosed                -0.0012863  0.0016000  -0.804 0.421760    
  # scale(tree_line_position) -0.1207165  0.0357055  -3.381 0.000769 ***
  # scale(elevation)           0.7300731  0.1163483   6.275 6.68e-10 ***
  # scale(elevation^2)        -0.7514253  0.1225348  -6.132 1.57e-09 ***
  
library(MuMIn)
options(na.action = "na.fail")
dd <- dredge(glm1)
subset(dd, delta < 2)


#corrected by transect length
summary(siteMeans$meanNu/siteMeans$meanTL)
hist(siteMeans$meanNu/siteMeans$meanTL)
glm1 <- lm(meanNu/meanTL ~ scale(bio1) + scale(bio5) + scale(bio6) + 
              Open + PrefOpen + PrefClosed +
              scale(tree_line_position) + 
              scale(elevation)+scale(elevation^2),
              data=siteMeans)
summary(glm1)
Coefficients:
  #                             Estimate Std. Error t value Pr(>|t|)    
  # (Intercept)                3.816e-03  5.429e-04   7.028 5.68e-12 ***
  # scale(bio1)                3.771e-04  2.319e-04   1.626  0.10451    
  # scale(bio5)               -6.277e-04  1.370e-04  -4.582 5.61e-06 ***
  # scale(bio6)                4.577e-04  2.142e-04   2.137  0.03301 *  
  # Open                      -3.156e-06  3.234e-06  -0.976  0.32953    
  # PrefOpen                   1.573e-06  2.280e-06   0.690  0.49043    
  # PrefClosed                -1.565e-05  6.038e-06  -2.592  0.00976 ** 
  # scale(tree_line_position) -2.049e-04  1.347e-04  -1.521  0.12881    
  # scale(elevation)           2.808e-03  4.391e-04   6.394 3.23e-10 ***
  # scale(elevation^2)        -2.440e-03  4.624e-04  -5.276 1.84e-07 ***

### brt #########################################################################

#Boosted regression tree
library(dismo)
library(gbm)

siteMeans$meanNu <- log(siteMeans$meanNu/siteMeans$meanTL+1)

brt1 <- gbm.step(data=siteMeans, 
                 gbm.x = c(4:13,16:22), 
                 gbm.y = 2,
                 family = 'gaussian')

summary(brt1)
#                                 var    rel.inf
# tree_line                   tree_line 25.0060644
# bio6                             bio6 20.8827893
# bio5                             bio5 12.3456494
# elevation                   elevation  9.8402625
# PrefOpen                     PrefOpen  9.1553879
# tree_line_position tree_line_position  4.2568394
# bio1                             bio1  4.0936455
# PrefClosed                 PrefClosed  3.1023779
# Open                             Open  3.0298239
# Forest                         Forest  2.3967398
# Top                               Top  1.9409646
# alpine_habitat3       alpine_habitat3  1.5355663
# alpine_habitat1       alpine_habitat1  0.8527290
# Bottom                         Bottom  0.6128414
# alpine_habitat2       alpine_habitat2  0.5938059
# Agriculture               Agriculture  0.3545129
# alpine_habitat4       alpine_habitat4  0.0000000

#plot main effects
gbm.plot(brt1, n.plots=8, write.title = TRUE)
gbm.plot.fits(brt1)
#non-linear plot for tree line position and bio1

#interactions?
find.int <- gbm.interactions(brt1)
find.int$interactions#none!
find.int$rank.list

### end ################################################################
