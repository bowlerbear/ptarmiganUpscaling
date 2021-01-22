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

siteMeans$adm <- mergeCounties(siteMeans$adm,further=TRUE)#needs to be done again

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

glm1 <- lm(log(meanNu+1) ~ scale(bio1) + scale(bio5) + scale(bio6) + 
            Open + PrefOpen + PrefClosed +
            scale(tree_line) + 
            scale(elevation)+
            scale(elevation^2),
            data=siteMeans)
summary(glm1)
#                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         1.9044860  0.1466637  12.985  < 2e-16 ***
#   scale(bio1)         0.1607169  0.0633338   2.538  0.01142 *  
#   scale(bio5)        -0.0525747  0.0364160  -1.444  0.14935    
#   scale(bio6)         0.0653411  0.0593270   1.101  0.27119    
#   Open                0.0014211  0.0007037   2.019  0.04389 *  
#   PrefOpen            0.0006222  0.0005126   1.214  0.22537    
#   PrefClosed         -0.0006118  0.0013092  -0.467  0.64044    
#   scale(tree_line)   -0.2987196  0.0979224  -3.051  0.00239 ** 
#   scale(elevation)    1.0294777  0.1732242   5.943 4.79e-09 ***
#   scale(elevation^2) -0.6955471  0.1222922  -5.688 2.03e-08 ***
  
library(MuMIn)
options(na.action = "na.fail")
dd <- dredge(glm1)
subset(dd, delta < 2)


#corrected by transect length
summary(siteMeans$meanNu/siteMeans$meanTL)
hist(siteMeans$meanNu/siteMeans$meanTL)
glm1 <- lm(log(meanNu/meanTL+1) ~ scale(bio1) + scale(bio5) + scale(bio6) + 
              Open + PrefOpen + PrefClosed +
              scale(tree_line) + 
              scale(elevation)+scale(elevation^2),
              data=siteMeans)
summary(glm1)
Coefficients:
  Coefficients:
  #                       Estimate Std. Error t value Pr(>|t|)    
  # (Intercept)         3.664e-03  5.397e-04   6.788 2.78e-11 ***
  # scale(bio1)         3.241e-04  2.331e-04   1.390   0.1649    
  # scale(bio5)        -5.421e-04  1.340e-04  -4.045 5.93e-05 ***
  # scale(bio6)         5.047e-04  2.183e-04   2.312   0.0211 *  
  # Open               -2.339e-06  2.590e-06  -0.903   0.3669    
  # PrefOpen            1.398e-06  1.887e-06   0.741   0.4590    
  # PrefClosed         -9.982e-06  4.818e-06  -2.072   0.0387 *  
  # scale(tree_line)   -5.929e-04  3.604e-04  -1.645   0.1005    
  # scale(elevation)    3.321e-03  6.375e-04   5.209 2.63e-07 ***
  # scale(elevation^2) -2.236e-03  4.500e-04  -4.969 8.82e-07 ***

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
#                                   var     rel.inf
# bio6                             bio6 23.59032517
# tree_line                   tree_line 22.84091205
# bio5                             bio5 12.76908910
# elevation                   elevation 10.83059745
# PrefOpen                     PrefOpen  7.75299072
# bio1                             bio1  3.76294838
# alpine_habitat2       alpine_habitat2  3.16039800
# tree_line_position tree_line_position  3.05851428
# Top                               Top  2.79262098
# Open                             Open  2.37575864
# alpine_habitat3       alpine_habitat3  2.08569589
# PrefClosed                 PrefClosed  1.86873326
# Bottom                         Bottom  1.45702810
# Forest                         Forest  1.41170978
# Agriculture               Agriculture  0.20984804
# alpine_habitat1       alpine_habitat1  0.03283015
# alpine_habitat4       alpine_habitat4  0.00000000

#plot main effects
gbm.plot(brt1, n.plots=8, write.title = TRUE)
gbm.plot.fits(brt1)
#non-linear plot for tree line position and bio1

#interactions?
find.int <- gbm.interactions(brt1)
find.int$interactions#none!
find.int$rank.list

### end ################################################################
