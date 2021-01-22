#script to analysis line transect data on the HPC
library(tidyverse)
library(sp)
library(rgeos)
library(raster)
library(maptools)

#specify top level folder
#myfolder <- "Data" #on local PC
myfolder <- "/data/idiv_ess/ptarmiganUpscaling" #HPC

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

### glm model ######################################################

hist(siteMeans$meanNu)
summary(siteMeans$meanNu)

glm1 <- glm(log(meanNu+1) ~ scale(bio1) + scale(bio5) + scale(bio6) + 
            Open + PrefOpen + PrefClosed +
            scale(tree_line_position) + 
            scale(elevation)+
            scale(elevation^2),
            data=siteMeans)
summary(glm1)

library(MuMIn)
options(na.action = "na.fail")
dd <- dredge(glm1)
subset(dd, delta < 2)


#corrected by transect length
summary(siteMeans$meanNu/siteMeans$meanTL)
hist(siteMeans$meanNu/siteMeans$meanTL)
glm1 <- glm(meanNu/meanTL ~ scale(bio1) + scale(bio5) + scale(bio6) + 
              Open + PrefOpen + PrefClosed +
              scale(tree_line_position) + 
              scale(elevation),
              data=siteMeans)
summary(glm1)

#look at predictions
siteInfo_ArtsDaten$preds <- exp(predict(glm1,newdata=siteInfo_ArtsDaten,type="response"))
mygrid[] <- NA
mygrid[siteInfo_ArtsDaten$grid] <- siteInfo_ArtsDaten$preds
plot(mygrid)

### brt #########################################################################

#Boosted regression tree
library(dismo)
library(gbm)

siteMeans$meanNu <- log(siteMeans$meanNu+1)
siteMeans$meanNu <- siteMeans$meanNu/siteMeans$meanTL

brt1 <- gbm.step(data=siteMeans, 
                 gbm.x = c(4:13,16:22), 
                 gbm.y = 2,
                 family = 'gaussian')

summary(brt1)
# var     rel.inf
# bio6                             bio6 18.54525269
# bio5                             bio5 17.75563279
# elevation                   elevation 11.92057815
# tree_line                   tree_line 11.44292757
# PrefOpen                     PrefOpen  7.35145713
# tree_line_position tree_line_position  7.28765948
# bio1                             bio1  6.06444292
# Top                               Top  3.71964277
# Bottom                         Bottom  3.28750231
# Open                             Open  3.07606849
# Forest                         Forest  2.87752486
# PrefClosed                 PrefClosed  2.81810993
# alpine_habitat3       alpine_habitat3  1.98746598
# alpine_habitat2       alpine_habitat2  1.06964121
# Agriculture               Agriculture  0.72161226
# alpine_habitat1       alpine_habitat1  0.07448146
# alpine_habitat4       alpine_habitat4  0.00000000

#plot main effects
gbm.plot(brt1, n.plots=11, write.title = TRUE)
gbm.plot.fits(brt1)
#non-linear plot for tree line position and bio1

#interactions?
find.int <- gbm.interactions(brt1)
find.int$interactions#none!
find.int$rank.list

### end ################################################################
