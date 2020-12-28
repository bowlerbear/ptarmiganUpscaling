
#get the occupancy data:

#########################################################################################

#run ArtenDatenBank_missing_data_5km_sparta.R script

#get the bugs.data object

bugs.data_ArtsDaten <- bugs.data

########################################################################################

#get BUGS functions
source('bugsFunctions.R')

#need to specify initial values
library(reshape2)
zst <- acast(listlengthDF, siteIndex~yearIndex, value.var="y",fun=max,na.rm=T)
#line below if you data subset - see section below
#zst <- acast(selected_listlengthDF, s_siteIndex~yearIndex, value.var="y",fun=max,na.rm=T)
zst [is.infinite(zst)] <- NA

#fill in the blanks more cleverly
replace_na_with_last<-function(x,a=!is.na(x)){
  x[which(a)[c(1,1:sum(a))][cumsum(a)+1]]
}

#inits <- function(){list(z = zst)}
for(i in 1:nrow(zst)){
  zst[i,] <- replace_na_with_last(zst[i,])
}

inits <- function(){list(z = zst)}

########################################################################################

#get the abundance data:

#run the linetransects_combinedModel.Rmd script up until the 2nd bugs objects
bugs.data_LineTransects <- bugs.data 

#######################################################################################

#get effective strip width data

#load("model-outputs/out1_linetransectModel_detection.RData")
#siteInfo<-siteInfo_LineTransects
#esw<-getBUGSFits(out1,param="predESW")
#bugs.data_LineTransects$ESW.mean <- acast(esw,lineIndex~Year,value.var="mean")
#bugs.data_LineTransects$ESW.sd <- acast(esw,lineIndex~Year,value.var="sd")

########################################################################################

#sort out each bugs data object
names(bugs.data_ArtsDaten)

#site is the 5 x 5 km grid

#the thing labelled site is really a grid index
#the thing labelled grid is grid number mapping to the Norway raster

#we dont use "grid" in the model so remove
bugs.data_ArtsDaten <- bugs.data_ArtsDaten[-which(names(bugs.data_ArtsDaten)=="grid")]

#first change the site label to grid
names(bugs.data_ArtsDaten) <- gsub("site","grid",names(bugs.data_ArtsDaten))

###############################################################################

#now for the line transects
names(bugs.data_LineTransects)

names(bugs.data_LineTransects)<-sapply(names(bugs.data_LineTransects),
                                       function(x)
                                         paste0(x,"LT"))

bugs.data<-c(bugs.data_ArtsDaten,bugs.data_LineTransects)
names(bugs.data)

###############################################################################

#check survey fractions:
summary(bugs.data$TransectLengthLT)
temp <- bugs.data$TransectLengthLT/1000 * bugs.data$ESW.mean/1000 * 2
summary(temp)

#############################################################################
#pull out elements of the bugs data object that are used so far

bugs.dataS <- bugs.data[c("ngrid",'nyear',"nvisit","grid","year",
                           "Effort","adm","n.adm","n.adm2",
                           "adm2","y","NuIndivsLT","TransectLengthLT",
                           "bio1")]

save(bugs.data,file="bugs.data_simplified.RData")

######################################################################################

#also cut to smaller number of grids?
siteInfo_ArtsDaten$bio1 <- listlengthDF_SiteCovs$bio1[match(siteInfo_ArtsDaten$grid,
                                                            listlengthDF_SiteCovs$grid)]

selectedData <- subset(siteInfo_ArtsDaten,adm %in% unique(siteInfo_ArtsDaten$adm)[2:5])
years <- 6:9

#grid includes in both datasets
oGrids <- unique(subset(listlengthDF,!is.na(y))$grid)
aGrids <- unique(allDataObs@data$grid)
commonGrids <- intersect(oGrids,aGrids)
commonGrids <- intersect(commonGrids,selectedData$grid)
selectedData <- subset(selectedData,grid %in% commonGrids)

selectedData$s_siteIndex <- as.numeric(as.factor(selectedData$siteIndex))
selectedData$s_admN <- as.numeric(as.factor(selectedData$admN))
selectedData$s_admN2 <- as.numeric(as.factor(selectedData$admN2))
selected_listlengthDF <- subset(listlengthDF,
                                grid %in% selectedData$grid & yearIndex %in% years)
selected_listlengthDF$s_siteIndex <-  selectedData$s_siteIndex[match(selected_listlengthDF$grid,selectedData$grid)]

bugs.dataSS = list(ngrid=length(unique(selectedData$grid)),
                   nyear=length(unique(selected_listlengthDF$yearIndex)),            
                   nvisit=nrow(selected_listlengthDF),           
                   grid=selected_listlengthDF$s_siteIndex,             
                   year=as.numeric(as.factor(selected_listlengthDF$year)),
                   Effort=selected_listlengthDF$Effort,
                   y=selected_listlengthDF$y, 
                   adm=selectedData$s_admN,
                   n.adm=length(unique(selectedData$s_admN)),
                   n.adm2=length(unique(selectedData$s_admN2)),
                   adm2=selectedData$s_admN2,
                   NuIndivsLT=bugs.data$NuIndivsLT[selectedData$siteIndex,years],       
                   TransectLengthLT=bugs.data$TransectLengthLT[selectedData$siteIndex,years], 
                   bio1=selectedData$bio1)

######################################################################################

# bugs.dataS$Effort2 <- bugs.dataS$Effort^2
# summary(bugs.data$Effort)
# summary(bugs.data$Effort2)

# bugs.data$occDM <- model.matrix(~ bugs.data$tree_line_position + 
#                                   bugs.data$tree_line_position2+
#                                   bugs.data$bio1 + 
#                                   bugs.data$bio1_2 +
#                                   bugs.data$prefopen +
#                                   bugs.data$prefclosed +
#                                   bugs.data$open)[,-1]
# 
# covariateOrder<-c("tree line","quad tree line","mean temp", "quad mean temp","preferred open",
#                   "preferred forest","open")
#bugs.data$n.covs <- ncol(bugs.data$occDM)

params <- c("totAbund","int.d","mean.psi","mean.p",
            "meanZ","realAbund")

#specify model structure
out1 <- jags(bugs.dataSS, inits=inits, params, "models/hurdleModel_v2.txt", n.thin=nt,
             n.chains=nc, n.burnin=20000,n.iter=40000,parallel=T)
print(out1,2)

save(out1,file="model-outputs/out1_hurdle_v2.RData")

#######################################################################################

#look at predicted relative abundances
load("model-outputs/out1_hurdle_varable_exd_occ_abund_random_effort2.RData")
siteInfo<-selectedData
siteInfo$fits<-out1$mean$realAbund
mygrid[]<-0
mygrid[as.numeric(siteInfo$grid)]<-log(siteInfo$fits)
par(mfrow=c(1,2))
plot(mygrid,main="abundance")

#look at predicted meanZ
siteInfo$fits<-out1$mean$meanZ
mygrid[]<-0
mygrid[as.numeric(siteInfo$grid)]<-siteInfo$fits
plot(mygrid,main="occupancy")

########################################################################################