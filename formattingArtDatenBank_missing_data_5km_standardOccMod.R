########################################################################

#load in general functions

source('C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/generalFunctions.R')

#get norway##############################################################

library(raster)
library(maptools)
library(rgeos)
data(wrld_simpl)
Norway<- subset(wrld_simpl,NAME=="Norway")
plot(Norway)
Norway <- gBuffer(Norway,width=1)

########################################################################

#choose modelling resolution:

#newres=1

newres=0.05

#########################################################################

#get willow ptarmigan data - occurrence data

#get GBIF data
tdir<-"C:/Users/diana.bowler/OneDrive - NINA/Alpine/GBIF"
lirype <- read.delim(paste(tdir,"willow_ptarmigan_GBIF.txt",sep="/"),as.is=T,row.names=NULL)
lirype$decimalLatitude<-as.numeric(lirype$decimalLatitude)
lirype$decimalLongitude<-as.numeric(lirype$decimalLongitude)
lirype<-subset(lirype,!is.na(lirype$decimalLatitude)|!is.na(lirype$decimalLongitude))
lirype<-lirype[,c("decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters",
                  "name","year","month","day","recordedBy")]

#get artsdatenbanken data
tdir<-"C:/Users/diana.bowler/OneDrive - NINA/Alpine/SpeciesMapServices/Ptarmigan"
lirype2 <- read.delim(paste(tdir,"DataExportFraArtskart.txt",sep="/"),skipNul=T,dec=",")
lirype2$Species<-apply(lirype2,1,function(x)paste(x["Genus"],x["Species"],sep=" "))
lirype2<-lirype2[,c("Latitude","Longitude","CoordinatePrecision",
                  "Species","YearCollected","MonthCollected","DayCollected","Collector")]
names(lirype2)<-names(lirype)

#combine them
lirype<-rbind(lirype,lirype2)
lirypeD<-lirype

#remove duplicates
lirype$date<-paste(lirype$year,lirype$month,lirype$day,sep="-")
lirype$id<-paste(lirype$date,round(lirype$decimalLatitude,digits=4),round(lirype$decimalLongitude,digits=4),sep="_")
lirype<-subset(lirype,!duplicated(id))

#remove dodgy records
summary(lirype$coordinateUncertaintyInMeters)
lirype <- subset(lirype, !coordinateUncertaintyInMeters>10000|is.na(coordinateUncertaintyInMeters))

#remove those without 2 decimal places
lirype$latDP<-getDecimalPlaces(lirype$decimalLatitude)
lirype$lonDP<-getDecimalPlaces(lirype$decimalLongitude)
lirype<-subset(lirype,latDP>2 & lonDP>2)

#subsetting
lirype<-subset(lirype,year>2006&year<2018)
#exclude Nov to March  - rock and willow both have white coats#
lirype<-subset(lirype,month >3 & month <11)
table(lirype$year)

#make spatial
coordinates(lirype)<-c("decimalLongitude","decimalLatitude")
plot(Norway)
plot(lirype,add=T)
#all points look pretty good!! few in the sea

#overlay to the grid
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
gridSummary<-ddply(lirype@data,.(grid),summarise,nuRecs=length(name))
mygrid[]<-0
mygrid[gridSummary$grid]<-gridSummary$nuRecs
length(gridSummary$grid)==length(gridSummary$nuRecs)
plot(mygrid)

#how many records per site are there
hist(gridSummary$nuRecs)
gridSummary$RepeatedVisits<-sapply(gridSummary$nuRecs,function(x)ifelse(x>1,1,0))
table(gridSummary$RepeatedVisits)
#0    1 
#1654 1662

#grid cell per year
gridSummary<-ddply(lirype@data,.(grid,year),summarise,nuRecs=length(name))
gridSummary$RepeatedVisits<-sapply(gridSummary$nuRecs,function(x)ifelse(x>1,1,0))
table(gridSummary$RepeatedVisits)
#0    1 
#4726 2435

##########################################################################

#get data for all birds - as an effort layer

#first time

#get GBIF data
tdir<-"C:/Users/diana.bowler/OneDrive - NINA/Alpine/GBIF"
speciesFiles<-list.files(tdir)[grepl("all_birds",list.files(tdir))]
library(plyr)
allbirds<-ldply(speciesFiles,function(x){
  read.delim(paste(tdir,x,sep="/"),as.is=T,row.names=NULL)
})

allbirds<-allbirds[,c("decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters",
                      "name","year","month","day","recordedBy")]
allbirds$decimalLatitude<-as.numeric(allbirds$decimalLatitude)
allbirds$decimalLongitude<-as.numeric(allbirds$decimalLongitude)
allbirds$year<-as.numeric(allbirds$year)
allbirds$month<-as.numeric(allbirds$month)
allbirds<-subset(allbirds,!is.na(allbirds$decimalLatitude))
allbirds<-subset(allbirds,!is.na(allbirds$decimalLongitude))
               
#get artsdatenbanken data
tdir<-"C:/Users/diana.bowler/OneDrive - NINA/Alpine/SpeciesMapServices/Birds"
allbirds2 <- read.delim(paste(tdir,"DataExportFraArtskart.txt",sep="/"),skipNul=T,dec=",",as.is=T)
allbirds2$Species<-apply(allbirds2,1,function(x)paste(x["Genus"],x["Species"],sep=" "))
allbirds2<-allbirds2[,c("Latitude","Longitude","CoordinatePrecision",
                    "Species","YearCollected","MonthCollected","DayCollected","Collector")]
names(allbirds2)<-names(allbirds)

#combine them
allbirds<-rbind(allbirds,allbirds2)

#also with ptarmigan data
allbirds<-rbind(allbirds,lirypeD)

#remove duplicates
allbirds$date<-paste(allbirds$year,allbirds$month,allbirds$day,sep="-")
allbirds$id<-paste(allbirds$name,allbirds$date,round(allbirds$decimalLatitude,digits=4),
                   round(allbirds$decimalLongitude,digits=4),sep="_")
allbirds<-subset(allbirds,!duplicated(id))

#remove dodgy records
summary(allbirds$coordinateUncertaintyInMeters)
allbirds <- subset(allbirds, !coordinateUncertaintyInMeters>10000|is.na(coordinateUncertaintyInMeters))
allbirds<-subset(allbirds,!is.na(allbirds$decimalLatitude)|!is.na(allbirds$decimalLongitude))

#remove those without 2 decimal places
allbirds$latDP<-getDecimalPlaces(allbirds$decimalLatitude)
allbirds$lonDP<-getDecimalPlaces(allbirds$decimalLongitude)
allbirds<-subset(allbirds,latDP>2 & lonDP>2)

#extract survyed grids within time period of interest
allbirds<-subset(allbirds,month >3 & month <11)
allbirds<-subset(allbirds,year>2006&year<2018)
table(allbirds$year)

#sort out observer ids
allRecorders<-names(table(allbirds$recordedBy))
#if comma separate the names
allRecorders<-sort(unique(trim(unlist(lapply(allRecorders,function(x)strsplit(x,",")[[1]])))))
#work out how many times each has contributed
effort<-as.numeric()
for(i in 1:length(allRecorders)){
  temp<-allbirds[grepl(allRecorders[i],allbirds$recordedBy),]
  effort[i]<-length(unique(interaction(temp$month,temp$year)))
}
bestRecorders<-data.frame(allRecorders,effort)  
bestRecorders<-bestRecorders[order(bestRecorders$effort),]
#remove those without spaces in the name (i.e., first name and surname)
bestRecorders<-bestRecorders[grepl("\\ ",bestRecorders$allRecorders),]
summary(bestRecorders$effort)
bestRecorders<-subset(bestRecorders,effort>52)#upper quartile  
allbirds$Best<-sapply(allbirds$recordedBy,function(x){
  all<-strsplit(x,",")[[1]]
  all<-trim(all)
  any(is.element(all,bestRecorders$allRecorders))
})
table(allbirds$Best)
allbirds<-subset(allbirds,Best=="TRUE")


save(allbirds,file="allbirds_uniqueRecords.RData")

#####################################################################################

#subsequent times
load("allbirds_uniqueRecords.RData")
names(allbirds)[which(names(allbirds)=="name")]<-"Species"

#####################################################################################

#make spatial
coordinates(allbirds)<-c("decimalLongitude","decimalLatitude")
proj4string(allbirds)<-crs(Norway)
#plot(Norway)
#plot(allbirds,add=T)
#all points look pretty good!! few in the sea...

#exclude points beyond the mask
out<-over(allbirds,Norway)
allbirds<-allbirds[!is.na(out),]

#overlay to the grid
library(raster)
mygrid<-raster(extent(Norway))
res(mygrid)<-newres
mygrid[]<-1:ncell(mygrid)
plot(mygrid)
#get grid number for each point
allbirds$grid<-extract(mygrid,allbirds)

#identify all grid cells within the buffer region
library(rgeos)
mygridMask<-mask(mygrid,Norway)
plot(mygridMask)
mygridMaskDF<-as.data.frame(mygridMask,xy=T)
mygridMaskDF<-subset(mygridMaskDF,!is.na(layer))

#recoganise
allbirds<-subset(allbirds,!is.na(Species))
allbirds_data<-data.frame(allbirds@data,allbirds@coords) 
allbirds_data <- allbirds_data[,c("year","month","day","grid","Species")]
all_data <- subset(allbirds_data,!is.na(grid))

#########################################################################

#get info on administrative names for the grid
library(rgdal)
NorwayADM<-readOGR(dsn="C:/Users/diana.bowler/OneDrive - NINA/Alpine/NOR_adm",layer="NOR_adm2")
crs(NorwayADM)
mygridDF<-as.data.frame(mygrid,xy=T)
mygridPoints<-mygridDF
coordinates(mygridPoints)<-c("x","y")
proj4string(mygridPoints)<-crs(NorwayADM)
myAdm<-over(mygridPoints,NorwayADM)
myAdm$grid<-mygridPoints@data$layer

myAdm$VARNAME_1<-as.character(myAdm$NAME_1)
myAdm$VARNAME_1[is.na(myAdm$VARNAME_1)]<-"outside"
myAdm$VARNAME_2<-as.character(myAdm$NAME_2)
myAdm$VARNAME_2[is.na(myAdm$VARNAME_2)]<-"outside"
table(myAdm$VARNAME_1)
table(myAdm$VARNAME_2)
mygrid[myAdm$grid]<-as.numeric(as.factor(myAdm$VARNAME_1))
plot(mygrid)#looks good!
mygrid[myAdm$grid]<-as.numeric(as.factor(myAdm$VARNAME_2))
plot(mygrid)#looks good!

#########################################################################

#get list length info

#add list length of all birds on each sampling visit
listlengthDF<-ddply(all_data,.(year,month,day,grid),summarise,
                    L = length(unique(Species[!is.na(Species)])), #number of species
                    L2 = length(Species[!is.na(Species)]), #number of records
                    L3 = L2/L) #records per species

summary(listlengthDF$L)
#cap list length
table(listlengthDF$L)
listlengthDF$L[listlengthDF$L>quantile(listlengthDF$L,0.975,na.rm=T)]<-quantile(listlengthDF$L,0.975,na.rm=T)

#other measures
hist(listlengthDF$L3)
listlengthDF$L2[listlengthDF$L2>quantile(listlengthDF$L2,0.975,na.rm=T)]<-quantile(listlengthDF$L2,0.975,na.rm=T)
listlengthDF$L3[listlengthDF$L3>quantile(listlengthDF$L3,0.975,na.rm=T)]<-quantile(listlengthDF$L3,0.975,na.rm=T)

#get total number of species ever seen per grid
gridRichness <- ddply(all_data,.(grid),summarise,nuSpecies=length(unique(Species[!is.na(Species)])),nuRecs=length(Species[!is.na(Species)]))
mygrid[]<-0
gridRichness$nuRecs[gridRichness$nuRecs>quantile(gridRichness$nuRecs,0.975)]<-quantile(gridRichness$nuRecs,0.975)
mygrid[gridRichness$grid]<-gridRichness$nuRecs/gridRichness$nuSpecies
plot(mygrid)
mygrid[gridRichness$grid]<-gridRichness$nuRecs
plot(mygrid)

#add species richnes to the listlength df
listlengthDF$nuSpecies <- gridRichness$nuSpecies[match(listlengthDF$grid,gridRichness$grid)]
summary(listlengthDF$nuSpecies)

##################################################################################

#sort out occupancy matrix of ptarmigan
listlengthDF$visit<-paste(listlengthDF$year,listlengthDF$month,
                          listlengthDF$day,listlengthDF$grid,sep="-")
lirype$visit<-paste(lirype$year,lirype$month,
                     lirype$day,lirype$grid,sep="-")
listlengthDF$y<-sapply(listlengthDF$visit,function(x)ifelse(x%in%lirype$visit,1,0))

table(listlengthDF$y)
#0      1 
#604912   6063
  
#add NAs for all years and grids
newgrid<-expand.grid(year=unique(listlengthDF$year),
                     grid=sort(unique(mygridMaskDF$layer)))
listlengthDF<-merge(listlengthDF,newgrid,by=c("year","grid"),all=T)
summary(listlengthDF)
table(listlengthDF$grid,listlengthDF$year)

#########################################
#adding non-detections??
#listlengthDF$L[is.na(listlengthDF$L)]<-0
#listlengthDF$y[listlengthDF$L==0]<-0
#########################################
#######################################################################

#order data by site and year
listlengthDF<-arrange(listlengthDF,year,grid)

#add indices
listlengthDF$siteIndex<-as.numeric(as.factor(listlengthDF$grid))
listlengthDF$yearIndex<-as.numeric(as.factor(listlengthDF$year))

#add adm
listlengthDF$adm<-myAdm$VARNAME_1[match(listlengthDF$grid,myAdm$grid)]
table(listlengthDF$adm)
listlengthDF$admN<-as.numeric(factor(listlengthDF$adm))

listlengthDF$adm2<-myAdm$VARNAME_2[match(listlengthDF$grid,myAdm$grid)]
table(listlengthDF$adm2)
listlengthDF$admN2<-as.numeric(factor(listlengthDF$adm2))

#remove data not falling into an admN or admN2
#listlengthDF<-subset(listlengthDF,adm2!="outside")

######################################################################################################

#combine environmental and pop data:

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling")
load("varDF_allEnvironData_5km.RData")
varDF<-subset(varDF,!duplicated(grid))
listlengthDF<-merge(listlengthDF,varDF,by="grid",sort=F)

#check where we have data
mygrid[]<-0
mygrid[listlengthDF$grid]<-listlengthDF$admN
plot(mygrid)

#if testing, reduce the number of sites - get number of times each grid was visited
#gridSummary <- ddply(listlengthDF,.(grid),summarize,nuVisits=length(y[!is.na(y)]))
#listlengthDF <- subset(listlengthDF,grid%in%gridSummary$grid[gridSummary$nuVisits>0])

#add indices
listlengthDF$site<-paste(listlengthDF$adm2,listlengthDF$grid,sep="-")
listlengthDF$siteIndex<-as.numeric(factor(listlengthDF$grid))
listlengthDF$yearIndex<-as.numeric(factor(listlengthDF$year))

#order data by site and year
listlengthDF<-arrange(listlengthDF,yearIndex,siteIndex)

#extract site data
listlengthDF_SiteCovs<-subset(listlengthDF,!duplicated(grid))

#adm data
listlengthDF$admN<-as.numeric(factor(listlengthDF$adm))
listlengthDF$admN2<-as.numeric(factor(listlengthDF$adm2))
siteInfo<-unique(listlengthDF[,c("siteIndex","admN","grid","admN2")])

####################################################################################################
#use for other model
####################################################################################################
#organise occupancy data into a matrix

#add repeat survey data
listlengthDF<-ddply(listlengthDF,.(siteIndex,yearIndex),function(x){
  x$visitNu<-sample(1:nrow(x))
  x<-x[order(x$visitNu),]
  return(subset(x,visitNu<=10))
})
summary(listlengthDF$visitNu)

#convert into a matrix
library(reshape2)
occupancyMatrix<-acast(listlengthDF,siteIndex~year~visitNu,value.var="y")
dim(occupancyMatrix)
#[1] 28081    11    10


#check row names makes site data
all(dimnames(occupancyMatrix)[[1]]==siteInfo$siteIndex)
all(dimnames(occupancyMatrix)[[1]]==listlengthDF_SiteCovs$siteIndex)

############################################################################

#also get number of unique records for each year/site as another measure of effort

siteDF<-ddply(all_data,.(year,grid),summarise,
                    L = length(unique(Species[!is.na(Species)])), #number of species
                    L2 = length(Species[!is.na(Species)]), #number of records
                    L3 = L2/L,
                    nuDays = ) #records per species

#############################################################################
#explore detection prob by hand

#per year - for each line, what fraction of surveys had a ptarmigan in it
#get rid of NAs

listlengthDF_NAfree<-subset(listlengthDF,!is.na(visit))
propSuccess<-ddply(listlengthDF_NAfree,.(grid,adm,adm2,year),
                   summarise,
                   propY=ifelse(length(y)>1,mean(y),NA),#want to calculate it only when there are repeat visits
                   meanL=mean(L),
                   meanL2=mean(L2),
                   meanL3=mean(L3))
propSuccess<-subset(propSuccess,!is.na(propY))

hist(propSuccess$propY)#very all or nothing
mean(propSuccess$propY[propSuccess$propY>0])
#0.3335538

#checking
propSuccess<-subset(propSuccess,propY!=0)
all_data<-arrange(all_data,year,month,day,grid)

subset(all_data, grid==11256 & year==2015)

subset(all_data, grid==11789 & year==2012)

subset(all_data, grid==11802 & year==2013)

subset(all_data, grid==12926 & year==2017)

subset(listlengthDF,grid==12926 & year==2017)

subset(allbirds@data,grid==12926 & year==2017)

####################################################################################################

#for BUGS

bugs.data <- list(nsite = dim(occupancyMatrix)[1],
                  nyear = dim(occupancyMatrix)[2],
                  nvisit = dim(occupancyMatrix)[3],
                  grid = siteInfo$grid,
                  site = siteInfo$siteIndex,
                  year = 1:length(dim(occupancyMatrix)[2]),
                  #L = listlengthDF$L,
                  #nuSpecies = listlengthDF$nuSpecies,
                  y = occupancyMatrix,
                  #add an adm effect
                  adm = siteInfo$admN,
                  n.adm = length(unique(siteInfo$admN)),
                  adm2 = siteInfo$admN2,
                  n.adm2 = length(unique(siteInfo$admN2)),
                  habitat = ifelse(listlengthDF_SiteCovs$habitat=="Forest",1,0),
                  bio1 = scale(listlengthDF_SiteCovs$bio1),
                  bio1_2 = listlengthDF_SiteCovs$bio1^2,
                  bio6 = listlengthDF_SiteCovs$bio6,
                  bio5 = listlengthDF_SiteCovs$bio5,
                  bio5_2 = listlengthDF_SiteCovs$bio5^2,
                  forest = listlengthDF_SiteCovs$Forest,
                  open = scale(listlengthDF_SiteCovs$Open),
                  top = log(listlengthDF_SiteCovs$Top+1),
                  alpine_habitat1 = listlengthDF_SiteCovs$alpine_habitat1,
                  alpine_habitat2 = log(listlengthDF_SiteCovs$alpine_habitat2+1),
                  alpine_habitat3 = log(listlengthDF_SiteCovs$alpine_habitat3+1),
                  alpine_habitat4 = log(listlengthDF_SiteCovs$alpine_habitat4+1),
                  tree_line_position = scale(listlengthDF_SiteCovs$tree_line_position),
                  tree_line_position2 = scale(listlengthDF_SiteCovs$tree_line_position^2))

#alpine_habitat:
#1= Open lowland, 
#2 = Low alpine zone, 
#3 = intermediate alpine zone, 
#4 = high alpine zone 

#need to specify initial values
library(reshape2)
zst <- acast(listlengthDF, siteIndex~year, value.var="y",fun=max,na.rm=T)
zst [is.infinite(zst)] <- 0
inits <- function(){list(z = zst)}
bugs.data$n.covs <- 1

#get BUGS functions
source('C:/Users/diana.bowler/OneDrive - NINA/methods/models/bugsFunctions.R')

##########################################################################################

#costant detection proability and occupancy

#specify parameters to monitor
params <- c("mean.p","mean.psi")

#specify model structure
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "standard_occModel.txt", n.thin=10,
             n.chains=3, n.burnin=600,n.iter=2000,parallel=T)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_OM_random_missing_5km_standard.RData")

##########################################################################################

#adm and adm2 effects on occupancy and detection

#specify parameters to monitor
params <- c("mean.p","mean.psi","random.adm.sd","random.adm2.sd","random.adm.p.sd","random.adm2.p.sd")

#specify model structure
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "standard_occModel.txt", n.thin=10,
             n.chains=3, n.burnin=300,n.iter=1000,parallel=T)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_OM_random_missing_5km_standard_adm.RData")

########################################################################################

#test models on a subset of the data

#######################################################################################

#constant detection, random adm effects

#specify parameters to monitor
params <- c("mean.p","mean.psi","random.adm")

#specify model structure
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "standard_occModel.txt", n.thin=10,
             n.chains=3, n.burnin=600,n.iter=2000,parallel=T)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_OM_random_missing_5km_standard1.RData")

#######################################################################################

#constant detection, random adm and ad2 effects

#specify parameters to monitor
params <- c("mean.p","mean.psi","random.adm","random.adm2")

#specify model structure
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "standard_occModel.txt", n.thin=10,
             n.chains=3, n.burnin=600,n.iter=2000,parallel=T)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_OM_random_missing_5km_standard2.RData")

#######################################################################################

#run model including explanatory variables on observation

#specify model structure
bugs.data$occDM <- model.matrix(~ bio1 + open)[,-1]
bugs.data$n.covs <- ncol(bugs.data$occDM)


params <- c("mean.p","mean.psi","beta","random.adm","a","eta")

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "standard_occModel.txt", n.thin=10,
             n.chains=3, n.burnin=600,n.iter=2000,parallel=T)

print(out1,2)
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_OM_variables_missing_5km_standard3.RData")
colnames(bugs.data$occDM)

##########################################################################################

#run model including explanatory variables on observation

#specify model structure
bugs.data$occDM <- model.matrix(~ bugs.data$tree_line_position + bugs.data$tree_line_position2+
                                  bugs.data$bio1 + bugs.data$open)[,-1]
bugs.data$n.covs <- ncol(bugs.data$occDM)


params <- c("mean.p","mean.psi","beta")

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "standard_occModel.txt", n.thin=10,
             n.chains=3, n.burnin=600,n.iter=2000,parallel=T)

print(out1,2)
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_OM_variables_missing_5km_standard4.RData")
colnames(bugs.data$occDM)

##########################################################################################

#use model to predict to whole coutry
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
load("out1_OM_variables_missing_5km_standard4.RData")

varDF$tree_line_position<-scale(varDF$tree_line_position)
varDF$bio1<-scale(varDF$bio1)
varDF$Open<-scale(varDF$Open)

varDF$psi<-apply(varDF,1,function(x){
  invlogit(logit(0.3127) + as.numeric(x["tree_line_position"]) * 0.3065 -
             as.numeric(x["tree_line_position"])^2 * 2.5616 -
             as.numeric(x["bio1"]) * 0.5410 + 
             as.numeric(x["Open"]) * 0.5622)
  
})

mygrid[]<-0
mygrid[varDF$grid]<-varDF$psi
plot(mygrid)

##########################################################################################

#plot occupancy model
out2<-update(out1,parameters.to.save=c("psi","z"),n.iter=100)
plotZ(out2,param="z")
plotZ(out2,param="psi")
plotZerror(out2)

##########################################################################################



