library(tidyverse)
library(sp)
library(rgeos)
library(raster)
library(maptools)

#local
#myfolder <- "data"

#HPC
myfolder <- "/data/idiv_ess/ptarmiganUpscaling" 

### get norway##############################################################

#using a m grid
equalM<-"+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
data(wrld_simpl)
Norway <- subset(wrld_simpl,NAME=="Norway")
NorwayOrig <- Norway
Norway <- gBuffer(Norway,width=1)
Norway <- spTransform(Norway,crs(equalM))
NorwayOrigProj <- spTransform(NorwayOrig,crs(equalM))

### ref grid ########################################################

#create grid
newres = 5000#5 km grid
mygrid <- raster(extent(projectExtent(Norway,equalM)))
res(mygrid) <- newres
mygrid[] <- 1:ncell(mygrid)
plot(mygrid)
gridTemp <- mygrid
myGridDF <- as.data.frame(mygrid,xy=T)

### listlength #####################################################

#source('formattingArtDatenBank_missing_data.R')

#read in list length object (made on the Rstudio server)
listlengthDF <- readRDS(paste(myfolder,"listlength_iDiv.rds",sep="/"))

#subset to May to September
listlengthDF <- subset(listlengthDF,is.na(month)|(month > 4 & month <10))

### subset ##########################################################

#subset to focal grids and those with environ covariate data

focusGrids <- readRDS(paste(myfolder,"focusGrids.rds",sep="/"))
varDF <- readRDS(paste(myfolder,"varDF_allEnvironData_5km_idiv.rds",sep="/"))
listlengthDF <- subset(listlengthDF,grid %in% focusGrids)
listlengthDF <- subset(listlengthDF,grid %in% varDF$grid)

### effort ######################################################################

#we will treat it as a categorical variable
table(listlengthDF$L)
listlengthDF$singleton <- ifelse(listlengthDF$L==1,1,0)

### Absences ###################################################################

#the dataset contains missing values

listlengthDF$L[is.na(listlengthDF$L)] <- 1 #set to nominal effort
#listlengthDF$y[listlengthDF$L==0] <- 0

### indices #####################################################################

#order data by site and year:
listlengthDF$siteIndex <- as.numeric(factor(listlengthDF$grid))
listlengthDF$yearIndex <- as.numeric(factor(listlengthDF$year))
listlengthDF <- arrange(listlengthDF,siteIndex,yearIndex)

#merge with environ data
names(listlengthDF)[which(names(listlengthDF)=="y")] <- "species"
listlengthDF <- merge(listlengthDF,varDF,by="grid",all.x=T)

#have missing singletons as 1
listlengthDF$singleton[is.na(listlengthDF$singleton)] <- 1

### adm inidices #############################################################

listlengthDF$admN <- as.numeric(factor(listlengthDF$adm))
listlengthDF$admN2 <- as.numeric(factor(listlengthDF$adm2))

#extract site data
siteInfo <- subset(listlengthDF,!duplicated(grid))
siteInfo_ArtsDaten <- siteInfo
#saveRDS(siteInfo,
#        file = "data/siteInfo_ArtsDaten.rds")

### BUGS object ################################################################

#for BUGS

bugs.data <- list(nsite = length(unique(listlengthDF$siteIndex)),
                  nyear = length(unique(listlengthDF$yearIndex)),
                  nvisit = nrow(listlengthDF),
                  site = listlengthDF$siteIndex,
                  year = listlengthDF$yearIndex,
                  Effort = listlengthDF$singleton,
                  y = listlengthDF$species,#includes NAs
                  #add an adm effect
                  adm = siteInfo$admN,
                  det.adm = listlengthDF$admN,
                  det.open = listlengthDF$Open,
                  n.adm = length(unique(siteInfo$admN)),#19
                  adm2 = siteInfo$admN2,
                  det.adm2 = listlengthDF$admN2,
                  n.adm2 = length(unique(siteInfo$admN2)))

#bugs.data_ArtsDaten <- bugs.data

#alpine_habitat:
#1= Open lowland, 
#2 = Low alpine zone, 
#3 = intermediate alpine zone, 
#4 = high alpine zone 

### initials ################################################

#get JAGS libraries
library(rjags)
library(jagsUI)

#need to specify initial values
zst <- reshape2::acast(listlengthDF, siteIndex~yearIndex, value.var="species",
                       fun=max,na.rm=T)
zst [is.infinite(zst)] <- 0
inits <- function(){list(z = zst)}

#saveRDS(zst,file="data/zst_ArtsDaten.rds")


### bpv index ##############################################

#for each i, sum into t
StrIdx <- array(data=0, dim = c(bugs.data$nvisit,bugs.data$nyear))
for(i in 1:bugs.data$nvisit){
  StrIdx[i,bugs.data$year[i]] <- 1
}

bugs.data$StrIdx <- StrIdx

### scale variables ##########################################

#total grids
nrow(varDF)#11788
#sampled grid
nrow(siteInfo)#11788

siteInfo[,-c(1:11)] <- plyr::numcolwise(scale)(siteInfo[,-c(1:11)])

### choose model ##############################################

modelTaskID <- read.delim(paste(myfolder,"modelTaskID_occuModel.txt",sep="/"),as.is=T)

#get task id
task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

#get model for this task
mymodel <- modelTaskID$Model[which(modelTaskID$TaskID==task.id)]

### choose covariates ##########################################################

#standard model
if(mymodel == "BUGS_occuModel_upscaling.txt"){
  
#specify model structure - a priori picking variables
bugs.data$occDM <- model.matrix(~ siteInfo$tree_line_position +
                                   I(siteInfo$tree_line_position^2) +
                                   siteInfo$y +
                                   siteInfo$bio1 +
                                   I(siteInfo$bio1^2) +
                                   siteInfo$Bog +
                                   siteInfo$Meadows +
                                   siteInfo$ODF +
                                   siteInfo$OSF +
                                   siteInfo$MountainBirchForest)[,-1]
#saveRDS(bugs.data,file="data/bugs.data_ArtsDaten.rds")

#lasso model or model selection model
}else{
  
bugs.data$occDM <- model.matrix(~ siteInfo$bio6 +
                                   siteInfo$bio5 +
                                   siteInfo$tree_line_position +
                                   siteInfo$MountainBirchForest +
                                   siteInfo$Bog +
                                   siteInfo$ODF + 
                                   siteInfo$Meadows +
                                   siteInfo$OSF +
                                   siteInfo$Mire +
                                   siteInfo$SnowBeds +
                                   siteInfo$y +
                                   siteInfo$distCoast +
                                   I(siteInfo$bio6^2) +
                                   I(siteInfo$bio5^2) +
                                   I(siteInfo$tree_line_position^2) +
                                   I(siteInfo$MountainBirchForest^2) +
                                   I(siteInfo$Bog^2) +
                                   I(siteInfo$ODF^2) + 
                                   I(siteInfo$Meadows^2) +
                                   I(siteInfo$OSF^2) +
                                   I(siteInfo$Mire^2) +
                                   I(siteInfo$SnowBeds^2) +
                                   I(siteInfo$y^2) +
                                   I(siteInfo$distCoast^2))[,-1]

}

bugs.data$n.covs <- ncol(bugs.data$occDM)

# myvars <- c("bio6","bio5","tree_line_position","MountainBirchForest","Bog","ODF","Meadows",
#             "OSF","Mire","SnowBeds","y","distCoast",
#             "bio6_2","bio5_2","tree_line_position_2",
#             "MountainBirchForest_2","Bog_2","ODF_2","Meadows_2","OSF_2","Mire_2",
#             "SnowBeds_2","y_2","distCoast_2")
# myvars <- c("tree_line_position","tree_line_position_2","geo_y","bio1","bio1_2","Bog",
#             "Meadows","ODF","OSF","MountainBirchForest")


### fit model ########################################################

params <- c("mean.p","beta","g","beta.effort",
            "beta.det.open","grid.z","grid.psi")

#chosen already earlier
#modelfile <- "/data/idiv_ess/ptarmiganUpscaling/BUGS_occuModel_upscaling.txt"
#modelfile <- "/data/idiv_ess/ptarmiganUpscaling/BUGS_occuModel_upscaling_ModelSelection.txt"
#modelfile <- "/data/idiv_ess/ptarmiganUpscaling/BUGS_occuModel_upscaling_LASSO.txt"
modelfile <- paste(myfolder,mymodel,sep="/")

#n.cores = as.integer(Sys.getenv("NSLOTS", "1")) 
n.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
#n.cores = 3

n.iterations = 20000

out1 <- jags(bugs.data, 
             inits = inits, 
             params, 
             modelfile, 
             n.thin = 10, 
             n.chains = n.cores, 
             n.burnin = round(n.iterations/2),
             n.iter = n.iterations,
             parallel = T)

saveRDS(out1$summary,file=paste0("outSummary_occModel_upscaling_",task.id,".rds"))

#update by a small number and get full model
out2 <- update(out1, parameters.to.save = c("mid.psi","mid.z","Py", 
                                            "e.count","sim.count",
                                            "fit","fit.new"),n.iter=1000)

saveRDS(out2,file=paste0("out_update_occModel_upscaling_",task.id,".rds"))

