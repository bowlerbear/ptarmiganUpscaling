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

#source('formattingArtDatenBank_missing_data_iDiv.R')

#read in list length object (made on the Rstudio server)
listlengthDF <- readRDS(paste(myfolder,"listlength_iDiv.rds",sep="/"))

#subset to May to September
listlengthDF <- subset(listlengthDF,is.na(month)|(month > 4 & month <10))

#2008 to 2017....
listlengthDF$Year <- listlengthDF$year
listlengthDF <- subset(listlengthDF, Year>2007 & Year <=2017) 

### subset ##########################################################

#subset to focal grids and those with environ covariate data

focusGrids <- readRDS(paste(myfolder,"focusGrids.rds",sep="/"))
varDF <- readRDS(paste(myfolder,"varDF_allEnvironData_5km_idiv.rds",sep="/"))
listlengthDF <- subset(listlengthDF,grid %in% focusGrids)
listlengthDF <- subset(listlengthDF,grid %in% varDF$grid)

### plot data #######################################################

# #presences
# 
# listlengthDF_Positive <- subset(listlengthDF, y==1)
# listlengthDF_Positive <- subset(listlengthDF_Positive, !duplicated(grid))
# 
# mygrid[] <- NA
# mygrid[listlengthDF_Positive $grid] <- listlengthDF_Positive$y
# crs(mygrid) <- equalM
# occ_tmap <- tm_shape(NorwayOrigProj) +
#   tm_borders() +
# tm_shape(mygrid)+
#   tm_raster(title="Occupancy", legend.show = FALSE)
# saveRDS(occ_tmap,"plots/occupancy_raw_presences.RDS")
# 
# #absences
# listlengthDF_Abs<- subset(listlengthDF, y==0)
# listlengthDF_Abs <- subset(listlengthDF_Abs, !duplicated(grid))
# 
# mygrid[] <- NA
# mygrid[listlengthDF_Abs$grid] <- listlengthDF_Abs$y
# crs(mygrid) <- equalM
# occ_tmap <- tm_shape(NorwayOrigProj) + 
#   tm_borders() +
#   tm_shape(mygrid)+
#   tm_raster(title="Occupancy", palette="YlGnBu", legend.show = FALSE)
# saveRDS(occ_tmap,"plots/occupancy_raw_absences.RDS")

### Absences ###################################################################

#the dataset contains missing values

listlengthDF$L[is.na(listlengthDF$L)] <- 1 #set to nominal effort
#listlengthDF$y[listlengthDF$L==0] <- 0
table(is.na(listlengthDF$y))

### effort ######################################################################

#we will treat it as a categorical variable
table(listlengthDF$L)
listlengthDF$singleton <- ifelse(listlengthDF$L==1,1,0)
listlengthDF$short <- ifelse(listlengthDF$L %in% c(2:4),1,0)

### indices #####################################################################

#order data by site and year:
listlengthDF$siteIndex <- as.numeric(factor(listlengthDF$grid))
listlengthDF$yearIndex <- as.numeric(factor(listlengthDF$year))
listlengthDF <- arrange(listlengthDF,siteIndex,yearIndex)

#merge with environ data
names(listlengthDF)[which(names(listlengthDF)=="y")] <- "species"
listlengthDF <- merge(listlengthDF,varDF,by="grid",all.x=T)

### adm inidices #############################################################

listlengthDF$admN <- as.numeric(factor(listlengthDF$adm))
listlengthDF$admN2 <- as.numeric(factor(listlengthDF$adm2))

#extract site data
siteInfo <- subset(listlengthDF,!duplicated(grid))
siteInfo_ArtsDaten <- siteInfo
#saveRDS(siteInfo,
#        file = "data/siteInfo_ArtsDaten.rds")


### data summary ##############################################################

length(unique(listlengthDF$grid))
# 
# listlengthDF_obs <- subset(listlengthDF, !is.na(species))
# length(unique(listlengthDF_obs$grid))
# 
# listlengthDF_posobs <- subset(listlengthDF, species==1)
# length(unique(listlengthDF_posobs$grid))
# 
# #when a ptarmigan was observed at a grid - how many times was it seen there
# ptarmiganGrids <- listlengthDF_posobs %>%
#                   group_by(grid) %>%
#                   summarise(nuVisits = length(unique(visit)),
#                             nuYears = length(unique(year)))
# summary(ptarmiganGrids$nuVisits)
# summary(ptarmiganGrids$nuYears)
# 
# #how often was a grid surveyed
# visitGrids <- listlengthDF_obs %>%
#   group_by(grid) %>%
#   summarise(nuVisits = length(unique(visit)),
#             nuYears = length(unique(year)))
# summary(visitGrids$nuVisits)
# summary(visitGrids$nuYears)

### BUGS object ################################################################

myScale <- function(x){
  
  as.numeric(scale(x))
  
}

#for BUGS

bugs.data <- list(nsite = length(unique(listlengthDF$siteIndex)),
                  nyear = length(unique(listlengthDF$yearIndex)),
                  nvisit = nrow(listlengthDF),
                  site = listlengthDF$siteIndex,
                  year = listlengthDF$yearIndex,
                  y = listlengthDF$species,#includes NAs
                  #detection covariates
                  Effort = listlengthDF$singleton,
                  Effort2 = listlengthDF$short,
                  det.tlp = myScale(listlengthDF$tree_line/100),
                  det.tlp2 = myScale(listlengthDF$tree_line^2/100000),
                  det.open = myScale(listlengthDF$Open),
                  det.bio5 = myScale(listlengthDF$bio5/100),
                  det.bio6 = myScale(listlengthDF$bio6/100),
                  #add an adm effect
                  adm = siteInfo$admN,
                  det.adm = listlengthDF$admN,
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
names(siteInfo)[1:12]

siteInfo[,-c(1:12)] <- plyr::numcolwise(scale)(siteInfo[,-c(1:12)])

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
                                   siteInfo$bio6 +
                                   siteInfo$Bog +
                                   siteInfo$Mire +
                                   siteInfo$Meadows +
                                   siteInfo$ODF +
                                   I(siteInfo$ODF^2) +
                                   siteInfo$OSF +
                                   I(siteInfo$OSF^2) +
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


### save objects ####################################################

# saveRDS(bugs.data, file="data/bugs.data_ArtsDaten.rds")
# saveRDS(zst, file="data/zst_ArtsDaten.rds")

### fit model ########################################################

params <- c("beta","g",
            "average.p","average.psi","propOcc",
            "beta.effort","beta.effort2",
            "beta.det.open","beta.det.bio5","beta.det.bio6",
            "beta.det.tlp","beta.det.tlp2",
            "bpv","grid.z","grid.psi")

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

### get AUC ##########################################################

#update by a small number and get full model
out2 <- update(out1, parameters.to.save = c("mid.psi","mid.z","Py","Py.pred"),n.iter=2000,
               n.thin = 10)

library(ggmcmc)
ggd2 <- ggs(out2$samples)

#Z against psi
Preds <- subset(ggd2,grepl("mid.psi",ggd2$Parameter))
Preds$Iteration <- as.numeric(interaction(Preds$Iteration,Preds$Chain))
Z_Preds <- subset(ggd2,grepl("mid.z",ggd2$Parameter))
Z_Preds$Iteration <- as.numeric(interaction(Z_Preds$Iteration,Z_Preds$Chain))
nu_Iteractions <- max(Preds$Iteration)
head(Preds)

#look through all iterations
AUC_psi <- rep(NA, nu_Iteractions)

for (i in 1:nu_Iteractions){
  
  psi.vals <- Preds$value[Preds$Iteration==i]
  z.vals <- Z_Preds$value[Z_Preds$Iteration==i]
  
  pred <- ROCR::prediction(psi.vals, z.vals)
  
  #get AUC
  perf <- ROCR::performance(pred, "auc")
  AUC_psi[i] <- perf@y.values[[1]]

}

summary(AUC_psi)
saveRDS(summary(AUC_psi),file=paste0("AUC_psi_occModel_upscaling_",task.id,".rds"))


#Y against Py - pull out one year
Py_preds <- subset(ggd2,grepl("Py",ggd2$Parameter))
Py_preds <- subset(Py_preds,!grepl("Py.pred",Py_preds$Parameter))
Py_preds$Iteration <- as.numeric(interaction(Py_preds$Iteration,Py_preds$Chain))
nu_Iteractions <- max(Py_preds$Iteration)
head(Py_preds)

#look through all iterations
AUC_py <- rep(NA, nu_Iteractions)
              
for (i in 1:nu_Iteractions){
  
  py.vals <- Py_preds$value[Py_preds$Iteration==i]
  
  #select data
  useData <- !is.na(bugs.data$y) & bugs.data$year==6
  
  pred <- ROCR::prediction(py.vals[useData],bugs.data$y[useData])
  
  #get AUC
  perf <- ROCR::performance(pred, "auc")
  AUC_py[i] <- perf@y.values[[1]]
  
}

summary(AUC_py)
saveRDS(summary(AUC_py),file=paste0("AUC_py_occModel_upscaling_",task.id,".rds"))

#Y against Py_pred (psi x p) - pull out one year
Py_preds <- subset(ggd2,grepl("Py.pred",ggd2$Parameter))
Py_preds$Iteration <- as.numeric(interaction(Py_preds$Iteration,Py_preds$Chain))
nu_Iteractions <- max(Py_preds$Iteration)
head(Py_preds)

#look through all iterations
AUC_py <- rep(NA, nu_Iteractions)

for (i in 1:nu_Iteractions){
  
  py.vals <- Py_preds$value[Py_preds$Iteration==i]
  
  #select data
  useData <- !is.na(bugs.data$y) & bugs.data$year==6
  
  pred <- ROCR::prediction(py.vals[useData],bugs.data$y[useData])
  
  #get AUC
  perf <- ROCR::performance(pred, "auc")
  AUC_py[i] <- perf@y.values[[1]]
  
}

summary(AUC_py)

saveRDS(summary(AUC_py),file=paste0("AUC_pypred_occModel_upscaling_",task.id,".rds"))

### get full z and psi ##########################################################

out2 <- update(out1, parameters.to.save = c("z","psi"),n.iter=2000, n.thin = 10)

#save summary
saveRDS(out2$summary,file=paste0("Z_summary_occModel_upscaling_",task.id,".rds"))

#save z and psi samples
library(ggmcmc)
ggd <- ggs(out2$samples)

ggd2 <- subset(ggd,grepl("z",ggd$Parameter))
saveRDS(ggd2,file=paste0("Z_occModel_upscaling_",task.id,".rds"))

ggd3 <- subset(ggd,grepl("psi",ggd$Parameter))
saveRDS(ggd3,file=paste0("psi_occModel_upscaling_",task.id,".rds"))
