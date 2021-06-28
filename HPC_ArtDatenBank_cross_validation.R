library(sp)
library(raster)
library(maptools)
library(ggplot2)
library(rgeos)
library(plyr)

#HPC
myfolder <- "/data/idiv_ess/ptarmiganUpscaling" 
#local
#myfolder <- "data"

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
newres=5000#5 km grid
mygrid<-raster(extent(projectExtent(Norway,equalM)))
res(mygrid) <- newres
mygrid[] <- 1:ncell(mygrid)
plot(mygrid)
gridTemp <- mygrid
myGridDF <- as.data.frame(mygrid,xy=T)

### listlength #####################################################

#source('formattingArtDatenBank_missing_data.R')

#read in list length object (made on the Rstudio server)
listlengthDF <- readRDS(paste(myfolder,"listlength_iDiv.rds",sep="/"))
names(listlengthDF)[which(names(listlengthDF)=="y")] <- "species"

#subset to May to September
listlengthDF <- subset(listlengthDF,is.na(month)|(month > 4 & month <10))

### subset ##########################################################

#subset to focal grids and those with environ covariate data

focusGrids <- readRDS(paste(myfolder,"focusGrids.rds",sep="/"))
varDF <- readRDS(paste(myfolder,"varDF_allEnvironData_5km_idiv.rds",sep="/"))
listlengthDF <- subset(listlengthDF,grid %in% focusGrids)
listlengthDF <- subset(listlengthDF,grid %in% varDF$grid)

#merge with environ data
listlengthDF <- merge(listlengthDF,varDF,by="grid",all.x=T)

#adm indices
listlengthDF$admN <- as.numeric(factor(listlengthDF$adm))
listlengthDF$admN2 <- as.numeric(factor(listlengthDF$adm2))

### effort ######################################################################

#we will treat it as a categorical variable
table(listlengthDF$L)
listlengthDF$singleton <- ifelse(listlengthDF$L==1,1,0)

### Absences ###################################################################

#remove missing observations
listlengthDF <- subset(listlengthDF,!is.na(species))
siteInfo <- subset(listlengthDF,!duplicated(grid))
siteInfo <- arrange(siteInfo,grid)

### folds #############################################################

folds <- readRDS(paste(myfolder,"folds_occModel_bands.rds",sep="/"))
listlengthDF$fold <- folds$fold[match(listlengthDF$grid,folds$grid)]
table(listlengthDF$fold)

#select fold of this task
fold.id = as.integer(Sys.getenv("SGE_TASK_ID", "1"))
#fold.id = 1

#split intp test and train
listlengthDF_test <- subset(listlengthDF,fold == fold.id)
listlengthDF_train<- subset(listlengthDF,fold != fold.id)

### indices #####################################################################

#order data by site and year:
listlengthDF_test$siteIndex <- as.numeric(factor(listlengthDF_test$grid))
listlengthDF_test$yearIndex <- as.numeric(factor(listlengthDF_test$year))
listlengthDF_test <- arrange(listlengthDF_test,siteIndex,yearIndex)

listlengthDF_train$siteIndex <- as.numeric(factor(listlengthDF_train$grid))
listlengthDF_train$yearIndex <- as.numeric(factor(listlengthDF_train$year))
listlengthDF_train <- arrange(listlengthDF_train,siteIndex,yearIndex)

### site info #################################################################

siteInfo_test <- subset(listlengthDF_test,!duplicated(grid))
siteInfo_train <- subset(listlengthDF_train,!duplicated(grid))

### BUGS object ################################################################

#for BUGS

bugs.data <- list(nsite = length(unique(listlengthDF$siteIndex)),
                  nsite_test = length(unique(listlengthDF_test$siteIndex)),
                  nsite_train = length(unique(listlengthDF_train$siteIndex)),
                  nyear_test = length(unique(listlengthDF_test$yearIndex)),
                  nyear_train = length(unique(listlengthDF_train$yearIndex)),
                  nvisit_test = nrow(listlengthDF_test),
                  nvisit_train = nrow(listlengthDF_train),
                  site_test = listlengthDF_test$siteIndex,
                  site_train = listlengthDF_train$siteIndex,
                  year_test = listlengthDF_test$yearIndex,
                  year_train = listlengthDF_train$yearIndex,
                  Effort_test = listlengthDF_test$singleton,
                  Effort_train = listlengthDF_train$singleton,
                  y_test = listlengthDF_test$species,
                  y_train = listlengthDF_train$species,
                  #add an adm effect
                  adm_train = siteInfo_train$admN,
                  det.adm_train = listlengthDF_train$admN,
                  det.open_train = scale(listlengthDF_train$Open),
                  n.adm_train = length(unique(siteInfo_train$admN)),
                  adm2 = siteInfo$admN2,
                  det.adm2_train = listlengthDF_train$admN2,
                  n.adm2 = length(unique(siteInfo$admN2)))

#bugs.data_ArtsDaten <- bugs.data

### initials ####################################################################

#get JAGS libraries
library(rjags)
library(jagsUI)

#need to specify initial values
zst <- reshape2::acast(listlengthDF_train, 
                       siteIndex~yearIndex, value.var="species",fun=max,na.rm=T)
zst [is.infinite(zst)] <- 0
inits <- function(){list(z = zst)}

### scale vars #################################################################

siteInfo$admGrouped <- as.numeric(as.factor(siteInfo$admGrouped))
siteInfo[,-c(1:10)] <- plyr::numcolwise(scale)(siteInfo[,-c(1:10)])

siteInfo_test <- subset(siteInfo,grid %in% siteInfo_test$grid)
siteInfo_train <- subset(siteInfo,grid %in% siteInfo_train$grid)

### fit model ########################################################

#specify model structure
bugs.data$occDM_train <- model.matrix(~ scale(siteInfo_train$tree_line_position) + 
                                  scale(siteInfo_train$tree_line_position^2) +
                                  scale(siteInfo_train$y) +
                                  scale(siteInfo_train$bio1) +
                                  scale(siteInfo_train$bio1^2) +
                                  scale(siteInfo_train$Bog) +
                                  scale(siteInfo_train$Meadows) +
                                  scale(siteInfo_train$ODF) +
                                  scale(siteInfo_train$OSF) +
                                  scale(siteInfo_train$MountainBirchForest))[,-1]

bugs.data$occDM_test <- model.matrix(~ scale(siteInfo_test$tree_line_position) + 
                                        scale(siteInfo_test$tree_line_position^2) +
                                        scale(siteInfo_test$y) +
                                        scale(siteInfo_test$bio1) +
                                        scale(siteInfo_test$bio1^2) +
                                        scale(siteInfo_test$Bog) +
                                        scale(siteInfo_test$Meadows) +
                                        scale(siteInfo_test$ODF) +
                                        scale(siteInfo_test$OSF) +
                                        scale(siteInfo_test$MountainBirchForest))[,-1]

bugs.data$n.covs <- ncol(bugs.data$occDM_test)

params <- c("mean.p","beta","beta.effort")

modelfile <- paste(myfolder,"BUGS_occuModel_upscaling_CV.txt",sep="/")

n.cores = as.integer(Sys.getenv("NSLOTS", "1")) 
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

summary(out1$Rhat$beta)

#once converged update to get z, py and psi
out2 <- update(out1,
               parameters.to.save = c("mid.z_train","mid.psi_train","Py_train",
                                      "mid.z_test","mid.psi_test","Py_test"),
               n.iter = 1000)

### AUC check ############################################################

library(ggmcmc)
ggd2 <- ggs(out2$samples)

#Y against Py 

#test 
Py_preds <- subset(ggd2,grepl("Py_test",ggd2$Parameter))
Py_preds$Iteration <- as.numeric(factor(paste(Py_preds$Iteration,Py_preds$Chain)))
nu_Iteractions <- max(Py_preds$Iteration)
head(Py_preds)

#look through all iterations
AUC_py_test <- rep(NA, nu_Iteractions)

for (i in 1:nu_Iteractions){
  
  py.vals <- Py_preds$value[Py_preds$Iteration==i]
  
  pred <- ROCR::prediction(py.vals, 
                           bugs.data$y_test)
  
  #get AUC
  perf <- ROCR::performance(pred, "auc")
  AUC_py_test[i] <- perf@y.values[[1]]
  
}

summary(AUC_py_test)

#train
Py_preds <- subset(ggd2,grepl("Py_train",ggd2$Parameter))
Py_preds$Iteration <- as.numeric(factor(paste(Py_preds$Iteration,Py_preds$Chain)))
nu_Iteractions <- max(Py_preds$Iteration)
head(Py_preds)

#look through all iterations
AUC_py_train <- rep(NA, nu_Iteractions)

for (i in 1:nu_Iteractions){
  
  py.vals <- Py_preds$value[Py_preds$Iteration==i]
  
  pred <- ROCR::prediction(py.vals, 
                     bugs.data$y_train)
  
  #get AUC
  perf <- ROCR::performance(pred, "auc")
  AUC_py_train[i] <- perf@y.values[[1]]

}

summary(AUC_py_train)


#Z against psi

#test
Preds <- subset(ggd2,grepl("mid.psi_test",ggd2$Parameter))
Z_Preds <- subset(ggd2,grepl("mid.z_test",ggd2$Parameter))
Preds$Iteration <- as.numeric(factor(paste(Preds$Iteration,Preds$Chain)))
nu_Iteractions <- max(Preds$Iteration)
head(Preds)

#look through all iterations
AUC_psi_test <- rep(NA, nu_Iteractions)

for (i in 1:nu_Iteractions){
  
  psi.vals <- Preds$value[Preds$Iteration==i]
  z.vals <- Z_Preds$value[Preds$Iteration==i]
  
  pred <- ROCR::prediction(psi.vals, z.vals)
  
  #get AUC
  perf <- ROCR::performance(pred, "auc")
  AUC_psi_test[i] <- perf@y.values[[1]]
  
}

summary(AUC_psi_test)


#train
Preds <- subset(ggd2,grepl("mid.psi_train",ggd2$Parameter))
Z_Preds <- subset(ggd2,grepl("mid.z_train",ggd2$Parameter))
Preds$Iteration <- as.numeric(factor(paste(Preds$Iteration,Preds$Chain)))
nu_Iteractions <- max(Preds$Iteration)
head(Preds)

#look through all iterations
AUC_psi_train <- rep(NA, nu_Iteractions)

for (i in 1:nu_Iteractions){
  
  psi.vals <- Preds$value[Preds$Iteration==i]
  z.vals <- Z_Preds$value[Preds$Iteration==i]
  
  pred <- ROCR::prediction(psi.vals, z.vals)
  
  #get AUC
  perf <- ROCR::performance(pred, "auc")
  AUC_psi_train[i] <- perf@y.values[[1]]
  
}

summary(AUC_psi_train)


#combine all together and save
all_AUCs <- data.frame(AUC_psi_test,AUC_psi_train,AUC_py_test,AUC_py_train)
saveRDS(all_AUCs,file=paste0("all_AUCs_fold.id_",fold.id,".rds"))

### end ##################################################################
