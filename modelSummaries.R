library(raster)
library(sp)
library(maptools)
library(rgeos)

### check siteInfo #######################################################

siteInfo_Occ <- readRDS("data/siteInfo_ArtsDaten.rds")
siteInfo_Abund <- readRDS("data/siteInfo_LineTransects.rds")

### common grid ###########################################################

#using a m grid
equalM<-"+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

#get Norway
data(wrld_simpl)
Norway <- subset(wrld_simpl,NAME=="Norway")
NorwayOrig <- Norway
Norway <- gBuffer(Norway,width=1)
Norway <- spTransform(Norway,crs(equalM))
NorwayOrigProj <- spTransform(NorwayOrig,crs(equalM))

#create grid
newres = 5000#5 km grid
mygrid<-raster(extent(projectExtent(Norway,equalM)))
res(mygrid) <- newres
mygrid[] <- 1:ncell(mygrid)
plot(mygrid)
gridTemp <- mygrid
myGridDF <- as.data.frame(mygrid,xy=T)

### OCCU models ##############################################

siteInfo <- readRDS("data/siteInfo_ArtsDaten.rds")

### plot map #################################################

#from full model
out1 <- readRDS("model-outputs/out_occModel_upscaling.rds")
print(out1$summary,3)
siteInfo$preds <- out1$mean$grid.muZ
mygrid[] <- NA
mygrid[siteInfo$grid] <- siteInfo$preds
plot(mygrid)

#from model summary
out1 <- readRDS("model-outputs/outSummary_occModel_upscaling.rds")
out1 <- data.frame(out1)
out1$Param <- row.names(out1)
preds <- subset(out1,grepl("grid.z",out1$Param))
siteInfo$preds <- preds$mean
mygrid[] <- NA
mygrid[siteInfo$grid] <- siteInfo$preds
plot(mygrid)

# using tmap package
crs(mygrid) <- equalM
#tmaptools::palette_explorer()
library(tmap)
tm_shape(mygrid)+
  tm_raster(title="Occupancy prob",palette="YlGnBu")

### plot coefficients #######################################

betas <- subset(out1,grepl("beta",out1$Param))
betas <- betas[1:8,]
betas$variables <- c("tree_line_position", "tree_line_position2","bio1","bio1_2", "bio6","elevation","prefopen", "open")

ggplot(betas)+
  geom_crossbar(aes(x=variables,y=mean,
                    ymin=X2.5.,ymax=X97.5.))+
  coord_flip()+
  theme_bw()+
  ylab("effect size on occupancy")+
  geom_hline(yintercept=0,color="red",
             linetype="dashed")


### auc ######################################################

out2 <- readRDS("model-outputs/out_update_occModel_upscaling.rds")

library(ggmcmc)
ggd2 <- ggs(out2$samples)

#Z against psi
Preds <- subset(ggd2,grepl("mid.psi",ggd2$Parameter))
Z_Preds <- subset(ggd2,grepl("mid.z",ggd2$Parameter))
Preds$Iteration <- as.numeric(factor(paste(Preds$Iteration,Preds$Chain)))
nu_Iteractions <- max(Preds$Iteration)
head(Preds)

#look through all iterations
AUC_psi <- rep(NA, nu_Iteractions)

for (i in 1:nu_Iteractions){
  
  psi.vals <- Preds$value[Preds$Iteration==i]
  z.vals <- Z_Preds$value[Preds$Iteration==i]
  
  pred <- ROCR::prediction(psi.vals, z.vals)
  
  #get AUC
  perf <- ROCR::performance(pred, "auc")
  AUC_psi[i] <- perf@y.values[[1]]
  
}

summary(AUC_psi)

#Y against Py 

Py_preds <- subset(ggd2,grepl("Py",ggd2$Parameter))
Py_preds$Iteration <- as.numeric(factor(paste(Py_pred$Iteration,Py_preds$Chain)))
nu_Iteractions <- max(Py_preds$Iteration)
head(Py_preds)

#look through all iterations
AUC_py <- rep(NA, nu_Iteractions)

for (i in 1:nu_Iteractions){
  
  py.vals <- Py_preds$value[Py_preds$Iteration==i]
  
  pred <- ROCR::prediction(py.vals, 
                           bugs.data$y_test)
  
  #get AUC
  perf <- ROCR::performance(pred, "auc")
  AUC_py[i] <- perf@y.values[[1]]
  
}

summary(AUC_py)


### DISTANCE models ##########################################

siteInfo <- readRDS("data/siteInfo_LineTransects.rds")
environData <- readRDS("data/environData.rds")

### plot map ##############################################

out1 <- readRDS("model-outputs/outSummary_linetransectModel_variables.rds")
out1[row.names(out1)=="totalPop",]

out1 <- data.frame(out1)
out1$Param <- row.names(out1)
preds <- subset(out1,grepl("Density",out1$Param))
environData$preds <- preds$mean
environData$predsSD <- preds$sd

mygrid[] <- NA
mygrid[environData$grid] <- environData$preds
plot(mygrid)

# using tmap package
crs(mygrid) <- equalM
#tmaptools::palette_explorer()
library(tmap)
tm_shape(NorwayOrigProj)+tm_polygons(col="white")+
  tm_shape(mygrid)+ tm_raster(title="Density",palette="YlGnBu",n=10)+
  tm_layout(legend.position = c("left","top"))

#over range of data
preds <- subset(out1,grepl("Dens_lt",out1$Param))
environData$preds <- NA
environData$preds[environData$surveys==1] <- preds$mean
summary(environData$preds)
mygrid[] <- NA
mygrid[environData$grid] <- environData$preds

tm_shape(NorwayOrigProj)+tm_polygons(col="white")+
  tm_shape(mygrid)+ tm_raster(title="Density",palette="YlGnBu")+
  tm_layout(legend.position = c("left","top"))

#uncertainty
environData$Uncertain <- environData$preds/environData$predsSD 
mygrid[] <- NA
mygrid[environData$grid] <- environData$Uncertain
plot(mygrid)

# using tmap package
crs(mygrid) <- equalM
#tmaptools::palette_explorer()
library(tmap)
tm_shape(NorwayOrigProj)+tm_polygons(col="white")+
  tm_shape(mygrid)+ tm_raster(title = "Density uncertainty",
                              style="pretty")+
  tm_layout(legend.position = c("left","top"))

### plot coefficients #####################################

betas <- subset(out1,grepl("beta",out1$Param))
betas$Param <- c("bio1","open","tree_line_position", "tree_line_position2")

ggplot(betas)+
  geom_crossbar(aes(x=Param,y=mean,
                    ymin=X2.5.,ymax=X97.5.))+
  coord_flip()+
  theme_bw()+
  ylab("effect size on occupancy")+
  geom_hline(yintercept=0,color="red",
             linetype="dashed")



### predictive fits #################

#full data model

out1 <- readRDS("model-outputs/outSummary_linetransectModel_variables_v3.rds")
out1 <- data.frame(out1)
out1$Param <- row.names(out1)

#mean across all years
preds <- subset(out1,grepl("meanExpNu",out1$Param))
dataMeans <- apply(bugs.data$NuIndivs,1,mean,na.rm=T)
preds$data <- as.numeric(dataMeans)
#main plot
qplot(data,mean,data=preds)+
  geom_abline(intercept=0,slope=1)
#correlated but noisey
#on log-scale
qplot(data,mean,data=preds)+
  geom_abline(intercept=0,slope=1)+
  scale_x_log10()+scale_y_log10()#not bad, but not great
cor.test(log(preds$data),log(preds$mean))
#v1 - correlation is 0.358
#v2 - correlation is 0.367
#v3 - correlation is 0.374
#BPV
hist(out1$mean[out1$Param=="bpv"])
summary(out1$mean[out1$Param=="bpv"])

#CV fold models
library(ggmcmc)

#get all files for fold 1
fold <- 1
#folder <- "null"
folder <- "environ_only"

myFiles <- list.files(paste0("model-outputs/linetransectModel_CV/",folder))
myFolds <- myFiles[grepl(paste0("_",fold,".rds"),myFiles)]

#read in main model
out1 <- readRDS(paste("model-outputs/linetransectModel_CV",folder,myFolds[1],sep="/"))
ggd <- ggs(out1$samples)
#we have mid.expNuIndivs_test, 
#mid.expNuIndivs_train

#get training data

out1_train <- subset(ggd,grepl("mean.expNuIndivs_train",ggd$Parameter))
out1_train$siteIndex <- sub(".*\\[([^][]+)].*", "\\1", as.character(out1_train$Parameter))
out1_train$index <- as.numeric(interaction(out1_train$Iteration,out1_train$Chain))

#get actual NuIndiv
train_fold <- myFolds[grepl("train",myFolds)]
totalsInfo <- readRDS(paste("model-outputs/linetransectModel_CV/",folder,train_fold,sep="/"))
totalsInfo_mean <- rowMeans(totalsInfo,na.rm=T)

#check they are consistent:
length(unique(out1_train$siteIndex)) == dim(totalsInfo)[1]

#plot correlation between mean values
meanVals <- plyr::ddply(out1_train,"siteIndex",summarise,pred=median(value))
meanVals$obs <- totalsInfo_mean
qplot(obs,pred,data=meanVals)#bad:) 

#get difference between this value and the simulated values
mad_train <- as.numeric()
rmse_train <- as.numeric()
n.index <- max(out1_train$index)
for(i in 1:n.index){
  mad_train[i] <- mean(abs(totalsInfo_year5[!is.na(totalsInfo_year5)] - 
                       out1_train$value[out1_train$index==i][!is.na(totalsInfo_year5)]))
  rmse_train[i] <- sqrt(mean((totalsInfo_year5[!is.na(totalsInfo_year5)] - 
                        out1_train$value[out1_train$index==i][!is.na(totalsInfo_year5)])^2))
}

summary(mad)
hist(mad)
hist(rmse)

#compare the density estimates and the raw density estimates
out1_train <- subset(ggd,grepl("mean.Density_train",ggd$Parameter))
out1_train$siteIndex <- sub(".*\\[([^][]+)].*", "\\1", as.character(out1_train$Parameter))
out1_train$index <- as.numeric(interaction(out1_train$Iteration,out1_train$Chain))

#get observed densities
train_fold <- myFolds[grepl("train",myFolds)]
totalsInfo <- readRDS(paste("model-outputs/linetransectModel_CV/",folder,train_fold,sep="/"))
densityInfo <- totalsInfo/(transectLengths/1000 * 106/1000 *2)
densityInfo_mean <- rowMeans(densityInfo,na.rm=T)

#check they are consistent:
length(unique(out1_train$siteIndex)) == length(densityInfo_mean)

#plot correlation between mean values
meanVals <- plyr::ddply(out1_train,"siteIndex",summarise,pred=median(value))
meanVals$obs <- densityInfo_mean
qplot(obs,pred,data=meanVals) 


#test data
#get training data
out1_test <- subset(ggd,grepl("mid.expNuIndivs_test",ggd$Parameter))
out1_test$siteIndex <- sub(".*\\[([^][]+)].*", "\\1", as.character(out1_test$Parameter))
out1_test$index <- as.numeric(interaction(out1_test$Iteration,out1_test$Chain))

#get actual NuIndiv
test_fold <- myFolds[grepl("test",myFolds)]
totalsInfo <- readRDS(paste("model-outputs/linetransectModel_CV",test_fold,sep="/"))
totalsInfo_year5 <- totalsInfo[,5]

#check they are consistent:
length(unique(out1_test$siteIndex)) == dim(totalsInfo)[1]

#plot correlation between mean values
meanVals <- plyr::ddply(out1_test,"siteIndex",summarise,pred=median(value))
meanVals$obs <- totalsInfo_year5
qplot(obs,pred,data=meanVals)#bad:) 

#get difference between this value and the simulated values
mad_test <- as.numeric()
rmse_test <- as.numeric()
n.index <- max(out1_test$index)
for(i in 1:n.index){
  mad_test[i] <- mean(abs(totalsInfo_year5[!is.na(totalsInfo_year5)] - 
                       out1_test$value[out1_test$index==i][!is.na(totalsInfo_year5)]))
  rmse_test[i] <- sqrt(mean((totalsInfo_year5[!is.na(totalsInfo_year5)] - 
                          out1_test$value[out1_test$index==i][!is.na(totalsInfo_year5)])^2))
}

summary(mad_test)
hist(mad_test)
hist(rmse_test)

#### COMBINED model #####################################

out1 <- readRDS("model-outputs/outSummary_simpleCombinedModel.rds")
out1[row.names(out1)=="totalAbund",]

out1 <- data.frame(out1)
out1$Param <- row.names(out1)
preds <- subset(out1,grepl("realAbund",out1$Param))
environData$preds <- preds$mean
environData$predsSD <- preds$sd

mygrid[] <- NA
mygrid[environData$grid] <- environData$preds
plot(mygrid)

# using tmap package
crs(mygrid) <- equalM
#tmaptools::palette_explorer()
library(tmap)
tm_shape(NorwayOrigProj)+tm_polygons(col="white")+
  tm_shape(mygrid)+ tm_raster(title="Density",palette="YlGnBu",n=10)+
  tm_layout(legend.position = c("left","top"))

#over range of data
preds <- subset(out1,grepl("Dens_lt",out1$Param))
environData$preds <- NA
environData$preds[environData$surveys==1] <- preds$mean
summary(environData$preds)
mygrid[] <- NA
mygrid[environData$grid] <- environData$preds

tm_shape(NorwayOrigProj)+tm_polygons(col="white")+
  tm_shape(mygrid)+ tm_raster(title="Density",palette="YlGnBu")+
  tm_layout(legend.position = c("left","top"))

#uncertainty
environData$predsSD <- preds$mean/preds$sd
mygrid[] <- NA
mygrid[environData$grid] <- environData$predsSD
plot(mygrid)

# using tmap package
crs(mygrid) <- equalM
#tmaptools::palette_explorer()
library(tmap)
tm_shape(NorwayOrigProj)+tm_polygons(col="white")+
  tm_shape(mygrid)+ tm_raster(title = "Density uncertainty")+
  tm_layout(legend.position = c("left","top"))

### end #################################################