library(raster)
library(sp)
library(sf)
library(maptools)
library(rgeos)
library(tmap)
library(ggplot2)
#tmaptools::palette_explorer()

### check siteInfo #######################################################

siteInfo_Occ <- readRDS("data/siteInfo_ArtsDaten.rds")

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

### plot map #################################################

#or slurm models
out1 <- readRDS("model-outputs/SLURM/occModel/outSummary_occModel_upscaling_1.rds")
out1 <- readRDS("model-outputs/SLURM/occModel/outSummary_occModel_upscaling_2.rds")
out1 <- readRDS("model-outputs/SLURM/occModel/outSummary_occModel_upscaling_3.rds")

out1 <- data.frame(out1)
out1$Param <- row.names(out1)
table(out1$Rhat<1.1)

#Z
preds <- subset(out1,grepl("grid.z",out1$Param))
siteInfo_Occ$preds <- preds$mean
mygrid[] <- NA
mygrid[siteInfo_Occ$grid] <- siteInfo_Occ$preds
plot(mygrid)

# using tmap package
crs(mygrid) <- equalM
occ_tmap <- tm_shape(mygrid)+
  tm_raster(title="Pr(Occupancy)",palette="YlGnBu")
occ_tmap

#psi
preds <- subset(out1,grepl("grid.psi",out1$Param))
siteInfo_Occ$preds <- preds$mean
mygrid[] <- NA
mygrid[siteInfo_Occ$grid] <- siteInfo_Occ$preds
plot(mygrid)

# using tmap package
crs(mygrid) <- equalM
tm_shape(mygrid)+
  tm_raster(title="Occupancy prob",palette="YlGnBu")

### model selection ##########################################

#model 1 - a priori selection
out1 <- readRDS("model-outputs/SLURM/occModel/outSummary_occModel_upscaling_1.rds")
out1 <- data.frame(out1)
out1$Param <- row.names(out1)
betas <- subset(out1,grepl("beta",out1$Param))
betas <- subset(betas,Param!="beta.det.open")
betas <- subset(betas,Param!="beta.effort")
betas$variables <- c("tree_line_position","tree_line_position_2","geo_y","bio1","bio1_2","Bog","Meadows","ODF","OSF","MountainBirchForest")

ggplot(betas)+
  geom_crossbar(aes(x=variables,y=mean, ymin=X2.5.,ymax=X97.5.))+
  coord_flip()+
  theme_bw()+
  ylab("effect size on occupancy")+
  geom_hline(yintercept=0,color="red",
             linetype="dashed")


#model 2 = lasso
out1 <- readRDS("model-outputs/SLURM/occModel/outSummary_occModel_upscaling_2.rds")
out1 <- data.frame(out1)
out1$Param <- row.names(out1)
betas <- subset(out1,grepl("beta",out1$Param))
betas <- subset(betas,Param!="beta.det.open")
betas <- subset(betas,Param!="beta.effort")
betas$variables <- c("bio6","bio5","tree_line_position","MountainBirchForest","Bog","ODF","Meadows",
                               "OSF","Mire","SnowBeds","y","distCoast",
                               "bio6_2","bio5_2","tree_line_position_2",
                               "MountainBirchForest_2","Bog_2","ODF_2","Meadows_2","OSF_2","Mire_2",
                               "SnowBeds_2","y_2","distCoast_2")
ggplot(betas)+
  geom_crossbar(aes(x=variables,y=mean, ymin=X2.5.,ymax=X97.5.))+
  coord_flip()+
  theme_bw()+
  ylab("effect size on occupancy")+
  geom_hline(yintercept=0,color="red",
             linetype="dashed")

#model 3 = variable indicator
out1 <- readRDS("model-outputs/SLURM/occModel/outSummary_occModel_upscaling_3.rds")
out1 <- data.frame(out1)
out1$Param <- row.names(out1)
betas <- subset(out1,grepl("beta",out1$Param))
betas <- subset(betas,Param!="beta.det.open")
betas <- subset(betas,Param!="beta.effort")
betas$variables <- c("bio6","bio5","tree_line_position","MountainBirchForest","Bog","ODF","Meadows",
                     "OSF","Mire","SnowBeds","y","distCoast",
                     "bio6_2","bio5_2","tree_line_position_2",
                     "MountainBirchForest_2","Bog_2","ODF_2","Meadows_2","OSF_2","Mire_2",
                     "SnowBeds_2","y_2","distCoast_2")
ggplot(betas)+
  geom_crossbar(aes(x=variables,y=mean, ymin=X2.5.,ymax=X97.5.))+
  coord_flip()+
  theme_bw()+
  ylab("effect size on occupancy")+
  geom_hline(yintercept=0,color="red",
             linetype="dashed")

### AUC ######################################################

out2 <- readRDS("model-outputs/out_update_occModel_upscaling_v1.rds")

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
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.8677  0.8767  0.8800  0.8799  0.8828  0.8938

#Y against Py 
Py_preds <- subset(ggd2,grepl("Py",ggd2$Parameter))
Py_preds$Iteration <- as.numeric(factor(paste(Py_preds$Iteration,Py_preds$Chain)))
nu_Iteractions <- max(Py_preds$Iteration)
head(Py_preds)

#look through all iterations
AUC_py <- rep(NA, nu_Iteractions)

for (i in 1:nu_Iteractions){
  
  py.vals <- Py_preds$value[Py_preds$Iteration==i]
  
  pred <- ROCR::prediction(py.vals, 
                           bugs.data$y)
  
  #get AUC
  perf <- ROCR::performance(pred, "auc")
  AUC_py[i] <- perf@y.values[[1]]
  
}

summary(AUC_py)
#doesnt run

### cross validation #######################################

myfolds <- list.files("model-outputs/occModel_CV") %>%
            str_subset("fold")

temp <- plyr::ldply(1:5,function(x){
  temp <- readRDS(paste0("model-outputs/occModel_CV/",myfolds[x]))
  temp$fold.id <- x
  return(temp)
})
  
temp %>% 
  group_by(fold.id) %>%
  summarise(across(everything(),mean))
# fold.id AUC_psi_test AUC_psi_train AUC_py_test AUC_py_train
# 1       1        0.932         0.967       0.614        0.960
# 2       2        0.970         0.962       0.799        0.956
# 3       3        0.950         0.966       0.590        0.975
# 4       4        0.974         0.960       0.894        0.947
# 5       5        0.946         0.965       0.742        0.952
                
temp %>% 
  group_by(fold.id) %>%
  summarise(across(everything(),mean)) %>%
  colMeans()

### DISTANCE models ##########################################

bufferData <- readRDS("data/varDF_allEnvironData_buffers_idiv.rds")

#slurm model

#negative binomial models
out1 <- readRDS("model-outputs/SLURM/distanceModel/outSummary_linetransectModel_variables_1.rds")

out1 <- readRDS("model-outputs/SLURM/distanceModel/outSummary_linetransectModel_variables_2.rds")

out1 <- readRDS("model-outputs/SLURM/distanceModel/outSummary_linetransectModel_variables_3.rds")


#poisson models
out1 <- readRDS("model-outputs/SLURM/distanceModel/poisson/outSummary_linetransectModel_variables_1.rds")

out1 <- readRDS("model-outputs/SLURM/distanceModel/poisson/outSummary_linetransectModel_variables_2.rds")

out1 <- readRDS("model-outputs/SLURM/distanceModel/poisson/outSummary_linetransectModel_variables_3.rds")

### plot map ##############################################

#density over range of line transects
out1 <- data.frame(out1)
out1$Param <- row.names(out1)
preds <- subset(out1,grepl("meanDensity",out1$Param))
bufferData$preds <- preds$mean
bufferData$predsSD <- preds$sd
summary(bufferData$preds)

#model 1
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#5.118  10.755  12.963  12.955  15.566  24.552
#model 2
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.882  10.927  12.688  12.920  14.880  24.797
#model 3
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#3.58   11.03   13.02   12.95   15.07   21.90

#tmap
bufferData_st <- st_as_sf(bufferData,coords=c("x","y"),crs=equalM)
density_tmap <- tm_shape(NorwayOrigProj)+
  tm_borders()+
tm_shape(bufferData_st)+
  tm_dots("preds",title="Est. Density",palette="YlGnBu",size=0.05)+
  tm_layout(legend.position=c("left","top"))

#plot side by side with occupancy
tmap_arrange(occ_tmap,density_tmap,nrow=1)

#across whole range
preds <- subset(out1,grepl("Density",out1$Param))
preds <- subset(preds,!grepl("meanDensity",preds$Param))
siteInfo_Occ$density_preds <- preds$mean
siteInfo_Occ$density_predsSD <- preds$sd

#tmap
mygrid[] <- NA
mygrid[siteInfo_Occ$grid] <- siteInfo_Occ$density_preds
tm_shape(mygrid)+
  tm_raster(title="Pred. Density",palette="YlGnBu")
#weird!!

#is it just bad prediction outside of covariate range??

#relationship between occupancy and abundance?
ggplot(siteInfo_Occ)+
  geom_point(aes(x=preds,
                 y=density_preds))

#what does realised density look like?
siteInfo_Occ$real_density <- siteInfo_Occ$preds * siteInfo_Occ$density_preds
summary(siteInfo_Occ$real_density)
mygrid[] <- NA
mygrid[siteInfo_Occ$grid] <- siteInfo_Occ$real_density
tm_shape(mygrid)+
  tm_raster(title="Real. Density",palette="YlGnBu")

#is it just outliers?
summary(siteInfo_Occ$real_density)
siteInfo_Occ$real_density[siteInfo_Occ$real_density>300] <- 300
mygrid[] <- NA
mygrid[siteInfo_Occ$grid] <- siteInfo_Occ$real_density
tm_shape(mygrid)+
  tm_raster(title="Real. Density",palette="YlGnBu")
#yes! it was just weird outliers!!
sum(siteInfo_Occ$real_density)

### predictive fits #################

#mean across all years
out1 <- data.frame(out1)
out1$Param <- row.names(out1)

preds <- subset(out1,grepl("meanExpNu",out1$Param))
dataMeans <- apply(bugs.data$NuIndivs,1,mean,na.rm=T)
preds$data <- as.numeric(dataMeans)
summary(preds$data)
summary(preds$mean)

#main plot
qplot(data,mean,data=preds)+
  geom_abline(intercept=0,slope=1)
#correlated but noisey
#on log-scale
qplot(data,mean,data=preds)+
  geom_abline(intercept=0,slope=1)+
  scale_x_log10()+scale_y_log10()
cor.test(log(preds$data),log(preds$mean))
#v1 - correlation is 0.358
#v2 - correlation is 0.367
#v3 - correlation is 0.374
#v4 - correlation is 0.506
#v5 - correlation is 0.524
#v6 - correlation is 0.546
#model 1 - 0.55
#model 2 - 0.60
#model 3 - 0.58

#poisson model
#model 1 - 0.7159
#model 2 - 0.7232
#model 3 - 0.7211

### BPV ###################################################

#compare expNuIndivs and expNuIndivs.new
exp <- subset(out1,grepl("expNuIndivs",out1$Param))
expNu <- subset(out1,grepl("expNuIndivs.new",out1$Param))

hist(exp$mean)
hist(expNu$mean)
#look very similar!!!
hist(exp$mean-expNu$mean)
summary(exp$mean-expNu$mean)#centered on zero...

#full models
#negative binomial
out1 <- readRDS("model-outputs/SLURM/distanceModel/out_linetransectModel_variables_1.rds")
out1 <- readRDS("model-outputs/SLURM/distanceModel/out_linetransectModel_variables_2.rds")
out1 <- readRDS("model-outputs/SLURM/distanceModel/out_linetransectModel_variables_3.rds")

#poisson
out1 <- readRDS("model-outputs/SLURM/distanceModel/poisson/out_linetransectModel_variables_1.rds")
out1 <- readRDS("model-outputs/SLURM/distanceModel/poisson/out_linetransectModel_variables_2.rds")
out1 <- readRDS("model-outputs/SLURM/distanceModel/poisson/out_linetransectModel_variables_3.rds")

#BPV
mean(out1$sims.list$fit.new > out1$sims.list$fit)

hist(out1$sims.list$fit.new)
summary(out1$sims.list$fit.new)
median(out1$sims$fit)
abline(v=median(out1$sims.list$fit),col="red")

### Dharma ###############################################

library(DHARMa)
#get simulated data
simulations = out1$sims.list$expNuIndivs.new

#change into a 2-D matrix
#https://stackoverflow.com/questions/37662433/r-3d-array-to-2d-matrix
#dim(simulations)
dim(simulations) <- c(dim(simulations)[1],599 * 11)
simsMean <- simulations
dim(simsMean) <- c(dim(simsMean)[1],599*11)
simsMean = apply(simsMean, 2, median)

#get model predictions
preds = out1$sims.list$expNuIndivs
dim(preds)
dim(preds) <- c(dim(simulations)[1],599*11)
preds = apply(preds, 2, median)
#length(preds)

#get observed data
obs <- bugs.data$NuIndivs
dim(obs)
dim(obs) <- c(599*11)
#length(obs)

#now need to remove missing values
notMiss <- !is.na(obs)
length(notMiss)

sim = createDHARMa(simulatedResponse = t(simulations[,notMiss]),
                   observedResponse = obs[notMiss],
                   fittedPredictedResponse = preds[notMiss],
                   integerResponse = T)

plot(sim)


#poisson model - residuals less good!!

### model selection #######################################

#model 1
out1 <- readRDS("model-outputs/SLURM/distanceModel/outSummary_linetransectModel_variables_1.rds")
out1 <- data.frame(out1)
out1$Param <- row.names(out1)
betas <- subset(out1,grepl("beta",out1$Param))
betas$variables <- c("y",'bio6',"bio6_2","distCoast","distCoast_2",
                     "bio5","bio5_2","tree_line","tree_line_2","OSF","SnowBeds")

ggplot(betas)+
  geom_crossbar(aes(x=variables,y=mean, ymin=X2.5.,ymax=X97.5.))+
  coord_flip()+
  theme_bw()+
  ylab("effect size on abundance")+
  geom_hline(yintercept=0,color="red",
             linetype="dashed")


#model 2 = lasso
out1 <- readRDS("model-outputs/SLURM/distanceModel/outSummary_linetransectModel_variables_2.rds")
out1 <- data.frame(out1)
out1$Param <- row.names(out1)
betas <- subset(out1,grepl("beta",out1$Param))
betas$variables <- c("bio6","bio5","y","distCoast","tree_line","MountainBirchForest",
                     "Bog","ODF","Meadows","OSF","Mire","SnowBeds",
                     "bio6_2","bio5_2","y_2","distCoast_2","tree_line_2","MountainBirchForest_2",
                     "Bog_2","ODF_2","Meadows_2","OSF_2","Mire_2","SnowBeds_2")
ggplot(betas)+
  geom_crossbar(aes(x=variables,y=mean, ymin=X2.5.,ymax=X97.5.))+
  coord_flip()+
  theme_bw()+
  ylab("effect size on abundance")+
  geom_hline(yintercept=0,color="red",
             linetype="dashed")

#model 3 = variable indicator
out1 <- readRDS("model-outputs/SLURM/distanceModel/outSummary_linetransectModel_variables_3.rds")
out1 <- data.frame(out1)
out1$Param <- row.names(out1)
betas <- subset(out1,grepl("beta",out1$Param))
betas$variables <- c("bio6","bio5","y","distCoast","tree_line","MountainBirchForest",
                     "Bog","ODF","Meadows","OSF","Mire","SnowBeds",
                     "bio6_2","bio5_2","y_2","distCoast_2","tree_line_2","MountainBirchForest_2",
                     "Bog_2","ODF_2","Meadows_2","OSF_2","Mire_2","SnowBeds_2")
ggplot(betas)+
  geom_crossbar(aes(x=variables,y=mean, ymin=X2.5.,ymax=X97.5.))+
  coord_flip()+
  theme_bw()+
  ylab("effect size on abundance")+
  geom_hline(yintercept=0,color="red",
             linetype="dashed")

#### cross validation ######################################

#CV fold models

library(ggmcmc)
#get all files for fold 1
fold <- 4
#folder <- "null"
folder <- "environ_v2"

myFiles <- list.files(paste0("model-outputs/linetransectModel_CV/",folder))
myFolds <- myFiles[grepl(paste0("_",fold,".rds"),myFiles)]

#read in main model
out1 <- readRDS(paste("model-outputs/linetransectModel_CV",folder,myFolds[1],sep="/"))
ggd <- ggs(out1$samples)

#extraxt training data
out1_train <- subset(ggd,grepl("mean.expNuIndivs_train",ggd$Parameter))
out1_train$siteIndex <- sub(".*\\[([^][]+)].*", "\\1", as.character(out1_train$Parameter))
out1_train$siteIndex <- as.numeric(as.character(out1_train$siteIndex))
out1_train$index <- as.numeric(interaction(out1_train$Iteration,out1_train$Chain))

#get actual NuIndiv
train_fold <- myFolds[grepl("train",myFolds)]
totalsInfo <- readRDS(paste("model-outputs/linetransectModel_CV/",folder,train_fold,sep="/"))
totalsInfo_mean <- rowMeans(totalsInfo,na.rm=T)

#check they are consistent:
length(unique(out1_train$siteIndex)) == dim(totalsInfo)[1]

#plot correlation between mean values
meanVals <- plyr::ddply(out1_train,"siteIndex",summarise,pred=mean(value))
meanVals$obs <- totalsInfo_mean
qplot(log(obs),log(pred),data=meanVals)
cor.test(log(meanVals$obs),log(meanVals$pred))#0.33

#get difference between this value and the simulated values
mad_train <- as.numeric()
rmse_train <- as.numeric()
n.index <- max(out1_train$index)

for(i in 1:n.index){
  mad_train[i] <- mean(abs(totalsInfo_mean[!is.na(totalsInfo_mean)] - 
                       out1_train$value[out1_train$index==i][!is.na(totalsInfo_mean)]))
  
  rmse_train[i] <- sqrt(mean((totalsInfo_mean[!is.na(totalsInfo_mean)] - 
                        out1_train$value[out1_train$index==i][!is.na(totalsInfo_mean)])^2))
  
}

summary(mad_train)
hist(mad_train)
hist(rmse_train)
summary(totalsInfo_mean)


#test data
out1_test <- subset(ggd,grepl("mean.expNuIndivs_test",ggd$Parameter))
out1_test$siteIndex <- sub(".*\\[([^][]+)].*", "\\1", as.character(out1_test$Parameter))
out1_test$siteIndex <- as.numeric(as.character(out1_test$siteIndex))
out1_test$index <- as.numeric(interaction(out1_test$Iteration,out1_test$Chain))

#get actual NuIndiv
test_fold <- myFolds[grepl("test",myFolds)]
totalsInfo <- readRDS(paste("model-outputs/linetransectModel_CV",folder,test_fold,sep="/"))
totalsInfo_mean <- rowMeans(totalsInfo,na.rm=T)

#check they are consistent:
length(unique(out1_test$siteIndex)) == dim(totalsInfo)[1]

#plot correlation between mean values
meanVals <- plyr::ddply(out1_test,"siteIndex",summarise,pred=median(value))
meanVals$obs <- totalsInfo_mean
qplot(log(obs),log(pred),data=meanVals)#bad:) 
cor.test(log(meanVals$obs),log(meanVals$pred))#0.42!! 

#get difference between this value and the simulated values
mad_test <- as.numeric()
rmse_test <- as.numeric()
n.index <- max(out1_test$index)
for(i in 1:n.index){
  mad_test[i] <- mean(abs(totalsInfo_mean[!is.na(totalsInfo_mean)] - 
                       out1_test$value[out1_test$index==i][!is.na(totalsInfo_mean)]))
  rmse_test[i] <- sqrt(mean((totalsInfo_mean[!is.na(totalsInfo_mean)] - 
                          out1_test$value[out1_test$index==i][!is.na(totalsInfo_mean)])^2))
}

summary(mad_test)
hist(mad_test)
hist(rmse_test)
summary(totalsInfo_mean)

#### COMBINED model #####################################

#from combined model:

out1 <- readRDS("model-outputs/outSummary_simpleCombinedModel.rds")
out1[row.names(out1)=="totalAbund",]

out1 <- data.frame(out1)
out1$Param <- row.names(out1)
preds <- subset(out1,grepl("realAbund",out1$Param))

#or get data from "HPC_simple_combined_analysis:

out2 <- data.frame(out1$summary)
out2$Param <- row.names(out2)
preds <- subset(out2,grepl("realAbund",out2$Param))
siteInfo_Occ$preds <- preds$mean
summary(siteInfo_Occ$preds)#minus numbers??? 

#bound to reasonable
siteInfo_Occ$preds[siteInfo_Occ$preds<0] <- 0
siteInfo_Occ$preds[siteInfo_Occ$preds>400] <- 400

mygrid[] <- NA
mygrid[siteInfo_Occ$grid] <- siteInfo_Occ$preds
crs(mygrid) <- equalM
total_tmap <- tm_shape(mygrid)+
  tm_raster(title="Abundance",palette="YlGnBu")
total_tmap

#plot uncertainty too
siteInfo_Occ$preds <- preds$sd
mygrid[] <- NA
mygrid[siteInfo_Occ$grid] <- siteInfo_Occ$preds
crs(mygrid) <- equalM
sd_tmap <- tm_shape(mygrid)+
  tm_raster(title="Abundance uncertainty",palette="YlGnBu",n=6)
sd_tmap

#relationship between mean and uncertainity
qplot(mean,sd^2,data=preds)+stat_smooth(method="lm")
#positive....

#get residuals of the relationship
preds$var <- preds$sd^2
preds$resid_sd <- preds$var/preds$mean
preds$resid_sd[preds$resid_sd<0] <- 0 #where we have more variation than expected given the mean
summary(preds$resid_sd)
preds$resid_sd[preds$resid_sd>60] <- 60

#plot it
siteInfo_Occ$preds <- preds$resid_sd
mygrid[] <- NA
mygrid[siteInfo_Occ$grid] <- siteInfo_Occ$preds
resid_tmap <- tm_shape(mygrid)+
  tm_raster(title="Overdispersion",palette="YlGnBu")
resid_tmap

#side by side
tmap_arrange(sd_tmap,resid_tmap,nrow=1)

#multiple them???
siteInfo_Occ$preds <- sqrt(preds$resid_sd * preds$sd)
siteInfo_Occ$preds[siteInfo_Occ$preds<=0] <- 0.001
summary(siteInfo_Occ$preds)
mygrid[] <- NA
mygrid[siteInfo_Occ$grid] <- siteInfo_Occ$preds
combined <- tm_shape(mygrid)+
  tm_raster(title="Combined",palette="YlGnBu",style="cont")

tmap_arrange(sd_tmap,resid_tmap,combined, nrow=1)
#model effect of number of CS data points of each type? not here, elsewhere