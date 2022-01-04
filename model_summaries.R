library(raster)
library(sp)
library(sf)
library(maptools)
library(rgeos)
library(tmap)
library(ggplot2)
library(ggmcmc)
library(tidyverse)
library(ggthemes)

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

#with detection covariates
out1 <- readRDS("model-outputs/SLURM/occModel/detection_covariates2/outSummary_occModel_upscaling_1.rds")
out1 <- readRDS("model-outputs/SLURM/occModel/detection_covariates2/outSummary_occModel_upscaling_2.rds")
out1 <- readRDS("model-outputs/SLURM/occModel/detection_covariates2/outSummary_occModel_upscaling_3.rds")


#null models
out1 <- readRDS("model-outputs/SLURM/occModel/null_models/outSummary_occModel_upscaling_1.rds")
out1 <- readRDS("model-outputs/SLURM/occModel/null_models/outSummary_occModel_upscaling_2.rds")
out1 <- readRDS("model-outputs/SLURM/occModel/null_models/outSummary_occModel_upscaling_3.rds")
out1 <- readRDS("model-outputs/SLURM/occModel/null_models/outSummary_occModel_upscaling_4.rds")
out1 <- readRDS("model-outputs/SLURM/occModel/null_models/outSummary_occModel_upscaling_5.rds")
out1 <- readRDS("model-outputs/SLURM/occModel/null_models/outSummary_occModel_upscaling_6.rds")


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
  tm_raster(title="Occupancy",palette="YlGnBu", style="cont")
occ_tmap

#and sd
siteInfo_Occ$preds <- preds$sd
mygrid[] <- NA
mygrid[siteInfo_Occ$grid] <- siteInfo_Occ$preds
plot(mygrid)

# using tmap package
crs(mygrid) <- equalM
occ_tmap_sd <- tm_shape(mygrid)+
  tm_raster(title="SD",palette="YlGnBu", style="cont")
occ_tmap_sd

tmap_arrange(occ_tmap,occ_tmap_sd,nrow=1)


#psi
preds <- subset(out1,grepl("grid.psi",out1$Param))
siteInfo_Occ$preds <- preds$mean
mygrid[] <- NA
mygrid[siteInfo_Occ$grid] <- siteInfo_Occ$preds
plot(mygrid)

# using tmap package
crs(mygrid) <- equalM
occ_tmap <- tm_shape(mygrid)+
  tm_raster(title="a) Occupancy prob",palette="YlGnBu", style="cont")
occ_tmap

siteInfo_Occ$preds <- preds$sd
mygrid[] <- NA
mygrid[siteInfo_Occ$grid] <- siteInfo_Occ$preds
plot(mygrid)

# using tmap package
crs(mygrid) <- equalM
occ_tmap_sd <- tm_shape(mygrid)+
  tm_raster(title="b) Occupancy SD",palette="YlGnBu", style="cont")
occ_tmap_sd

temp <- tmap_arrange(occ_tmap,occ_tmap_sd,nrow=1)
temp

tmap_save(temp, "plots/Fig_2.png",width = 6, height = 4)

### summary stats ###########################################

subset(out1,grepl("average",out1$Param))
subset(out1,grepl("propOcc",out1$Param))
subset(out1,grepl("beta.",out1$Param))

### model selection ##########################################

#model 1 - a priori selection
betas <- subset(out1,grepl("beta",out1$Param))
betas <- subset(betas,!grepl("beta.det",betas$Param))
betas <- subset(betas,!grepl("beta.effort",betas$Param))
betas$variables <- c("tree_line_position","tree_line_position_2","y","bio6",
                     "Bog","Mire","Meadows","ODF","ODF2",
                     "OSF","OSF2","MountainBirchForest")

ggplot(betas)+
  geom_crossbar(aes(x=variables,y=mean, ymin=X2.5.,ymax=X97.5.))+
  coord_flip()+
  theme_bw()+
  ylab("effect size on occupancy")+
  geom_hline(yintercept=0,color="red",
             linetype="dashed")


#model 2 = lasso
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
betas <- subset(out1,grepl("beta",out1$Param))
betas <- subset(betas,!grepl("effort",betas$Param))
betas <- subset(betas,!grepl("det",betas$Param))
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
             linetype="dashed")+
  theme_few()


#inclusion 
betas <- subset(out1,grepl("g",out1$Param))
betas <- subset(betas,!grepl("grid",betas$Param))
betas <- subset(betas,!grepl("average",betas$Param))
betas$variables <- c("bio6","bio5","tree_line_position","MountainBirchForest","Bog","ODF","Meadows",
                     "OSF","Mire","SnowBeds","y","distCoast",
                     "bio6_2","bio5_2","tree_line_position_2",
                     "MountainBirchForest_2","Bog_2","ODF_2","Meadows_2","OSF_2","Mire_2",
                     "SnowBeds_2","y_2","distCoast_2")
ggplot(betas)+
  geom_crossbar(aes(x=variables,y=mean, ymin=X2.5.,ymax=X97.5.))+
  coord_flip()+
  theme_bw()+
  ylab("Model inclusion")+
  geom_hline(yintercept=0,color="red",
             linetype="dashed")

subset(betas,mean>0.25)

### BPV #####################################################

#or slurm models
out1 <- readRDS("model-outputs/SLURM/occModel/out_update_occModel_upscaling_1.rds")
out1 <- readRDS("model-outputs/SLURM/occModel/out_update_occModel_upscaling_2.rds")
out1 <- readRDS("model-outputs/SLURM/occModel/out_update_occModel_upscaling_3.rds")

#with detection covariates
out1 <- readRDS("model-outputs/SLURM/occModel/detection_covariates2/out_update_occModel_upscaling_1.rds")
out1 <- readRDS("model-outputs/SLURM/occModel/detection_covariates2/out_update_occModel_upscaling_2.rds")
out1 <- readRDS("model-outputs/SLURM/occModel/detection_covariates2/out_update_occModel_upscaling_3.rds")

#with null models
out1 <- readRDS("model-outputs/SLURM/occModel/null_models/out_update_occModel_upscaling_1.rds")
out1 <- readRDS("model-outputs/SLURM/occModel/null_models/out_update_occModel_upscaling_2.rds")
out1 <- readRDS("model-outputs/SLURM/occModel/null_models/out_update_occModel_upscaling_3.rds")
out1 <- readRDS("model-outputs/SLURM/occModel/null_models/out_update_occModel_upscaling_4.rds")
out1 <- readRDS("model-outputs/SLURM/occModel/null_models/out_update_occModel_upscaling_5.rds")
out1 <- readRDS("model-outputs/SLURM/occModel/null_models/out_update_occModel_upscaling_6.rds")


#fit vs fit new
hist(out1$sims.list$fit.new)
summary(out1$sims.list$fit.new)
median(out1$sims$fit)
abline(v=median(out1$sims.list$fit),col="red")
mean(out1$sims.list$fit.new > out1$sims.list$fit)
#good!!

#with new models in which BPV is calculated within the model
subset(out1,grepl("bpv",out1$Param))

### AUC ######################################################

#now included direct in the HPC code

library(ggmcmc)
ggd2 <- ggs(out1$samples)

#Z against psi
Preds <- subset(ggd2,grepl("mid.psi",ggd2$Parameter))
Z_Preds <- subset(ggd2,grepl("mid.z",ggd2$Parameter))
Preds$Iteration <- as.numeric(factor(paste(Preds$Iteration,Preds$Chain)))
nu_Iteractions <- max(Preds$Iteration)
head(Preds)

#look through all iterations
AUC_psi <- rep(NA, nu_Iteractions)
TSS_psi <- rep(NA, nu_Iteractions)

for (i in 1:nu_Iteractions){

  psi.vals <- Preds$value[Preds$Iteration==i]
  z.vals <- Z_Preds$value[Preds$Iteration==i]

  pred <- ROCR::prediction(psi.vals, z.vals)

  #get AUC
  perf <- ROCR::performance(pred, "auc")
  AUC_psi[i] <- perf@y.values[[1]]
  
  #get TSS
  tp_perf <- ROCR::performance(pred, "tpr")
  tn_perf <- ROCR::performance(pred, "tnr")
  TSS_psi[i] <- tp_perf@y.values[[1]] + tn_perf@y.values[[1]] - 1
  
}

summary(AUC_psi)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.8677  0.8767  0.8800  0.8799  0.8828  0.8938

#Y against Py
Py_preds <- subset(ggd2,grepl("Py",ggd2$Parameter))
Py_preds$Iteration <- as.numeric(paste(factor(Py_preds$Iteration,Py_preds$Chain)))
nu_Iteractions <- max(Py_preds$Iteration)
head(Py_preds)

#look through all iterations
AUC_py <- rep(NA, nu_Iteractions)

for (i in 1:nu_Iteractions){

  py.vals <- Py_preds$value[Py_preds$Iteration==i]

  #remove the NAs
  missing <- is.na(bugs.data$y)

  pred <- ROCR::prediction(py.vals[!missing],bugs.data$y[!missing])

  #get AUC
  perf <- ROCR::performance(pred, "auc")
  AUC_py[i] <- perf@y.values[[1]]

}

summary(AUC_py)


### model AUC ##############################################

#in new models, AUC is output as part of model running on HPC

#full model - 3
readRDS("model-outputs/SLURM/occModel/null_models/AUC_psi_occModel_upscaling_3.rds")
readRDS("model-outputs/SLURM/occModel/null_models/AUC_py_occModel_upscaling_3.rds")

#null model - 4
readRDS("model-outputs/SLURM/occModel/null_models/AUC_psi_occModel_upscaling_4.rds")
readRDS("model-outputs/SLURM/occModel/null_models/AUC_py_occModel_upscaling_4.rds")

#null detection - 5
readRDS("model-outputs/SLURM/occModel/null_models/AUC_psi_occModel_upscaling_5.rds")
readRDS("model-outputs/SLURM/occModel/null_models/AUC_py_occModel_upscaling_5.rds")


#null state - 6
readRDS("model-outputs/SLURM/occModel/null_models/AUC_psi_occModel_upscaling_6.rds")
readRDS("model-outputs/SLURM/occModel/null_models/AUC_py_occModel_upscaling_6.rds")

### cross validation #######################################

myfolds <- list.files("model-outputs/SLURM/occModel/CV") %>%
            str_subset("BUGS_occuModel_upscaling_CV.txt")#model 1

myfolds <- list.files("model-outputs/SLURM/occModel/CV") %>%
  str_subset("BUGS_occuModel_upscaling_ModelSelection_CV")#model 3

temp <- plyr::ldply(1:5,function(x){
  temp <- readRDS(paste0("model-outputs/SLURM/occModel/CV/",myfolds[x]))
  temp$fold.id <- x
  return(temp)
})
  
model_3 <- temp %>% 
  group_by(fold.id) %>%
  summarise(across(everything(),mean))
# A tibble: 5 x 5
#fold.id AUC_psi_test AUC_psi_train AUC_py_test AUC_py_train
#  1       1        0.864         0.861       0.736        0.964
#2       2        0.835         0.869       0.658        0.964
#3       3        0.871         0.868       0.725        0.965
#4       4        0.868         0.849       0.743        0.962
#5       5        0.878         0.849       0.772        0.959

#satisfactory....

temp %>% 
  group_by(fold.id) %>%
  summarise(across(everything(),mean)) %>%
  colMeans()


# #plot all
# model_1$model <- "Gaussian priors" 
# model_2$model <- "LASSO priors" 
# model_3$model <- "Variable indicator" 
# allCV <- rbind(model_1,
#                model_2,
#                model_3)
# 
# allCV <- allCV %>%
#   pivot_longer(!c("fold.id","model"),
#                names_to = "variable", values_to="value")
# allCV$dataset <- sapply(allCV$variable,function(x)
#   strsplit(x,"_")[[1]][3])
# allCV$parameter <- sapply(allCV$variable,function(x)
#   strsplit(x,"_")[[1]][2])
# allCV$Parameter <- ifelse(allCV$parameter=="psi","Occupancy","Detection")
# 
# allCV$model <- factor(allCV$model, levels = c("LASSO priors", "Variable indicator","Gaussian priors"))
# 
# ggplot(allCV) +
#   geom_point(aes(x=fold.id,y=value, colour=model))+
#   facet_grid(Parameter~dataset)+
#   theme_bw()+
#   ylab("AUC") + xlab("Fold")

### DISTANCE models ##########################################

bufferData <- readRDS("data/varDF_allEnvironData_buffers_idiv.rds")
bufferData <- subset(bufferData, ! LinjeID %in%
                       c(935,874,876,882,884,936,2317,2328,2338,878,886,1250,1569,2331,2339,1925))

#slurm model
out1 <- readRDS("model-outputs/SLURM/distanceModel/outSummary_linetransectModel_variables_1.rds")

out1 <- readRDS("model-outputs/SLURM/distanceModel/outSummary_linetransectModel_variables_2.rds")

out1 <- readRDS("model-outputs/SLURM/distanceModel/outSummary_linetransectModel_variables_3.rds")

#make into a df
out1 <- data.frame(out1)
out1$Param <- row.names(out1)

### compare models #######################################

#get preds from section below

# #compare with limited line transect model
# out1 <- readRDS("model-outputs/SLURM/distanceModel/limited_line_transects_random/outSummary_linetransectModel_variables_1.rds")
# 
# out1 <- data.frame(out1)
# out1$Param <- row.names(out1)
# preds_limited <- subset(out1,grepl("meanDensity",out1$Param))
# 
# plot(preds$mean,preds_limited$mean)
# abline(0,1)
# 
# #compare with limited gridyear model
# out1 <- readRDS("model-outputs/SLURM/distanceModel/limited_gridyear_random/outSummary_linetransectModel_variables_1.rds")
# 
# out1 <- data.frame(out1)
# out1$Param <- row.names(out1)
# preds_limited <- subset(out1,grepl("meanDensity",out1$Param))
# 
# plot(preds$mean,preds_limited$mean)
# abline(0,1)

#compare different response formulations
# out1 <- readRDS("model-outputs/SLURM/distanceModel/limited_grid_transects_random/outSummary_linetransectModel_variables_1.rds")
# out1 <- data.frame(out1)
# out1$Param <- row.names(out1)
# preds_grid <- subset(out1,grepl("meanDensity",out1$Param))
# 
# #with direct response (and surveyArea as effort term)
# out1 <- readRDS("model-outputs/SLURM/distanceModel/limited_grid_transects_random/direct_response/outSummary_linetransectModel_variables_1.rds")
# out1 <- data.frame(out1)
# out1$Param <- row.names(out1)
# preds_grid_response <- subset(out1,grepl("meanDensity",out1$Param))
# plot(preds_grid$mean,preds_grid_response$mean)
# abline(0,1)
#same!!!!

### plot map ##############################################

#density over range of line transects
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
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.112   8.620  11.773  12.782  15.871  30.875

#including line and grid random effect

#model1
##Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.032   8.535  11.549  12.803  15.975  34.741

#final - model 3
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.007   9.222  12.304  13.502  16.726  32.257

#tmap

bufferData_st <- st_as_sf(bufferData,coords=c("x","y"),crs=equalM)

#mean preds
density_tmap <- tm_shape(NorwayOrigProj)+
  tm_borders()+
tm_shape(bufferData_st)+
  tm_dots("preds",title="Est. Density",
          palette="YlGnBu",size=0.05, style="cont")+
  tm_layout(legend.position=c("left","top"))

#sd of preds
density_tmap_sd <- tm_shape(NorwayOrigProj)+
  tm_borders()+
  tm_shape(bufferData_st)+
  tm_dots("predsSD",title="SD",
          palette="YlGnBu",size=0.05, style="cont")+
  tm_layout(legend.position=c("left","top"))

#saving
temp <- tmap_arrange(density_tmap,density_tmap_sd,nrow=1)

temp

tmap_save(temp, "plots/Fig_3.png",width = 6, height = 4)

### coefficients ################

#average strip width
summary(subset(out1,grepl("meanESW",out1$Param))[,"mean"])

#group size effect
subset(out1,grepl("b.group.size",out1$Param))

#mean density per km
summary(subset(out1,grepl("meanDensity",out1$Param))[,"mean"])

### correlations #################

#mean across all years
preds <- subset(out1,grepl("exp.j",out1$Param))
dataMeans <- apply(bugs.data$NuIndivs,1,mean,na.rm=T)
preds$data <- as.numeric(dataMeans)
summary(preds$data)
summary(preds$mean)

#main plot
qplot(data,mean,data=preds)+
  geom_abline(intercept=0,slope=1)+
  theme_bw()+
  xlab("Observed data")+ylab("Model prediction")

#correlated but noisey
#on log-scale
qplot(data,mean,data=preds)+
  geom_abline(intercept=0,slope=1)+
  scale_x_log10()+scale_y_log10()+
  theme_bw()+
  xlab("Observed data")+ylab("Model prediction")


cor.test(preds$data,preds$mean)#0.87
cor.test(log(preds$data),log(preds$mean))
#model 3 - 0.85

### BPV ###################################################

# #compare data versus preds
# exp <- subset(out1,grepl("exp.j",out1$Param))
# obs <- subset(out1,grepl("NuIndivs.j",out1$Param))
# expNu <- subset(out1,grepl("NuIndivs.new.j",out1$Param))
# 
# par(mfrow=c(2,2))
# hist(dataMeans)
# hist(obs$mean)
# hist(exp$mean)
# hist(expNu$mean)
# #look very similar!!!
# 
# summary(dataMeans)#highest values
# summary(obs$mean)
# summary(exp$mean)
# summary(expNu$mean)
# 
# #fit and fit new
# subset(out1,grepl("fit",out1$Param))
# #fit new is bigger
# 
# #full models
# out1 <- readRDS("model-outputs/SLURM/distanceModel/out_linetransectModel_variables_3.rds")
# out1 <- readRDS("model-outputs/SLURM/distanceModel/out_linetransectModel_variables_1.rds")
# 
# 
# #BPV
# mean(out1$sims.list$fit.new > out1$sims.list$fit)
# 
# par(mfrow=c(2,1))
# hist(out1$sims.list$fit)
# hist(out1$sims.list$fit.new)#larger
# summary(out1$sims.list$fit.new)
# summary(out1$sims.list$fit)
# median(out1$sims$fit)
# abline(v=median(out1$sims.list$fit),col="red")

#limited line transect model
#model 3 - 0.53

#limited grid transect model
#model 1 - 0.17

#with limited gridyear
#model 1 - -0.27

#with adm grid year
#model 1 - 0.213

#with line and grid random effect
#model 1 - 0.72
#model 3 - 0.73


#bpv calculated in the model
#subset(out1,grepl("bpv",out1$Param))
#                 mean        sd X2.5. X25. X50. X75. X97.5.     Rhat n.eff overlap0 f Param
#bpv          0.8291667 0.3765202     0    1    1    1      1 1.001043  1200        1 1   bpv

### Dharma ###############################################

#for each site mean

# library(DHARMa)
# #get simulated data
# simulations = out1$sims.list$NuIndivs.new.j
# #change into a 2-D matrix
# #https://stackoverflow.com/questions/37662433/r-3d-array-to-2d-matrix
# #dim(simulations)
# #dim(simulations) <- c(dim(simulations)[1],599 * 11)
# 
# #get model predictions
# preds = out1$mean$exp.j
# #dim(preds)
# #dim(preds) <- c(dim(simulations)[1],599*11)
# #preds = apply(preds, 2, median)
# #length(preds)
# 
# #get mean of observed data
# obs <- out1$mean$NuIndivs.j
# #obs <- apply(bugs.data$NuIndivs,1,mean,na.rm=T)
# #dim(obs)
# #dim(obs) <- c(599*11)
# #length(obs)
# 
# #now need to remove missing values
# #notMiss <- !is.na(obs)
# #length(notMiss)
# 
# #sim = createDHARMa(simulatedResponse = t(simulations[,notMiss]),
# #                   observedResponse = obs[notMiss],
# #                   fittedPredictedResponse = preds[notMiss],
# #                   integerResponse = T)
# 
# sim = createDHARMa(simulatedResponse = t(simulations),
#                    observedResponse = obs,
#                    fittedPredictedResponse = preds,
#                    integerResponse = T)
# 
# plot(sim)
# 
# #DHARMA tests
# testZeroInflation(sim)
# testDispersion(sim)

### check model calc #####################################

# #mean of the data, calculated by the model
# out1$mean$NuIndivs.j
# 
# #identify sites with complete data
# nuYears <- apply(bugs.data$NuIndivs,1,function(x)sum(!is.na(x)))
# bugs.data$NuIndivs[1:10,]
# #how many sampled in each year
# table(nuYears)
# 
# #45 was sampled in all years
# bugs.data$NuIndivs[45,]
# mean(bugs.data$NuIndivs[45,])
# out1$mean$NuIndivs.j[45]#works!
# 
# #test another
# #543
# bugs.data$NuIndivs[543,]
# mean(bugs.data$NuIndivs[543,])
# out1$mean$NuIndivs.j[543]#works!
# 
# #few were visited in 2007 and 2008

### MAD ###################################################

#quick check
# out1 <- readRDS("model-outputs/SLURM/distanceModel/admyear_random/outSummary_linetransectModel_variables_4.rds")
# 
# qplot(out1$mean[grepl("mean.expNuIndivs",row.names(out1))],
#       out1$mean[grepl("exp.j",row.names(out1))])
# #the same!!!

#null model
out1 <- readRDS("model-outputs/SLURM/distanceModel/out_linetransectModel_variables_4.rds")

#extract dataset data
ggd <- ggs(out1$samples)
out1_dataset <- subset(ggd,grepl("expNuIndivs",ggd$Parameter))
out1_dataset <- subset(out1_dataset,!grepl("mean.expNuIndivs",out1_dataset$Parameter))
out1_dataset$index <- as.numeric(interaction(out1_dataset$Iteration,out1_dataset$Chain))
subset(out1_dataset,Iteration==1 & Chain==1)

#get actual NuIndiv
totalsInfo_long <- as.numeric(gdata::unmatrix(totalsInfo,byrow=T))

#get difference between this value and the simulated values
mad_dataset <- as.numeric()
rmse_dataset <- as.numeric()
n.index <- max(out1_dataset$index)

for(i in 1:n.index){
  mad_dataset[i] <- mean(abs(totalsInfo_long[!is.na(totalsInfo_long)] - 
                               out1_dataset$value[out1_dataset$index==i][!is.na(totalsInfo_long)]))
  
  rmse_dataset[i] <- sqrt(mean((totalsInfo_long[!is.na(totalsInfo_long)] - 
                                  out1_dataset$value[out1_dataset$index==i][!is.na(totalsInfo_long)])^2))
  
}

summary(mad_dataset)
summary(rmse_dataset)

#calculation within code

#MAD
madFiles <- list.files("model-outputs/SLURM/distanceModel", full.names=TRUE) %>%
  str_subset("MAD")  %>%
  map(~ readRDS(.x)) %>%
  reduce(rbind)

#RMSE
rmseFiles <- list.files("model-outputs/SLURM/distanceModel", full.names=TRUE) %>%
  str_subset("RMSE") %>%
  map(~ readRDS(.x)) %>%
  reduce(rbind)

### model selection #######################################

#model 1
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


#or plot gs
gs <- subset(out1,grepl("g",out1$Param))
gs <- subset(gs,!grepl("b.group.size",gs$Param))
gs$variables <- c("bio6","bio5","y","distCoast","tree_line","MountainBirchForest",
                     "Bog","ODF","Meadows","OSF","Mire","SnowBeds",
                     "bio6_2","bio5_2","y_2","distCoast_2","tree_line_2","MountainBirchForest_2",
                     "Bog_2","ODF_2","Meadows_2","OSF_2","Mire_2","SnowBeds_2")
ggplot(gs)+
  geom_col(aes(x=variables,y=mean))+
  coord_flip()+
  theme_bw()+
  ylab("Proportion of model inclusion")+
  geom_hline(yintercept=0,color="red",
             linetype="dashed")+
  xlab("Predictor")

#ones included in more than 25% of models:
#tree line - additive and squared
#bio5
#bio6

### cross validation 1 #####################################

#CV fold models
allFiles <- list.files("model-outputs/SLURM/distanceModel/CV")

#write a function to test CV for each fold and dataset (train or test)
testAbundanceCV <- function(dataset = "train",
                            mymodel = "linetransectModel_variables_CV.txt",
                            fold = 1){

#get files for model and fold  
myFiles <- allFiles[grepl(mymodel,allFiles)]
myFolds <- myFiles[grepl(fold,myFiles)]

#read in main model
out1 <- readRDS(paste("model-outputs/SLURM/distanceModel/CV",myFolds[1],sep="/"))
ggd <- ggs(out1$samples)

#extract dataset data
out1_dataset <- subset(ggd,grepl(paste0("mean.expNuIndivs_",dataset),ggd$Parameter))
out1_dataset$siteIndex <- sub(".*\\[([^][]+)].*", "\\1", as.character(out1_dataset$Parameter))
out1_dataset$siteIndex <- as.numeric(as.character(out1_dataset$siteIndex))
out1_dataset$index <- as.numeric(interaction(out1_dataset$Iteration,out1_dataset$Chain))

#get actual NuIndiv
dataset_fold <- myFolds[grepl(dataset,myFolds)]
totalsInfo <- readRDS(paste("model-outputs/SLURM/distanceModel/CV",dataset_fold,sep="/"))
totalsInfo_mean <- rowMeans(totalsInfo,na.rm=T)

#check they are consistent:
length(unique(out1_dataset$siteIndex)) == dim(totalsInfo)[1]

#plot correlation between mean values
meanVals <- plyr::ddply(out1_dataset,"siteIndex",summarise,pred=mean(value))
meanVals$obs <- totalsInfo_mean

#get difference between this value and the simulated values
mad_dataset <- as.numeric()
rmse_dataset <- as.numeric()
n.index <- max(out1_dataset$index)

for(i in 1:n.index){
  mad_dataset[i] <- mean(abs(totalsInfo_mean[!is.na(totalsInfo_mean)] - 
                       out1_dataset$value[out1_dataset$index==i][!is.na(totalsInfo_mean)]))
  
  rmse_dataset[i] <- sqrt(mean((totalsInfo_mean[!is.na(totalsInfo_mean)] - 
                        out1_dataset$value[out1_dataset$index==i][!is.na(totalsInfo_mean)])^2))
  
}

#package results into a dataset
data.frame(model = mymodel,
           fold = fold,
           dataset = dataset,
           cor = cor(log(meanVals$obs),log(meanVals$pred)),
           mad_mean = mean(mad_dataset),
           mad_sd = sd(mad_dataset),
           mad_median = quantile(mad_dataset, 0.5),
           mad_lower = quantile(mad_dataset, 0.25),
           mad_upper = quantile(mad_dataset, 0.75),
           rmse_mean = mean(rmse_dataset),
           rmse_sd = sd(rmse_dataset),
           rmse_median = quantile(rmse_dataset, 0.5),
           rmse_lower = quantile(rmse_dataset, 0.25),
           rmse_upper = quantile(rmse_dataset, 0.75))
}


#apply function to the model outputs

#standard model
train_1 <- 1:5 %>%
map_dfr(function(i)
  testAbundanceCV(dataset = "train", mymodel = "linetransectModel_variables_CV.txt",fold = i)) 
test_1 <- 1:5 %>%
  map_dfr(function(i)
    testAbundanceCV(dataset = "test", mymodel = "linetransectModel_variables_CV.txt",fold = i))

#lasso model
train_2 <- 1:5 %>%
  map_dfr(function(i)
    testAbundanceCV(dataset = "train", mymodel = "linetransectModel_variables_LASSO_CV.txt",fold = i)) 
test_2 <- 1:5 %>%
  map_dfr(function(i)
    testAbundanceCV(dataset = "test", mymodel = "linetransectModel_variables_LASSO_CV.txt",fold = i))

#model selection model
train_3 <- 1:5 %>%
  map_dfr(function(i)
    testAbundanceCV(dataset = "train", mymodel = "linetransectModel_variables_ModelSelection_CV.txt",fold = i)) 
test_3 <- 1:5 %>%
  map_dfr(function(i)
    testAbundanceCV(dataset = "test", mymodel = "linetransectModel_variables_ModelSelection_CV.txt",fold = i))

#combine
allCV <- bind_rows(train_1,test_1,train_2,test_2,train_3,test_3) 

allCV$model <- as.factor(allCV$model)
levels(allCV$model) <- c("Gaussian priors", "LASSO priors", "Variable indicator")
allCV$model <- factor(allCV$model, levels = c("LASSO priors", "Variable indicator","Gaussian priors"))

#plots
ggplot(allCV) + 
  geom_point(aes(x=fold,y=cor,colour=model))+
  facet_wrap(~dataset)+
  theme_bw() + ylab("Correlation coefficient")

ggplot(allCV) + 
  geom_pointrange(aes(x = fold, y = mad_median, ymin = mad_lower, ymax = mad_upper,
                      colour=model))+
  facet_wrap(~dataset)

ggplot(allCV) + 
  geom_pointrange(aes(x = fold, y = rmse_median, ymin = rmse_lower, ymax = rmse_upper,
                      colour=model))+
  facet_wrap(~dataset)
  
colMeans(subset(allCV, model=="Variable indicator" & dataset=="train")[,-c(1:3)])

colMeans(subset(allCV, model=="Variable indicator" & dataset=="test")[,-c(1:3)])


### cross validation 2 ###################################

#for those that were processed in the HPC script

modelTaskID <- read.delim(paste("data","modelTaskID_occuModel_CV.txt",sep="/"),as.is=T)


#MAD
madFiles <- list.files("model-outputs/SLURM/distanceModel/CV", full.names=TRUE) %>%
                  str_subset("MAD") 

#read in and combine
madData <- madFiles %>%
            map(~ readRDS(.x)) %>%
            reduce(rbind) %>%
            as_tibble()%>%
            add_column(file = madFiles) %>%
            dplyr::mutate(tmp = map_chr(file, ~strsplit(.x, "_")[[1]][5])) %>%
            dplyr::mutate(TaskID = parse_number(tmp)) %>%
            dplyr::mutate(Type = map_chr(file, ~strsplit(.x, "_")[[1]][2])) %>% 
            left_join(.,modelTaskID)
            

#get mean and range for each Model
madData %>%
    dplyr::group_by(Model, Type) %>%
    dplyr::summarise(median = median(Mean),min=min(Mean),max(max(Mean)))
#Model                                          Type  median   min `max(max(Mean))`
#<chr>                                          <chr>  <dbl> <dbl>            <dbl>
#  1 BUGS_occuModel_upscaling_CV.txt                test    7.49  4.95             8.14
#2 BUGS_occuModel_upscaling_CV.txt                train   7.49  6.81             7.69
#3 BUGS_occuModel_upscaling_LASSO_CV.txt          test    7.45  4.98            13.7 
#4 BUGS_occuModel_upscaling_LASSO_CV.txt          train   7.41  6.79             7.64
#5 BUGS_occuModel_upscaling_ModelSelection_CV.txt test    7.43  4.91             8.52
#6 BUGS_occuModel_upscaling_ModelSelection_CV.txt train   7.47  6.84             7.68

#RMSE        
rmseFiles <- list.files("model-outputs/SLURM/distanceModel/CV", full.names=TRUE) %>%
  str_subset("RMSE") 

#read in and combine
rmseData <- rmseFiles %>%
  map(~ readRDS(.x)) %>%
  reduce(rbind) %>%
  as_tibble()%>%
  add_column(file = madFiles) %>%
  dplyr::mutate(tmp = map_chr(file, ~strsplit(.x, "_")[[1]][5])) %>%
  dplyr::mutate(TaskID = parse_number(tmp)) %>%
  dplyr::mutate(Type = map_chr(file, ~strsplit(.x, "_")[[1]][2])) %>% 
  left_join(.,modelTaskID)


#get mean and range for each Model
rmseData %>%
  dplyr::group_by(Model, Type) %>%
  dplyr::summarise(median = median(Mean),min=min(Mean),max(max(Mean)))

### COMBINED model #####################################
#### full models ####

#from combined model:
out1 <- readRDS("model-outputs/SLURM/combinedModel/outSummary_fullCombinedModel.rds")

out1 <- as.data.frame(out1)
out1$Param <- row.names(out1)
subset(out1,grepl("total.pop",out1$Param))
subset(out1,grepl("annual.pop",out1$Param))#infinity
subset(out1,grepl("grid.pop",out1$Param))

#all huge

#### simple models ####
#or get data from "HPC_simple_combined_analysis:
out2 <- data.frame(out1$summary)
out2$Param <- row.names(out2)
preds <- subset(out2,grepl("realAbund",out2$Param))
siteInfo_Occ$preds <- preds$mean
summary(siteInfo_Occ$preds)#minus numbers??? 

#bound to reasonable
siteInfo_Occ$preds[siteInfo_Occ$preds<0] <- 0
siteInfo_Occ$preds[siteInfo_Occ$preds>500] <- 500

mygrid[] <- NA
mygrid[siteInfo_Occ$grid] <- siteInfo_Occ$preds
crs(mygrid) <- equalM
total_tmap <- tm_shape(mygrid)+
  tm_raster(title="Abundance",palette="YlGnBu",style="cont")
total_tmap

#plot uncertainty too
siteInfo_Occ$preds <- preds$sd
mygrid[] <- NA
mygrid[siteInfo_Occ$grid] <- siteInfo_Occ$preds
crs(mygrid) <- equalM
sd_tmap <- tm_shape(mygrid)+
  tm_raster(title="SD",palette="YlGnBu",n=6,style="cont")
sd_tmap

tmap_arrange(total_tmap, sd_tmap, nrow = 1)

### other ##################################################

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