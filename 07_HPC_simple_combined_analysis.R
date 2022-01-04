#choose folder:

#myfolder <- "data"
myfolder <- "/data/idiv_ess/ptarmiganUpscaling" #HPC

### get data #############################################

# get occupancy model predictions
out1 <- readRDS("model-outputs/SLURM/occModel/outSummary_occModel_upscaling_1.rds")
out1 <- readRDS("model-outputs/SLURM/occModel/outSummary_occModel_upscaling_2.rds")
out1 <- readRDS("model-outputs/SLURM/occModel/outSummary_occModel_upscaling_3.rds")

out1 <- readRDS("model-outputs/SLURM/occModel/detection_covariates2/outSummary_occModel_upscaling_1.rds")
out1 <- readRDS("model-outputs/SLURM/occModel/detection_covariates2/outSummary_occModel_upscaling_2.rds")
out1 <- readRDS("model-outputs/SLURM/occModel/detection_covariates2/outSummary_occModel_upscaling_3.rds")


#extract part of output we want
out1 <- data.frame(out1)
out1$Param <- row.names(out1)
predsOcc <- subset(out1,grepl("grid.z",out1$Param))
head(predsOcc)

# get abundance predictions

#negative binomial
out1 <- readRDS("model-outputs/SLURM/distanceModel/outSummary_linetransectModel_variables_1.rds")
out1 <- readRDS("model-outputs/SLURM/distanceModel/outSummary_linetransectModel_variables_2.rds")
out1 <- readRDS("model-outputs/SLURM/distanceModel/outSummary_linetransectModel_variables_3.rds")

out1 <- readRDS("model-outputs/SLURM/distanceModel/admyear_random/outSummary_linetransectModel_variables_1.rds")
out1 <- readRDS("model-outputs/SLURM/distanceModel/admyear_random/outSummary_linetransectModel_variables_2.rds")
out1 <- readRDS("model-outputs/SLURM/distanceModel/admyear_random/outSummary_linetransectModel_variables_3.rds")

out1 <- readRDS("model-outputs/SLURM/distanceModel/full_linetransect_random/outSummary_linetransectModel_variables_1.rds")
out1 <- readRDS("model-outputs/SLURM/distanceModel/full_linetransect_random/outSummary_linetransectModel_variables_2.rds")
out1 <- readRDS("model-outputs/SLURM/distanceModel/full_linetransect_random/outSummary_linetransectModel_variables_3.rds")

#extract part of output we want
out1 <- data.frame(out1)
out1$Param <- row.names(out1)
preds <- subset(out1,grepl("Density",out1$Param))
predsDensity <- subset(preds,!grepl("meanDensity",preds$Param))
head(predsDensity)

predsDensity$grid <- siteInfo_ArtsDaten$grid
predsDensity <- arrange(predsDensity,grid)

### combine data #########################################

modelPredictions <- data.frame(occ_mean = predsOcc$mean, occ_sd = predsOcc$sd,
                       density_mean = predsDensity$mean, density_sd = predsDensity$sd)
#saveRDS(modelPredictions,file="data/modelPredictions.rds")

summary(modelPredictions)

#model 2 - supper high density-mean prediction!
upperQ <- as.numeric(quantile(modelPredictions$density_mean,0.95))
modelPredictions$density_mean[modelPredictions$density_mean>upperQ] <- upperQ

upperQ <- as.numeric(quantile(modelPredictions$density_sd,0.95))
modelPredictions$density_sd[modelPredictions$density_sd>upperQ] <- upperQ

### make bugs object ######################################

#modelPredictions <- readRDS(paste(myfolder,"modelPredictions.rds",sep="/"))

bugs.data <- list(occ_mean = modelPredictions$occ_mean,
                  occ_sd = modelPredictions$occ_sd,
                  density_mean = modelPredictions$density_mean,
                  density_sd = modelPredictions$density_sd,
                  nsite = nrow(modelPredictions),
                  dummy = rnorm(10))

### run model #############################################

#just run locally

library(rjags)
library(jagsUI)

params <- c("totalAbund","realAbund")

#modelfile <- paste(myfolder,"simpleCombinedModel.txt",sep="/")
modelfile <- paste("models","CombinedModel_simple.txt",sep="/")

out1 <- jags(bugs.data, 
             inits=NULL, 
             params, 
             modelfile, 
             n.thin=10,
             n.chains=3, 
             n.burnin=100,
             n.iter=1000,
             parallel=T)

#### plot predictions ############################

#see model summaries - COMBINED ANALYSIS section

### draw posterior of total abund #################

library(ggmcmc)
library(ggthemes)

ggd <- ggs(out1$samples)
ggd <- subset(ggd, grepl("totalAbund",ggd$Parameter))

ggplot(ggd)+
  geom_density(aes(value),fill="red",colour="red",
               alpha=0.2)+
  theme_few()+
  ylab("Probability density") + xlab("Total population size estimate")

#### geom ridge ########################################

library(ggridges)

ggd1 <- ggd
ggd1$Model <- "Gaussian priors"
ggd2 <- ggd 
ggd2$Model <- "LASSO priors"
ggd3 <- ggd
ggd3$Model <- "Variable indicator priors"


all_ggd <- rbind(ggd1,ggd2,ggd3)

ggplot(all_ggd, aes(x = value, y = Model)) + geom_density_ridges2()+
  theme_minimal()+xlab("Total population size")

quantile(ggd$value,c(0.025,0.5,0.975))

### compare #############################################

summary(ggd$value)
quantile(ggd$value,c(0.025,0.5,0.975))

#slurm 1
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1171631 1180391 1183465 1183565 1186924 1197616


#slurm 2
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-87971 1118461 1448297 1450196 1786439 2897063

#2.5%       50%     97.5% 
#472991.1 1448296.8 2601972.

#cutting extremes at 95%
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1269736 1280204 1283506 1283720 1287227 1299197

#  2.5%     50%   97.5% 
#1272175 1283506 1296064

#slurm 3
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1254265 1264654 1267603 1267962 1271253 1285433

#2.5%     50%   97.5% 
#1256830 1267603 1278905

#with a Poisson model
#c 1 million

#with full random adm model

#model 1
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1048663 1073020 1080751 1081351 1089259 1117990 
#quantile(ggd$value,c(0.025,0.5,0.975))
#2.5%     50%   97.5% 
#1059141 1080751 1105082

#model 2
#summary(ggd$value)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1046279 1066540 1076615 1076669 1086662 1115896 
#quantile(ggd$value,c(0.025,0.5,0.975))
#2.5%     50%   97.5% 
#1052550 1076615 1103255

#model 3
#summary(ggd$value)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1015647 1045521 1053592 1053198 1061058 1091031 
#quantile(ggd$value,c(0.025,0.5,0.975))
#2.5%     50%   97.5% 
#1029577 1053592 1073368 

#### occu vs abund #####

library(ggplot2)

summary(modelPredictions$density_mean)

pairs(modelPredictions)

ggplot(modelPredictions,aes(x = occ_mean, y = density_mean))+
  geom_point(method="gam")+
  stat_smooth()

cor.test(modelPredictions$occ_mean, modelPredictions$density_mean)
#slightly negative???

#limit to sites with line transects
modelPredictions$grid <- siteInfo_ArtsDaten$grid
modelPredictions_subset <- subset(modelPredictions, grid %in% siteIndex_linetransects$grid)
cor.test(modelPredictions_subset$occ_mean, modelPredictions_subset$density_mean)
ggplot(modelPredictions_subset,aes(x = occ_mean, y = density_mean))+
  geom_point(method="gam")+
  stat_smooth()

#density in areas predicted to be occupied
modelPredictions_high <- subset(modelPredictions, occ_mean>0.9)
nrow(modelPredictions_high)
summary(modelPredictions_high$density_mean)
#or now
modelPredictions_low <- subset(modelPredictions, occ_mean<0.1)
nrow(modelPredictions_low)
summary(modelPredictions_low$density_mean)


#relationship after correction
modelPredictions$estAbund <- modelPredictions$occ_mean * modelPredictions$density_mean
ggplot(modelPredictions,aes(x = occ_mean, y = estAbund))+
  geom_point(method="gam")+
  stat_smooth()

#long-way to get total density
#number of occupied grids
sum(modelPredictions$occ_mean) * mean(modelPredictions$density_mean)
