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

out1 <- readRDS("model-outputs/SLURM/occModel/outSummary_occModel_upscaling_1.rds")
out1 <- readRDS("model-outputs/SLURM/occModel/outSummary_occModel_upscaling_2.rds")
out1 <- readRDS("model-outputs/SLURM/occModel/outSummary_occModel_upscaling_3.rds")
out1 <- readRDS("model-outputs/SLURM/occModel/outSummary_occModel_upscaling_4.rds")
out1 <- readRDS("model-outputs/SLURM/occModel/outSummary_occModel_upscaling_5.rds")
out1 <- readRDS("model-outputs/SLURM/occModel/outSummary_occModel_upscaling_6.rds")

### plot map #################################################

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

subset(out1,grepl("bpv",out1$Param))

### AUC ##############################################

list.files("model-outputs/SLURM/occModel", full.names=TRUE) %>%
  str_subset("AUC_") %>%
  set_names() %>%
  map_dfr(readRDS, .id = "File")

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

out1 <- readRDS("model-outputs/SLURM/distanceModel/outSummary_linetransectModel_variables_1.rds")

out1 <- readRDS("model-outputs/SLURM/distanceModel/outSummary_linetransectModel_variables_2.rds")

out1 <- readRDS("model-outputs/SLURM/distanceModel/outSummary_linetransectModel_variables_3.rds")

out1 <- readRDS("model-outputs/SLURM/distanceModel/outSummary_linetransectModel_variables_4.rds")

out1 <- readRDS("model-outputs/SLURM/distanceModel/outSummary_linetransectModel_variables_5.rds")

out1 <- readRDS("model-outputs/SLURM/distanceModel/outSummary_linetransectModel_variables_6.rds")

#make into a df
out1 <- data.frame(out1)
out1$Param <- row.names(out1)
table(out1$Rhat<1.1)

### plot map ##############################################

#density over range of line transects
preds <- subset(out1,grepl("meanDensity",out1$Param))
bufferData$preds <- preds$mean
bufferData$predsSD <- preds$sd
summary(bufferData$preds)

#final - model 3
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.186   8.755  11.731  12.836  16.028  32.253

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

#overdispersion?
summary(subset(out1,grepl("r",out1$Param))[,"mean"])

### correlations #################

#mean across all years
preds <- subset(out1,grepl("exp.j",out1$Param))
preds$data <- apply(bugs.data$NuIndivs,1,mean,na.rm=T)
summary(preds$data)
summary(preds$mean)
summary(abs(preds$data-preds$mean))

#main plot
qplot(data,mean,data=preds)+
  geom_abline(intercept=0,slope=1)+
  theme_bw()+
  xlab("Observed data")+ylab("Model prediction")

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

subset(out1,grepl("bpv",out1$Param))

### MAD ###################################################

#MAD
list.files("model-outputs/SLURM/distanceModel", full.names=TRUE) %>%
  str_subset("MAD")  %>%
  map(~ readRDS(.x)) %>%
  reduce(rbind)

#RMSE
list.files("model-outputs/SLURM/distanceModel", full.names=TRUE) %>%
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


#tree line position not important for density

### cross validation ###################################

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