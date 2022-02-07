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
preds <- z_preds <- subset(out1,grepl("grid.z",out1$Param))
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
#fold.id AUC_psi_test AUC_psi_train AUC_py_test AUC_py_train AUC_pypred_test AUC_pypred_train
#<int>        <dbl>         <dbl>       <dbl>        <dbl>           <dbl>            <dbl>
#1       1        0.860         0.859       0.740        0.967           0.914            0.922
#2       2        0.829         0.865       0.661        0.970           0.879            0.927
#3       3        0.870         0.866       0.719        0.969           0.883            0.925
#4       4        0.865         0.846       0.734        0.966           0.942            0.913
#5       5        0.877         0.850       0.737        0.965           0.930            0.917

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

subset(out1,grepl("MAD",out1$Param))

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
    dplyr::summarise(median = median(Mean),min=min(Mean),max=max(Mean))
# Groups:   Model [3]
#<chr>                                          <chr>  <dbl> <dbl> <dbl>
#  1 BUGS_occuModel_upscaling_CV.txt                test    5.20  3.39  6.65
#2 BUGS_occuModel_upscaling_CV.txt                train   5.01  3.81  6.24
#3 BUGS_occuModel_upscaling_LASSO_CV.txt          test    6.08  3.67 11.6 
#4 BUGS_occuModel_upscaling_LASSO_CV.txt          train   4.83  3.73  6.18
#5 BUGS_occuModel_upscaling_ModelSelection_CV.txt test    5.21  3.53  7.12
#6 BUGS_occuModel_upscaling_ModelSelection_CV.txt train   4.93  3.84  6.22

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

modelFolder <- "model-outputs/SLURM/combinedModel/"

#choose model
models = c(" _1 ")
models = c(" _2 ")
models = c(" _3 ")

siteInfo_Occ <- readRDS("data/siteInfo_ArtsDaten.rds")
siteInfo_Abund <- readRDS("data/siteInfo_AbundanceModels.rds")

#plot occupancy preds
predsOcc_summary <- readRDS(paste0(modelFolder,"predsOcc_summary", models,".rds"))
siteInfo_Occ$preds <- predsOcc_summary$myMean
mygrid[] <- NA
mygrid[siteInfo_Occ$grid] <-siteInfo_Occ$preds
plot(mygrid)

#plot density preds
predsDensity_summary <- readRDS(paste0(modelFolder,"predsDensity_summary", models,".rds"))
siteInfo_Abund$preds <- predsDensity_summary$myMean
mygrid[] <- NA
mygrid[siteInfo_Abund$grid] <-siteInfo_Abund$preds
plot(mygrid)

#simple combination (just to check)
siteInfo_Occ$abund <- siteInfo_Abund$preds[match(siteInfo_Occ$grid,
                                                 siteInfo_Abund$grid)]
siteInfo_Occ$preds <- siteInfo_Occ$preds * siteInfo_Occ$abund
mygrid[] <- NA
mygrid[siteInfo_Occ$grid] <- siteInfo_Occ$preds
plot(mygrid)

#look at annual preds
(populationAnnual <- readRDS(paste0(modelFolder,"population_Annual", 
                                       models,".rds")))

g1 <- populationAnnual %>%
  mutate(Year = as.numeric(Year)+2007) %>%
  ggplot()+
  geom_pointrange(aes(x=Year, y=medianPop, 
                      ymin=lowerPop, ymax=upperPop))+
  ylab("Predicted national population size")+
  scale_x_continuous(breaks=c(2008,2012,2016),
                     labels=c(2008,2012,2016))+
  theme_few()

#look at mean preds
totalSummary <- readRDS(paste0(modelFolder,"totalSummary", models,".rds"))

g2 <- ggplot(totalSummary, aes(totalPop))+
  geom_density(fill="red",colour="black",alpha=0.2)+
  geom_density(colour="black", size=2)+
  theme_few()+
  ylab("Probability density") + xlab("Total population size estimate")

cowplot::plot_grid(g2,g1,labels=c("a","b"),nrow=1,ncol=2)

ggsave("plots/Fig.5.png",width=9.5,height=4)

summary(totalSummary$totalPop)
quantile(totalSummary$totalPop,c(0.025,0.5,0.975))

#short cut
sum(predsOcc_summary$myMean)
mean(predsDensity_summary$myMean)
sum(predsOcc_summary$myMean)*mean(predsDensity_summary$myMean)
#with model 3 - 1207997

#plot grid predictions
(populationGrid <- readRDS(paste0(modelFolder,"population_grid", 
                                    models,".rds")))

mygrid[] <- NA
mygrid[populationGrid$grid] <- populationGrid$medianPop
crs(mygrid) <- equalM
occ_tmap <- tm_shape(mygrid)+
  tm_raster(title="a) Abundance",palette="YlGnBu",
            breaks=c(0,50,100,200,400,800))
occ_tmap

mygrid[] <- NA
mygrid[populationGrid$grid] <- populationGrid$sdPop
plot(mygrid)
crs(mygrid) <- equalM
occ_tmap_sd <- tm_shape(mygrid)+
  tm_raster(title="b) SD   ",palette="YlGnBu",
            breaks=c(0,50,100,200,400))
occ_tmap_sd

temp <- tmap_arrange(occ_tmap,occ_tmap_sd,nrow=1)
temp

tmap_save(temp, "plots/Fig_4.png",width = 6, height = 4)

#compare priors
library(ggridges)
all_ggd <- list.files(paste(modelFolder,"12198176",sep="/"), full.names=TRUE) %>%
            str_subset("totalSummary") %>%
            set_names() %>%
            map_dfr(readRDS, .id="source") %>%
            as_tibble()

all_ggd$Model <- NA
all_ggd$Model[grepl("_1", all_ggd$source)] <- "Only important vars"
all_ggd$Model[grepl("_2", all_ggd$source)] <- "LASSO priors"
all_ggd$Model[grepl("_3", all_ggd$source)] <- "Variable indicator"

ggplot(all_ggd, aes(x = totalPop, y = Model)) + 
  geom_density_ridges2()+
   theme_minimal()+xlab("Total population size")+
   xlim(900000,1800000)
 
ggplot(all_ggd, aes(totalPop))+
  geom_density(fill="red",colour="black",alpha=0.2)+
  geom_density(colour="black", size=2)+
  theme_few()+
  facet_wrap(~Model,scales="free")+
  ylab("Probability density") + xlab("Total population size estimate")

### make model predictions ####

#for use in the uncertainity analysis

varDF <- readRDS("data/varDF_allEnvironData_5km_idiv.rds")

#run sections above to get z_preds
out_occ <- z_preds
out_occ$grid <- siteInfo_Occ$grid
out_occ <- arrange(out_occ,grid)
all(varDF$grid==out_occ$grid)

#abundance predictions
out1 <- readRDS("model-outputs/SLURM/distanceModel/outSummary_linetransectModel_variables_3.rds")
out1 <- data.frame(out1)
out1$Param <- row.names(out1)
out1 <- subset(out1,grepl("Density.p",out1$Param))
out1$grid <- siteInfo_Abund$grid
out_dens <- arrange(out1, grid)

#check align
all(out_occ$grid, out_dens$grid)

#make modelPredictions
modelPredictions <- data.frame(density_mean=out_dens$mean, density_sd = out_dens$sd,
                               occ_mean=out_occ$mean, occ_sd = out_occ$sd)

#add x and y
modelPredictions$x <- varDF$x
modelPredictions$y <- varDF$y

saveRDS(modelPredictions, file="data/modelPredictions.rds")

### end #####
  