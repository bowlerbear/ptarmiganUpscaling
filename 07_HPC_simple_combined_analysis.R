library(tidyverse)
library(sp)
library(sf)
library(raster)
library(rgeos)
library(maptools)
library(ggthemes)


#local
#myfolder <- "data"
#on HPC
myfolder <- "/data/idiv_ess/ptarmiganUpscaling"

#models 
modelfolderOccu <- "/work/bowler/=ptarmigan_occuModel/11670197"
modelfolderAbund <- "/work/bowler/=ptarmigan_distanceModel/12041062"
list.files(modelfolderOccu)
list.files(modelfolderAbund)

### check siteInfo #######################################################

siteInfo_Occ <- readRDS(paste(myfolder,"siteInfo_ArtsDaten.rds",sep="/"))
siteInfo_Abund <- readRDS(paste(myfolder,"siteInfo_AbundanceModels.rds",sep="/"))

### model choice ########################################

modelList <- c("_1","_2","_3")

task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

mymodel <- modelList[task.id]

### get data #############################################

# get occupancy model predictions
occModel <- ifelse(mymodel=="_2","_3",mymodel)

predsOcc <- readRDS(paste0(modelfolderOccu,"/",
                          "Z_occModel_upscaling",occModel,".rds"))
names(predsOcc)
head(predsOcc)

# get abundance predictions

predsDensity <- readRDS(paste0(modelfolderAbund, "/",
                               "Density.pt_linetransectModel_variables",mymodel,".rds")) 

names(predsDensity)
head(predsDensity)

### sort indices #######################################

#predsDensity <- as.data.frame(predsDensity)
#predsDensity$Parameter <- row.names(predsDensity)
predsDensity$ParamNu <-  as.character(sub(".*\\[([^][]+)].*", "\\1", predsDensity$Parameter))
predsDensity <- subset(predsDensity,Parameter!="deviance")
predsDensity <- predsDensity %>%
                  separate(ParamNu, c("Site", "Year"))
head(predsDensity)

#predsOcc <- as.data.frame(predsOcc)
#predsOcc$Parameter <- row.names(predsOcc)
predsOcc$ParamNu <-  as.character(sub(".*\\[([^][]+)].*", "\\1", predsOcc$Parameter))
predsOcc <- subset(predsOcc,Parameter!="deviance")
predsOcc <- predsOcc %>%
                  separate(ParamNu, c("Site", "Year"))
head(predsOcc)
### mean occ ############################################

#occupancy
predsOcc_summary <- predsOcc %>%
                    group_by(Site) %>%
                    summarise(myMean = mean(value),
                              mySD = sd(value)) %>%
                    mutate(Site = as.numeric(Site)) %>%
                    arrange(Site)

saveRDS(predsOcc_summary, file=paste("predsOcc_summary",mymodel,".rds"))

### mean density ##########################################

#density - average across all years

predsDensity_summary <- predsDensity %>%
  group_by(Site) %>%
  summarise(myMean = mean(value),
            mySD = sd(value)) %>%
  mutate(Site = as.numeric(Site)) %>%
  arrange(Site)

saveRDS(predsDensity_summary, file=paste("predsDensity_summary",mymodel,".rds"))

### cap density ########################################

#cap high values
#get maximum elevation at which one is seen at
#1658 m - but few above 1500 m
predsDensity_summary$myMean[siteInfo_Abund$elevation>1660] <- 0

#what is the maximum observed density? cap by this
#high mean density was 34 per km2
34*25#850
summary(predsDensity_summary$myMean)
mean(predsDensity_summary$myMean>850)#only 0.1% of grids

predsDensity_summary$myMean[predsDensity_summary$myMean>850] <- 850
siteInfo_Abund$preds <- predsDensity_summary$myMean

#save again
saveRDS(predsDensity_summary, file=paste("predsDensity_summary",mymodel,".rds"))

### mean multiplication ####################################

siteInfo_Abund$meanDensity <- predsDensity_summary$myMean
summary(siteInfo_Abund$meanDensity)

siteInfo_Occ$meanOcc <- predsOcc_summary$myMean
summary(siteInfo_Occ$meanOcc)

siteInfo <- inner_join(siteInfo_Abund, siteInfo_Occ, by="grid")
siteInfo$realAbund <- siteInfo$meanDensity * siteInfo$meanOcc
sum(siteInfo$realAbund)
#still around 1 million!!!

### iteration multiplication ###############################

head(predsOcc)
head(predsDensity)

#check they align
all(predsOcc$ParamNu==predsDensity$ParamNu)

#add on density data to occupancy data frame
predsOcc$Density <- predsDensity$value
predsOcc$realDensity <- predsOcc$value * predsOcc$Density

### annual total population size #######################

population_Annual <- predsOcc %>%
                      group_by(Iteration,Chain,Year) %>%
                      summarise(totalPop = sum(realDensity)) %>%
                      group_by(Year) %>%
                      summarise(medianPop = median(totalPop),
                        lowerPop = quantile(totalPop, 0.025),
                        upperPop = quantile(totalPop, 0.975))

saveRDS(population_Annual, file=paste("population_Annual",mymodel,".rds"))
                  
### mean total abundance #################################

nuYears <- n_distinct(predsOcc$Year)

totalSummary <- predsOcc %>%
                  group_by(Iteration,Chain) %>%
                  summarise(totalPop = sum(realDensity)/nuYears) 

saveRDS(totalSummary, file=paste("totalSummary",mymodel,".rds"))

#### end ###################################



