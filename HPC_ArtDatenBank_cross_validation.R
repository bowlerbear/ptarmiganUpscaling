library(sp)
library(raster)
library(maptools)
library(ggplot2)
library(rgeos)
library(plyr)

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
listlengthDF <- readRDS("/data/idiv_ess/ptarmiganUpscaling/listlength_iDiv.rds")

### subset ##########################################################

#subset to focal grids and those with environ covariate data

focusGrids <- readRDS("/data/idiv_ess/ptarmiganUpscaling/focusGrids.rds")
varDF <- readRDS("/data/idiv_ess/ptarmiganUpscaling/varDF_allEnvironData_5km_idiv.rds")
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
listlengthDF <- subset(listlengthDF,!is.na(y))

### folds #############################################################

folds <- readRDS("data/folds_occuModels.rds")
listlengthDF$fold <- folds$fold[match(listlengthDF$grid,folds$grid)]

#select fold of this task
fold.id = as.integer(Sys.getenv("SGE_TASK_ID", "1"))

#split intp test and train
listlengthDF_test <- subset(listlength,fold == fold.id)
listlengthDF_train<- subset(listlength,fold != fold.id)

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

bugs.data <- list(nsite_test = length(unique(listlengthDF_test$siteIndex)),
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
                  y_test = listlengthDF_test$y,
                  y_train = listlengthDF_train$y,
                  #add an adm effect
                  adm = siteInfo$admN,
                  det.adm = listlengthDF$admN,
                  det.open = scale(listlengthDF$Open),
                  n.adm = length(unique(siteInfo$admN)),
                  adm2 = siteInfo$admN2,
                  det.adm2 = listlengthDF$admN2,
                  n.adm2 = length(unique(siteInfo$admN2)),
                  #environcovs
                  bio1_test = scale(siteInfo_test$bio1),
                  bio1_train = scale(siteInfo_train$bio1),
                  bio1_2_test = scale(siteInfo_test$bio1^2),
                  bio1_2_train = scale(siteInfo_train$bio1^2),
                  bio6_test = scale(siteInfo_test$bio6),
                  bio6_train = scale(siteInfo_train$bio6),
                  bio5_test = scale(siteInfo_test$bio5),
                  bio5_train = scale(siteInfo_train$bio5),
                  forest_test = scale(siteInfo_test$Forest),
                  forest_train = scale(siteInfo_train$Forest),
                  open_test = scale(siteInfo_test$Open),
                  open_train = scale(siteInfo_train$Open),
                  prefopen_test = scale(log(siteInfo_test$PrefOpen+1)),
                  prefopen_train = scale(log(siteInfo_train$PrefOpen+1)),
                  prefclosed_test = scale(log(siteInfo_test$PrefClosed+1)),
                  prefclosed_train = scale(log(siteInfo_train$PrefClosed+1)),
                  top_test = scale(log(siteInfo_test$Top+1)),
                  top_train = scale(log(siteInfo_train$Top+1)),
                  elevation_test = scale(siteInfo_test$elevation),
                  elevation_train = scale(siteInfo_train$elevation),
                  tree_line_position_test = scale(siteInfo_test$tree_line_position),
                  tree_line_position_train = scale(siteInfo_train$tree_line_position),
                  tree_line_position2_test = scale(siteInfo_test$tree_line_position^2),
                  tree_line_position2_train = scale(siteInfo_train$tree_line_position^2))

#bugs.data_ArtsDaten <- bugs.data

#alpine_habitat:
#1= Open lowland, 
#2 = Low alpine zone, 
#3 = intermediate alpine zone, 
#4 = high alpine zone 

### initials ####################################################################

#get JAGS libraries
library(rjags)
library(jagsUI)

#need to specify initial values
zst <- reshape2::acast(listlengthDF_train, 
                       siteIndex~yearIndex, value.var="y",fun=max,na.rm=T)
zst [is.infinite(zst)] <- 0
inits <- function(){list(z = zst)}

### fit model ########################################################

#specify model structure
bugs.data$occDM_train <- model.matrix(~ bugs.data$tree_line_position_train + 
                                  bugs.data$tree_line_position2_train +
                                  bugs.data$bio1_train + 
                                  bugs.data$bio1_2_train + 
                                  bugs.data$bio6_train +
                                  bugs.data$elevation_train +
                                  bugs.data$prefopen_train + 
                                  bugs.data$open_train)[,-1]

bugs.data$occDM_test <- model.matrix(~ bugs.data$tree_line_position_test + 
                                        bugs.data$tree_line_position2_test +
                                        bugs.data$bio1_test + 
                                        bugs.data$bio1_2_test + 
                                        bugs.data$bio6_test +
                                        bugs.data$elevation_test +
                                        bugs.data$prefopen_test + 
                                        bugs.data$open_test)[,-1]

bugs.data$n.covs <- ncol(bugs.data$occDM_test)

params <- c("mean.p","beta","beta.effort","beta.det.open","grid.z")

modelfile <- "/data/idiv_ess/ptarmiganUpscaling/BUGS_occuModel_upscaling_crossvalid.txt"

n.cores = as.integer(Sys.getenv("NSLOTS", "1")) 

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

saveRDS(out1$summary,file="outSummary_occModel_upscaling_crossvalid.rds")
