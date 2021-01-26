library(sp)
library(raster)
library(maptools)
library(ggplot2)
library(rgeos)
library(plyr)

#HPC
myfolder <- "/data/idiv_ess/ptarmiganUpscaling" 
#local
myfolder <- "data"

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
listlengthDF <- subset(listlengthDF,!is.na(y))
siteInfo <- subset(listlengthDF,!duplicated(grid))

### folds #############################################################

folds <- readRDS(paste(myfolder,"folds_occModel.rds",sep="/"))
listlengthDF$fold <- folds$fold[match(listlengthDF$grid,folds$grid)]

#select fold of this task
fold.id = as.integer(Sys.getenv("SGE_TASK_ID", "1"))

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
                  y_test = listlengthDF_test$y,
                  y_train = listlengthDF_train$y,
                  #add an adm effect
                  adm_train = siteInfo_train$admN,
                  det.adm_train = listlengthDF_train$admN,
                  det.open_train = scale(listlengthDF_train$Open),
                  n.adm_train = length(unique(siteInfo_train$admN)),
                  adm2 = siteInfo$admN2,
                  det.adm2_train = listlengthDF_train$admN2,
                  n.adm2 = length(unique(siteInfo$admN2)))

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

### scale vars #################################################################

siteInfo$tree_line_position <- scale(siteInfo$tree_line_position)
siteInfo$bio1 <- scale(siteInfo$bio1)
siteInfo$bio5 <- scale(siteInfo$bio5)
siteInfo$bio6 <- scale(siteInfo$bio6)
siteInfo$elevation <- scale(siteInfo$elevation)
siteInfo$PrefOpen <- scale(siteInfo$PrefOpen)
siteInfo$Open <- scale(siteInfo$Open)
siteInfo$Top <- scale(siteInfo$Top)

siteInfo_test <- subset(listlengthDF_test,!duplicated(grid))
siteInfo_train <- subset(listlengthDF_train,!duplicated(grid))


### fit model ########################################################

#specify model structure
bugs.data$occDM_train <- model.matrix(~ siteInfo_train$tree_line_position + 
                                  siteInfo_train$tree_line_position^2 +
                                  siteInfo_train$bio1 + 
                                  siteInfo_train$bio5 + 
                                  siteInfo_train$bio6 +
                                  siteInfo_train$elevation +
                                  siteInfo_train$Open + 
                                  siteInfo_train$Top)[,-1]

bugs.data$occDM_test <- model.matrix(~ siteInfo_test$tree_line_position + 
                                        siteInfo_test$tree_line_position^2 +
                                        siteInfo_test$bio1 + 
                                        siteInfo_test$bio5 + 
                                        siteInfo_test$bio6 +
                                        siteInfo_test$elevation +
                                        siteInfo_test$Open + 
                                        siteInfo_test$Top)[,-1]

bugs.data$n.covs <- ncol(bugs.data$occDM_test)

params <- c("mean.p","beta","beta.effort","beta.det.open","pred.muZ")

modelfile <- paste(myfolder,"BUGS_occuModel_upscaling_CV.txt",sep="/")

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

saveRDS(out1$summary,file="outSummary_occModel_upscaling_CV.rds")
