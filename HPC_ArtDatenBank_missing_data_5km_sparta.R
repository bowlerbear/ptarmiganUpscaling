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

### effort ######################################################################

#we will treat it as a categorical variable
table(listlengthDF$L)
listlengthDF$singleton <- ifelse(listlengthDF$L==1,1,0)

### Absences ###################################################################

listlengthDF$L[is.na(listlengthDF$L)] <- 1 #set to nominal effort
#listlengthDF$y[listlengthDF$L==0] <- 0

### indices #####################################################################

#order data by site and year:
listlengthDF$siteIndex <- as.numeric(factor(listlengthDF$grid))
listlengthDF$yearIndex <- as.numeric(factor(listlengthDF$year))
listlengthDF <- arrange(listlengthDF,siteIndex,yearIndex)

#merge with environ data
listlengthDF <- merge(listlengthDF,varDF,by="grid",all.x=T)

#have missing singletons as 1
listlengthDF$singleton[is.na(listlengthDF$singleton)] <- 1

### adm inidices #############################################################

listlengthDF$admN <- as.numeric(factor(listlengthDF$adm))
listlengthDF$admN2 <- as.numeric(factor(listlengthDF$adm2))

#extract site data
siteInfo <- subset(listlengthDF,!duplicated(grid))

### BUGS object ################################################################

#for BUGS

bugs.data <- list(nsite = length(unique(listlengthDF$siteIndex)),
                  nyear = length(unique(listlengthDF$yearIndex)),
                  nvisit = nrow(listlengthDF),
                  site = listlengthDF$siteIndex,
                  year = listlengthDF$yearIndex,
                  Effort = listlengthDF$singleton,
                  y = listlengthDF$y,
                  #add an adm effect
                  adm = siteInfo$admN,
                  det.adm = listlengthDF$admN,
                  n.adm = length(unique(siteInfo$admN)),
                  adm2 = siteInfo$admN2,
                  det.adm2 = listlengthDF$admN2,
                  n.adm2 = length(unique(siteInfo$admN2)),
                  #environcovs
                  bio1 = scale(siteInfo$bio1),
                  bio1_2 = scale(siteInfo$bio1^2),
                  bio6 = scale(siteInfo$bio6),
                  bio5 = scale(siteInfo$bio5),
                  forest = scale(siteInfo$Forest),
                  open = scale(siteInfo$Open),
                  prefopen = scale(log(siteInfo$PrefOpen+1)),
                  prefclosed = scale(log(siteInfo$PrefClosed+1)),
                  top = scale(log(siteInfo$Top+1)),
                  alpine_habitat1 = siteInfo$alpine_habitat1,
                  alpine_habitat2 = log(siteInfo$alpine_habitat2+1),
                  alpine_habitat3 = log(siteInfo$alpine_habitat3+1),
                  alpine_habitat4 = log(siteInfo$alpine_habitat4+1),
                  elevation = scale(siteInfo$elevation),
                  tree_line_position = scale(siteInfo$tree_line_position),
                  tree_line_position2 = scale(siteInfo$tree_line_position^2))

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
zst <- reshape2::acast(listlengthDF, siteIndex~yearIndex, value.var="y",fun=max,na.rm=T)
zst [is.infinite(zst)] <- 0
inits <- function(){list(z = zst)}

### fit model ########################################################

#specify model structure
bugs.data$occDM <- model.matrix(~ bugs.data$tree_line_position + 
                                  bugs.data$tree_line_position2 +
                                  bugs.data$bio1 + 
                                  bugs.data$bio1_2 + 
                                  bugs.data$bio6 +
                                  bugs.data$elevation +
                                  bugs.data$prefopen + 
                                  bugs.data$open)[,-1]

bugs.data$n.covs <- ncol(bugs.data$occDM)

params <- c("mean.p","beta","beta.effort","beta.det.open","grid.muZ")

modelfile <- "/data/idiv_ess/ptarmiganUpscaling/BUGS_occuModel_upscaling.txt"

n.cores = as.integer(Sys.getenv("NSLOTS", "1")) 

n.iterations = 100

out1 <- jags(bugs.data, 
             inits = inits, 
             params, 
             modelfile, 
             n.thin = 10, 
             n.chains = n.cores, 
             n.burnin = round(n.iterations/3),
             n.iter = n.iterations,
             parallel = T)

saveRDS(out1,file="out_occModel_upscaling.rds")
