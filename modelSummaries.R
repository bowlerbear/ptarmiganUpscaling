library(raster)
library(sp)
library(maptools)
library(rgeos)

### check siteInfo #######################################################

siteInfo_Occ <- readRDS("data/siteInfo_ArtsDaten.rds")
siteInfo_Abund <- readRDS("data/siteInfo_LineTransects.rds")

#check all are the same
all(siteInfo_Occ$grid==siteInfo_Abund$grid)
all(siteInfo_Occ$siteIndex==siteInfo_Abund$siteIndex)
#TRUE

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

### occu models ##############################################

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

### distance models ##########################################

siteInfo <- readRDS("data/siteInfo_LineTransects.rds")
environData <- readRDS("data/environData.rds")

### plot map ##############################################

out1 <- readRDS("model-outputs/outSummary_linetransectModel_variables.rds")
out1[row.names(out1)=="totalPop",]

out1 <- data.frame(out1)
out1$Param <- row.names(out1)
preds <- subset(out1,grepl("Density",out1$Param))
environData$preds <- preds$mean
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

### end #################################################


