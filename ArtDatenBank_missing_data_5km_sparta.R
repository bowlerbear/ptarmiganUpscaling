setwd("C:/Users/db40fysa/Dropbox/ptarmigan Upscaling")

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
listlengthDF <- readRDS("data/listlength_iDiv.rds")

#subset to May to September
listlengthDF <- subset(listlengthDF,is.na(month)|(month > 4 & month <10))

### subset ##########################################################

#subset to focal grids and those with environ covariate data

focusGrids <- readRDS("data/focusGrids.rds")
varDF <- readRDS("data/varDF_allEnvironData_5km_idiv.rds")
listlengthDF <- subset(listlengthDF,grid %in% focusGrids)
listlengthDF <- subset(listlengthDF,grid %in% varDF$grid)

### effort ######################################################################

#we will treat it as a categorical variable
table(listlengthDF$L)
listlengthDF$singleton <- ifelse(listlengthDF$L==1,1,0)

### Absences ###################################################################

listlengthDF$L[is.na(listlengthDF$L)] <- 1 #set to nominal effort
#listlengthDF$y[listlengthDF$L==0] <- 0

### GLM - binary ###############################################################

#fit as glm with explanatory variables - binary response

#get y per grid
occupancyGrid <- ddply(listlengthDF,.(grid),summarise,species=max(y,na.rm=T))
occupancyGrid <- subset(occupancyGrid,!is.infinite(species))
mean(occupancyGrid$species)

#get environ data
tempDF <- merge(varDF,occupancyGrid,by="grid",all.y=T)

#simple glms

#temp vars
summary(glm(species~scale(bio1),data=tempDF,family="binomial"))#large
summary(glm(species~scale(bio5),data=tempDF,family="binomial"))
summary(glm(species~scale(bio6),data=tempDF,family="binomial"))#large

#tree line or elevation
summary(glm(species~scale(elevation),data=tempDF,family="binomial"))#large 
summary(glm(species~scale(tree_line),data=tempDF,family="binomial"))
summary(glm(species~scale(tree_line_position),data=tempDF,family="binomial"))#best
qplot(elevation,tree_line_position,data=tempDF)

#habitat
summary(glm(species~MountainBirchForest,data=tempDF,family="binomial"))#small positive
summary(glm(species~Bog,data=tempDF,family="binomial"))#large positive
summary(glm(species~Forest,data=tempDF,family="binomial"))#negative large
summary(glm(species~ODF,data=tempDF,family="binomial"))#large positive
summary(glm(species~Meadows,data=tempDF,family="binomial"))#large positive
summary(glm(species~OSF,data=tempDF,family="binomial"))#large positive
summary(glm(species~SnowBeds,data=tempDF,family="binomial"))#no effect
summary(glm(species~Mire,data=tempDF,family="binomial"))#small positive
summary(glm(species~Human,data=tempDF,family="binomial"))#large negative


#all together
glm1 <- glm(species~bio1 + MountainBirchForest+ Bog + Forest + ODF + Meadows + OSF + Mire + Human + SnowBeds + tree_line_position + I(tree_line_position^2),
            data=tempDF,family="binomial")
#meadows, bog and mountainbirch forest and ODF
DescTools::PseudoR2(glm1)

#large range of variables to test
library(MuMIn)
options(na.action = "na.fail")

glm1<-glm(species ~ scale(bio1) + 
            scale(bio5) + #mostly
            scale(bio6) + #always
            scale(y) + #always
            scale(distCoast) + #mostly
            scale(MountainBirchForest) + #always
            scale(Bog) + #always 
            scale(Forest) + 
            scale(Mire) +
            scale(Human) +
            scale(ODF) + #always
            scale(Meadows) + #always
            scale(OSF) + 
            scale(SnowBeds) +
            scale(tree_line_position) + 
            scale(I(tree_line_position^2))+ #always
            scale(elevation)+ #mostly
            scale(elevation^2),
          family="binomial",data=tempDF)
summary(glm1)
DescTools::PseudoR2(glm1)

dd <- dredge(glm1)
subset(dd, delta < 2)

### plot predictions #############################################

tempDF$x <- myGridDF$x[match(tempDF$grid,myGridDF$layer)]
tempDF$y <- myGridDF$y[match(tempDF$grid,myGridDF$layer)]

tempDF$fits<-predict(glm1,type="response",newdata=tempDF)
ggplot(tempDF)+
  geom_point(aes(x,y,colour=fits),shape=15,size=rel(1))+
  scale_colour_viridis_c()

tempDF$fits.se<-predict(glm1,type="response",newdata=tempDF,se.fit=TRUE)$se.fit
ggplot(tempDF)+
  geom_point(aes(x,y,colour=fits.se),shape=15,size=rel(1))+
  scale_colour_gradient(low="steelblue",high="red")

### brt #########################################################################

#Boosted regression tree
library(dismo)
library(gbm)
brt1 <- gbm.step(data=tempDF, gbm.x = c(2:14,20:21,29:30), gbm.y = 33, family = "bernoulli")

summary(brt1)
#                                   var    rel.inf
# tree_line_position   tree_line_position 31.1882397#complex
# Bog                                 Bog 21.1275937
# ODF                                 ODF 13.8494494
# y                                     y  6.2408870
# tree_line                     tree_line  5.0624717
# bio1                               bio1  3.7947544
# Meadows                         Meadows  3.5494534
# bio6                               bio6  3.0857581
# MountainBirchForest MountainBirchForest  2.7713510
# elevation                     elevation  1.9931910
# OSF                                 OSF  1.8197470
# distCoast                     distCoast  1.3255056
# SnowBeds                       SnowBeds  1.0535460

#plot main effects
gbm.plot(brt1, n.plots=12, write.title = TRUE)
gbm.plot.fits(brt1)
#non-linear plot for tree line position and bio1

#### quadratic effects #########################################################

glm1<-glm(species ~ tree_line_position +
            I(tree_line_position^2)+
            Bog +
            ODF + 
            OSF +
            y + 
            bio1 +
            I(bio1^2)+
            Meadows +
            MountainBirchForest +
            OSF,
          family="binomial",data=tempDF)
summary(glm1)
DescTools::PseudoR2(glm1)

### ignore from here on ########################################################

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

#add on observations 
occupancyGrid <- ddply(listlengthDF,.(grid),summarise,
                       species = max(y,na.rm=T))
occupancyGrid$species[is.infinite(occupancyGrid$species)] <- NA
table(occupancyGrid$species)
siteInfo$species <- occupancyGrid$species[match(siteInfo$grid,occupancyGrid$grid)]

siteInfo$x <- myGridDF$x[match(siteInfo$grid,myGridDF$layer)]
siteInfo$y <- myGridDF$y[match(siteInfo$grid,myGridDF$layer)]

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

### constant p #####################################################

#costant detection proability
params <- c("mean.p","mean.psi")

#specify model structure
out1 <- jags(bugs.data, inits=inits, params, "models/BUGS_sparta_constant.txt", n.thin=10,
             n.chains=3, n.burnin=300,n.iter=1000,parallel=T)

print(out1,2)

save(out1,file="model-outputs/out1_OM_random_missing_5km_basic0.RData")

#test sample of 5000
            mean   sd     2.5%      50%    97.5% overlap0 f Rhat n.eff
mean.p       0.15  0.0     0.15     0.15     0.16    FALSE 1 1.00   210
mean.psi     0.17  0.0     0.16     0.17     0.18    FALSE 1 1.01   175
deviance 10716.33 62.6 10600.99 10716.29 10840.94    FALSE 1 1.02   119

### RE###############################################################

#run model with random effects on adm

#specify parameters to monitor
params <- c("mean.p","mean.psi","random.adm.sd","random.adm2.sd")

#specify model structure
out1 <- jags(bugs.data, inits=inits, params, "models/BUGS_sparta_random_missing.txt", n.thin=5,
             n.chains=3, n.burnin=1000,n.iter=2500,parallel=T)
print(out1,2)

#test sample of 1000
                    mean    sd      2.5%       50%     97.5% overlap0 f Rhat n.eff
mean.p              0.09  0.01      0.08      0.09      0.10    FALSE 1 1.07    35
mean.psi            0.04  0.03      0.01      0.03      0.12    FALSE 1 1.57     7
random.adm.sd       3.71  1.16      1.93      3.55      6.56    FALSE 1 1.40     9
random.adm2.sd      3.93  0.53      3.04      3.88      5.11    FALSE 1 1.73     6
deviance       290193.98 31.06 290132.65 290193.92 290254.42    FALSE 1 1.14    19

### EXP VARS ########################################################

#run model including explanatory variables on observation

#specify model structure
bugs.data$occDM <- model.matrix(~ bugs.data$tree_line_position + bugs.data$tree_line_position2+
                                  bugs.data$bio1 + bugs.data$open)[,-1]
bugs.data$n.covs <- ncol(bugs.data$occDM)

params <- c("mean.p","mean.psi","beta")

out1 <- jags(bugs.data, inits=inits, params, "models/BUGS_sparta_variables_missing.txt", n.thin=10,
             n.chains=3, n.burnin=1000,n.iter=3000,parallel=T)

print(out1,2)
colnames(bugs.data$occDM)

#test sample of 1000

              mean    sd      2.5%       50%     97.5% overlap0    f Rhat n.eff
mean.p        0.08  0.01      0.07      0.08      0.09    FALSE 1.00 1.01   163
mean.psi      0.36  0.08      0.22      0.35      0.52    FALSE 1.00 1.24    13
beta[1]       0.87  0.61     -0.33      0.90      2.00     TRUE 0.92 1.24    13
beta[2]      -3.17  0.64     -4.18     -3.27     -1.88    FALSE 1.00 1.90     5
beta[3]      -1.40  0.22     -1.85     -1.40     -1.00    FALSE 1.00 1.04    64
beta[4]       1.51  0.32      0.88      1.49      2.17    FALSE 1.00 1.02   278
deviance 290302.36 29.10 290249.40 290301.02 290360.29    FALSE 1.00 1.03    79

### EFFORT ##########################################################

#run model including explanatory variable on observation and detection
#on detection also put - number of sampling days and records per date

#specify model structure
bugs.data$occDM <- model.matrix(~ bugs.data$tree_line_position + bugs.data$tree_line_position2+
                                  bugs.data$bio1 + bugs.data$open)[,-1]
bugs.data$n.covs <- ncol(bugs.data$occDM)

params <- c("mean.p","mean.psi","beta","beta.det","beta.sd","beta.e")

out1 <- jags(bugs.data, inits=inits, params, "models/BUGS_sparta_variables_missing.txt", n.thin=3,
             n.chains=3, n.burnin=2500,n.iter=5000,parallel=T)

save(out1,file="model-outputs/out1_OM_random_missing_5km.RData")

print(out1,2)
                  mean     sd       2.5%        50%      97.5% overlap0    f Rhat n.eff
mean.p            0.06   0.00       0.06       0.06       0.07    FALSE 1.00 1.05    47
mean.psi          0.12   0.00       0.12       0.12       0.13    FALSE 1.00 1.11    26
beta[1]           0.68   0.06       0.56       0.68       0.78    FALSE 1.00 1.00   480
beta[2]          -1.96   0.07      -2.12      -1.96      -1.83    FALSE 1.00 1.20    17
beta[3]          -0.58   0.02      -0.63      -0.58      -0.54    FALSE 1.00 1.02   104
beta[4]           0.39   0.03       0.33       0.39       0.45    FALSE 1.00 1.01   158
beta.det[1]       0.08   0.05      -0.02       0.08       0.19     TRUE 0.94 1.01   284
beta.det[2]      -0.49   0.08      -0.64      -0.49      -0.34    FALSE 1.00 1.09    29
beta.det[3]      -0.11   0.02      -0.16      -0.11      -0.07    FALSE 1.00 1.01   153
beta.det[4]       0.20   0.02       0.15       0.20       0.25    FALSE 1.00 1.01   351
beta.sd          -0.67   0.04      -0.74      -0.67      -0.60    FALSE 1.00 1.01   329 ## negative effect!!
beta.e            0.66   0.02       0.62       0.66       0.69    FALSE 1.00 1.00  2502
deviance    8384730.56 250.63 8384259.12 8384725.70 8385231.56    FALSE 1.00 1.01   371

### PLOT ###########################################################

#plot occupancy model
out2<-update(out1,parameters.to.save=c("muZ","z"),n.iter=100)
plotZ(out2,param="z")
plotZ(out2,param="muZ")
plotZerror(out2)

### JAGAM ###########################################################

#also using jagam

#fit as gam
library(mgcv)
gam1 <- gam(species~ 1 + s(x, y,k=10), 
                    data=varDF, 
                    family="binomial")

#plot it
varDF$fits<-gam1$fitted.values
library(ggplot2)
ggplot(varDF)+
  geom_point(aes(x,y,colour=fits),shape=15,size=rel(4))+
  scale_colour_gradient(low="blue",high="red")+
  geom_point(data=subset(varDF,species==1),aes(x,y))


#Use the JAGAM function to get the BUGS code for the spline
#in the BUGS model
jags.ready <- jagam(species ~ 1 + s(x, y,k=10), 
                    data = varDF, 
                    family="binomial", 
                    sp.prior="log.uniform",
                    file="jagam.txt")

#get the data bits we need from jags data
bugs.data$X = jags.ready$jags.data$X
bugs.data$S1 = jags.ready$jags.data$S1
bugs.data$zero = jags.ready$jags.data$zero

#specify parameters to monitor
params <- c("dtype.p","mu.lp","rho")

#run model
out1 <- jags(bugs.data, inits=inits, params, "models/BUGS_sparta_JAGAM.txt", n.thin=nt,
             n.chains=3, n.burnin=10000,n.iter=50000)

traceplot(out1)
print(out1,2)

### DETECTION PROB ######################################################################################

#explore detection prob by hand

#per year - for each line, what fraction of surveys had a ptarmigan in it
#get rid of NAs

listlengthDF_NAfree<-subset(listlengthDF,!is.na(visit))
propSuccess<-ddply(listlengthDF_NAfree,.(grid,adm,adm2,YearCollected),
                   summarise,
                   propY=ifelse(length(y)>1,mean(y),NA),#want to calculate it only when there are repeat visits
                   meanL=mean(L),
                   meanL2=mean(L2),
                   meanL3=mean(L3))
propSuccess<-subset(propSuccess,!is.na(propY))

hist(propSuccess$propY)#very all or nothing
mean(propSuccess$propY[propSuccess$propY>0])#0.86

#look at relationship
library(ggplot2)
library(cowplot)
q1<-qplot(meanL,propY,data=propSuccess)#negative relationship
q2<-qplot(meanL2,propY,data=propSuccess)#negative relationship
q3<-qplot(meanL3,propY,data=propSuccess)#more positive
plot_grid(q1,q2,q3)

#plot in relation to site covariates
propSuccess<-merge(propSuccess,listlengthDF_SiteCovs,by=c("grid","adm","adm2","site"))
q1<-qplot(tree_line_position,propY,data=propSuccess)#quadratic relationship
q2<-qplot(access,propY,data=propSuccess)#positive??
q3<-qplot(bio1,propY,data=propSuccess)#negative relationship
plot_grid(q1,q2,q3)



