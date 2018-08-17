
#get the occupancy data:

#########################################################################################

#run ArtenDatenBank_missing_data_5km_sparta.R script

#get the bugs.data object

bugs.data_ArtsDaten <- bugs.data

########################################################################################

#get BUGS functions
source('C:/Users/diana.bowler/OneDrive - NINA/methods/models/bugsFunctions.R')

#need to specify initial values
library(reshape2)
zst <- acast(listlengthDF, siteIndex~yearIndex, value.var="y",fun=max,na.rm=T)
zst [is.infinite(zst)] <- 0
inits <- function(){list(z = zst)}

########################################################################################

#get the abundance data:

#run the linetransects.Rmd script unp until the 2nd bugs objects

bugs.data_LineTransects<-bugs.data 

#######################################################################################
#get effective strip width data

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
load("out1_linetransectModel_detection.RData")
siteInfo<-siteInfo_LineTransects
esw<-getBUGSFits(out1,param="predESW")

bugs.data_LineTransects$ESW.mean <- acast(esw,lineIndex~Year,value.var="mean")
bugs.data_LineTransects$ESW.sd <- acast(esw,lineIndex~Year,value.var="sd")

########################################################################################

#Add the line transect info to the bugs.data object

names(bugs.data_ArtsDaten)
names(bugs.data_LineTransects)

names(bugs.data_LineTransects)<-sapply(names(bugs.data_LineTransects),
                                       function(x)
                                         paste0(x,"LT"))

bugs.data<-c(bugs.data_ArtsDaten,bugs.data_LineTransects)

#check survey fractions:
temp <- bugs.data$TransectLengthLT/1000 * bugs.data$ESW.mean/1000 * 2
summary(temp)

########################################################################################

#run hurdle model - using adm and adm2 as random effects
#effort on detection probablity

params <- c("totAbund","int.d","beta.e","random.adm.sd",
            "random.adm2.sd","meanZ","realAbund")

#specify model structure
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "hurdleModel.txt", n.thin=nt,
             n.chains=nc, n.burnin=300,n.iter=1000,parallel=T)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_hurdle_random.RData")

hist(out1$mean$meanZ)
hist(out1$mean$realAbund)
print(out1,2)

##########################################################################################

#with random effects on occupancy and abundance:

bugs.data$Effort2 <- bugs.data$Effort^2

params <- c("int.p","beta.e","beta.e2",
            "random.adm.sd","random.adm2.sd","random.site.sd",
            "random.line.adm.sd","random.line.adm2.sd","random.line.site.sd")

#specify model structure
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "hurdleModel.txt", n.thin=nt,
             n.chains=nc, n.burnin=2000,n.iter=5000,parallel=T)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_hurdle_random2.RData")

#random.line.adm.sd, not really converging...
out2<-update(out1,"totAbund",n.iter=1000)

##########################################################################################

bugs.data$Effort2 <- bugs.data$Effort^2

params <- c("int.p","beta.e","beta.e2","beta.d","totAbund",
            "random.o.site","random.o.adm","random.o.adm2","random.o.year",
            "random.a.site","random.a.adm","random.a.adm2","random.a.year")

#add NAs
NuIndivs

#specify model structure
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "hurdleModel_v2.txt", n.thin=nt,
             n.chains=nc, n.burnin=2000,n.iter=5000,parallel=T)



##########################################################################################

#with variables on occupancy data

#specify model structure
bugs.data$occDM <- model.matrix(~ bugs.data$tree_line_position + 
                                  bugs.data$tree_line_position2+
                                  bugs.data$bio1 + 
                                  bugs.data$open)[,-1]

bugs.data$n.covs <- ncol(bugs.data$occDM)

params <- c("totAbund","int.d","beta.e","beta","random.adm.sd","meanZ","realAbund")

#specify model structure
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "hurdleModel.txt", n.thin=nt,
             n.chains=nc, n.burnin=300,n.iter=1000,parallel=T)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_hurdle_varableOcc.RData")

########################################################################################

#with variables on occupancy data and detection

#specify model structure
bugs.data$occDM <- model.matrix(~ bugs.data$tree_line_position + 
                                  bugs.data$tree_line_position2+
                                  bugs.data$bio1 + 
                                  bugs.data$open)[,-1]

bugs.data$n.covs <- ncol(bugs.data$occDM)

params <- c("totAbund","int.d","beta.e","beta","beta.det","random.adm.sd","meanZ","realAbund")

#specify model structure
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "hurdleModel.txt", n.thin=nt,
             n.chains=nc, n.burnin=300,n.iter=1000,parallel=T)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_hurdle_varable.RData")

#######################################################################################

#including effort as a squared term
#along with covariates on detection and occupancy

bugs.data$Effort2 <- bugs.data$Effort^2
summary(bugs.data$Effort2)

bugs.data$occDM <- model.matrix(~ bugs.data$tree_line_position + 
                                  bugs.data$tree_line_position2+
                                  bugs.data$bio1 + 
                                  bugs.data$open)[,-1]

bugs.data$n.covs <- ncol(bugs.data$occDM)

params <- c("totAbund","int.d","beta.e","beta.e2","beta","beta.det",
            "random.adm.sd","meanZ","realAbund")

#specify model structure
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "hurdleModel.txt", n.thin=nt,
             n.chains=nc, n.burnin=300,n.iter=1000,parallel=T)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_hurdle_varable_effort2.RData")

######################################################################################

#including effort as a squared term
#along with covariates on detection and occupancy
#full list of covariates

bugs.data$Effort2 <- bugs.data$Effort^2
summary(bugs.data$Effort2)

bugs.data$occDM <- model.matrix(~ bugs.data$tree_line_position + 
                                  bugs.data$tree_line_position2+
                                  bugs.data$bio1 + 
                                  bugs.data$bio1_2 +
                                  bugs.data$prefopen +
                                  bugs.data$prefclosed +
                                  bugs.data$open)[,-1]

bugs.data$n.covs <- ncol(bugs.data$occDM)

params <- c("totAbund","int.d","beta.e","beta.e2","beta","beta.det",
            "random.adm.sd","meanZ","realAbund")

#specify model structure
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "hurdleModel.txt", n.thin=nt,
             n.chains=nc, n.burnin=300,n.iter=1000,parallel=T)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_hurdle_varable_exd_effort2.RData")

#######################################################################################

#including effort as a squared term
#along with covariates on detection and occupancy
#full list of covariates
#add in random effects 1 and 2 on detection as well 

bugs.data$Effort2 <- bugs.data$Effort^2
summary(bugs.data$Effort2)

bugs.data$occDM <- model.matrix(~ bugs.data$tree_line_position + 
                                  bugs.data$tree_line_position2+
                                  bugs.data$bio1 + 
                                  bugs.data$bio1_2 +
                                  bugs.data$prefopen +
                                  bugs.data$prefclosed +
                                  bugs.data$open)[,-1]

bugs.data$n.covs <- ncol(bugs.data$occDM)

params <- c("totAbund","int.d","beta.e","beta.e2","beta","beta.det",
            "random.adm.sd","random.adm2.sd",
            "random.det.adm.sd","random.det.adm2.sd",
            "meanZ","realAbund")

#specify model structure
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "hurdleModel.txt", n.thin=nt,
             n.chains=nc, n.burnin=300,n.iter=1000,parallel=T)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_hurdle_varable_exd_random_effort2.RData")

######################################################################################

#including effort as a squared term
#along with covariates on occupancy
#full list of covariates
#add in random effects 1 and 2 on detection as well 

bugs.data$Effort2 <- bugs.data$Effort^2
summary(bugs.data$Effort2)

bugs.data$occDM <- model.matrix(~ bugs.data$tree_line_position + 
                                  bugs.data$tree_line_position2+
                                  bugs.data$bio1 + 
                                  bugs.data$bio1_2 +
                                  bugs.data$prefopen +
                                  bugs.data$prefclosed +
                                  bugs.data$open)[,-1]

bugs.data$n.covs <- ncol(bugs.data$occDM)

params <- c("totAbund","int.d","beta.e","beta.e2","beta",
            "random.adm.sd","random.adm2.sd",
            "random.det.adm.sd","random.det.adm2.sd",
            "meanZ","realAbund")

#specify model structure
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "hurdleModel.txt", n.thin=nt,
             n.chains=nc, n.burnin=300,n.iter=1000,parallel=T)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_hurdle_varable_exd_occ_random_effort2.RData")

######################################################################################

#including effort as a squared term
#along with covariates on occupancy
#full list of covariates
#add in random effects 1 and 2 on detection as well 
#covariates also on abundance

bugs.data$Effort2 <- bugs.data$Effort^2
summary(bugs.data$Effort2)

#which line transect grids match the occDM
matchGrids<-sapply(siteInfo_LineTransects$grid,function(x)which(listlengthDF_SiteCovs$grid==x))

bugs.data$occDM <- model.matrix(~ bugs.data$tree_line_position + 
                                  bugs.data$tree_line_position2+
                                  bugs.data$bio1 + 
                                  bugs.data$bio1_2 +
                                  bugs.data$prefopen +
                                  bugs.data$prefclosed +
                                  bugs.data$open)[,-1]

bugs.data$occDM_Occupancy<-bugs.data$occDM[matchGrids,]

bugs.data$n.covs <- ncol(bugs.data$occDM)

params <- c("totAbund","int.d","beta.e","beta.e2","beta","beta.occupancy",
            "random.adm.sd","random.adm2.sd",
            "random.det.adm.sd","random.det.adm2.sd",
            "meanZ","realAbund")

#specify model structure
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "hurdleModel.txt", n.thin=nt,
             n.chains=nc, n.burnin=300,n.iter=1000,parallel=T)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_hurdle_varable_exd_occ_abund_random_effort2.RData")

#######################################################################################

#including effort as a squared term
#along with covariates on occupancy
#full list of covariates
#add in random effects 1 and 2 on detection as well as ext covariates
#covariates also on abundance

bugs.data$Effort2 <- bugs.data$Effort^2
summary(bugs.data$Effort2)

#which line transect grids match the occDM
matchGrids<-sapply(siteInfo_LineTransects$grid,function(x)which(listlengthDF_SiteCovs$grid==x))

bugs.data$occDM <- model.matrix(~ bugs.data$tree_line_position + 
                                  bugs.data$tree_line_position2+
                                  bugs.data$bio1 + 
                                  bugs.data$bio1_2 +
                                  bugs.data$prefopen +
                                  bugs.data$prefclosed +
                                  bugs.data$open)[,-1]

bugs.data$occDM_Occupancy<-bugs.data$occDM[matchGrids,]

bugs.data$n.covs <- ncol(bugs.data$occDM)

params <- c("totAbund","int.d","beta.e","beta.e2","beta","beta.det",
            "beta.occupancy",
            "random.adm.sd","random.adm2.sd",
            "random.det.adm.sd","random.det.adm2.sd",
            "meanZ","realAbund")

#specify model structure
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "hurdleModel.txt", n.thin=nt,
             n.chains=nc, n.burnin=300,n.iter=1000,parallel=T)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_hurdle_varable_exd_occ_abund_det_effort2.RData")

#######################################################################################

#look at predicted relative abundances
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
load("out1_hurdle_varable_exd_occ_abund_random_effort2.RData")
siteInfo<-siteInfo_ArtsDaten
siteInfo$fits<-out1$mean$realAbund
mygrid[]<-0
mygrid[as.numeric(siteInfo$grid)]<-siteInfo$fits
par(mfrow=c(1,2))
plot(mygrid,main="abundance")

#look at predicted meanZ
siteInfo$fits<-out1$mean$meanZ
mygrid[]<-0
mygrid[as.numeric(siteInfo$grid)]<-siteInfo$fits
plot(mygrid,main="occupancy")

########################################################################################