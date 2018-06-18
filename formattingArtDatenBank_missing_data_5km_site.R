####################################################################################################
source('C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/formattingArtDatenBank_missing_data.R')
###################################################################################################




###################################################################################################
#fit as glm with explanatory variables

#add to the dataset, the coordinates of the grid
gridDF<-as.data.frame(gridTemp,xy=T)
load("varDF_allEnvironData_5km.RData")
varDF<-merge(varDF,gridDF,by.x="grid",by.y="layer",all.x=T)

#get y
occupancyGrid<-ddply(listlengthDF,.(grid),summarise,species=max(y,na.rm=T))
varDF<-merge(varDF,occupancyGrid,by="grid",all.x=T)
varDF$species[is.infinite(varDF$species)]<-NA

#remove nonsurvyed data
varDF2<-subset(varDF,!is.na(species))
mean(varDF2$species)

summary(glm(species~alpine_habitat2,data=varDF2))
summary(glm(species~bio1,data=varDF2))
summary(glm(species~bio5,data=varDF2))
summary(glm(species~tree_line_position + I(tree_line_position^2),data=varDF2))
summary(glm(species~Top,data=varDF2))

#all together
summary(glm(species~bio1 + Open + tree_line_position + I(tree_line_position^2),data=varDF))
#all significant...
glm1<-glm(species~bio1 + Open + tree_line_position + I(tree_line_position^2),data=varDF2)


varDF2$fits<-predict(glm1,type="response",newdata=varDF)

#need to get x and y coords below
library(ggplot2)
ggplot(varDF2)+
  geom_point(aes(x,y,colour=fits),shape=15,size=rel(1),alpha=0.5)+
  scale_colour_gradient(low="steelblue",high="red")+
  geom_point(data=subset(varDF2,species==1),aes(x,y))

###################################################################################################

#Boosted regression tree

library(dismo)
library(gbm)

str(varDF)

brt1 <- gbm.step(data=varDF2, gbm.x = 4:15, gbm.y = 18,family = "bernoulli")

summary(brt1)
#                                var    rel.inf
#                                tree_line_position tree_line_position 49.5344301
#                                bio1                             bio1 13.0001113
#                                Open                             Open  8.2227746
#                                Top                               Top  8.0722788
#                                bio6                             bio6  7.0045898
#                                bio5                             bio5  3.8682442
#                                Forest                         Forest  3.5676963
#                                alpine_habitat3       alpine_habitat3  3.0023880
#                                alpine_habitat1       alpine_habitat1  2.3505843
#                                alpine_habitat4       alpine_habitat4  0.6835894
#                                Agriculture               Agriculture  0.3634804
#                                alpine_habitat2       alpine_habitat2  0.3298327
                                


#plot main effects
gbm.plot(brt1, n.plots=12, write.title = TRUE)
gbm.plot.fits(brt1)

#interactions?
find.int <- gbm.interactions(brt1)
find.int$interactions#none!
find.int$rank.list

####################################################################################################

#restrict the data??

set.seed(3)

listlengthDF<-ddply(listlengthDF,.(siteIndex,yearIndex),function(x){
  x$visitNu<-sample(1:nrow(x))
  x<-x[order(x$visitNu),]
})
oldlistlengthDF<-listlengthDF
  
summary(listlengthDF$visitNu[!is.na(listlengthDF$y)])
table(listlengthDF$visitNu)
#subset to 13 (median)
listlengthDF<-subset(listlengthDF,visitNu<=13)

###################################################################################################

#for BUGS

bugs.data <- list(nsite = length(unique(listlengthDF$siteIndex)),
                  nyear = length(unique(listlengthDF$yearIndex)),
                  nvisit = nrow(listlengthDF),
                  grid = listlengthDF$grid,
                  site = listlengthDF$siteIndex,
                  year = listlengthDF$yearIndex,
                  L = listlengthDF$L,
                  nuSpecies = listlengthDF$nuSpecies,
                  y = listlengthDF$y,
                  #add an adm effect
                  adm = siteInfo$admN,
                  n.adm = length(unique(siteInfo$admN)),
                  adm2 = siteInfo$admN2,
                  n.adm2 = length(unique(siteInfo$admN2)),
                  habitat = ifelse(listlengthDF_SiteCovs$habitat=="Forest",1,0),
                  bio1 = scale(listlengthDF_SiteCovs$bio1),
                  bio1_2 = listlengthDF_SiteCovs$bio1^2,
                  bio6 = listlengthDF_SiteCovs$bio6,
                  bio5 = listlengthDF_SiteCovs$bio5,
                  bio5_2 = listlengthDF_SiteCovs$bio5^2,
                  forest = listlengthDF_SiteCovs$Forest,
                  open = scale(listlengthDF_SiteCovs$Open),
                  top = log(listlengthDF_SiteCovs$Top+1),
                  alpine_habitat1 = listlengthDF_SiteCovs$alpine_habitat1,
                  alpine_habitat2 = log(listlengthDF_SiteCovs$alpine_habitat2+1),
                  alpine_habitat3 = log(listlengthDF_SiteCovs$alpine_habitat3+1),
                  alpine_habitat4 = log(listlengthDF_SiteCovs$alpine_habitat4+1),
                  tree_line_position = scale(listlengthDF_SiteCovs$tree_line_position),
                  tree_line_position2 = scale(listlengthDF_SiteCovs$tree_line_position^2))

#alpine_habitat:
#1= Open lowland, 
#2 = Low alpine zone, 
#3 = intermediate alpine zone, 
#4 = high alpine zone 

#need to specify initial values
library(reshape2)
zst <- acast(listlengthDF, siteIndex~yearIndex, value.var="y",fun=max,na.rm=T)
zst [is.infinite(zst)] <- 0
inits <- function(){list(z = zst)}

#get BUGS functions
source('C:/Users/diana.bowler/OneDrive - NINA/methods/models/bugsFunctions.R')

##########################################################################################

#costant detection proability
params <- c("mean.p","mean.psi")

#specify model structure
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "BUGS_sparta_constant.txt", n.thin=10,
             n.chains=3, n.burnin=300,n.iter=1000,parallel=T)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_OM_random_missing_5km_basic0.RData")


############################################################################################
#run model with random effects on adm

#specify parameters to monitor
params <- c("psi.fs","sd.y","sd.s","mean.p","mean.psi","random.adm.sd","random.adm")

#specify model structure
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "BUGS_sparta_random_missing.txt", n.thin=5,
             n.chains=3, n.burnin=300,n.iter=1000,parallel=T)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_OM_random_missing_5km_basic.RData")

############################################################################################

#run model with random effects on adm and year and site

#specify parameters to monitor
params <- c("psi.fs","sd.y","sd.s","mean.p","mean.psi","random.adm.sd","random.adm")

#specify model structure
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "BUGS_sparta_random_missing.txt", n.thin=10,
             n.chains=3, n.burnin=600,n.iter=2000,parallel=T)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_OM_random_missing_5km_basic2.RData")

#include year on the detection probability

#run model with random effects on adm

#specify parameters to monitor
params <- c("int","psi.fs","dtype.p","mu.lp","a","random.adm.sd","alpha.p")

#specify model structure
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "BUGS_sparta_random_missing.txt", n.thin=10,
             n.chains=3, n.burnin=300,n.iter=1000,parallel=T)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_OM_random_missing_5km1.RData")

#include other covariates on detection

#specify parameters to monitor
params <- c("int","psi.fs","dtype.p","mu.lp","beta","a","d.tree","d.tree2")

#specify model structure
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "BUGS_sparta_random_missing.txt", n.thin=10,
             n.chains=3, n.burnin=300,n.iter=1000,parallel=T)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_OM_random_missing_5km_det.RData")
#d.type effect is negative!!

#just elevation effect on detection (wasnt scaled)

#specify parameters to monitor
params <- c("int","psi.fs","mu.lp","beta","a","int.d","d.tree","d.tree2")

#specify model structure
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "BUGS_sparta_random_missing.txt", n.thin=10,
             n.chains=3, n.burnin=300,n.iter=1000,parallel=T)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_OM_random_missing_5km_det2.RData")

#elevation effect and L2/L

#specify parameters to monitor
params <- c("int","psi.fs","mu.lp","beta","a","int.d","dtype.p","d.tree","d.tree2","random.adm.sd")

#specify model structure
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "BUGS_sparta_random_missing.txt", n.thin=10,
             n.chains=3, n.burnin=300,n.iter=1000,parallel=T)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_OM_random_missing_5km_det3.RData")

########################################################################################

#run model including explanatory variables on observation

#specify model structure
bugs.data$occDM <- model.matrix(~ bugs.data$tree_line_position + bugs.data$tree_line_position2)[,-1]
bugs.data$n.covs <- ncol(bugs.data$occDM)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "BUGS_sparta_variables_missing.txt", n.thin=10,
             n.chains=3, n.burnin=300,n.iter=1000,parallel=T)

print(out1,2)
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_OM_variables_missing_5km.RData")
colnames(bugs.data$occDM)

##########################################################################################

#plot occupancy model
out2<-update(out1,parameters.to.save=c("muZ","z"),n.iter=100)
plotZ(out2,param="z")
plotZ(out2,param="muZ")
plotZerror(out2)

##########################################################################################
 
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
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "BUGS_sparta_JAGAM.txt", n.thin=nt,
             n.chains=3, n.burnin=10000,n.iter=50000)

traceplot(out1)
print(out1,2)

#########################################################################################

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
q2<-qplot(access,propY,data=propSuccess)#positivee??
q3<-qplot(bio1,propY,data=propSuccess)#negative relationship
plot_grid(q1,q2,q3)

########################################################################################

