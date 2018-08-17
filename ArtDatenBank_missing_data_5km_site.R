####################################################################################################
source('C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/formattingArtDatenBank_missing_data.R')
###################################################################################################

oldlistlengthDF<-listlengthDF

#collapse data to site-level rather than site/year level (too many NAs and we arent interested in year..)
listlengthDF <- subset(listlengthDF,!is.na(y))
listlengthDF <- merge(listlengthDF,listlengthDF_SiteCovs[,-c(2:10,12)],#remove variables that are not site-specific
                    by=names(listlengthDF_SiteCovs)[-c(2:10,12)],all=T)
#for missing year data give random values
listlengthDF$yearIndex[is.na(listlengthDF$y)]<-
  sample(listlengthDF$yearIndex[!is.na(listlengthDF$y)],length(listlengthDF$yearIndex[is.na(listlengthDF$y)]))

#how much missing data is there?
length(listlengthDF$y[listlengthDF$y==0])#786743
length(listlengthDF$y[listlengthDF$y==1])#20541
length(listlengthDF$y[is.na(listlengthDF$y)])#9775

#where do we have not data???
siteSummary<-ddply(listlengthDF,.(grid,adm),summarise,
                   Notsurveyed=ifelse(all(is.na(y)),1,0))
mygrid[]<-0
mygrid[siteSummary$grid]<-siteSummary$Notsurveyed
plot(mygrid,main="notsurveyed")          
nrow(subset(siteSummary,Notsurveyed==1))#9775
nrow(subset(siteSummary,Notsurveyed==1 & adm!="outside"))#8483

#where do we have data
siteSummary<-ddply(listlengthDF,.(grid,adm),summarise,
                   Surveyed=ifelse(any(!is.na(y)),1,0))
mygrid[]<-0
mygrid[siteSummary$grid]<-siteSummary$Surveyed
plot(mygrid,main="surveyed")
nrow(subset(siteSummary,Surveyed==1))#18007
nrow(subset(siteSummary,Surveyed==1 & adm!="outside"))#15527

#out<-unique(listlengthDF$grid[!is.na(listlengthDF$y)])
#out2<-unique(listlengthDF$grid[is.na(listlengthDF$y)])
#sum(out%in%out2)

#####################################################################
#restrict the data?? No
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

###################################################################"

#sample a random sample of the grid

sampledSites<-sample(unique(listlengthDF$grid),1000)
listlengthDF<-subset(listlengthDF,grid%in%sampledSites)

#####################################################################

#order data by site and year
listlengthDF$siteIndex<-as.numeric(factor(listlengthDF$grid))
listlengthDF$yearIndex<-as.numeric(factor(listlengthDF$year))
listlengthDF<-arrange(listlengthDF,siteIndex,yearIndex)

#extract site data
listlengthDF_SiteCovs<-subset(listlengthDF,!duplicated(grid))

#adm data
listlengthDF$admN<-as.numeric(factor(listlengthDF$adm))
listlengthDF$admN2<-as.numeric(factor(listlengthDF$adm2))
siteInfo<-unique(listlengthDF[,c("siteIndex","admN","grid","admN2")])

#add random siye/year effect
listlengthDF$siteyearIndex<-as.numeric(factor(interaction(listlengthDF$siteIndex,listlengthDF$yearIndex)))

###################################################################################################

#estimate detection probability by hand:

detProb<-ddply(listlengthDF,.(grid),summarise,nuSurveys=length(y[!is.na(y)]),
                                              nuObs = length(y[y==1&!is.na(y)]))
detProb<-subset(detProb,nuObs>0)
hist(detProb$nuObs/detProb$nuSurveys)
mean(detProb$nuObs/detProb$nuSurveys)#0.3

#plot the distribution of detections
mygrid[]<-0
mygrid[detProb$grid]<-detProb$nuObs/detProb$nuSurveys
plot(mygrid,main="Detection probs")

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

###################################################################################################

#for BUGS

bugs.data <- list(nsite = length(unique(listlengthDF$siteIndex)),
                  nyear = length(unique(listlengthDF$yearIndex)),
                  nvisit = nrow(listlengthDF),
                  grid = listlengthDF$grid,
                  site = listlengthDF$siteIndex,
                  year = listlengthDF$yearIndex,
                  nsiteyear = length(unique(listlengthDF$siteyearIndex)),
                  siteyear = listlengthDF$siteyearIndex,
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
zst <- as.numeric(tapply(listlengthDF$y, listlengthDF$siteIndex,max,na.rm=T))
zst [is.infinite(zst)] <- 0
inits <- function(){list(z = zst)}

#get BUGS functions
source('C:/Users/diana.bowler/OneDrive - NINA/methods/models/bugsFunctions.R')

##########################################################################################

#costant detection proability
params <- c("mean.p","mean.psi")

#specify model structure
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "BUGS_sparta_constant_site.txt", n.thin=2,
             n.chains=3, n.burnin=300,n.iter=1000,parallel=T)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_OM_missing_5km_basic0_site.RData")

############################################################################################

#run model with random effects on adm

#specify parameters to monitor
params <- c("mean.p","mean.psi","random.adm.sd","random.adm2.sd")

#specify model structure
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "BUGS_sparta_constant_site.txt", n.thin=5,
             n.chains=3, n.burnin=600,n.iter=2000,parallel=T)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_OM_random_missing_5km_basic_test.RData")

############################################################################################

#include other covariates on detection

#number of sampling days
#listlength per date

#specify parameters to monitor
params <- c("int","psi.fs","dtype.p","mu.lp","beta","a","d.tree","d.tree2")

#specify model structure
setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/models")
out1 <- jags(bugs.data, inits=inits, params, "BUGS_sparta_random_missing.txt", n.thin=10,
             n.chains=3, n.burnin=300,n.iter=1000,parallel=T)

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
save(out1,file="out1_OM_random_missing_5km_det.RData")

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
 
