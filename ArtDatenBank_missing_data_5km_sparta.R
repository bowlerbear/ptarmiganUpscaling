####################################################################################################

setwd("C:/Users/db40fysa/Dropbox/ptarmigan Upscaling")

source('formattingArtDatenBank_missing_data.R')
#remember only a subset of the data is currently used
#need to fix "outside" adms still

###################################################################################################

oldlistlengthDF<-listlengthDF

nrow(listlengthDF_SiteCovs)#15530

#sample a random sample of the grid?
#sampledSites<-sample(unique(listlengthDF$grid),5000)
#listlengthDF<-subset(listlengthDF,grid%in%sampledSites)

###################################################################################################

#order data by site and year:
listlengthDF$siteIndex <- as.numeric(factor(listlengthDF$grid))
listlengthDF$yearIndex <- as.numeric(factor(listlengthDF$year))
listlengthDF <- arrange(listlengthDF,siteIndex,yearIndex)

#extract site data
listlengthDF_SiteCovs <- subset(listlengthDF,!duplicated(grid))

#adm indicies
listlengthDF$admN <- as.numeric(factor(listlengthDF$adm))
listlengthDF$admN2 <- as.numeric(factor(listlengthDF$adm2))
siteInfo <- unique(listlengthDF[,c("siteIndex","admN","grid","admN2","adm","adm2")])
siteInfo_ArtsDaten <- siteInfo

#add random siye/year effect
listlengthDF$siteyearIndex <- as.numeric(factor(interaction(listlengthDF$siteIndex,listlengthDF$yearIndex)))

###################################################################################################

#Add sampling days per year/site
samplingDays<-ddply(listlengthDF,.(siteyearIndex),summarise,
                   nuDays=length(unique(month[!is.na(y)],day[!is.na(y)])))

listlengthDF$samplingDays<-samplingDays$nuDays[match(listlengthDF$siteyearIndex,samplingDays$siteyearIndex)]

listlengthDF$Effort <- listlengthDF$L2/listlengthDF$nuSpecies
listlengthDF$Effort[is.na(listlengthDF$Effort)] <- 0

###################################################################################################

#fit as glm with explanatory variables - binary response

tempDF <- subset(listlengthDF,!duplicated(grid))

#get y
occupancyGrid <- ddply(listlengthDF,.(grid),summarise,species=max(y,na.rm=T))
tempDF <- merge(tempDF,occupancyGrid,by="grid",all.x=T)
tempDF$species[is.infinite(tempDF$species)]<-NA

#remove nonsurvyed grids
tempDF <- subset(tempDF,!is.na(species))
mean(tempDF$species)

summary(glm(species~alpine_habitat2,data=tempDF))
summary(glm(species~bio1,data=tempDF))
summary(glm(species~bio5,data=tempDF))
summary(glm(species~tree_line_position + I(tree_line_position^2),data=tempDF))
summary(glm(species~Top,data=tempDF))

#all together
summary(glm(species~bio1 + Open + tree_line_position + I(tree_line_position^2),
            data=tempDF,family="binomial"))
#all significant...
glm1<-glm(species ~ scale(bio1) + Open + scale(tree_line_position) + 
            scale(I(tree_line_position^2)),family="binomial",data=tempDF)
summary(glm1)

#plot predictions
tempDF$x <- mygridMaskDF$x[match(tempDF$grid,mygridMaskDF$layer)]
tempDF$y <- mygridMaskDF$y[match(tempDF$grid,mygridMaskDF$layer)]

tempDF$fits<-predict(glm1,type="response",newdata=tempDF)
ggplot(tempDF)+
  geom_point(aes(x,y,colour=fits),shape=15,size=rel(1))+
  scale_colour_gradient(low="steelblue",high="red")
  #geom_point(data=subset(varDF2,species==1),aes(x,y),alpha=0.1)

tempDF$fits.se<-predict(glm1,type="response",newdata=tempDF,se.fit=TRUE)$se.fit
ggplot(tempDF)+
  geom_point(aes(x,y,colour=fits.se),shape=15,size=rel(1))+
  scale_colour_gradient(low="steelblue",high="red")

###############################################################################

#glm - number of visits as the response

visitGrid <- ddply(listlengthDF,.(grid),summarise,
                   nuVisits = length(y[!is.na(y)]),
                   nuObs = sum(y[!is.na(y)]))

tempDF <- merge(tempDF,visitGrid,by="grid")

#remove nonsurvyed data
tempDF <- subset(tempDF,nuVisits>0)
table(tempDF$nuObs)

#all together
glm1 <- glm(cbind(nuObs,nuVisits-nuObs)~
              scale(bio1) + scale(Open) + 
              scale(tree_line_position) + scale(I(tree_line_position^2)),
            data=tempDF,family="binomial")

#put probability scale on the same as last one

#plot predictions
tempDF$fits<-predict(glm1,type="response",newdata=tempDF)
ggplot(tempDF)+
  geom_point(aes(x,y,colour=fits),shape=15,size=rel(1))+
  scale_colour_gradient(low="steelblue",high="red")

tempDF$fits.se<-predict(glm1,type="response",newdata=tempDF,se.fit=TRUE)$se.fit
ggplot(tempDF)+
  geom_point(aes(x,y,colour=fits.se),shape=15,size=rel(1))+
  scale_colour_gradient(low="steelblue",high="red")

qplot(fits,fits.se,data=tempDF)

##############################################################################

tempDF$fitsL <- predict(glm1,type="link",newdata=tempDF)
tempDF$fitsL.se <- predict(glm1,type="link",newdata=tempDF,se.fit=TRUE)$se.fit

#sensitivity analysis
allPreds <- tempDF[,c("grid","x","y","nuVisits",
                      "fits","fits.se","fitsL","fitsL.se")]

#total predicted number of grid cells
sum(allPreds$fits)
sum(allPreds$fits.se)

#for each grid
#all current occupancies for all other grids
#take new samples for focus grid
#get new predicted total sample size

myGridSurveys <- function(allPreds,mygrid){
  
  currentFits <- allPreds$fits[allPreds$grid!=mygrid]
  
  totalNu <- as.numeric()
  require(boot)
  for(i in 1:100){
  newSample <- invlogit(rnorm(1,allPreds$fitsL[allPreds$grid==mygrid],
                     allPreds$fitsL.se[allPreds$grid==mygrid]))
  
  #new predicted total number of grids
  totalNu[i] <- sum(currentFits,newSample)
  }
  
  data.frame(grid=mygrid,mean(totalNu),sd(totalNu))
}

#apply function to each grid
sensitivityOutput <- ldply(allPreds$grid,function(x){
  myGridSurveys(allPreds,mygrid=x)
})

allPreds <- merge(allPreds,sensitivityOutput,by="grid")
ggplot(allPreds)+
  geom_point(aes(x,y,colour=sd.totalNu.),shape=15,size=rel(1))+
  scale_colour_gradient2(low="grey",mid="pink",high="purple",
                         midpoint=0.004)

qplot(fits,sd.totalNu.,data=allPreds)
qplot(sd.totalNu.,nuVisits,data=allPreds)

summary(lm(sqrt(sd.totalNu.) ~ scale(fits) + scale(nuVisits),data=allPreds))

#simulate change in binomial variance
#with different values of n

#variation only partly caused by number of visits
#pull out the predicted effect of another visit?
summary(lm(sqrt(sd.totalNu.) ~ scale(fits)*nuVisits,data=allPreds))

allPreds$fitsB <- ifelse(allPreds$fits < median(allPreds$fits),"L","H")
summary(lm(sqrt(sd.totalNu.) ~ fitsB + fitsB:nuVisits,data=allPreds))
#effect of visit more important in higher density areas


n <- 1:100

p1 <- 0.1
p2 <- 0.5
p3 <- 0.9

varBin <- function(n,p){
  n*p*(1-p)
}

varBin(n=n,p=p1)

df1 <- data.frame(n,varBin(n,p1))
df1$P <- p1
names(df1)<-c("n","var","p")

df2 <- data.frame(n,varBin(n,p2))
df2$P <- p2
names(df2)<-c("n","var","p")

df3 <- data.frame(n,varBin(n,p3))
df3$P <- p3
names(df3)<-c("n","var","p")

allDF <- rbind(df1,df2,df3)
qplot(n,var,data=allDF,facets=~factor(p))

###influence measures#######################################################################

#https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/influence.measures

#Bar Plot of Cook's distance to detect observations that strongly influence fitted values of the mode
#cooks.distance

tempDF$cooks <- cooks.distance(glm1)
ggplot(tempDF)+
  geom_point(aes(x,y,colour=cooks),shape=15,size=rel(1))+
  scale_colour_gradient(low="steelblue",high="red")

#DFBETA measures the difference in each parameter estimate with and without the influential point
#dfbeta - for the overall mean
temp <- dfbeta(glm1)

tempDF$dfbeta <- temp[,1]
ggplot(tempDF)+
  geom_point(aes(x,y,colour=abs(dfbeta)),shape=15,size=rel(1))+
  scale_colour_gradient2(low="grey",mid="pink",high="purple")

###################################################################################################

#Boosted regression tree

library(dismo)
library(gbm)

str(varDF)

brt1 <- gbm.step(data=tempDF, gbm.x = 4:15, gbm.y = 18,family = "bernoulli")

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
                  siteyearIndex = listlengthDF$siteyearIndex,
                  L = listlengthDF$L,
                  nuSpecies = listlengthDF$nuSpecies,
                  Effort = listlengthDF$Effort,
                  #nuSamplingDays = scale(listlengthDF$samplingDays),
                  y = listlengthDF$y,
                  #add an adm effect
                  adm = siteInfo$admN,
                  det.adm = listlengthDF$admN,
                  n.adm = length(unique(siteInfo$admN)),
                  adm2 = siteInfo$admN2,
                  det.adm2 = listlengthDF$admN2,
                  n.adm2 = length(unique(siteInfo$admN2)),
                  #habitat = ifelse(listlengthDF_SiteCovs$habitat=="Forest",1,0),
                  bio1 = scale(listlengthDF_SiteCovs$bio1),
                  bio1_2 = scale(listlengthDF_SiteCovs$bio1^2),
                  bio6 = listlengthDF_SiteCovs$bio6,
                  bio5 = listlengthDF_SiteCovs$bio5,
                  bio5_2 = listlengthDF_SiteCovs$bio5^2,
                  forest = listlengthDF_SiteCovs$Forest,
                  open = scale(listlengthDF_SiteCovs$Open),
                  prefopen = scale(listlengthDF_SiteCovs$PrefOpen),
                  prefclosed = scale(listlengthDF_SiteCovs$PrefClosed),
                  top = log(listlengthDF_SiteCovs$Top+1),
                  alpine_habitat1 = listlengthDF_SiteCovs$alpine_habitat1,
                  alpine_habitat2 = log(listlengthDF_SiteCovs$alpine_habitat2+1),
                  alpine_habitat3 = log(listlengthDF_SiteCovs$alpine_habitat3+1),
                  alpine_habitat4 = log(listlengthDF_SiteCovs$alpine_habitat4+1),
                  tree_line_position = scale(listlengthDF_SiteCovs$tree_line_position),
                  tree_line_position2 = scale(listlengthDF_SiteCovs$tree_line_position^2))

bugs.data_ArtsDaten <- bugs.data

#alpine_habitat:
#1= Open lowland, 
#2 = Low alpine zone, 
#3 = intermediate alpine zone, 
#4 = high alpine zone 

########################################################################################

#get BUGS functions
source('bugsFunctions.R')

#need to specify initial values
library(reshape2)
zst <- acast(listlengthDF, siteIndex~yearIndex, value.var="y",fun=max,na.rm=T)
zst [is.infinite(zst)] <- 0
inits <- function(){list(z = zst)}

##########################################################################################

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

############################################################################################

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

############################################################################################

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

##########################################################################################

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

##########################################################################################

#plot occupancy model
out2<-update(out1,parameters.to.save=c("muZ","z"),n.iter=100)
plotZ(out2,param="z")
plotZ(out2,param="muZ")
plotZerror(out2)

##########################################################################################

#idenfity areas with high probability of occurence but no observations 





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
out1 <- jags(bugs.data, inits=inits, params, "models/BUGS_sparta_JAGAM.txt", n.thin=nt,
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
q2<-qplot(access,propY,data=propSuccess)#positive??
q3<-qplot(bio1,propY,data=propSuccess)#negative relationship
plot_grid(q1,q2,q3)

########################################################################################

