setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
load("out1_hurdle_random3.RData")
load("siteInfo_ArtsDaten.RData")
siteInfo<-siteInfo_ArtsDaten

#print(out1$summary,2)

#get random effects
randomSummary <- data.frame(out1$summary[grepl("random",row.names(out1$summary)),])
randomSummary$Param <- row.names(randomSummary)

#separate the random effects:

#get number
randomSummary$ParamNu <- as.numeric(sub(".*\\[([^][]+)].*", "\\1", randomSummary$Param))
randomSummary$ParamNu <- as.numeric(sub(".*\\[([^][]+)].*", "\\1", randomSummary$Param))

#get o or a
randomSummary$Type <- sapply(randomSummary$Param,function(x)
  strsplit(as.character(x),"\\.")[[1]][2])

unique(randomSummary$Type)

#get spatial level
randomSummary$Level <- sapply(randomSummary$Param,function(x)
  strsplit(as.character(x),"\\.")[[1]][3])
unique(randomSummary$Level)
randomSummary$Level <- sapply(randomSummary$Level,function(x)
  strsplit(as.character(x),"\\[")[[1]][1])
unique(randomSummary$Level)

####################################################################################################

#look at correlation of random effects between a and o
library(reshape2)
randomCast <- dcast(randomSummary,ParamNu+Level~Type,value.var="mean")
library(ggplot2)
ggplot(randomCast)+
  geom_point(aes(x=a,y=o))+
  facet_wrap(~Level,scales="free")

#####################################################################################################

#get estimated population size for each regional and look at its sd

#plot the random effects in space!!

#get info on administrative names for the grid

#make both spatial objects in the same crs
library(rgdal)
equalM<-"+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
NorwayADM<-readOGR(dsn="C:/Users/diana.bowler/OneDrive - NINA/Alpine/NOR_adm",layer="NOR_adm2")

#get the adm level random effects:

#for occupancy
library(ggplot2)
randomSummary_adm <- subset(randomSummary,Level=="adm" & Type=="o")
randomSummary_adm$adm <- siteInfo$adm[match(randomSummary_adm$ParamNu,siteInfo$admN)]
AG <- fortify(NorwayADM,region="NAME_1")
AG$re_mean <- randomSummary_adm$mean[match(AG$id,randomSummary_adm$adm)]
AG$re_sd <- randomSummary_adm$sd[match(AG$id,randomSummary_adm$adm)]

g1<-ggplot(AG) + 
  aes(long,lat,group=group,fill=re_mean) + 
  geom_polygon() +
  geom_path(color="black")+
  scale_fill_gradient2(low="blue",high="red",mid="white",midpoint=median(AG$re_mean))+
  theme_minimal()+
  ggtitle("occupancy")

g2<-ggplot(AG) + 
  aes(long,lat,group=group,fill=re_sd) + 
  geom_polygon() +
  geom_path(color="black")+
  scale_fill_gradient2(low="blue",high="red",mid="white",midpoint=median(AG$re_sd))+
  theme_minimal()+
  ggtitle("occupancy")

#for abundance
randomSummary_adm <- subset(randomSummary,Level=="adm" & Type=="a")
randomSummary_adm$adm <- siteInfo$adm[match(randomSummary_adm$ParamNu,siteInfo$admN)]
AG <- fortify(NorwayADM,region="NAME_2")
AG$re_mean <- randomSummary_adm$mean[match(AG$id,randomSummary_adm$adm)]
AG$re_sd <- randomSummary_adm$sd[match(AG$id,randomSummary_adm$adm)]

g3<-ggplot(AG) + 
  aes(long,lat,group=group,fill=re_mean) + 
  geom_polygon() +
  geom_path(color="black")+
  scale_fill_gradient2(low="blue",high="red",mid="white",midpoint=median(AG$re_mean))+
  theme_minimal()+
  ggtitle("abundance")

g4<-ggplot(AG) + 
  aes(long,lat,group=group,fill=re_sd) + 
  geom_polygon() +
  geom_path(color="black")+
  scale_fill_gradient2(low="blue",high="red",mid="white",midpoint=median(AG$re_sd))+
  theme_minimal()+
  ggtitle("abundance")

#plot together
library(cowplot)
plot_grid(g1,g2,g3,g4,ncol=2)

#get the adm2 level random effects

randomSummary_adm <- subset(randomSummary,Level=="adm2" & Type=="o")
randomSummary_adm$adm <- siteInfo$adm2[match(randomSummary_adm$ParamNu,siteInfo$admN2)]
AG <- fortify(NorwayADM,region="NAME_2")
AG$re_mean <- randomSummary_adm$mean[match(AG$id,randomSummary_adm$adm)]
AG$re_sd <- randomSummary_adm$sd[match(AG$id,randomSummary_adm$adm)]

g1<-ggplot(AG) + 
  aes(long,lat,group=group,fill=re_mean) + 
  geom_polygon() +
  scale_fill_gradient2(low="blue",high="red",mid="white",midpoint=median(AG$re_mean,na.rm=T))+
  theme_minimal()+
  ggtitle("occupancy")

g2<-ggplot(AG) + 
  aes(long,lat,group=group,fill=re_sd) + 
  geom_polygon() +
  scale_fill_gradient2(low="blue",high="red",mid="white",midpoint=median(AG$re_sd,na.rm=T))+
  theme_minimal()+
  ggtitle("occupancy")

#for abundance
randomSummary_adm <- subset(randomSummary,Level=="adm2" & Type=="a")
randomSummary_adm$adm <- siteInfo$adm2[match(randomSummary_adm$ParamNu,siteInfo$admN2)]
AG <- fortify(NorwayADM,region="NAME_2")
AG$re_mean <- randomSummary_adm$mean[match(AG$id,randomSummary_adm$adm)]
AG$re_sd <- randomSummary_adm$sd[match(AG$id,randomSummary_adm$adm)]

g3<-ggplot(AG) + 
  aes(long,lat,group=group,fill=re_mean) + 
  geom_polygon() +
  scale_fill_gradient2(low="blue",high="red",mid="white",midpoint=median(AG$re_mean,na.rm=T))+
  theme_minimal()+
  ggtitle("abundance")

g4<-ggplot(AG) + 
  aes(long,lat,group=group,fill=re_sd) + 
  geom_polygon() +
  scale_fill_gradient2(low="blue",high="red",mid="white",midpoint=median(AG$re_sd,na.rm=T))+
  theme_minimal()+
  ggtitle("abundance")

#plot together
library(cowplot)
plot_grid(g1,g2,g3,g4,ncol=2)

######################################################################################################

#use simulations to explore the sensitivity of the total population size estimate to 
#reducing the variation in sd








#######################################################################################################

#how does effort relate to sd?





#######################################################################################################
