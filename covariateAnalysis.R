#covariateAnalysis

setwd("C:/Users/diana.bowler/OneDrive - NINA/Alpine/ptarmiganUpscaling/model-outputs")
load("out1_hurdle_varable_exd_occ_abund_det_effort2.RData")

#plot relationship between detectability and effort
possSampleRecs<-seq(1,7,by=0.1)
predDetection <- as.numeric()
for(i in 1:length(possSampleRecs)){
  predDetection[i] <- out1$mean$int.p + possSampleRecs[i]*out1$mean$beta.e + possSampleRecs[i]^2*out1$mean$beta.e2
}
plot(possSampleRecs,predDetection,type="line")

#plot environmental covariates
covariateSummary <- data.frame(out1$summary[grepl("beta",row.names(out1$summary)),])
covariateSummary$Param <- row.names(covariateSummary)
covariateSummary<-subset(covariateSummary,!Param %in% c("beta.e","beta.e2"))

#get number
covariateSummary$ParamNu <- as.numeric(sub(".*\\[([^][]+)].*", "\\1", covariateSummary$Param))

#get parameter type
covariateSummary$Type <- sapply(as.character(covariateSummary$Param),function(x)
                              strsplit(x,"\\[")[[1]][1])

#match with covariate names
covariateOrder<-c("tree line","quad tree line","mean temp", "quad mean temp","preferred open",
                  "preferred forest","open")
covariateSummary$covariate <- sapply(covariateSummary$ParamNu,function(x)
  covariateOrder[as.numeric(x)])
covariateSummary$covariate <- as.factor(covariateSummary$covariate)
covariateSummary$covariate <- factor(covariateSummary$covariate, levels=rev(c("mean temp","quad mean temp",
                                                                          "tree line","quad tree line",
                                                                          "preferred open","preferred forest",
                                                                          "open")))

#rename levels of Type
covariateSummary$Type<-as.factor(covariateSummary$Type)
levels(covariateSummary$Type)<-c("occupancy","detection","abundance")

#plot data
library(ggplot2)
ggplot(covariateSummary)+
  geom_crossbar(aes(x=covariate,y=mean,ymin=X2.5.,ymax=X97.5.,colour=Type))+
  coord_flip()+
  geom_hline(yintercept=0,colour="red",linetype="dashed")+
  theme_bw()+
  facet_grid(~Type)+
  theme(legend.position="none")

