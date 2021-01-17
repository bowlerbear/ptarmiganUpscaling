myfolder <- "data"
myfolder <- "/data/idiv_ess/ptarmiganUpscaling" #HPC

### get data ################################################

#get the occupancy data:
bugs.data_ArtsDaten <- readRDS(paste(myfolder,"bugs.data_ArtsDaten.rds",sep="/"))
zst <- readRDS("zst_ArtsDaten.rds")
inits <- function(){list(z = zst)}

#get the abundance data:
bugs.data_LineTransects <- readRDS(paste(myfolder,"bugs.data_linetransects.rds",sep="/")) 

### relabel ##################################################

names(bugs.data_LineTransects)<-sapply(names(bugs.data_LineTransects),
                                       function(x)
                                         paste0(x,"LT"))

bugs.data<-c(bugs.data_ArtsDaten,bugs.data_LineTransects)
names(bugs.data)

### define covariates ##########################################

#just for subset with line surveys
bugs.data$occDM <- model.matrix(~ visitedData$tree_line +
                                  I(visitedData$tree_line^2) +
                                  visitedData$bio1 +
                                  visitedData$bio5 +
                                  visitedData$elevation +
                                  I(visitedData$elevation^2) +
                                  visitedData$Open +
                                  visitedData$Top)[,-1]

bugs.data$n.covs <- ncol(bugs.data$occDM)

#for whole grid
bugs.data$predDM <- model.matrix(~ environData$tree_line +
                                   I(environData$tree_line^2) +
                                   environData$bio1 +
                                   environData$bio5 +
                                   I(environData$elevation^2) +
                                   environData$Open +
                                   environData$Top)[,-1]

bugs.data$npreds <- nrow(environData)

### fit model ##################################################

library(rjags)
library(jagsUI)

params <- c("totAbund","int.d","mean.psi","mean.p",
            "meanZ","realAbund")

modelfile <- paste(myfolder,"hurdleModel_v2.txt",sep="/")

out1 <- jags(bugs.data, 
             inits=inits, 
             params, 
             modelfile, 
             n.thin=20,
             n.chains=20, 
             n.burnin=20000,
             n.iter=40000,
             parallel=T)

save(out1$summary,file="out1_hurdle.RData")

#######################################################################################