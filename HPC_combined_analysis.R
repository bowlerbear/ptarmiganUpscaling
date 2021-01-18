#myfolder <- "data"
myfolder <- "/data/idiv_ess/ptarmiganUpscaling" #HPC

### get data ################################################

#get the occupancy data:
bugs.data_ArtsDaten <- readRDS(paste(myfolder,"bugs.data_ArtsDaten.rds",sep="/"))
zst <- readRDS(paste(myfolder,"zst_ArtsDaten.rds",sep="/"))
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
# bugs.data$occDM <- model.matrix(~ visitedData$tree_line +
#                                   I(visitedData$tree_line^2) +
#                                   visitedData$bio1 +
#                                   visitedData$bio5 +
#                                   visitedData$elevation +
#                                   I(visitedData$elevation^2) +
#                                   visitedData$Open +
#                                   visitedData$Top)[,-1]
# 
# bugs.data$n.covs <- ncol(bugs.data$occDM)
# 
# #for whole grid
# bugs.data$predDM <- model.matrix(~ environData$tree_line +
#                                    I(environData$tree_line^2) +
#                                    environData$bio1 +
#                                    environData$bio5 +
#                                    I(environData$elevation^2) +
#                                    environData$Open +
#                                    environData$Top)[,-1]
# 
# bugs.data$npreds <- nrow(environData)

### fit model ##################################################

library(rjags)
library(jagsUI)

n.cores = as.integer(Sys.getenv("NSLOTS", "1")) 

params <- c("totalAbund","realAbund")

modelfile <- paste(myfolder,"simpleCombinedModel.txt",sep="/")

out1 <- jags(bugs.data, 
             inits=inits, 
             params, 
             modelfile, 
             n.thin=20,
             n.chains=n.cores, 
             n.burnin=20000,
             n.iter=40000,
             parallel=T)

saveRDS(out1$summary,file="outSummary_simpleCombinedModel.rds")

#######################################################################################
#https://stats.stackexchange.com/questions/234887/predicting-future-values-with-hurdle-poisson-model
#https://stats.stackexchange.com/questions/15967/how-can-i-set-up-a-zero-inflated-poisson-in-jags