#retriving the bird data

library(rgbif)

#get willow ptarmigan data
name_backbone(name='Lagopus lagopus')

library(plyr)
wp<-ldply(2007:2017,function(x){
  out<-occ_data(taxonKey=2473421, country="NO",hasCoordinate=TRUE,year=x,basisOfRecord="HUMAN_OBSERVATION",hasGeospatialIssue=FALSE,limit=200000)
  out$data
})
write.table(wp,file="willow_ptarmigan_GBIF.txt",sep="\t")


#get data for all birds
name_backbone(name='Aves')
library(plyr)

myyears<-2010:2017
mymonths<-4:10

#for months 4 to 10
for(y in 1:length(myyears)){
  for(m in 1:length(mymonths)){
  out<-occ_data(classKey=212, country="NO",hasCoordinate=TRUE,year=myyears[y],
                month=mymonths[m],hasGeospatialIssue=FALSE,basisOfRecord="HUMAN_OBSERVATION",
                limit=200000)
  write.table(out$data,file=paste0("all_birds",myyears[y],"_",mymonths[m],"_GBIF.txt"),sep="\t",row.names=FALSE)
  }
}
