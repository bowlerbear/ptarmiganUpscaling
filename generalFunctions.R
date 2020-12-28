#general functions:

#get function to extract raster into for these grids:

getEnvironData<-function(myraster,mygridTemp){
  require(maptools)
  require(plyr)
  
  #crop raster to Norway extent
  rasterCRS<-crs(myraster)
  NorwayB<-spTransform(Norway,rasterCRS)
  myraster<-crop(myraster,extent(NorwayB))
  
  #convert raster into points and convert
  myrasterDF<-as.data.frame(myraster,xy=T)
  names(myrasterDF)[3]<-"myraster"
  coordinates(myrasterDF)<-c("x","y")
  proj4string(myrasterDF)<-rasterCRS
  myrasterDF<-spTransform(myrasterDF,crs(equalM))
  
  #get general grid
  mygrid<-gridTemp
  projection(mygrid)<- equalM
  
  #get mean myraster values per grid cell
  mygrid[]<-1:ncell(mygrid)
  variable<-extract(mygrid,myrasterDF)
  myrasterDF<-data.frame(myrasterDF@data)
  myrasterDF$grid<-variable
  myrasterDF<-ddply(myrasterDF,.(grid),summarise,myraster=mean(myraster,na.rm=T))
  return(myrasterDF)
}

plotBUGSData<-function(myvar){
  temp<-listlengthDF
  temp<-subset(temp,siteIndex %in% bugs.data$site)
  temp$variable<-as.numeric(bugs.data[myvar][[1]])[match(temp$siteIndex,bugs.data$site)]
  temp$variablePA<-sapply(temp$variable,function(x)ifelse(is.na(x),0,1))
  temp<-subset(temp,!duplicated(siteIndex))
  mygrid<-gridTemp
  par(mfrow=c(1,2))
  mygrid[]<-0
  mygrid[temp$grid]<-temp$variable
  plot(mygrid)
  mygrid[]<-0
  mygrid[temp$grid]<-temp$variablePA
  plot(mygrid)
}

plotZ<-function(model,param="z"){
  zSummary<-data.frame(model$summary[grepl(param,row.names(model$summary)),])
  zSummary$ParamNu <- as.character(sub(".*\\[([^][]+)].*", "\\1", row.names(zSummary)))
  zSummary$Site<-sapply(zSummary$ParamNu,function(x)strsplit(x,",")[[1]][1])
  zSummary$Year<-sapply(zSummary$ParamNu,function(x)strsplit(x,",")[[1]][2])
  zSummary$grid<-listlengthDF$grid[match(zSummary$Site,listlengthDF$siteIndex)]
  #zSummary_year<-subset(zSummary,Year==11)
  
  #predicted occupancy across all years
  zSummary_year<-ddply(zSummary,.(grid),summarise,prop=mean(mean))
  #plot mean prop occupancy
  par(mfrow=c(1,1))
  mygrid[]<-0
  mygrid[zSummary_year$grid]<-zSummary_year$prop
  plot(mygrid)
  plot(Norway,add=T)
}

plotZerror<-function(model){
  zSummary<-data.frame(model$summary[grepl("z",row.names(model$summary)),])
  zSummary$ParamNu <- as.character(sub(".*\\[([^][]+)].*", "\\1", row.names(zSummary)))
  zSummary$Site<-sapply(zSummary$ParamNu,function(x)strsplit(x,",")[[1]][1])
  zSummary$Year<-sapply(zSummary$ParamNu,function(x)strsplit(x,",")[[1]][2])
  zSummary$grid<-listlengthDF$grid[match(zSummary$Site,listlengthDF$siteIndex)]
  
  #predicted occupancy across all years
  zSummary_year<-ddply(zSummary,.(grid),summarise,prop=mean(mean),prop_sd=mean(sd),prop_cov=prop_sd/prop)
  
  #plot sd
  par(mfrow=c(2,1))
  mygrid[]<-0
  mygrid[zSummary_year$grid]<-zSummary_year$prop_sd
  plot(mygrid)
  plot(Norway,add=T)
  mygrid[]<-0
  #plot cov
  mygrid[zSummary_year$grid]<-zSummary_year$prop_cov
  plot(mygrid)
  plot(Norway,add=T)
}

getParam<-function(model,param="z"){
  zSummary<-data.frame(model$summary[grepl(param,row.names(model$summary)),])
  zSummary$ParamNu <- as.character(sub(".*\\[([^][]+)].*", "\\1", row.names(zSummary)))
  zSummary$Site<-sapply(zSummary$ParamNu,function(x)strsplit(x,",")[[1]][1])
  zSummary$Year<-sapply(zSummary$ParamNu,function(x)strsplit(x,",")[[1]][2])
  return(zSummary)
}

logit<-function(x) log(x/(1-x))

invlogit<- function(x) exp(x)/(1+exp(x))

getDecimalPlaces<-function(x){
  nchar(gsub("(.*\\.)|([0]*$)", "", as.character(x)))
}

#get the mode Name for each mtbq
Mode <- function(x) {
  ux <- unique(x)
  x <- x[!is.na(x)]
  if(length(x)>0){
    return(ux[which.max(tabulate(match(x, ux)))])
  }
  else{
    return(NA)
  }
}
