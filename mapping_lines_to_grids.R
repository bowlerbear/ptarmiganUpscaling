#run half of the HPC_linetransect_combinedModel_Buffer.R script

head(siteInfo_ArtsDaten)

#connect lines and grids - based on line centroids
lines_to_grids <- readRDS("data/lines_to_grids.rds")
subset(lines_to_grids,LinjeID %in% c("92","93"))#in the same grid

#do we have the grid data for all the lines in this analysis
siteInfo$LinjeID[!siteInfo$LinjeID %in% lines_to_grids$LinjeID]#yes!!!
siteInfo$grid <- lines_to_grids$grid[match(siteInfo$LinjeID,lines_to_grids$LinjeID)]
siteInfo$grid #all present

#and are all these grids in the siteInfo_Arts_daten
siteInfo$grid %in% siteInfo_ArtsDaten$grid
siteInfo$grid[!siteInfo$grid %in% siteInfo_ArtsDaten$grid]
#there are NAs...
#probably lines at the edge

#find an alternative grid that overlaps for those missing
siteInfo_ArtsDaten$siteIndex <- as.numeric(as.factor(siteInfo_ArtsDaten$grid))
siteInfo$siteIndex_All <- siteInfo_ArtsDaten$siteIndex[match(siteInfo$grid,siteInfo_ArtsDaten$grid)]
missing <- siteInfo$LinjeID[is.na(siteInfo$siteIndex_All)]
length(missing)#23

#get other overlapping grids within 5km
lineBuffers_to_grids <- readRDS("data/lineBuffers_to_grids_5km.rds")
lineBuffers_to_grids <- subset(lineBuffers_to_grids, LinjeID %in% missing)
lineBuffers_to_grids <- subset(lineBuffers_to_grids, grid %in% siteInfo_ArtsDaten$grid)
lineBuffers_to_grids5km <- subset(lineBuffers_to_grids,!duplicated(LinjeID))
nrow(lineBuffers_to_grids5km)#13
siteInfo$grid_5km <- lineBuffers_to_grids5km$grid[match(siteInfo$LinjeID,lineBuffers_to_grids5km$LinjeID)]

#other overlapping gruds within 15km _these are in islands in the north and in Sweden (see linetransect_Buffers.R)
missing <- missing[!missing %in% lineBuffers_to_grids5km$LinjeID]
lineBuffers_to_grids <- readRDS("data/lineBuffers_to_grids_15km.rds")
lineBuffers_to_grids <- subset(lineBuffers_to_grids, LinjeID %in% missing)
lineBuffers_to_grids <- subset(lineBuffers_to_grids, grid %in% siteInfo_ArtsDaten$grid)
lineBuffers_to_grids15km <- subset(lineBuffers_to_grids,!duplicated(LinjeID))
nrow(lineBuffers_to_grids15km)#10
siteInfo$grid_15km <- lineBuffers_to_grids15km$grid[match(siteInfo$LinjeID,lineBuffers_to_grids15km$LinjeID)]

#add to siteInfo object

#switch original grid to those with data
#first at 5km
siteInfo$grid[!is.na(siteInfo$grid_5km)] <- siteInfo$grid_5km[!is.na(siteInfo$grid_5km)]
#then at 15km
siteInfo$grid[!is.na(siteInfo$grid_15km)] <- siteInfo$grid_15km[!is.na(siteInfo$grid_15km)]

#now check we have data for all
siteInfo$grid[!siteInfo$grid %in% siteInfo_ArtsDaten$grid]
siteInfo$siteIndex_All <- siteInfo_ArtsDaten$siteIndex[match(siteInfo$grid,siteInfo_ArtsDaten$grid)]
siteInfo$LinjeID[is.na(siteInfo$siteIndex_All)]

saveRDS(siteInfo[,c("LinjeID","grid","siteIndex_All")], file="data/siteIndex_linetransects.rds")

