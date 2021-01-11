### occu models ##############################################

### plotting map #################################################

#from full model
out1 <- readRDS("model-outputs/out_occModel_upscaling.rds")
print(out1$summary,3)
siteInfo$preds <- out1$mean$grid.muZ
mygrid[] <- NA
mygrid[siteInfo$grid] <- siteInfo$preds
plot(mygrid)

#from model summary
out1 <- readRDS("model-outputs/outSummary_occModel_upscaling.rds")
out1 <- data.frame(out1)
out1$Param <- row.names(out1)
preds <- subset(out1,grepl("grid.z",out1$Param))
siteInfo$preds <- preds$mean
mygrid[] <- NA
mygrid[siteInfo$grid] <- siteInfo$preds
plot(mygrid)

# using tmap package
crs(mygrid) <- equalM
#tmaptools::palette_explorer()
library(tmap)
tm_shape(mygrid)+
  tm_raster(title="Occupancy prob",palette="YlGnBu")

### plotting coefficients #######################################

betas <- subset(out1,grepl("beta",out1$Param))
betas <- betas[1:8,]
betas$variables <- c("tree_line_position", "tree_line_position2","bio1","bio1_2", "bio6","elevation","prefopen", "open")

ggplot(betas)+
  geom_crossbar(aes(x=variables,y=mean,
                    ymin=X2.5.,ymax=X97.5.))+
  coord_flip()+
  theme_bw()+
  ylab("effect size on occupancy")+
  geom_hline(yintercept=0,color="red",
             linetype="dashed")

### distance models ##########################################







