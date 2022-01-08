library(ade4)
library(ape)
library(gbm)
library(sp)
library(ecospat)
library(rgeos)
library(raster)
library(rworldmap)
library(parallel)
library(reshape)
library(ggplot2)
library(rgdal)

setwd("I:/cao/")######找到文件夹
DataSpecies1 <- read.csv("csv_file/All_data.csv") #####本土和入侵放一起的点
head(DataSpecies1)
myRespName1 <- 'Eucalyptus.globulus.Labill.' 
myResp1 <- as.numeric(DataSpecies1[,myRespName1])
myRespXY1 <- DataSpecies1[,c("X_WGS84","Y_WGS84")]


myExpl = stack(
               raster( "current/bio1.tif"),
               raster( "current/bio2.tif"),
               raster( "current/bio3.tif"),
               raster( "current/bio4.tif"),
               raster( "current/bio5.tif"),
               raster( "current/bio6.tif"),
               raster( "current/bio7.tif"),
               raster( "current/bio8.tif"),
               raster( "current/bio9.tif"),
               raster( "current/bio10.tif"),
               raster( "current/bio11.tif"),
               raster( "current/bio12.tif"),
               raster( "current/bio13.tif"),
               raster( "current/bio14.tif"),
               raster( "current/bio15.tif"),
               raster( "current/bio16.tif"),
               raster( "current/bio17.tif"),
               raster( "current/bio18.tif"),
               raster( "current/bio19.tif"))

xy.sp1<-subset(myRespXY1)#Bromus_erectus
env.sp1<-extract(myExpl,xy.sp1)
pca.cal <- dudi.pca(env.sp1, center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)

pca.cal$co
ecospat.plot.contrib(contrib=pca.cal$co,eigen=pca.cal$eig)######出图
