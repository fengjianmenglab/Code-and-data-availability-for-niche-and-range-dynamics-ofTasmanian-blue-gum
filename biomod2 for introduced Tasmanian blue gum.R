
library(sp)
library(raster)
library(parallel)
library(reshape)
library(ggplot2)
library(biomod2)
library(rgdal)


setwd("I:/cao/")
DataSpecies <- read.csv("csv_file/Introduce.csv")
head(DataSpecies)
myRespName <-'Eucalyptus.globulus.Labill.' 
myResp <- as.numeric(DataSpecies[,myRespName])
myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]




myExpl = stack( 
                raster( "current/bio1.tif"),
                raster( "current/bio2.tif"),
                
                raster( "current/bio4.tif"),
                
                raster( "current/bio9.tif"),
               
                raster( "current/bio12.tif"),
                
                raster( "current/bio17.tif"),
                raster( "current/bio18.tif"),
                raster( "current/bio19.tif"))
  
setwd("I:/cao/biomod2/Introduce/")
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp, 
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName,
                                     PA.nb.rep = 5,
                                     PA.nb.absences = 2090,
                                     PA.strategy = 'random')

write.csv(myBiomodData@coord,".PA_Introduce.CSV")

myBiomodData

plot(myBiomodData)
myBiomodOption <- BIOMOD_ModelingOptions()


getwd()
myBiomodModelOut <- BIOMOD_Modeling( 
  myBiomodData, 
  models = c("GLM", "GBM",  "CTA", "ANN", "FDA", "RF","MAXENT.Phillips"),
  models.options = myBiomodOption, 
  NbRunEval=5,
  DataSplit=70,
  Prevalence=0.5,
  VarImport=3, 
  models.eval.meth = c('TSS','ROC','KAPPA'), 
  SaveObj = TRUE, 
  rescal.all.models = TRUE,
  do.full.models = FALSE,
  modeling.id = paste(myRespName,"FirstModeling",sep=""))



myBiomodModelOut
myBiomodModelEval <- get_evaluations(myBiomodModelOut)
dimnames(myBiomodModelEval)
write.csv(myBiomodModelEval,file = "./myBiomodModelEval_Introduce.csv")
KAPPA=myBiomodModelEval["KAPPA","Testing.data",,,]
write.csv(KAPPA,"./KAPPA_Introduce.CSV")
ROC=myBiomodModelEval["ROC","Testing.data",,,]
write.csv(ROC,"./ROC_Introduce.CSV")
TSS=myBiomodModelEval["TSS","Testing.data",,,]
write.csv(TSS,"./TSS_Introduce.CSV")
varimportan=get_variables_importance(myBiomodModelOut)
write.csv(varimportan,"./varimportan_Introduce.csv")


myBiomodEM <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut, 
                                       chosen.models = 'all', 
                                       em.by='all', 
                                       eval.metric = c('TSS','ROC'), 
                                       eval.metric.quality.threshold = c(0.6,0.8),
                                       prob.mean = T, 
                                       prob.cv = T, 
                                       prob.ci = T,
                                       prob.ci.alpha = 0.05, 
                                       prob.median = T, 
                                       committee.averaging = T, 
                                       prob.mean.weight = T,
                                       prob.mean.weight.decay = 'proportional' )

myBiomodEM




EM_evaluations=get_evaluations(myBiomodEM)
write.csv(EM_evaluations,"/.EM_evaluations_Introduce.csv")
myBiomodProj <- BIOMOD_Projection( modeling.output = myBiomodModelOut, 
                                   new.env = myExpl, 
                                   proj.name = 'current', 
                                   selected.models = 'all', 
                                   binary.meth = 'TSS', 
                                   compress = 'xz', 
                                   build.clamping.mask = FALSE,
                                   output.format='.img')



myBiomodProj

plot(myBiomodProj)
myCurrentProj <- get_predictions(myBiomodProj)
myCurrentProj


myBiomodEF <- BIOMOD_EnsembleForecasting(projection.output = myBiomodProj,
                                         EM.output = myBiomodEM)


myBiomodEF
plot(myBiomodEF)








setwd("I:/cao/biomod2//Introduce/Eucalyptus.globulus.Labill./proj_current")
raster_<- stack('proj_current_Eucalyptus.globulus.Labill._ensemble.grd')
writeRaster(raster_, file="proj_current.asc", format="ascii", overwrite=TRUE, bylayer=TRUE, suffix=names(raster_))


setwd("I:/cao/biomod2/Introduce")