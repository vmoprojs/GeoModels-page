rm(list=ls())
require(devtools)
#install_github("vmoprojs/GeoModels")
library(GeoModels)
library(mapproj)
library(globe)
library(fields)
library(sphereplot)
set.seed(1891)
################################################################
###
### Example 8. Simulation of a Gaussian RF 
###  with a Wend0 correlation on the planet earth
###
###############################################################



radius=1
NN=700 ## number of location sites on the earth 
coords=pointsphere(NN,c(-180,180),c(-90,90),c(radius,radius))[,1:2]


globeearth(eye=place("newyorkcity"))
globepoints(loc=coords,pch=20,cex=0.4)



globeearth(eye=place("everest"))
globepoints(loc=coords,pch=20,cex=0.4)


corrmodel = "Smoke"    ## correlation model and parameters
scale=radius*0.2
smooth=0.5
sill=1
nugget=0
mean=0



###  simulation using geodesic distance
param=list(mean=mean, sill=sill, nugget=nugget, smooth=smooth,scale=scale)
data = GeoSim(coordx=coords, corrmodel=corrmodel, param=param,
             distance="Geod",radius=radius)$data



fixed<-list(nugget=nugget)
start<-list(mean=mean,scale=scale,sill=sill,smooth=smooth)


optimizer="nlminb"
I=Inf
upper<-list(mean=I,scale=I,sill=I,smooth=I)
lower<-list(mean=-I,scale=0,sill=0,smooth=0)


## maximum likelihood estimation  using geodesic distance: can be time consuming...
fit_geo_ml <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel, 
              likelihood="Full",type="Standard",
              optimizer=optimizer,lower=lower,upper=upper,
              start=start,fixed=fixed,distance="Geod",radius=radius)


geod_ds=rdist.earth(coords,miles=F,R=radius)
max_geod=max(geod_ds)
## maximum pairwise likelihood estimation using geodesic distance
fit_geo_pl <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel, 
             likelihood="Conditional",type="Pairwise", maxdist=max_geod/15,
             optimizer=optimizer,lower=lower,upper=upper,
             start=start,fixed=fixed,distance="Geod",radius=radius)


fit_geo_ml$param

fit_geo_pl$param
#---------------------------------------------------



chor_ds=2*sin(geod_ds/2)
max_chor=max(chor_ds)
## maximum pairwise likelihood estimation using chordal  distance
fit_chor_pl <- GeoFit(data=data,coordx=coords,corrmodel="Matern", 
                 likelihood="Marginal",type="Pairwise",maxdist=max_chor/15,
                    optimizer=optimizer,lower=lower,upper=upper,
                    start=start,fixed=fixed,distance="Chor",radius=radius)
print(fit_chor_pl)
################################



## sinusoidal   projection
prj=mapproject(coords[,1], coords[,2], projection="sinusoidal",orientation=c(90,0,0)) 
coords_prj=cbind(prj$x,prj$y)


sinusoidal.proj = map(database= "world", ylim=c(-90,90), xlim=c(-180,180),
 col="grey80", fill=TRUE, plot=FALSE, projection="sinusoidal")
 map(sinusoidal.proj)
points(coords_prj,pch=20,cex=0.4)






eucl_ds=dist(coords_prj)
max_eucl=max(eucl_ds)
## maximum pairwise likelihood estimation using euclidean  distance
distance="Eucl"
fit_eucl_pl <- GeoFit(data=data,coordx=coords_prj,corrmodel="Matern", 
                    maxdist=max_eucl/15,likelihood="Marginal",type="Pairwise",
                    optimizer=optimizer,lower=lower,upper=upper,
                    start=start,fixed=fixed,distance=distance,radius=radius)


fit_chor_pl$param
fit_eucl_pl$param


boxplot(c(geod_ds),c(chor_ds),c(eucl_ds))





# computing residuals
res=GeoResiduals(fit_geo_pl)

### checking model assumptions: marginal distribution
GeoQQ(res)
### checking model assumptions: ST variogram model
vario = GeoVariogram(data=res$data,coordx=coords,maxdist=max_geod/2,distance="Geod",radius=radius)
GeoCovariogram(res,vario=vario,show.vario=TRUE,pch=20)





########### predictiob using geodesic distance

loc_to_pred=as.matrix(expand.grid(seq(-180,180,2),seq(60,90,1)))
dim(loc_to_pred)
param_est<-as.list(c(fixed,fit_geo_pl$param))
pr_geo=GeoKrig(data=data,loc=loc_to_pred,coordx=coords,corrmodel=corrmodel,
            radius=radius,distance="Geod",param=param_est,mse=TRUE)


predictions=cbind(loc_to_pred,pr_geo$pred,pr_geo$mse)
head(predictions)



globeearth(eye=place("newyorkcity"))
globepoints(loc=loc_to_pred,pch=20,col =  heat.colors(length(pr_geo$pred),alpha=0.1)[rank(pr_geo$pred)])



globeearth(eye=place("newyorkcity"))
globepoints(loc=loc_to_pred,pch=20,col =  cm.colors(length(pr_geo$mse),alpha=0.1)[rank(pr_geo$mse)])



)





