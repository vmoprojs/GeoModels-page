rm(list=ls())
require(devtools)
#install_github("vmoprojs/GeoModels")
require(GeoModels)
require(fields)
require(maps)
require(maptools)
require(mapdata)
require(geoR)
require(sn)
library(mapproj)



#############
data(anomalies)
head(anomalies)
loc=cbind(anomalies[,1],anomalies[,2])
z=cbind(anomalies[,3])


### map data
quilt.plot(loc,z,xlab="long",ylab="lat")
map("usa", add = TRUE)

########### sinuisoidal projection  ##############
P.sinusoidal <- mapproject(loc[,1],loc[,2],projection="sinusoidal")
loc<-cbind(P.sinusoidal$x,P.sinusoidal$y)*6371
maxdist=max(dist(loc))
 

### histogram and scatterplot
summary(z)
hist(z,main="Anomalies histogram")
GeoScatterplot(data=z,coordx=loc,maxdist=50,numbins=4)



#### semivariogram 

evariog=GeoVariogram(data=z, coordx=loc,maxdist=maxdist/4)
plot(evariog$centers,evariog$variograms,ylim=c(0,1),pch=20,xlab="Km",ylab="Semi-variogram")





CorrParam("GenWend_Matern") 
NuisParam("Gaussian")
NuisParam("SkewGaussian")


corrmodel="GenWend_Matern"
###########
model="Gaussian"
start=list(mean=mean(z),sill=var(z),nugget=0.10,scale=200) 
fixed=list(smooth=0,power2=1/3.5)
pcl1=GeoFit(coordx=loc,corrmodel=corrmodel,data=z,
  model="Gaussian",
 maxpoints=10,optimizer="BFGS",
      start=start,fixed=fixed) 
pcl1
###########
model="SkewGaussian"
start=list(mean=mean(z),sill=var(z),nugget=0.1,skew=0.7,scale=200) 
fixed=list(smooth=0,power2=1/3.5)

pcl2=GeoFit(coordx=loc,corrmodel=corrmodel,data=z,
  likelihood="Marginal",type="Pairwise",model=model,
  maxpoints=10,optimizer="BFGS",start=start,fixed=fixed)
pcl2




#######################################################
############# residuals ###############################
res1=GeoResiduals(pcl1)

GeoQQ(res1)


evariog=GeoVariogram(data=res1$data,coordx=loc,maxdist=maxdist/4)

GeoCovariogram(res1,show.vario=TRUE, vario=evariog,pch=20)

#######################################################
res2=GeoResiduals(pcl2)

GeoQQ(res2)

evariog=GeoVariogram(data=res2$data,coordx=loc,maxdist=maxdist/4)

GeoCovariogram(res2,show.vario=TRUE,vario=evariog,pch=20)

#######################################################





################ covariance matrix ################################################
matrix1 = GeoCovmatrix(coordx=loc,corrmodel=corrmodel,model="Gaussian",sparse=TRUE,
  param=as.list(c(pcl1$param,pcl1$fixed)))
matrix2 = GeoCovmatrix(coordx=loc,corrmodel=corrmodel,model="SkewGaussian",sparse=TRUE,
  param=as.list(c(pcl2$param,pcl2$fixed)))
matrix1$nozero
matrix2$nozero




###### cross  validation
KK=100
seed=9
a1=GeoCV(pcl1, K=KK, n.fold=0.25,seed=seed,local=TRUE,maxpoints=100)
a2=GeoCV(pcl2, K=KK, n.fold=0.25,seed=seed,local=TRUE,maxpoints=100)
mean(a1$rmse);mean(a2$rmse);
mean(a1$mae);mean(a2$mae);

############# prediction map ####################
Sr1 = Polygon(loc)
Srs1 = Polygons(list(Sr1), "s1")
SpP = SpatialPolygons(list(Srs1))

long1=min(loc[,1])-5;long2=max(loc[,1])+5
lat1=min(loc[,2])-5;lat2=max(loc[,2])+5
lat_seq=seq(lat1,lat2,23)
lon_seq=seq(long1,long2,23)
coords_tot=as.matrix(expand.grid(lon_seq,lat_seq))
gr.in <- locations.inside(coords_tot, SpP)

pr1<-GeoKrigloc(loc=gr.in,coordx=loc,corrmodel=corrmodel,mse=TRUE,
  model="Gaussian",sparse=TRUE,maxpoints=100,
  param=as.list(c(pcl1$param,pcl1$fixed)),data=z)

pr2<-GeoKrigloc(loc=gr.in,coordx=loc,corrmodel=corrmodel,mse=TRUE,
  model="SkewGaussian",sparse=TRUE,maxpoints=100,
  param=as.list(c(pcl2$param,pcl2$fixed)),data=z)

par(mfrow=c(2,2))

quilt.plot(gr.in,pr1$pred)

quilt.plot(gr.in,pr1$mse) 

quilt.plot(gr.in,pr2$pred)

quilt.plot(gr.in,pr2$mse)


#######################################

