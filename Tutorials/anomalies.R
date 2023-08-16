rm(list=ls())
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
plot(evariog,ylim=c(0,1),pch=20,xlab="Km",ylab="Semi-variogram")





CorrParam("Matern") 
NuisParam("Gaussian")
NuisParam("SkewGaussian")


corrmodel="Matern"
###########
I=Inf
start=list(mean=mean(z),sill=var(z),nugget=0.10,scale=140)
lower=list(mean=-I,sill=0,nugget=0,scale=0)
upper=list(mean=I,sill=I,nugget=1,scale=I)

fixed=list(smooth=0.5)
pcl1=GeoFit2(coordx=loc,corrmodel=corrmodel,data=z,
likelihood="Conditional",type="Pairwise",model="Gaussian",
 optimizer="nlminb",lower=lower,upper=upper,sensitivity=TRUE,
neighb=3, start=start,fixed=fixed)
pcl1
###########




start=list(nugget=0,scale=100,mean=mean(z),sill=var(z),skew=0.4)
fixed=list(smooth=.5)

lower2=list(mean=-I,sill=0,skew=-I,scale=0,nugget=0)
upper2=list(mean=I,sill=I,skew=I,scale=I,nugget=1)
pcl2=GeoFit2(coordx=loc,corrmodel=corrmodel,data=z,
  likelihood="Conditional",type="Pairwise",model="SkewGaussian",sensitivity=TRUE,
optimizer="nlminb",lower=lower2,upper=upper2,
  neighb=3,start=start,fixed=fixed)
pcl2

pcl1_est=GeoVarestbootstrap(pcl1,K=100,
                             optimizer="nlminb",lower=lower, upper=upper, 
                             method="Vecchia",M=30,seed=5)



pcl2_est=GeoVarestbootstrap(pcl2,K=100,
                             optimizer="nlminb",lower=lower2, upper=upper2, 
                             method="Vecchia",M=30,seed=5)

pcl1_est$stderr;pcl2_est$stderr;
pcl1_est$conf.int;pcl2_est$conf.int
pcl1_est$claic;pcl2_est$claic
#######################################################
############# residuals ###############################
res1=GeoResiduals(pcl1)
res2=GeoResiduals(pcl2)


GeoQQ(res1);GeoQQ(res2)


evariog=GeoVariogram(data=res1$data,coordx=loc,maxdist=maxdist/4)
GeoCovariogram(res1,show.vario=TRUE, vario=evariog,pch=20,ylim=c(0,1.3))
evariog=GeoVariogram(data=res2$data,coordx=loc,maxdist=maxdist/4)
GeoCovariogram(res2,show.vario=TRUE, vario=evariog,pch=20,ylim=c(0,1.7))

#######################################################


###### cross  validation
a1=GeoCV(pcl1, K=100, estimation=TRUE, n.fold=0.25,seed=9,local=TRUE,neighb=100)
a2=GeoCV(pcl2, K=100,estimation=TRUE, n.fold=0.25,seed=9,local=TRUE,neighb=100)
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
  model="Gaussian",neighb=100,
  param=as.list(c(pcl1$param,pcl1$fixed)),data=z)

pr2<-GeoKrigloc(loc=gr.in,coordx=loc,corrmodel=corrmodel,mse=TRUE,
  model="SkewGaussian",neighb=100,
  param=as.list(c(pcl2$param,pcl2$fixed)),data=z)

par(mfrow=c(2,2))

quilt.plot(gr.in,pr1$pred)

quilt.plot(gr.in,pr1$mse) 

quilt.plot(gr.in,pr2$pred)

quilt.plot(gr.in,pr2$mse)


#######################################

