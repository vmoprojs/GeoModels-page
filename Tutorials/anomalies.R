################################################################
# GeoModels Tutorial: analysis of spatial precipitation anomalies
# Gaussian and skew-Gaussian random fields
################################################################
rm(list=ls())
require(GeoModels)
require(fields)
require(mapproj)
require(maps)              
require(sp)
require(geoR)


# Preliminary data analysis
data("anomalies", package="GeoModels")
head(anomalies)
loc_lonlat <- cbind(anomalies[,1],anomalies[,2])
z <- as.numeric(anomalies[,3])

# Data map (longitude/latitude coordinates)
quilt.plot(loc_lonlat,z,xlab="long",ylab="lat")
map("usa",add=TRUE)

# Sinusoidal projection. Save the orientation for the map boundary.
P.sinusoidal <- mapproj::mapproject(loc_lonlat[,1],loc_lonlat[,2],
                                    projection="sinusoidal")
proj_orientation <- get(".Last.projection",envir=.GlobalEnv)$orientation
loc <- cbind(P.sinusoidal$x,P.sinusoidal$y)*6371
# Exact diameter of a planar point cloud is attained on its convex hull.
maxdist <- max(stats::dist(loc[grDevices::chull(loc),,drop=FALSE]))

summary(z)
hist(z,main="Anomalies histogram")
GeoScatterplot(data=z,coordx=loc,neighb=c(1,2,3,4))

# Empirical semivariogram; same cutoff as used in the original PDF.
evariog <- GeoVariogram(data=z,coordx=loc,
                                    maxdist=maxdist/4)
plot(evariog,ylim=c(0,1),pch=20,xlab="Km",ylab="Semi-variogram")

# Correlation and nuisance parameters
CorrParam("Matern")
NuisParam("Gaussian")
NuisParam("SkewGaussian")
corrmodel <- "Matern"
I <- Inf

# Gaussian: conditional composite likelihood
start <- list(mean=mean(z),sill=var(z),nugget=0.10,scale=200)
lower <- list(mean=-I,sill=1e-6,nugget=0,scale=1e-6)
upper <- list(mean=I,sill=I,nugget=1-1e-6,scale=I)
fixed <- list(smooth=0.5)
pcl1 <- GeoFit2(coordx=loc,corrmodel=corrmodel,data=z,
                 likelihood="Conditional",type="Pairwise",model="Gaussian",
                 optimizer="nlminb",lower=lower,upper=upper,
                 sensitivity=TRUE,neighb=3,start=start,fixed=fixed)
pcl1

# Skew-Gaussian: conditional composite likelihood
start <- list(nugget=0.1,scale=200,mean=mean(z),
              sill=var(z),skew=0.4)
lower2 <- list(mean=-I,sill=1e-6,skew=-I,scale=1e-6,nugget=0)
upper2 <- list(mean=I,sill=I,skew=I,scale=I,nugget=1-1e-6)
pcl2 <- GeoFit2(coordx=loc,corrmodel=corrmodel,data=z,
                 likelihood="Conditional",type="Pairwise",model="SkewGaussian",
                 sensitivity=TRUE,optimizer="nlminb",lower=lower2,
                 upper=upper2,neighb=3,start=start,fixed=fixed)
pcl2

# Parametric bootstrap (computationally intensive)
set.seed(5)
pcl1_est <- GeoVarestbootstrap(pcl1,K=100,
              optimizer="nlminb",lower=lower,upper=upper,parallel=TRUE)
set.seed(5)
pcl2_est <- GeoVarestbootstrap(pcl2,K=100,
              optimizer="nlminb",lower=lower2,upper=upper2,parallel=TRUE)
pcl1_est$stderr; pcl2_est$stderr
pcl1_est$conf.int; pcl2_est$conf.int
pcl1_est$claic; pcl2_est$claic
GeoTests(pcl2_est,pcl1_est,statistic="Wald")

# Residual diagnostics
res1 <- GeoResiduals(pcl1)
res2 <- GeoResiduals(pcl2)
GeoQQ(res1,type="D",breaks=50,ylim=c(0,0.45))
GeoQQ(res2,type="D",breaks=50,ylim=c(0,0.45))
vario1 <- GeoVariogram(data=res1$data,coordx=loc,
                                   maxdist=maxdist/4)
GeoCovariogram(res1,show.vario=TRUE,vario=vario1,
                           pch=20,ylim=c(0,1.3))
vario2 <- GeoVariogram(data=res2$data,coordx=loc,
                                   maxdist=maxdist/4)
GeoCovariogram(res2,show.vario=TRUE,vario=vario2,
                           pch=20,ylim=c(0,1.7))

# Predictive validation (computationally intensive)
set.seed(9)
a1 <- GeoCV(pcl1,K=100,estimation=TRUE,n.fold=0.25,
                       local=TRUE,neighb=100,parallel=TRUE)
set.seed(9)
a2 <- GeoCV(pcl2,K=100,estimation=TRUE,n.fold=0.25,
                       local=TRUE,neighb=100,parallel=TRUE)
mean(a1$rmse);mean(a2$rmse)
mean(a1$mae);mean(a2$mae)

# Prediction grid within the USA outline 

long1 <- min(loc[,1])-5;long2 <- max(loc[,1])+5
lat1 <- min(loc[,2])-5;lat2 <- max(loc[,2])+5
lat_seq <- seq(lat1,lat2,23)
lon_seq <- seq(long1,long2,23)
coords_tot <- as.matrix(expand.grid(lon_seq,lat_seq))

Sr1 = Polygon(loc)
Srs1 = Polygons(list(Sr1), "s1")
SpP = SpatialPolygons(list(Srs1))

long1=min(loc[,1])-5;long2=max(loc[,1])+5
lat1=min(loc[,2])-5;lat2=max(loc[,2])+5
lat_seq=seq(lat1,lat2,23)
lon_seq=seq(long1,long2,23)
coords_tot=as.matrix(expand.grid(lon_seq,lat_seq))
gr.in <- locations.inside(coords_tot, SpP)


pr1 <- GeoKrigloc(pcl1,loc=gr.in,mse=TRUE,neighb=100)
pr2 <- GeoKrigloc(pcl2,loc=gr.in,mse=TRUE,neighb=100)

par(mfrow=c(2,2))
quilt.plot(gr.in,pr1$pred)
quilt.plot(gr.in,pr1$mse)
quilt.plot(gr.in,pr2$pred)
quilt.plot(gr.in,pr2$mse)

par(mfrow=c(2,2))

quilt.plot(gr.in,pr1$pred)
quilt.plot(gr.in,pr1$mse) 
quilt.plot(gr.in,pr2$pred)
quilt.plot(gr.in,pr2$mse)
