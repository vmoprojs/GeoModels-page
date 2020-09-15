rm(list=ls())
require(devtools)
#install_github("vmoprojs/GeoModels")
library(GeoModels)
library(fields)
require(limma)
require(oz)
library(maps)
require(maptools)
require(mapdata)
require(geoR)
#########################


data(austemp)
radius=6371
distance="Geod"


############# preliminary aanalysis  #####
coords=cbind(austemp[,1],austemp[,2])
temp=austemp[,3]
quilt.plot(coords,temp,xlab="long",ylab="lat")
oz(states=FALSE,add=T, lwd=2) #


hist(temp,main="Histogram of Maximum temperature",nclass=13)
maxdist=max(rdist.earth(coords,miles=F,R=radius))
fit <- GeoVariogram(coordx=coords,data=temp,distance="Geod",maxdist=maxdist/5,radius=radius)

GeoScatterplot(data=temp,coordx=coords,distance="Geod",maxdist=500,numbins=6)

# Results:
plot(fit$centers, fit$variograms, xlab='GC distance', ylab=expression(gamma(GC)),cex=1.5,cex.lab=1.5,
ylim=c(0, max(fit$variograms)), pch=20,
main="Semi-variogram")

########################################



corrmodel="Matern"
CorrParam(corrmodel)
NuisParam("Gaussian",num_betas=2)
NuisParam("StudentT",num_betas=2)
NN=nrow(coords)
X=cbind(rep(1,NN),austemp[,4])






############################################################################
##### starting values #########
mean=7.5
mean1=1
sill=8; nugget=0

smooth=0.5;scale=60; 

###############################

optimizer="nlminb"     # 

fixed1=list(nugget=0,smooth=smooth)

I=Inf

start1=list(mean=mean,mean1=mean1,scale=scale,sill=sill)

lower1=list(mean=-I,mean1=-I,scale=0,sill=0)
upper1= list(mean=I,mean1=I,scale=I,sill=I)
##############################################
########### Gaussian case ####################
##############################################
##### ml estimation

#fitml_m = GeoFit(data=temp,coordx=coords,corrmodel=corrmodel,X=X,model="Gaussian",optimizer=optimizer,
#                    distance=distance,radius=radius,
#                       lower=lower1,upper=upper1,
#                    likelihood="Full",type="Standard",#varest=TRUE,
#                    start=start1,fixed=fixed1)
#print(fitml_m)



dd=150
##### pl estimation
fit2 = GeoFit(data=temp,coordx=coords,corrmodel=corrmodel,X=X,model="Gaussian",
                    maxdist=dd,distance=distance,radius=radius,sensitivity=TRUE,
                       optimizer=optimizer, lower=lower1,upper=upper1,
                    start=start1,fixed=fixed1)
print(fit2)
##############################################
########### Student case  ####################
##############################################
df=4
start1=list(mean=mean,mean1=mean1, scale=scale,df=1/df,sill=sill)

lower2=list(mean=-I,mean1=-I,scale=0,sill=0,df=0)
upper2= list(mean=I,mean1=I,scale=I,sill=I,df=0.5)

fit3 = GeoFit(data=temp,coordx=coords,corrmodel=corrmodel,model="StudentT",
                    maxdist=dd,distance=distance,radius=radius,X=X,sensitivity=TRUE,
                    optimizer=optimizer,lower=lower2,upper=upper2,
                    start=start1,fixed=fixed1)       
print(fit3)
DF=as.numeric(round(1/fit3$param['df']))
if(DF==2) DF=3
print(DF)


fixed=list(nugget=nugget,df=1/DF,smooth=smooth)
start=list(mean=mean,mean1=mean1, scale=scale,sill=sill)
fit4 = GeoFit(data=temp,coordx=coords,corrmodel=corrmodel,X=X,optimizer=optimizer,
                    maxdist=dd,distance=distance,radius=radius, sensitivity=TRUE,
                       lower=lower1,upper=upper1,
                    start=start,fixed=fixed, model="StudentT")
                   
print(fit4)

fit2 #Gaussian
fit4#T


##############################################
### computing stderr  estimation with psarametric  bootstrap ############
##############################################
KK=100
v1=GeoVarestbootstrap(fit2,K=KK,optimizer=optimizer,seed=9)#Gaussian
v2=GeoVarestbootstrap(fit4,K=KK,,optimizer=optimizer,seed=9)#T


v1;v2




##############################################
### computing estimated residulas ############
##############################################
res_g=GeoResiduals(fit2)
res_t=GeoResiduals(fit4)

GeoQQ(res_g)
GeoQQ(res_t)

########### some diagnostic graphics ############################################
varionorm = GeoVariogram(data=res_t$data,coordx=coords,maxdist=maxdist/5,distance=distance,radius=radius)
GeoCovariogram(res_t,show.vario=TRUE, vario=varionorm,pch=20)





#####  prediction at one location point
coords_to_pred=matrix(c(135,-25),ncol=2)
## prediction of residuals
pr_student2<-GeoKrig(data=res_t$data, coordx=coords,loc=coords_to_pred,corrmodel=corrmodel,
    distance=distance,radius=radius,mse=TRUE,
   model="StudentT",param= as.list(c(res_t$param,res_t$fixed)))

pr_student2$pred
pr_student2$mse


####  prediction performance trough cross validation
KK=200
d=GeoCV(fit4,K=KK,n.fold=0.2,seed=9) #  T
e=GeoCV(fit2,K=KK,n.fold=0.2,seed=9) #gaussian

mean(d$rmse);mean(d$mae); #  T
mean(e$rmse);mean(e$mae); #  gaussian


################################################################################
####################  map prediction #############################
################################################################################
coord <- maps:::map.poly("worldHires", "Australia", exact=FALSE,
                         xlim=c(110,160),ylim=c(-45,-5),
                         boundary=FALSE,
                         interior=TRUE, fill=FALSE, as.polygon=TRUE)
coord.sp <- map2SpatialPolygons(coord,coord$names)

#par(mfrow=c(1,3))

#boundary=cbind(coord$x,coord$y)
#boundary1=boundary[!is.na(boundary[,1]),]

#plot(boundary1)

long1=110;long2=154
lat1=-39;lat2=-12

lat_seq=seq(lat1,lat2,0.35)
lon_seq=seq(long1,long2,0.35)
coords_tot=as.matrix(expand.grid(lon_seq,lat_seq))


gr.in <- locations.inside(coords_tot, coord.sp)

plot(gr.in)


boundary=cbind(coord$x,coord$y)
boundary1=boundary[!is.na(boundary[,1]),]

plot(boundary1)


pr_gaussian<-GeoKrig(data=res_g$data, coordx=coords,loc=gr.in,
  corrmodel=corrmodel,distance=distance,radius=radius,
  mse=TRUE,model="Gaussian",param= as.list(c(res_g$param,res_g$fixed)))

pr_student2<-GeoKrig(data=res_t$data, coordx=coords,loc=gr.in,
  corrmodel=corrmodel,distance=distance,radius=radius,mse=TRUE,
  model="StudentT",param= as.list(c(res_t$param,res_t$fixed)))

quilt.plot(gr.in,pr_gaussian$pred)
oz(states=FALSE,add=T, lwd=2)
quilt.plot(gr.in,pr_gaussian$mse)
oz(states=FALSE,add=T, lwd=2)
quilt.plot(gr.in,pr_student2$pred)
oz(states=FALSE,add=T, lwd=2)
quilt.plot(gr.in,pr_student2$mse)
oz(states=FALSE,add=T, lwd=2)
