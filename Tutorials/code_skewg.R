
################################################################
###
### Analysis of asymmetric spatial data with a skewGaussian RF
###
###############################################################
rm(list=ls())
require(devtools)
install_github("vmoprojs/GeoModels")
require(GeoModels)
require(fields)
require(sn)
model="SkewGaussian"
set.seed(8)



N=1200
# Define the spatial-coordinates of the points:
x = runif(N)
y = runif(N)
coords=cbind(x,y)
X=cbind(rep(1,N),runif(N))

NuisParam(model)
#### mean and covariance parameters ###
mean= 0.5
mean1=-0.5
sill=1.5
nugget=0
skew=-3



corrmodel = "GenWend"    ## correlation model and parameters
scale = 0.2
smooth=1
power2=5

#################
CorrParam("GenWend")


#############################
### Simulation
#############################
param=list(mean=mean,mean1=mean1,sill=sill,
	     nugget=nugget,scale=scale,skew=skew,power2=power2,smooth=smooth)
data = GeoSim(coordx=coords, corrmodel=corrmodel,model=model,X=X,
              param=param)$data
cc= GeoCovmatrix(coordx=coords, corrmodel=corrmodel,model=model,
             sparse=TRUE,X=X,
              param=param)


cc$nozero


#############################
## estimation with pairwise likelihood
#############################
start=list(sill=sill,mean=mean,mean1=mean1,scale=scale,skew=skew)
fixed=list(power2=power2,nugget=nugget,smooth=smooth)
fit=GeoFit(data=data,coordx=coords,corrmodel=corrmodel,X=X,
                    maxdist=0.04,model=model,
                    start=start,fixed=fixed)

fit


res=GeoResiduals(fit)
GeoQQ(res)
####
vario = GeoVariogram(data=res$data,coordx=coords,maxdist=0.5)
GeoCovariogram(res,show.vario=TRUE, vario=vario,pch=20)

#############################
## prediction
#############################                    
# locations to predict
xx=seq(0,1,0.012)
loc_to_pred=as.matrix(expand.grid(xx,xx))
Nloc=nrow(loc_to_pred)
Xloc=cbind(rep(1,Nloc),runif(Nloc))
param_est=as.list(c(fit$param,fixed))

pr=GeoKrig(data=data, coordx=coords,loc=loc_to_pred,corrmodel=corrmodel,model=model,mse=TRUE,X=X,Xloc=Xloc,
       sparse=TRUE,param= param_est)

colour = rainbow(100)
par(mfrow=c(1,3))
#### map of simulated data
quilt.plot(x, y, data,col=colour,main="Data")
map=matrix(pr$pred,ncol=length(xx))
## prediction map
image.plot(xx, xx, map,col=colour,xlab="",ylab="",main="Simple Kriging ")
## mse prediction map
map_mse=matrix(pr$mse,ncol=length(xx))
image.plot(xx, xx, map_mse,col=colour,xlab="",ylab="",main="mse ")
