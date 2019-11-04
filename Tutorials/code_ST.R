rm(list=ls())

require(devtools)
install_github("vmoprojs/GeoModels -OCL") 
require(GeoModels)
require(fields)

model="Gaussian" # model name in the GeoModels package 
set.seed(12)


#################################
coordt =1:40 # number of temporal instants 
T=length(coordt)
NN=400 # number of spatial locations 
x = runif(NN, 0, 1); 
y = runif(NN, 0, 1) 
coords=cbind(x,y) 
X=cbind(rep(1,NN*T),runif(NN*T)) #  matrix of covariates
################################

#### regression parameters
mean = 0.5; mean1= -0.25
sill =2; nugget =0 # variance and sill


####################################################
corrmodel="Wend0_Wend0"   # model correlation
####################################################
scale_s=0.2; scale_t=2     # scales parameters 
####################################################
param= list(nugget=0,mean=mean,mean1=mean1,scale_t=scale_t,scale_s=scale_s,sill=sill,power2_s=4,power2_t=4)


## simulation
ss1 = GeoSim(coordx=coords, coordt=coordt, corrmodel=corrmodel,X=X,sparse=TRUE,
	          model=model,param=param)$data
## covariance matrix
cc = GeoCovmatrix(coordx=coords, coordt=coordt, corrmodel=corrmodel,X=X,sparse=TRUE,
	          model=model,param=param)

cc$nozero



## estimation with pairwise likelihood
start=list(mean=mean, mean1=mean1,scale_s=scale_s,scale_t=scale_t,sill=sill)
fixed=list(nugget=nugget,power2_s=4,power2_t=4)
fit = GeoFit(data=ss1,coordx=coords, coordt=coordt, corrmodel=corrmodel,
            maxdist=0.04,maxtime=1,X=X,
            optimizer="BFGS",
            start=start,fixed=fixed)
fit

# computing residuals
res=GeoResiduals(fit)
### checking model assumptions: marginal distribution
qqnorm(unlist(res$data))
abline(0,1)
### checking model assumptions: Spacetime variogram model
vario = GeoVariogram(data=res$data,coordx=coords, coordt=coordt, maxdist=0.6,maxtime=4)
GeoCovariogram(res,vario=vario,fix.lagt=1,fix.lags=1,show.vario=TRUE,pch=20)

###############################################
##### kriging on a regular grid at time 2###### 
###############################################
xx=seq(0,1,0.02)
loc_to_pred=as.matrix(expand.grid(xx,xx))
n_loc=nrow(loc_to_pred)
times=c(2)
n_tim=length(times)
Xloc=cbind(rep(1,n_loc*n_tim),runif(n_loc*n_tim))

param_est=as.list(c(fit$param,fixed))
pr = GeoKrig(data=ss1,coordx=coords, coordt=coordt, corrmodel=corrmodel,
	         X=X,Xloc=Xloc,sparse=TRUE,
	          model=model,mse=TRUE,loc=loc_to_pred,time=times,param=param_est)

par(mfrow=c(1,3))
colour <- rainbow(100)
# data map
quilt.plot(coords[,1],coords[,2],ss1[2,],col=colour,main ="Time=2") 
# kriging map
image.plot(xx, xx, matrix(pr$pred,ncol=length(xx)),col=colour, main = paste(" Kriging Time=" , 2),ylab="")
#MSE  kriging map
image.plot(xx, xx, matrix(pr$mse,ncol=length(xx)),col=colour,
           main = paste("MSE Time=" , 2),ylab="")


