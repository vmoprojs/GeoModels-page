rm(list=ls())
require(GeoModels)
require(fields)
model="Gaussian" # model name in the GeoModels package 
set.seed(12)


coordt =1:30 # number of temporal instants 
T=length(coordt)
NN=400 # number of spatial locations 
x = runif(NN, 0, 1); 
y = runif(NN, 0, 1) 
coords=cbind(x,y) 
X=cbind(rep(1,NN*T),runif(NN*T))



mean = 0.5; mean1= -0.25
sill =2; nugget =0


####################################################

corrmodel="Wend0_Wend0";
####################################################
scale_s=0.2; scale_t=2     # scale 
####################################################
param= list(nugget=0,mean=mean,mean1=mean1,scale_t=scale_t,scale_s=scale_s,sill=sill,power2_s=4,power2_t=4)





## simulation
ss1 = GeoSim(coordx=coords, coordt=coordt, corrmodel=corrmodel,X=X,sparse=TRUE,
	          model=model,param=param)$data
cc = GeoCovmatrix(coordx=coords, coordt=coordt, corrmodel=corrmodel,X=X,sparse=TRUE,
	          model=model,param=param)


#U=chol.spam(cc$covmatrix)
#a=chol2inv.spam(U)

cc$nozero

## estimation with pairwise likelihood
start=list(mean=mean, mean1=mean1,scale_s=scale_s,scale_t=scale_t,sill=sill)
fixed=list(nugget=nugget,power2_s=4,power2_t=4)
fit = GeoFit(data=ss1,coordx=coords, coordt=coordt, corrmodel=corrmodel,
                neighb =5,maxtime=1,X=X,
            optimizer="BFGS",start=start,fixed=fixed)

fit
# computing residuals
res=GeoResiduals(fit)

### checking model assumptions: marginal distribution
GeoQQ(res)
### checking model assumptions: ST variogram model
vario = GeoVariogram(data=res$data,coordx=coords, coordt=coordt, maxdist=0.6,maxtime=4)
GeoCovariogram(res,vario=vario,fix.lagt=1,fix.lags=1,show.vario=TRUE,pch=20)


####################
##### kriging ###### 
####################
xx=seq(0,1,0.02)
loc_to_pred=as.matrix(expand.grid(xx,xx))
n_loc=nrow(loc_to_pred)
times=c(2)
n_tim=length(times)
Xloc=cbind(rep(1,n_loc*n_tim),runif(n_loc*n_tim))

pr = GeoKrig(fit,loc=loc_to_pred,time=times,Xloc=Xloc,sparse=TRUE,mse=TRUE)

par(mfrow=c(1,3))
colour <- rainbow(100)
quilt.plot(coords[,1],coords[,2],ss1[2,],col=colour,main ="Time=2")
image.plot(xx, xx, matrix(pr$pred,ncol=length(xx)),col=colour, main = paste(" Kriging Time=" , 2),ylab="")
image.plot(xx, xx, matrix(pr$mse,ncol=length(xx)),col=colour,
           main = paste("MSE Time=" , 2),ylab="")


prloc = GeoKrigloc(fit,loc=loc_to_pred,time=times,Xloc=Xloc,sparse=TRUE,mse=TRUE,neighb=150,maxtime=2,parallel=TRUE)

summary(c(pr$pred-prloc$pred))
