rm(list=ls())
require(GeoModels)
require(fields)
set.seed(89)


NN=500          # number of spatial locations
x = runif(NN, 0, 1); y = runif(NN, 0, 1)
coords=cbind(x,y)
X=cbind(rep(1,NN),runif(NN))


mean = 5; sill=1; nugget=0

NuisParam ("Gaussian")
corrmodel="GenWend_Matern"
CorrParam (corrmodel)

mu=3
scale=0.05; power2=1/mu;  smooth=0

scale*(1/power2)


param= list(nugget=nugget,mean=mean,
     scale=scale, sill=sill, power2=power2,smooth=smooth)
sim = GeoSim(coordx=coords,  corrmodel=corrmodel,
              model="Gaussian",param=param)$data


cc = GeoCovmatrix(coordx=coords, corrmodel=corrmodel,  sparse=TRUE, 
model="Gaussian", param=param)
is.spam(cc$covmatrix)
cc$nozero

GeoCovDisplay(cc)



optimizer="nlminb"
start=list(mean=mean,scale=scale, sill=sill, power2=power2)
I=Inf
lower=list(mean=-I,scale=0, sill=0, power2=0)
upper=list(mean=I,scale=I, sill=I, power2=1/(1.5))

fixed=list(nugget=nugget,smooth=smooth)
fitML <- GeoFit(data=sim,coordx=coords,corrmodel=corrmodel, 
                    model="Gaussian",varest=TRUE,
                    optimizer=optimizer,lower=lower,upper=upper,
                    likelihood="Full",type="Standard",
                    start=start,fixed=fixed)

(1/fitML$param$power2)*fitML$param$scale


fitML

res=GeoResiduals(fitML) 



vario = GeoVariogram(data=res$data,coordx=coords,maxdist=0.6)

GeoQQ(res)
GeoCovariogram(res,vario=vario,show.vario=TRUE,pch=20)


 GeoQQ(res,type="D",ylim=c(0,0.5),breaks=20)


xx=seq(0,1,0.015)
loc_to_pred=as.matrix(expand.grid(xx,xx))    # locations to predict


param_est=as.list(c(fitML$param,fixed))
pr = GeoKrig(data=sim,coordx=coords,  corrmodel=corrmodel,
 sparse=TRUE,model="Gaussian",mse=TRUE,loc=loc_to_pred,
	         param=param_est)




quilt.plot(coords[,1],coords[,2],sim,main ="Observed")

## kriging map

image.plot(xx, xx, matrix(pr$pred,ncol=length(xx)),
           main = paste("Prediction" ),ylab="")

## MSE kriging map 

image.plot(xx, xx, matrix(pr$mse,ncol=length(xx)),
           main = paste("MSE"),ylab="")



