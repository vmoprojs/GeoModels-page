require(GeoModels)
require(fields)
require(spam)
set.seed(89)


NN=500 # number of spatial locations
x=runif(NN,0,1); y=runif(NN,0,1)
coords=cbind(x,y)
X=cbind(rep(1,NN),runif(NN)) # unused; preserves original R seed sequence


GeoModels::NuisParam("Gaussian")
mean=5; sill=1; nugget=0



corrmodel="GenWend_Matern"
GeoModels::CorrParam(corrmodel)
mu=3
scale=0.05; power2=1/mu; smooth=0


scale*(1/power2)


param=list(nugget=nugget,mean=mean,scale=scale,
sill=sill,power2=power2,smooth=smooth)
sim=GeoModels::GeoSim(coordx=coords,corrmodel=corrmodel,
sparse=FALSE,model="Gaussian",param=param)$data


cc=GeoModels::GeoCovmatrix(coordx=coords,corrmodel=corrmodel,
sparse=TRUE,model="Gaussian",param=param)
is.spam(cc$covmatrix)
cc$nozero

GeoCovDisplay(cc)


optimizer="nlminb"
start=list(mean=mean,scale=scale,sill=sill,power2=power2)
I=Inf
lower=list(mean=-I,scale=1e-5,sill=1e-5,power2=1e-4)
upper=list(mean=I,scale=I,sill=I,power2=1/1.5)
fixed=list(nugget=nugget,smooth=smooth)

fitML <- GeoModels::GeoFit(data=sim,coordx=coords,
    corrmodel=corrmodel,model="Gaussian",varest=TRUE,
    optimizer=optimizer,lower=lower,upper=upper,
    likelihood="Full",type="Standard",start=start,
    fixed=fixed)


fitML$param[["scale"]] / fitML$param[["power2"]]




res=GeoModels::GeoResiduals(fitML)

GeoModels::GeoQQ(res)

vario=GeoModels::GeoVariogram(data=res$data,
                              coordx=coords,maxdist=0.6)

GeoModels::GeoCovariogram(res,vario=vario,
                           show.vario=TRUE,pch=20)


xx=seq(0,1,0.015)
loc_to_pred=as.matrix(expand.grid(xx,xx))


pr=GeoModels::GeoKrig(fitML,loc=loc_to_pred,
                       sparse=TRUE,mse=TRUE)


quilt.plot(coords[,1],coords[,2],sim,main="Observed")

## kriging map

image.plot(xx,xx,matrix(pr$pred,ncol=length(xx)),
                   main="Prediction",ylab="")

## MSE kriging map 
image.plot(xx,xx,matrix(pr$mse,ncol=length(xx)),
                   main="MSE",ylab="")
