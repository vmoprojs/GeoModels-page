#########################################################
#########################################################
# code fot the tutorial on Weibull  random fields analysis
#########################################################
#########################################################
rm(list=ls())
require(devtools)
install_github("vmoprojs/GeoModels")
require(GeoModels)
require(fields)
require(hypergeo)
model="Weibull" # model name in the GeoModels package



N=1000 # number of location sites 
set.seed(24)
x = runif(N, 0, 1)
y = runif(N, 0, 1)
coords=cbind(x,y) # spatial coordinates plot(coords,pch=20)
plot(coords,pch=20,xlab="",ylab="")


X=cbind(rep(1,N),runif(N)) # matrix covariates 

NuisParam(model,num_betas=2)

mean = -0.3; mean1 =0.5 # regression parameters
shape =2 # shape of the weibull RF 
nugget=0 # nugget parameter



corrmodel = "Wend0" ## correlation model and parameters 
scale = 0.2
power2 =4



param=list(mean=mean,mean1=mean1,sill=1-nugget, nugget=nugget, scale=scale ,power2=power2 ,shape=shape)
set.seed(312)
data = GeoSim(coordx=coords , corrmodel=corrmodel , model=model , param=param ,
X=X,sparse=TRUE)$data


start=list(mean=mean,mean1=mean1,scale=scale,shape=shape)
fixed=list(sill=1-nugget ,nugget=nugget ,power2=power2)
# Maximum composite-likelihood fitting of the Weibull random field:
fit = GeoFit(data=data,coordx=coords, corrmodel=corrmodel,model=model,X=X,
optimizer="BFGS",start=start,fixed=fixed,maxdist=0.02)




res=GeoResiduals(fit) # computing residuals

#### checking marginal assumptions
shape=fit$param["shape"]
probabilities = (1:N)/(N+1)
weibull.quantiles = qweibull(probabilities , shape=shape , scale = 1/(gamma(1+1/shape)))
plot(sort(weibull.quantiles), sort(c(res$data)), xlab="",ylab="",main="Weibull qq-plot") 
abline(0,1)
#### checking dependence assumptions: variogram
vario = GeoVariogram(data=res$data, coordx=coords,maxdist=0.3) # empirical variogram 
GeoCovariogram(res,show.vario=TRUE, vario=vario,pch=20)




##################
##### prediction
##################
xx=seq(0,1,0.013)
loc_to_pred=as.matrix(expand.grid(xx,xx))
Nloc=nrow(loc_to_pred)
Xloc=cbind(rep(1,Nloc),runif(Nloc))

param_est=as.list(c(fit$param,fixed))
pr=GeoKrig(data=data, coordx=coords,loc=loc_to_pred, X=X,Xloc=Xloc,
	corrmodel=corrmodel,model=model,mse=TRUE,
       sparse=TRUE,param=param_est)



colour = rainbow (100)

par(mfrow=c(1,3))
quilt.plot(x, y, data,col=colour,main="Data") ## map of simulated data 
map=matrix(pr$pred,ncol=length(xx))

map=matrix(pr$pred,ncol=length(xx))
image.plot(xx, xx, map,col=colour,xlab="",ylab="",main="Kriging") ## kriging map 

map_mse=matrix(pr$mse,ncol=length(xx)) ## associated MSE kriging map 
image.plot(xx, xx, map_mse,col=colour,xlab="",ylab="",main="MSE")