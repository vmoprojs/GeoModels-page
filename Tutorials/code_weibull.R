#########################################################
#########################################################
# code fot the tutorial on Weibull  random fields analysis
#########################################################
#########################################################
rm(list=ls())
require(GeoModels)
require(fields)
require(hypergeo)
set.seed(24)
model="Weibull" # model name in the GeoModels package




N=1000 # number of location sites 

x = runif(N, 0, 1)
y = runif(N, 0, 1)
coords=cbind(x,y) 

plot(coords,pch=20,xlab="",ylab="")

X=cbind(rep(1,N),runif(N)) # matrix covariates 

NuisParam(model,num_betas=2)

mean = -0.3; mean1 =0.5 # regression parameters
shape =2 # shape of the weibull RF 
nugget=0 # nugget parameter



corrmodel = "Wend0" ## correlation model and parameters 
CorrParam("Wend0")
scale = 0.2
power2 =4



param=list(mean=mean,mean1=mean1, nugget=nugget, scale=scale ,power2=power2 ,shape=shape)
set.seed(312)
data = GeoSim(coordx=coords , corrmodel=corrmodel , model=model , param=param ,
X=X)$data





start=list(mean=mean,mean1=mean1,scale=scale,shape=shape)
fixed=list(nugget=nugget ,power2=power2)
# Maximum composite-likelihood fitting of the Weibull random field:
fit = GeoFit(data=data,coordx=coords, corrmodel=corrmodel,model=model,X=X,
optimizer="BFGS",start=start,fixed=fixed,neighb=3)

fit

# computing residuals
res=GeoResiduals(fit) 
#### checking marginal assumptions
GeoQQ(res); GeoQQ(res,type="D")
#### checking dependence assumptions: variogram
vario = GeoVariogram(data=res$data, coordx=coords,maxdist=0.3) # empirical variogram 
GeoCovariogram(res, show.vario=TRUE, vario=vario,pch=20)



##################
##### prediction
##################
xx=seq(0,1,0.013)
loc_to_pred=as.matrix(expand.grid(xx,xx))
Nloc=nrow(loc_to_pred)
Xloc=cbind(rep(1,Nloc),runif(Nloc))


pr=GeoKrig(fit,loc=loc_to_pred,Xloc=Xloc,mse=TRUE,
       sparse=TRUE)



colour = rainbow (100)

par(mfrow=c(1,3))
quilt.plot(coords, data,col=colour,main="Data") ## map of simulated data 
quilt.plot(loc_to_pred, pr$pred,col=colour,xlab="",ylab="",main="Kriging") ## kriging map 
quilt.plot(loc_to_pred, pr$mse,col=colour,xlab="",ylab="",main="MSE") ## associated MSE kriging map