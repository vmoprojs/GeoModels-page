rm(list=ls())
require(devtools)
install_github("vmoprojs/GeoModels")
library(GeoModels)
library(fields)
set.seed(881)

model="Gaussian";

# Define the temporal coordinates 
coordt=seq(0,15,1)

###  list for spatial coordinates changing over time
coordx_dyn=list() 
X=list()
# Define the dynamical spatial coordinates
minN=150;maxN=250

# defining coordinates and matrix regression
for(k in 1:length(coordt))
{
NN=sample(minN:maxN,size=1)
x = runif(NN, 0, 1)
y = runif(NN, 0, 1)
coordx_dyn[[k]]=cbind(x,y)
X[[k]]=cbind(rep(1,NN),runif(NN))
}

#number of location sites for each temporal instant
unlist(lapply(coordx_dyn,nrow))

#par(mfrow=c(1,2))

plot(coordx_dyn[[1]],pch=20)
plot(coordx_dyn[[2]],pch=20)

# Set the  model parameters:
mean = 0.2
mean1= -0.3
sill = 1
nugget = 0

## correlation  model 
corrmodel="Exp_Exp"
scale_s = 0.3/3
scale_t = 4/3
param=list(mean=mean,mean1=mean1, sill=sill,nugget=nugget,scale_s=scale_s,scale_t=scale_t)

## simulation
ss1 = GeoSim(coordx_dyn=coordx_dyn, coordt=coordt, corrmodel=corrmodel,X=X,
	          model=model,param=param)$data


#cc = GeoCovmatrix(coordx_dyn=coordx_dyn, coordt=coordt, corrmodel=corrmodel,X=X,
#            model=model,param=param)$covmatrix
#min(eigen(cc)$values)

## estimation with pairwise likelihood
fixed=list(nugget=nugget)
start=list(mean=mean,mean1=mean1, sill=sill,scale_s=scale_s,scale_t=scale_t)
fit = GeoFit(data=ss1,coordx_dyn=coordx_dyn, coordt=coordt, corrmodel=corrmodel,
	          maxdist=0.1,maxtime=1,X=X,optimizer="BFGS",
	          start=start,fixed=fixed,model=model)


# computing residuals
res=GeoResiduals(fit)

### checking model assumptions: marginal distribution
qqnorm(unlist(res$data))
### checking model assumptions: ST variogram model
vario = GeoVariogram(data=res$data,coordx_dyn=coordx_dyn, coordt=coordt, maxdist=0.6,maxtime=7)
GeoCovariogram(res,vario=vario,fix.lagt=1,fix.lags=1,show.vario=TRUE,pch=20)





####################
##### kriging ###### 
####################
xx=seq(0,1,0.02)

## regular grid 
loc_to_pred=as.matrix(expand.grid(xx,xx))
## times  to predict
times=c(coordt[1],coordt[2])


Nloc=nrow(loc_to_pred)*length(times)
Xloc=cbind(rep(1,Nloc),runif(Nloc))


param_est=as.list(c(fit$param,fixed))
pr = GeoKrig(data=ss1,coordx_dyn=coordx_dyn, coordt=coordt, corrmodel=corrmodel,
	         X=X,Xloc=Xloc,model=model,loc=loc_to_pred,time=times,param=param_est)


par(mfrow=c(length(times),2))
colour = rainbow(100)

i=1
for(i in 1:length(times)) {
quilt.plot(coordx_dyn[[i]],ss1[[i]],col=colour)
image.plot(xx, xx, matrix(pr$pred[i,],ncol=length(xx)),col=colour,
            main = paste(" Kriging Time=" , times[i]),ylab="")
}
