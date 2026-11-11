############################################################
# GeoModels tutorial: Student-t spatial random fields
############################################################

rm(list=ls())
require(GeoModels)
require(fields)
require(hypergeo)
model="StudentT"  # model name in the GeoModels package
set.seed(199)

N=600
coords=cbind(runif(N),runif(N))

plot(coords,pch=20,xlab="",ylab="")

corrmodel="Matern"
scale=0.2/3
smooth=0.5
nugget=0

df=5    # degrees of freedom (GeoModels uses its reciprocal)
sill=1 # sigma^2; marginal variance is sill*df/(df-2)
mean=0.5; mean1=-1

a0=rep(1,N); a1=runif(N,-1,1)
X=cbind(a0,a1)
param=list(nugget=nugget,mean=mean,mean1=mean1,
           scale=scale,smooth=smooth,sill=sill,df=1/df)
data <-GeoSim(coordx=coords,corrmodel=corrmodel,
                           param=param,model=model,X=X)$data

optimizer="nlminb"
fixed1<-list(nugget=nugget,smooth=smooth)
start1<-list(mean=mean,mean1=mean1,scale=scale,sill=sill,df=1/df)
I=Inf
lower1<-list(mean=-I,mean1=-I,scale=1e-6,sill=1e-6,df=1e-4)
upper1<-list(mean=I,mean1=I,scale=I,sill=I,df=0.5-1e-4)
fit1 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel,
              optimizer=optimizer,lower=lower1,upper=upper1,
              neighb=3,X=X,start=start1,fixed=fixed1,model=model)


DF=max(3L,as.integer(round(1/fit1$param$df)))
print(DF)

start<-list(mean=mean,mean1=mean1,scale=scale,sill=sill)
fixed<-list(nugget=nugget,df=1/DF,smooth=smooth)
lower<-list(mean=-I,mean1=-I,scale=1e-6,sill=1e-6)
upper<-list(mean=I,mean1=I,scale=I,sill=I)
fit2 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel,
              optimizer=optimizer,lower=lower,upper=upper,
              neighb=3,X=X,start=start,fixed=fixed,model=model)
fit2

fit3 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel,
              optimizer=optimizer,lower=lower,upper=upper,
              neighb=3,X=X,start=start,fixed=fixed,
              model="Gaussian_misp_StudentT")
unlist(fit2$param)
unlist(fit3$param)

res=GeoResiduals(fit2)

GeoQQ(res)

vario <- GeoModels::GeoVariogram(data=res$data,coordx=coords,maxdist=0.4)

GeoCovariogram(res,show.vario=TRUE,vario=vario,pch=20,ylim=c(0,2))

xx=seq(0,1,0.012)
loc_to_pred=as.matrix(expand.grid(xx,xx))
Nloc=nrow(loc_to_pred)
Xloc=cbind(rep(1,Nloc),runif(Nloc))

# predict 
param_est=as.list(c(fit2$param,fixed))
pr=GeoKrig(data=data,coordx=coords,loc=loc_to_pred,
     X=X,Xloc=Xloc,corrmodel=corrmodel,model=model,mse=TRUE,param=param_est)


par(mfrow=c(1,3))
quilt.plot(coords[,1],coords[,2],data,main="Data")
map=matrix(pr$pred,ncol=length(xx))
image.plot(xx,xx,map,xlab="",ylab="",main="SimpleKriging")
map_mse=matrix(pr$mse,ncol=length(xx))
image.plot(xx,xx,map_mse,xlab="",ylab="",main="MSE")
