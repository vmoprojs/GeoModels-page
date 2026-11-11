###############################################################
### Analysis of asymmetric spatial data with a sinh-arcsinh RF
###############################################################
rm(list=ls())
require(GeoModels)
require(fields)
model="SinhAsinh"  # model name in the GeoModels package
set.seed(91)

N=1500
coords=cbind(runif(N),runif(N))
plot(coords,pch=20,xlab="",ylab="")

corrmodel="Matern"        # correlation model
scale=0.15/3             # scale parameter
smooth=0.5              # smoothness parameter
nugget=0                # nugget parameter

skew=3.5                # skewness parameter
tail=2                 # tail-weight parameter (>1: lighter than Gaussian)
sill=1.5               # scale parameter sigma^2, not the marginal variance

mean=1; mean1=-1         # regression LOCATION parameters
a0=rep(1,N); a1=runif(N,-1,1)
X=cbind(a0,a1)           # regression matrix

##### simulation
param=list(nugget=nugget,mean=mean,mean1=mean1,scale=scale,
           smooth=smooth,sill=sill,skew=skew,tail=tail)
data <- GeoSim(coordx=coords,corrmodel=corrmodel,
                          param=param,model=model,X=X)$data
hist(data,prob=TRUE,breaks=20)

##### estimation
optimizer="nlminb"
fixed<-list(nugget=nugget,smooth=smooth)
start<-list(mean=mean,mean1=mean1,scale=scale,sill=sill,
            skew=skew,tail=tail)
I=Inf
# Strictly positive parameters must not have zero as an optimization boundary.
lower<-list(mean=-I,mean1=-I,scale=1e-6,sill=1e-6,
            skew=-I,tail=1e-6)
upper<-list(mean=I,mean1=I,scale=I,sill=I,skew=I,tail=I)

fit1 <- GeoFit2(data=data,coordx=coords,corrmodel=corrmodel,
                           X=X,model=model,optimizer=optimizer,
                           lower=lower,upper=upper,
                           likelihood="Full",type="Standard",varest=TRUE,
                           start=start,fixed=fixed)
fit1

fit2 <- GeoFit2(data=data,coordx=coords,corrmodel=corrmodel,
                           X=X,model=model,optimizer=optimizer,
                           lower=lower,upper=upper,neighb=2,
                           likelihood="Marginal",type="Pairwise",
                           sensitivity=TRUE,start=start,fixed=fixed)
fit2

# Bootstrap standard errors.
fit2=GeoVarestbootstrap(fit2,K=500,parallel=TRUE)
round(fit1$stderr,5);round(fit2$stderr,5)


# Report Wald intervals computed from estimates and SEs.
wald_ci <- function(fit,level=0.95) {
  estimates <- unlist(fit$param)
  errors <- unlist(fit$stderr)
  common <- intersect(names(estimates),names(errors))
  if (!length(common)) stop("No matching estimates and standard errors")
  z <- stats::qnorm((1+level)/2)
  cbind(lower=estimates[common]-z*errors[common],
        upper=estimates[common]+z*errors[common])
}
round(wald_ci(fit1),5);round(wald_ci(fit2),5)

##### checking model assumptions using the full ML fit
res=GeoResiduals(fit1)
GeoQQ(res)
vario <- GeoVariogram(data=res$data,coordx=coords,maxdist=0.4)
GeoCovariogram(res,show.vario=TRUE,vario=vario,pch=20)

##### prediction
xx=seq(0,1,0.012)
loc_to_pred=as.matrix(expand.grid(xx,xx))
Nloc=nrow(loc_to_pred)
Xloc=cbind(rep(1,Nloc),runif(Nloc))

pr=GeoKrig(fit1,loc=loc_to_pred,Xloc=Xloc,mse=TRUE)


par(mfrow=c(1,3))
quilt.plot(coords[,1],coords[,2],data,main="Data")
map=matrix(pr$pred,ncol=length(xx))
image.plot(xx,xx,map,xlab="",ylab="",main="SimpleKriging")
map_mse=matrix(pr$mse,ncol=length(xx))
image.plot(xx,xx,map_mse,xlab="",ylab="",main="MSE")
