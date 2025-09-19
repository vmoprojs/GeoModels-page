
rm(list=ls());
require(GeoModels);
require(fields);
model="SinhAsinh"; # model name in the GeoModels package set.seed(989);
set.seed(91);

N=1500;
coords=cbind(runif(N),runif(N));


plot(coords ,pch=20,xlab="",ylab="");


corrmodel = "Matern";        ## correlation model
scale = 0.15/3;               ## scale parameter
smooth=0.5;                  ## smooth parameter
nugget=0;                    ## nugget parameter


skew=3.5;  ## skew parameter
tail=2;    ## tail parameter
sill= 1.5; ## variance parameter


mean = 1; mean1= -1 # regression paramteres 
a0=rep(1,N);a1=runif(N,-1,1)
X=cbind(a0,a1); ## regression matrix


##### simulation
param=list(nugget=nugget,mean=mean,mean1=mean1, scale=scale, smooth=smooth, sill=sill,skew=skew,tail=tail);
data <- GeoSim(coordx=coords ,corrmodel=corrmodel , param=param ,model=model ,X=X)$data;
hist(data,prob=TRUE,breaks=20)




##### estimation
optimizer="nlminb";

fixed<-list(nugget=nugget,smooth=smooth);
start<-list(mean=mean, mean1=mean1,scale=scale,sill=sill,skew=skew,tail=tail);
I=Inf;
lower<-list(mean=-I, mean1=-I,scale=0,sill=0,skew=-I,tail=0);
upper<-list(mean=I, mean1=I,scale=I,sill=I,skew=I,tail=I);


fit1 <- GeoFit2(data=data,coordx=coords,corrmodel=corrmodel,X=X,model=model,
               optimizer=optimizer,lower=lower,upper=upper,
               likelihood="Full",type="Standard", varest=TRUE,
               start=start,fixed=fixed)

fit1

fit2 <- GeoFit2(data=data,coordx=coords,corrmodel=corrmodel,X=X,model=model,
               optimizer=optimizer,lower=lower,upper=upper,neighb=2,
               likelihood="Marginal",type="Pairwise",sensitivity=TRUE,
               start=start,fixed=fixed)

fit2 

fit2=GeoVarestbootstrap(fit2,K=500,parallel=TRUE)

round(fit1$stderr,5);round(fit2$stderr,5)
round (fit1$conf.int,5); round (fit2$conf.int,5)




##### using ML estimation 
res=GeoResiduals(fit1) # computing residuals
GeoQQ(res)

vario <- GeoVariogram(data=res$data,coordx=coords,maxdist=0.4) 
GeoCovariogram(res,show.vario=TRUE, vario=vario,pch=20)
 
xx=seq(0,1,0.012) 
loc_to_pred=as.matrix(expand.grid(xx,xx)) 
Nloc=nrow(loc_to_pred)
Xloc=cbind(rep(1,Nloc),runif(Nloc))


pr=GeoKrig(fit1,loc=loc_to_pred,Xloc=Xloc,mse=TRUE)


#### map of  data
quilt.plot(coords[,1], coords[,2], data,main="Data")  

# linear kriging
map=matrix(pr$pred,ncol=length(xx))
image.plot(xx,xx,map,xlab="",ylab="",main="SimpleKriging")

#associated mean squared error
map_mse=matrix(pr$mse,ncol=length(xx))
image.plot(xx,xx,map_mse,xlab="",ylab="",main="MSE")




