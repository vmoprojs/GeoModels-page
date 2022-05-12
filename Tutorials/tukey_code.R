
rm(list=ls());
require(GeoModels);
require(fields);
model="Tukeyh"; # model name in the GeoModels package set.seed(989);
set.seed(818);

N=1500;
coords=cbind(runif(N),runif(N));




 plot(coords ,pch=20,xlab="",ylab="");


corrmodel = "Matern";
scale = 0.2/3;
smooth =0.5;
nugget=0;


tail=0.1; # tail parameter
sill= 1; # variance parameter


mean = 0.5; mean1= -1 # regression paramteres 
a0=rep(1,N);a1=runif(N,-1,1)
X=cbind(a0,a1); ## regression matrix


# simulation
param=list(nugget=nugget,mean=mean,mean1=mean1, scale=scale, smooth=smooth, sill=sill,tail=tail);

data <- GeoSim(coordx=coords ,corrmodel=corrmodel , param=param ,model=model ,X=X)$data;



### estimation
optimizer="nlminb";

fixed1<-list(nugget=nugget,smooth=smooth);
start1<-list(mean=mean, mean1=mean1,scale=scale,sill=sill,tail=tail);
I=Inf;
lower1<-list(mean=-I, mean1=-I,scale=0,sill=0,tail=0);
upper1<-list(mean=I, mean1=I,scale=I,sill=I,tail=0.5);



fit2 <- GeoFit2(data=data,coordx=coords,corrmodel=corrmodel,
optimizer=optimizer,lower=lower1,upper=upper1,
type="Pairwise",likelihood="Conditional",
neighb=4,X=X,start=start1,fixed=fixed1, model = model);

#### some graphics ######
fit2 

res=GeoResiduals(fit2) # computing residuals


GeoQQ(res)
GeoQQ(res,type="D",ylim=c(0,0.5),breaks=20)



 vario <- GeoVariogram(data=res$data,coordx=coords,maxdist=0.4) 

 GeoCovariogram(res,show.vario=TRUE, vario=vario,pch=20)

##########  kriging ###############
 xx=seq(0,1,0.012) 
 loc_to_pred=as.matrix(expand.grid(xx,xx)) 
 Nloc=nrow(loc_to_pred)
 Xloc=cbind(rep(1,Nloc),runif(Nloc))


param_est=as.list(c(fit2$param,fixed1))
pr=GeoKrig(data=data, coordx=coords,loc=loc_to_pred, X=X,Xloc=Xloc,
     corrmodel=corrmodel,model=model,mse=TRUE,param= param_est)

colour = rainbow(100)
#### map of  data
quilt.plot(coords[,1], coords[,2], data,col=colour,main="Data")  

# linear kriging
map=matrix(pr$pred,ncol=length(xx))

image.plot(xx,xx,map,col=colour,xlab="",ylab="",main="SimpleKriging")

#associated mean squared error
map_mse=matrix(pr$mse,ncol=length(xx))

image.plot(xx,xx,map_mse,col=colour,xlab="",ylab="",main="MSE")

