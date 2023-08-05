rm(list=ls());
require(devtools);
#install_github("vmoprojs/GeoModels");
require(GeoModels);
require(fields);

 ###################################
#  stationary case  ####
###################################


N=1000;
set.seed(269);
coords=cbind(runif(N),runif(N));
#plot(coords ,pch=20,xlab="",ylab="");



# correlation parameters
corrmodel = "Matern"; 
scale = 0.25/3; 
smooth =0.5;
nugget =0;

n=10 


# mean parameter
mean = 0.5 # 
param=list(nugget=nugget,mean=mean, scale=scale, smooth=smooth);

data_s <- GeoSim(coordx=coords ,corrmodel=corrmodel , 
        param=param ,model="Binomial",n=n)$data;




mean(data_s);var (data_s)
p=pnorm(mean)
n*p;n*p*(1-p)


plot(table(data_s),ylab = "Frequency")



quilt.plot(coords,data_s,nlevel=n+1,zlim=c(0,n))



fit = GeoVariogram(coordx=coords,data=data_s,maxdist=0.6)
plot(fit, xlab='h', ylab=expression(gamma(h)),
     ylim=c(0, max(fit$variograms)), pch=20,
     main="Semi-variogram")





###################################
#  nonstationary case  #
###################################
mean = 0.3 # regression paramteres 
mean1= -0.25

set.seed(29);
a0=rep(1,N);a1=runif(N)
X=cbind(a0,a1); ## regression matrix

n=sample(10:20,nrow(coords),replace=TRUE)
head(n)


param=list(nugget=nugget,mean=mean,mean1=mean1, scale=scale, 
          smooth=smooth)

data_ns<- GeoSim(coordx=coords ,corrmodel=corrmodel , param=param ,n=n,
               X=X,model="Binomial")$data;





###################################
# estimation    #
###################################


fixed2<-list(nugget=0,smooth=smooth);
start2<-list(mean=mean,mean1=mean1,scale=scale);
lower<-list(mean=-5,mean1=-5,scale=0);
upper<-list(mean=5,mean1=5,scale=10);



fit1_ns<- GeoFit2(data=data_ns,coordx=coords,corrmodel=corrmodel,
likelihood="Conditional",type="Pairwise",
n=n, X=X, neighb=4,sensitivity=TRUE,
optimizer="nlminb",lower=lower,upper=upper,
start=start2,fixed=fixed2, model="Binomial");

fit1_ns

fit1_ns_sd=GeoVarestbootstrap(fit1_ns,K=100,optimizer="nlminb",lower=lower, upper=upper)

fit1_ns_sd$stderr
###########################################################################################

###################################
# prediction  nonstationary case  #
###################################

xx=seq(0,1,0.02) 
loc_to_pred=as.matrix(expand.grid(xx,xx)) 
NN=nrow(loc_to_pred)

a0=rep(1,NN);a1=runif(NN)
Xloc=cbind(a0,a1); 

nloc=sample(10:20,nrow(loc_to_pred),replace=TRUE)
##### binomial 
param_est=as.list(c(fit1_ns$param,fixed2))
pr=GeoKrig(data=data_ns, coordx=coords,loc=loc_to_pred,X=X,Xloc=Xloc,
     n=n,nloc=nloc,corrmodel=corrmodel,model="Binomial",mse=TRUE,param= param_est)


par(mfrow=c(1,3))
#### map of  data
nn=max(n)

quilt.plot(coords, data_ns,nlevel=nn+1,zlim=c(0,nn)) 

# map predictions
map=matrix(pr$pred,ncol=length(xx))

image.plot(xx,xx,map,xlab="",zlim=c(0,nn),ylab="",main="Kriging")

#map MSE
map_mse=matrix(pr$mse,ncol=length(xx))

image.plot(xx,xx,map_mse,xlab="",ylab="",main="MSE")

