rm(list=ls());
require(GeoModels);
require(fields);
model="Poisson";  # model name in the  GeoModels package


###################################
# simulation  stationary case  ####
###################################

 set.seed(1989);
 N=500;
 coords=cbind(runif(N),runif(N));
 plot(coords ,pch=20,xlab="",ylab="");

# correlation parameters
corrmodel = "Matern"; 
scale = 0.25/3; 
smooth =0.5;
nugget =0;

# mean parameter
mean = 1.5 # regression paramteres 
param=list(nugget=nugget,mean=mean, scale=scale, smooth=smooth, sill=1);
data_s <- GeoSim(coordx=coords ,corrmodel=corrmodel , param=param ,model=model)$data;

mean(data_s);var (data_s)
exp(mean)


plot(table(data_s),ylab = "Frequency")


###################################
# simulation  nonstationary case  #
###################################
corrmodel = "Wend0";        ## correlation model 
scale = 0.2;               ## scale parameter
power2=4 ;                  ## power parameter
nugget=0;                    ## nugget parameter


mean = 1.5 # regression paramteres 
mean1= -0.25
a0=rep(1,N);a1=runif(N)
X=cbind(a0,a1); ## regression matrix


param=list(nugget=nugget,mean=mean,mean1=mean1, scale=scale, 
		power2=power2, sill=1);
data_ns <- GeoSim(coordx=coords ,corrmodel=corrmodel , param=param ,
			X=X,model=model)$data;


###################################
# estimation  stationary case  #
###################################
optimizer="nlminb";

fixed1<-list(sill=1,nugget=0,smooth=0.5);
start1<-list(mean=1.5,scale=0.25/3);

lower<-list(mean=-5,scale=0);
upper<-list(mean=5,scale=2);

maxdist=0.03;
corrmodel = "Matern"; 
## pairwise poisson likelihhod
fit1 <- GeoFit(data=data_s,coordx=coords,corrmodel=corrmodel,
optimizer=optimizer,lower=lower,upper=upper,
maxdist=maxdist,start=start1,fixed=fixed1, model = model);

## misspecified gaussian likelihhod
fit2 <- GeoFit(data=data_s,coordx=coords,corrmodel=corrmodel,
 optimizer=optimizer, lower=lower,upper=upper,maxdist=maxdist,
start=start1,fixed=fixed1, model = "Gaussian_misp_Poisson")

fit1$param
fit2$param
###########################################################################################
# Empirical estimation of the variogram:
vario <- GeoVariogram(data=data_s,coordx=coords,maxdist=0.4)
# comparing empirical and estimated semi-variograms 
GeoCovariogram(fit1,show.vario=TRUE, vario=vario,pch=20)




###################################
# estimation  nonstationary case  #
###################################
optimizer="nlminb";
corrmodel = "Wend0"; 
fixed2<-list(sill=1,nugget=0,power2=4);
start2<-list(mean=1.5,mean1=-0.25,scale=0.2);

lower<-list(mean=-5,mean1=-5,scale=0);
upper<-list(mean=5,mean1=5,scale=2);


fit1_ns <- GeoFit(data=data_ns,coordx=coords,corrmodel=corrmodel,
optimizer=optimizer,
lower=lower,upper=upper,X=X,
maxdist=maxdist,start=start2,fixed=fixed2, model = model);

###########################################################################################

fit2_ns <- GeoFit(data=data_ns,coordx=coords,corrmodel=corrmodel,
optimizer=optimizer,
lower=lower,upper=upper,X=X,
maxdist=maxdist,start=start2,fixed=fixed2, model = "Gaussian_misp_Poisson")

fit1_ns$param
fit2_ns$param





###################################
# prediction  stationary case  #
###################################
xx=seq(0,1,0.015) 
loc_to_pred=as.matrix(expand.grid(xx,xx)) 
param_est=as.list(c(fit1$param,fixed1))
corrmodel = "Matern"; 
pr=GeoKrig(data=data_s, coordx=coords,loc=loc_to_pred,
     corrmodel=corrmodel,model=model,mse=TRUE,param= param_est)
colour = rainbow(100)
#### map of  data
quilt.plot(coords[,1], coords[,2], data_s,col=colour,main="Data")  
#### map prediction
map=matrix(pr$pred,ncol=length(xx))
image.plot(xx,xx,map,col=colour,xlab="",ylab="",main="Kriging")
#map mean squared error
map_mse=matrix(pr$mse,ncol=length(xx))
image.plot(xx,xx,map_mse,col=colour,xlab="",ylab="",main="MSE")




###################################
# prediction  nonstationary case  #
###################################
set.seed(609)
NN=nrow(loc_to_pred)
a0=rep(1,NN);a1=runif(NN)
Xloc=cbind(a0,a1); ## 
corrmodel = "Wend0"; 
param_est=as.list(c(fit1_ns$param,fixed2))
pr=GeoKrig(data=data_ns, coordx=coords,loc=loc_to_pred,X=X,Xloc=Xloc,
     corrmodel=corrmodel,model=model,mse=TRUE,param= param_est)

colour = rainbow(100)
#### map of  data
quilt.plot(coords[,1], coords[,2], data_ns,col=colour,main="Data")  
# map predictions
map=matrix(pr$pred,ncol=length(xx))
image.plot(xx,xx,map,col=colour,xlab="",ylab="",main="Kriging")
#map MSE
map_mse=matrix(pr$mse,ncol=length(xx))
image.plot(xx,xx,map_mse,col=colour,xlab="",ylab="",main="MSE")

