rm(list=ls());
require(GeoModels);
require(fields);
model="Poisson";  # model name in the  GeoModels package


###################################
# simulation  stationary case  ####
###################################

set.seed(1989);
N =500;
coords = cbind ( runif ( N ) , runif ( N ));
plot ( coords , pch =20 , xlab ="", ylab ="");


# correlation parameters
corrmodel = "Matern"; ## correlation model
scale = 0.25 /3; ## scale parameter
smooth =0.5; ## smooth parameter
nugget =0; ## nugget parameter

# mean parameter
mean = 1.5 # regression paramteres 
param=list(nugget=nugget,mean=mean, scale=scale, smooth=smooth);
data_s <- GeoSim(coordx=coords ,corrmodel=corrmodel , param=param ,model=model)$data;

mean(data_s);var (data_s)
exp(mean)


plot(table(data_s),ylab = "Frequency")


###################################
# simulation  nonstationary case  #
###################################
corrmodel = "Wend0"; ## correlation model
scale = 0.2; ## scale parameter
power2 =4 ; ## power parameter
nugget =0; ## nugget parameter

mean = 1.5 # regression parameter beta _0
mean1 = -0.25 # regression parameter beta _1
a0 = rep (1 , N ); a1 = runif ( N )
X = cbind ( a0 , a1 ); ## regression matrix


param=list(nugget=nugget,mean=mean,mean1=mean1, scale=scale, 
		power2=power2);
data_ns <- GeoSim(coordx=coords ,corrmodel=corrmodel , param=param ,
			X=X,model=model)$data;


###################################
# estimation  stationary case  #
###################################
optimizer="nlminb";

fixed1<-list(nugget=0,smooth=0.5);
start1<-list(mean=1.5,scale=0.25/3);

lower<-list(mean=-5,scale=0);
upper<-list(mean=5,scale=2);

neighb=3;
corrmodel = "Matern"; 
## pairwise poisson likelihhod
fit1 <- GeoFit(data=data_s,coordx=coords,corrmodel=corrmodel,
optimizer=optimizer,lower=lower,upper=upper,
neighb=neighb,start=start1,fixed=fixed1, model = model);

## misspecified gaussian likelihhod
fit2 <- GeoFit(data=data_s,coordx=coords,corrmodel=corrmodel,
 optimizer=optimizer, lower=lower,upper=upper,neighb=neighb,
start=start1,fixed=fixed1, model = "Gaussian_misp_Poisson")

unlist(fit1$param)
unlist(fit2$param)
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
fixed2<-list(nugget=0,power2=4);
start2<-list(mean=1.5,mean1=-0.25,scale=0.2);

lower<-list(mean=-5,mean1=-5,scale=0);
upper<-list(mean=5,mean1=5,scale=2);


fit1_ns <- GeoFit(data=data_ns,coordx=coords,corrmodel=corrmodel,
optimizer=optimizer,
lower=lower,upper=upper,X=X,
neighb=neighb,start=start2,fixed=fixed2, model = model);

###########################################################################################

fit2_ns <- GeoFit(data=data_ns,coordx=coords,corrmodel=corrmodel,
optimizer=optimizer,
lower=lower,upper=upper,X=X,
neighb=neighb,start=start2,fixed=fixed2, model = "Gaussian_misp_Poisson")

unlist(fit1_ns$param)
unlist(fit2_ns$param)





###################################
# prediction  stationary case  #
###################################
xx=seq(0,1,0.015) 
loc_to_pred=as.matrix(expand.grid(xx,xx)) 
corrmodel = "Matern"; 
pr=GeoKrig(fit1,loc=loc_to_pred,mse=TRUE)
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
pr=GeoKrig(fit1_ns,loc=loc_to_pred,Xloc=Xloc,mse=TRUE)

colour = rainbow(100)
#### map of  data
quilt.plot(coords[,1], coords[,2], data_ns,col=colour,main="Data")  
# map predictions
map=matrix(pr$pred,ncol=length(xx))
image.plot(xx,xx,map,col=colour,xlab="",ylab="",main="Kriging")
#map MSE
map_mse=matrix(pr$mse,ncol=length(xx))
image.plot(xx,xx,map_mse,col=colour,xlab="",ylab="",main="MSE")