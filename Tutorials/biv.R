rm(list=ls())
require(devtools)
#install_github("vmoprojs/GeoModels") 
require(GeoModels)
require(fields)
model="Gaussian" # model name in the GeoModels package 
set.seed(24)


NN=500 # number of spatial locations 
x = runif(NN, 0, 1); 
y = runif(NN, 0, 1) 
coords=cbind(x,y) 



## setting parameters
mean_1 = 2; mean_2= -1
nugget_1 =0;nugget_2=0
sill_1 =0.5; sill_2 =1; 

### correlation parameters
corrmodel="Bi_Matern"
CorrParam("Bi_Matern")

scale_1=0.2/3; scale_2=0.15/3; scale_12=0.5*(scale_2+scale_1) 
smooth_1=smooth_2=smooth_12=0.5
pcol = -0.4

param= list(nugget_1=nugget_1,nugget_2=nugget_2,
	        sill_1=sill_1,sill_2=sill_2,
	        mean_1=mean_1,mean_2=mean_2,
	        smooth_1=smooth_1, smooth_2=smooth_2,smooth_12=smooth_12,
	        scale_1=scale_1, scale_2=scale_2,scale_12=scale_12,
	        pcol=pcol)





## simulation
ss1 = GeoSim(coordx=coords, corrmodel=corrmodel,model=model,param=param)$data
dim(ss1)
cc = GeoCovmatrix(coordx=coords, corrmodel=corrmodel,
	          model=model,param=param)
cc$covmatrix[1:4,1:4]





fixed=list(nugget_1=nugget_1,nugget_2=nugget_2, 
	     smooth_1=smooth_1, smooth_2=smooth_2,smooth_12=smooth_12)

start=list( mean_1=mean_1,mean_2=mean_2, sill_1=sill_1,sill_2=sill_2,
	      scale_1=scale_1,scale_2=scale_2,scale_12=scale_12, pcol=pcol)

## estimation with maximum likelihood 
fit = GeoFit(data=ss1,coordx=coords, corrmodel=corrmodel,
	likelihood="Full",type="Standard",optimizer="BFGS",
	start=start,fixed=fixed)

fit



fit_pl = GeoFit(data=ss1,coordx=coords, corrmodel=corrmodel,
 maxdist=c(0.1,0.1,0.1),likelihood="Marginal",type="Pairwise",
 optimizer="BFGS", start=start, fixed=fixed)

fit_pl

res=GeoResiduals(fit_pl)

GeoQQ(res)


### checking model assumptions: ST variogram model
vario = GeoVariogram(data=ss1,coordx=coords, bivariate=TRUE, maxdist=c(0.4,0.4,0.4))
GeoCovariogram(fit_pl,vario=vario,show.vario=TRUE,pch=20)



####################
##### kriging ###### 
####################
xx=seq(0,1,0.012)
loc_to_pred=as.matrix(expand.grid(xx,xx))


param_est=as.list(c(fit_pl$param,fixed))
pr1 = GeoKrig(data=ss1,coordx=coords,  corrmodel=corrmodel,which=1,
	          model=model,mse=TRUE,loc=loc_to_pred,param=param_est)

pr2 = GeoKrig(data=ss1,coordx=coords,  corrmodel=corrmodel,which=2,
	          model=model,mse=TRUE,loc=loc_to_pred,param=param_est)


par(mfrow=c(2,3))
colour <- rainbow(100)
quilt.plot(coords[,1],coords[,2],ss1[1,],col=colour,main ="Observed data: First variable")
image.plot(xx, xx, matrix(pr1$pred,ncol=length(xx)),col=colour, main = paste(" Kriging "),ylab="")
image.plot(xx, xx, matrix(pr1$mse,ncol=length(xx)),col=colour,
           main = paste("MSE"),ylab="")


quilt.plot(coords[,1],coords[,2],ss1[2,],col=colour,main ="Observed data: Second variable")
image.plot(xx, xx, matrix(pr2$pred,ncol=length(xx)),col=colour, main = paste(" Kriging "),ylab="")
image.plot(xx, xx, matrix(pr2$mse,ncol=length(xx)),col=colour,
           main = paste("MSE"),ylab="")




