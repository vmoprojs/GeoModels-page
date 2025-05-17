rm(list=ls())
library(GeoModels)
library(fields)
require(limma)
require(oz)
library(maps)
require(mapdata)
require(geoR)
require(sf)
library(dplyr)
#########################


data(austemp)
head(austemp)


coords <- cbind(austemp[,1],austemp[,2])
temp <- austemp[,3]



plot(temp,austemp[,4],pch=20,xlab="Geom. temperature",ylab="Max. temperature")



radius <- 6371
distance <- "Geod"


############# preliminary aanalysis  #####

quilt.plot(coords,temp,xlab="long",ylab="lat")
oz(states=FALSE,add=T, lwd=2) #




hist(temp,main="Histogram of Maximum temperature",nclass=13)

GeoScatterplot(data=temp,coordx=coords,distance="Geod",maxdist=500,numbins=6)

maxdist <- max(rdist.earth(coords,miles=F,R=radius))

vario <- GeoVariogram(coordx=coords,data=temp,distance="Geod",maxdist=maxdist/5,radius=radius)
plot(vario, xlab='GC distance', ylab=expression(gamma(GC)),cex=1.5,cex.lab=1.5,
ylim=c(0, max(vario$variograms)), pch=20,
main="Semi-variogram")

########################################



corrmodel <- "Matern"
CorrParam(corrmodel)
NuisParam("Gaussian",num_betas=2)
NuisParam("StudentT",num_betas=2)
NN <- nrow(coords)
X <- cbind(rep(1,NN),austemp[,4])






############################################################################
##### starting values #########
mean <- 7.5
mean1 <- 1
sill <- 8; nugget <- 0

smooth <- 0.5;scale <- 60; 

###############################

optimizer <- "nlminb"     # 

fixed1 <- list(nugget=0,smooth=smooth)

I <- Inf

start1 <- list(mean=mean,mean1=mean1,scale=scale,sill=sill)

lower1 <- list(mean=-I,mean1=-I,scale=0,sill=0)
upper1 <- list(mean=I,mean1=I,scale=I,sill=I)
##############################################
########### Gaussian case ####################
##############################################

##### pl estimation
fit2 <- GeoFit2(data=temp,coordx=coords,corrmodel=corrmodel,X=X,model="Gaussian",
                    neighb=5,distance=distance,radius=radius,sensitivity=TRUE,
                       optimizer=optimizer, lower=lower1,upper=upper1,
                       likelihood="Marginal",type='Pairwise',
                    start=start1,fixed=fixed1)
print(fit2)
##############################################
########### Student case  ####################
##############################################
df <- 4
start1 <- list(mean=mean,mean1=mean1, scale=scale,df=1/df,sill=sill)

lower2 <- list(mean=-I,mean1=-I,scale=0,sill=0,df=0)
upper2 <- list(mean=I,mean1=I,scale=I,sill=I,df=0.5)

fit3 <- GeoFit2(data=temp,coordx=coords,corrmodel=corrmodel,model="StudentT",
                   neighb=5,distance=distance,radius=radius,X=X,sensitivity=TRUE,
                    optimizer=optimizer,lower=lower2,upper=upper2,
                likelihood="Marginal",type='Pairwise',
                    start=start1,fixed=fixed1)       
print(fit3)
DF <- as.numeric(round(1/unlist(fit3$param)['df']))
if(DF==2) {DF <- 3}
print(DF)


fixed <- list(nugget=nugget,df=1/DF,smooth=smooth)
start <- list(mean=mean,mean1=mean1, scale=scale,sill=sill)
fit4 <- GeoFit2(data=temp,coordx=coords,corrmodel=corrmodel,X=X,optimizer=optimizer,
                    neighb=5,distance=distance,radius=radius, sensitivity=TRUE,
                       lower=lower1,upper=upper1,
                     likelihood="Marginal",type='Pairwise',
                    start=start,fixed=fixed, model="StudentT")
                   
print(fit4)

fit2 #Gaussian
fit4#T


##############################################
### computing stderr  estimation with psarametric  bootstrap ############
##############################################
KK <- 100
set.seed(9)
v1 <- GeoVarestbootstrap(fit2,K=KK,optimizer=optimizer,lower=lower1,upper=upper1)#Gaussian
set.seed(9)
v2 <- GeoVarestbootstrap(fit4,K=KK,optimizer=optimizer,lower=lower1,upper=upper1)#T


v1$stderr;v1$claic;v1$clbic
v2$stderr;v2$claic;v2$clbic




##############################################
### computing estimated residulas ############
##############################################
res_g <- GeoResiduals(fit2)
res_t <- GeoResiduals(fit4)

GeoQQ(res_g)
GeoQQ(res_t)

########### some diagnostic graphics ############################################

varionorm <- GeoVariogram(data=res_g$data,coordx=coords,radius=radius, maxdist=maxdist/5,distance=distance);
GeoCovariogram(res_g,show.vario=TRUE, vario=varionorm,pch=20,ylim=c(0,1.5))

### semi-variogram residuals of t Random fields
variot <- GeoVariogram(data=res_t$data,coordx=coords,maxdist=maxdist/5, distance=distance ,radius=radius);
GeoCovariogram(res_t,show.vario=TRUE, vario=variot,pch=20,ylim=c(0,2.5));




#####  prediction at one location point
coords_to_pred <- matrix(c(135,-25),ncol=2)
## prediction of residuals
pr_student2 <- GeoKrig(data=res_t$data, coordx=coords,loc=coords_to_pred,corrmodel=corrmodel,
    distance=distance,radius=radius,mse=TRUE,
   model="StudentT",param= append(res_t$param,res_t$fixed))

pr_student2$pred
pr_student2$mse


####  prediction performance trough cross validation
KK <- 100
e <- GeoCV(fit2,K=KK,n.fold=0.2,estimation=TRUE)  # Gaussian model
d <- GeoCV(fit4,K=KK,n.fold=0.2,estimation=TRUE)  # t model


mean(e$rmse);mean(e$mae); #  gaussian
mean(d$rmse);mean(d$mae); #  T

################################################################################
####################  map prediction #############################
################################################################################
coord <- maps:::map.poly("worldHires", "Australia", exact=FALSE,
                         xlim=c(110,160),ylim=c(-45,-5),
                         boundary=FALSE,
                         interior=TRUE, fill=FALSE, as.polygon=TRUE)
coord.sp <- data.frame(x = coord$x,y = coord$y,names = coord$names)

coord.sp <- cc %>%
  st_as_sf(coords = c("x", "y")) %>%
  group_by(names) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")

long1 <- 110;long2 <-154
lat1 <- -39;lat2 <- -12

lat_seq <- seq(lat1,lat2,0.35)
lon_seq <- seq(long1,long2,0.35)
coords_tot <- as.matrix(expand.grid(lon_seq,lat_seq))
coords_tot <- sf::st_as_sf(data.frame(coords_tot), coords=c("Var1","Var2"))

gr.in <- st_intersection(coords_tot, coord.sp)
gr.in <- st_coordinates(gr.in)

plot(gr.in)




pr_gaussian <- GeoKrig(data=res_g$data, coordx=coords,loc=gr.in,
  corrmodel=corrmodel,distance=distance,radius=radius,
  mse=TRUE,model="Gaussian",param= as.list(c(res_g$param,res_g$fixed)))

pr_student2 <- GeoKrig(data=res_t$data, coordx=coords,loc=gr.in,
  corrmodel=corrmodel,distance=distance,radius=radius,mse=TRUE,
  model="StudentT",param= as.list(c(res_t$param,res_t$fixed)))

quilt.plot(gr.in,pr_gaussian$pred)
oz(states=FALSE,add=T, lwd=2)
quilt.plot(gr.in,pr_gaussian$mse)
oz(states=FALSE,add=T, lwd=2)
quilt.plot(gr.in,pr_student2$pred)
oz(states=FALSE,add=T, lwd=2)
quilt.plot(gr.in,pr_student2$mse)
oz(states=FALSE,add=T, lwd=2)
