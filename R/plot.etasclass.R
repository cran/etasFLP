plot.etasclass <-
function(x,pdf=FALSE,file ="etasplot",  ngrid=201,tfixed=0,flag.3D=FALSE,flag.log=FALSE,...){
#  function(x,pdf=FALSE,file ="etasplot",  ngrid=201,flag.3D=FALSE,flag.log=FALSE,ellipse=FALSE,...){
	if(class(x)!="etasclass")stop("argument must be an etasclass object")
 
 ellipse=FALSE
# kern.var=x$kern.var
 kern.var=FALSE
# which=1:4
# maxgraph	=7
# ind.graph	=array(FALSE,maxgraph)
# ind.graph[which]=TRUE
 

#tfixed=0 add in the rd file
cat("Computation of the ETAS space intensity integrated intensity on a grid","\n")
cat("for large catalogs computation can take some minutes. Please wait ...","\n")
#####################################################
# plot time intensity # with points
typegraph=1

file=paste(file,".pdf",sep="")
if (pdf){
pdf(file=file,onefile=TRUE)}
else
{
dev.new()
}
#####################################################
### computation of back intensity and triggered intensity on a given x-y grid
### back.grid computed with the same kernel used for back.dens on observed  points
## back.grid computed on km coordinates

if(!x$is.backconstant){
xcat.km	=x$cat$long
ycat.km	=x$cat$lat
w	=x$rho.weights
hdef	=x$hdef
n	=length(xcat.km)
params	=x$params
        mu= params[1]
        k0= params[2]
        c= params[3]
        p= params[4]
        a= params[5]
        gamma= params[6]
        d= params[7]
        q= params[8]
#
rangex		=range(xcat.km)
rangey		=range(ycat.km)
ranget=diff(range(x$cat$time))

alpha	=0
if (kern.var) alpha=hdef[3]

space.grid	=xy.grid(rangex,rangey,ngrid)
ngridtot	=length(space.grid[,1])
#kde=kde2dnew.fortran(xcat.km,ycat.km,space.grid[,1],space.grid[,2],w=w,factor.xy=1,h=hdef,kern.var=kern.var,alpha=alpha)
kde=kde2dnew.fortran(xcat.km,ycat.km,space.grid[,1],space.grid[,2],w=w,factor.xy=1,h=hdef)
wmat1		=matrix(kde$wmat,n,4)
back.grid	=ranget*mu*kde$z
	
#### computation MUST be in km IN THIS VERSION because parameters are evaluated in km
ris=.Fortran("etasfull8tintegrated" ,NAOK=TRUE,
			n=as.integer(n),
			mu=as.double(mu),k=as.double(k0),
			c=as.double(c),p=as.double(p),
			a=as.double(a),g=as.double(gamma),
			d=as.double(d),q=as.double(q),
			x=as.double(xcat.km),y=as.double(ycat.km), t=as.double(x$cat$time),m=as.double(x$cat$magn1),l=back.grid,
			ngridtot=as.integer(ngridtot), xgrid=as.double(space.grid[,1]),
			ygrid=as.double(space.grid[,2]),
			tmax=as.double(max(x$cat$time)))
			trig.grid	=ris$l
				
### trig.grid intensity on a x-y grid   by time integration of ETAS intensity function
###########################################################################################
}
tot.grid=back.grid+trig.grid

# maps of triggered intensity for a single day

if(tfixed>0){
ind=x$cat$time<tfixed
ris=.Fortran("etasfull8tfixed" ,NAOK=TRUE,
			n=as.integer(sum(ind)),
			mu=as.double(mu),k=as.double(k0),
			c=as.double(c),p=as.double(p),
			a=as.double(a),g=as.double(gamma),
			d=as.double(d),q=as.double(q),
			x=as.double(xcat.km[ind]),y=as.double(ycat.km[ind]), t=as.double(x$cat$time[ind]),m=as.double(x$cat$magn1[ind]),l=back.grid,
			ngridtot=as.integer(ngridtot), xgrid=as.double(space.grid[,1]), 
			ygrid=as.double(space.grid[,2]),tfixed=as.double(tfixed))
		        totfixed.grid	=ris$l+back.grid/ranget

if(flag.log)totfixed.grid=log(totfixed.grid)
}

if(flag.log) {
    back.grid=log(back.grid)
    trig.grid=log(trig.grid)
     tot.grid=log( tot.grid)
}
## change grid to degrees
rangex		=range(x$cat.longlat$long)
rangey		=range(x$cat.longlat$lat)
space.grid	=xy.grid(rangex,rangey,ngrid)
ngridtot	=length(space.grid[,1])

			x.grid		=seq(rangex[1],rangex[2],length=ngrid)
			y.grid		=seq(rangey[1],rangey[2],length=ngrid)

if(tfixed>0){
### start triggered intensity plotting for a tfixed ########################################
if(!pdf) dev.new()
ind1=x$cat.longlat$time>=tfixed & x$cat.longlat$time<(tfixed+1)

mapxy=map("worldHires",xlim=range(x$cat.longlat$long),ylim=range(x$cat.longlat$lat),plot=FALSE)

image.plot(x.grid,y.grid,(matrix(totfixed.grid,ngrid,ngrid))
,col=gray.colors(128, start = 0., end = 1., gamma =2 )
,xlab="x-longitude",ylab="y-latitude"
,main=paste("Triggered intensity and observed points at day ",tfixed)
)

grid(col="grey")

map("worldHires",add=TRUE,xlab="x-longitude",ylab="y-latitude",
xlim=range(x$cat.longlat$long),ylim=range(x$cat.longlat$lat),col="green"
)
contour(x.grid,y.grid,(matrix(totfixed.grid,ngrid,ngrid)),col="red",add=TRUE)

points(x$cat.longlat$long[ind1],x$cat.longlat$lat[ind1],cex=sqrt(exp(x$cat.longlat$magn1[ind1]))/8,col=4,pch=19)


### end triggered intensity plotting for a tfixed ########################################
}

if(!pdf) dev.new()


### start triggered intensity plotting

mapxy=map("worldHires",xlim=range(x$cat.longlat$long),ylim=range(x$cat.longlat$lat),plot=FALSE)

image.plot(x.grid,y.grid,(matrix(trig.grid,ngrid,ngrid))
,col=gray.colors(128, start = 0., end = 1., gamma =2 )
,xlab="x-longitude",ylab="y-latitude"
,main="Triggered Intensity"
)

grid(col="grey")

map("worldHires",add=TRUE,xlab="x-longitude",ylab="y-latitude",
xlim=range(x$cat.longlat$long),ylim=range(x$cat.longlat$lat),col="green"
)
contour(x.grid,y.grid,(matrix(trig.grid,ngrid,ngrid)),col="red",add=TRUE)

### end triggered intensity plotting

### start background intensity plotting

box()

if(!pdf) dev.new()

image.plot(x.grid,y.grid,(matrix(back.grid,ngrid,ngrid))
,col=gray.colors(128, start = 0., end = 1., gamma =2 ),
xlab="x-longitude",ylab="y-latitude",
main="Background Intensity"
)
      
grid(col="grey")
map("worldHires",add=TRUE,xlab="x-longitude",ylab="y-latitude",
xlim=range(x$cat.longlat$long),ylim=range(x$cat.longlat$lat),col="green")

contour(x.grid,y.grid,(matrix(back.grid,ngrid,ngrid)),col="red",add=TRUE)

box()

### start total intensity plotting

if(!pdf) dev.new()

image.plot(x.grid,y.grid,(matrix(tot.grid,ngrid,ngrid)),
col=gray.colors(128, start = 0., end = 1., gamma =2 ),
xlab="x-longitude",ylab="y-latitude",main="Total Intensity"
)
      
grid(col="grey")
map("worldHires",add=TRUE,xlab="x-longitude",ylab="y-latitude",
xlim=range(x$cat.longlat$long),ylim=range(x$cat.longlat$lat),col="green")
contour(x.grid,y.grid,(matrix(tot.grid,ngrid,ngrid)),col="red",add=TRUE)

box()
### start total intensity plotting with observed points
ts=(x$cat$time-min(x$cat$time))/diff(range(x$cat$time))

if(!pdf) dev.new()

image.plot(x.grid,y.grid,
(matrix(tot.grid,ngrid,ngrid)),
col=gray.colors(128, start = 0., end = 1., gamma =2 ),
xlab="x-longitude",ylab="y-latitude",main="Total Intensity with observed points \n Circles area proportional to magnitude; red: recent, blu:older"
)
      
grid(col="grey")
map("worldHires",add=TRUE,xlab="x-longitude",ylab="y-latitude",
xlim=range(x$cat.longlat$long),ylim=range(x$cat.longlat$lat),col="green")
points(x$cat.longlat$long,x$cat.longlat$lat,cex=sqrt(exp(x$cat.longlat$magn1))/8,col=rgb(ts,0,1-ts),pch=19)
contour(x.grid,y.grid,(matrix(tot.grid,ngrid,ngrid)),col="yellow",add=TRUE)

box()



#####################################################
if((as.numeric(kern.var)*as.numeric(ellipse))==1){

if(!pdf) dev.new()


wx=wmat1[,2]^2
wy=wmat1[,3]^2
wxy=wmat1[,4]*wmat1[,2]*wmat1[,3]
x=xcat.km
y=ycat.km
 plot(x,y,type="p",pch=22,col=1,main="adaptive kernel")
  for (j in 1:n)
  cat(j)
  # internal ellipses
  lines(ellipse(matrix(c(wx[j],wxy[j],wxy[j],wy[j]),2,2),centre=c(x[j],y[j]),level=c(0.500)),col=2)
  # external ellipses
  for (j in 1:n) lines(ellipse(matrix(c(wx[j],wxy[j],wxy[j],wy[j]),2,2),centre=c(x[j],y[j]),level=c(0.999)),col=3)
  }

 dev.off()

etas.l=x$l

if(flag.log) etas.l=log(etas.l)

if(flag.3D){

typegraph=2
plot3d(x$cat.longlat$long,x$cat.longlat$lat,etas.l,type="n",zlab=paste("estimated intensity   ","lambda(x,y)"),
xlab="x-longitude",ylab="y-latitude",main="Estimated intensities in observed points")
lines3d(cbind(mapxy$x,mapxy$y,min(etas.l)),col="red")
plot3d(x$cat.longlat$long,x$cat.longlat$lat,etas.l,add=TRUE, type="h")

}

return(list(
x.grid=x.grid,
y.grid=y.grid,
back.grid=back.grid,
trig.grid=trig.grid,
tot.grid=tot.grid,
#back.obs=back.obs,
#trig.obs=trig.obs,
#tot.obs=tot.obs,
ranget=ranget
#,wmat=wmat1
))
}


