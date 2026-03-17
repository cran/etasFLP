cat.select <- function(etas.obj,longlat.to.km,sectoday){
    
    
    time.work	=	etas.obj$cat.orig$time
    if(sectoday)  time.work	=	time.work/86400.
    
    #    etas.obj$tmax		=	etas.obj$tmax/unit
    
    
    
    unit=1.
    if (longlat.to.km) unit=6371.3*pi/180
    
    ycat.work      =  etas.obj$cat.orig$lat*unit
    xcat.work      =  etas.obj$cat.orig$long*unit
    
    magnitude      =etas.obj$cat.orig$magn1-etas.obj$magn.threshold
    
    
    
    ind.magn1=(etas.obj$cat.orig$magn1>=etas.obj$magn.threshold)
    
    if (is.na(etas.obj$tmax))etas.obj$tmax=max(etas.obj$cat.orig$time[ind.magn1])
    if (is.na(etas.obj$long.range[1]))etas.obj$long.range=range(etas.obj$cat.orig$long[ind.magn1])
    if (is.na(etas.obj$lat.range[1]))etas.obj$lat.range=range(etas.obj$cat.orig$lat[ind.magn1])
    
    ind.time=	(etas.obj$cat.orig$time<=etas.obj$tmax)
    ind.long=	(etas.obj$cat.orig$long<=etas.obj$long.range[2])&(etas.obj$cat.orig$long>=etas.obj$long.range[1])
    ind.lat=	(etas.obj$cat.orig$lat<=etas.obj$lat.range[2])&(etas.obj$cat.orig$lat>=etas.obj$lat.range[1])
    
    ind=as.logical(ind.magn1*ind.time*ind.long*ind.lat)
    ord	<-	order(time.work)
    cat =  data.frame( etas.obj$cat.orig,time.work,xcat.work,ycat.work,magnitude,ind,ord)

    etas.obj$cat    =   cat
    return(etas.obj)
    
}