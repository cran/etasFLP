cat.select <- function(etas.obj,longlat.to.km,sectoday){
 
# cat 	=   	subset(etas.obj$cat.orig,(etas.obj$cat.orig$magn1>=etas.obj$magn.threshold)&(etas.obj$cat.orig$time<=etas.obj$tmax))
 
 cat 	=   	subset(etas.obj$cat.orig,(etas.obj$cat.orig$magn1>=etas.obj$magn.threshold))

 if (is.na(etas.obj$tmax))etas.obj$tmax=max(cat$time)

 cat 	=   	subset(cat,cat$time<=etas.obj$tmax)

        unit=1.
  	if(sectoday) unit=86400
        time.work	=	cat$time/unit
	etas.obj$tmax		=	etas.obj$tmax/unit
	

  		
		if (longlat.to.km){
		
		radius=6371.3
		

		ycat.work      =   radius*cat$lat*pi/180
		xcat.work      =   radius*cat$long*pi/180

		
		  }
	        else
	        {
 
		ycat.work      =  cat$lat
		xcat.work      =   cat$long
}
magn1.work      =cat$magn1-etas.obj$magn.threshold
cat             =  as.data.frame( cbind(cat,time.work,xcat.work,ycat.work,magn1.work))
ord	<-	order(cat$time.work)
cat 	<-	cat[ord,]

etas.obj$cat    =   cat
return(etas.obj)

}