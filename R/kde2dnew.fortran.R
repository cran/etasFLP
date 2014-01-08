kde2dnew.fortran <-
function(xkern,ykern,gx,gy,h,factor.xy=1,w=replicate(length(xkern),1)){
#
#
# xkern, ykern    x-y coordinates of the   points used in the kernel estimate with  w weights  
# gx, gy    x-y coordinates of the   points where the estimate must be computed
# 
# interface for FORTRAN subroutinedensity2	
#
            nx  <-  length(xkern)
            if(length(ykern)!=nx)   stop("Data vectors must be the same length")
            if(missing(h))      h<-c(bwd.nrd(xkern,w),bwd.nrd(ykern,w))

            h   <-  h*factor.xy
            h   <-  ifelse(h==0,max(h),h)
            n   <-  length(gx)

ris2d<-.Fortran("density2",x=as.double(gx),y=as.double(gy),m=as.integer(n),xkern=as.double(xkern),ykern=as.double(ykern),nkern=as.integer(nx),h=as.double(h),w=as.double(w),dens=as.double(gx))

	return(list(x=gx,y=gy,z=ris2d$dens,h=h))
                        }
