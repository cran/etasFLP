kde2dnew.fortran <-
function(xkern,ykern,gx,gy,h,factor.xy=1,eps=0,w=replicate(length(xkern),1),
hvarx=replicate(length(xkern),1),hvary=replicate(length(xkern),1)
)
{
#
#
# xkern, ykern    x-y coordinates of the   points used in the kernel estimate with  w weights  
# gx, gy    x-y coordinates of the   points where the estimate must be computed
# 
# interface for FORTRAN subroutinedensity2	
#
           n   <-  length(gx)
           nkern  <-  length(xkern)
            if(length(ykern)!=nkern)   stop("Data vectors must be the same length")
            if(missing(h))      h<-c(bwd.nrd(xkern,w),bwd.nrd(ykern,w))

            h   <-  h*factor.xy
            h   <-  ifelse(h==0,max(h),h)

        ris2d<-.Fortran("density2serial",x=as.double(gx),y=as.double(gy),m=as.integer(n),xkern=as.double(xkern),ykern=as.double(ykern),nkern=as.integer(nkern),h=as.double(h),w=as.double(w),
        hvarx=as.double(hvarx),
        hvary=as.double(hvary),

dens=as.double(gx))


#integral=kde2d.integral(xkern,ykern,gx,gy,factor.xy=factor.xy,eps=eps,w=w,h=h,kern.var=kern.var,wmat=wmat)
integral=kde2d.integral(xkern,ykern,gx,gy,factor.xy=factor.xy,eps=eps,w=w,h=h)


#	print(c("kde2dnew.fortran; dens ",ris2d$dens ))        
#	print(c("kde2dnew.fortran; integral ",integral))        

	return(list(x=gx,y=gy,z=ris2d$dens,h=h,integral=integral))
                        }
                
