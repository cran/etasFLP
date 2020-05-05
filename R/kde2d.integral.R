kde2d.integral <-
function(xkern,ykern,gx=xkern,gy=ykern,eps=0,factor.xy=1,
           h  = c(   bwd.nrd(xkern,w),bwd.nrd(ykern,w)),w=replicate(length(xkern),1),
         hvarx=replicate(length(xkern),1),hvary=replicate(length(xkern),1)
           ){
            eps.x   =diff(range(gx))*eps/2
            eps.y   =diff(range(gy))*eps/2
            rx  =c(min(gx)-eps.x,max(gx)+eps.x)
            ry  =c(min(gy)-eps.y,max(gy)+eps.y)
                hx  <- factor.xy*h[1]
                hy  <- factor.xy*h[2]
            nkern=length(xkern)
      ix=pnorm(rx[2],xkern,hx*hvarx)-pnorm(rx[1],xkern,hx*hvarx)
        iy=pnorm(ry[2],ykern,hy*hvary)-pnorm(ry[1],ykern,hy*hvary)
        integral    =   sum(w*ix*iy)/sum(w)
        
            return(integral)
         }
