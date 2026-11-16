
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        bivariate isotropic density estimation routines
!
!        xkern,ykern, vectors of lenght nkern (centers for  kernels)
!         x,y vectors of lenght  m (where we want to compute the density)
!         dens vector of lenght  m with the densities in x,y
!         h vector of 2 bandwidths
!         w        vector of weights
             subroutine density2serial(x,y,m,xkern,ykern,nkern,h,w,hvarx,hvary,dens)
        INTEGER nkern,m
!        input variables declaration:
      double precision xkern(nkern),ykern(nkern),x(m),y(m),w(nkern),ww,h(2)
      double precision hvarx(nkern),hvary(nkern)
!        output variables declaration:     
      double precision dens(m)
!        work variables declaration:
      double precision stdx(nkern),stdy(nkern),fac,hx(nkern),hy(nkern)
      DATA               fac/0.3989422804014327D0/ !        fac is 1/sqrt(2*PI)
        ww=sum(w)
        hx=h(1)*hvarx
        hy=h(2)*hvary
        dens=0
         do j=1,m
                stdx        =((x(j)-xkern))/hx        
                stdy        =((y(j)-ykern))/hy        
                dens(j)=sum(w*exp(-(stdx*stdx+stdy*stdy)/2))
            dens(j)=dens(j)*fac*fac/(ww*h(1)*h(2))
        end do
       end

       
     
       
       
!        etasfull8newserial computes only the triggered intensity of the ETAS model
!       with 7+ncov parametrs.
!       input magnitude vector must be  m(j)-m0
!
!       called by etas.mod2NEW.R
!
      subroutine etasfull8newserial(tflag,n,mu,k,c,p,g,d,q,x,y,t,m,predictor,l)
        INTEGER n,tflag
      double precision mu,k,c,p,g,d,q,x(n),y(n),t(n),m(n),predictor(n),l(n)
      double precision dx,dy,ds,dt,xi,yi,ti,etas,inc,dum
      double precision, allocatable :: ep(:), eg(:)
      allocate(ep(n),eg(n))
      dum=mu
      do j=1,n
        ep(j)        =exp(predictor(j))
        eg(j)        =exp(g*m(j))
      end do
      do i=2,n
        inc        =0
        ti        =t(i)
        xi        =x(i)
        yi        =y(i)
      do j=1,i-1
      dt        =ti-t(j)
        etas        =0
        if (dt > 0)  then
        dx        =xi-x(j)
        dy        =yi-y(j)
        ds        =dx*dx+dy*dy
        if (tflag>0) then
        etas        =exp(-p*log(dt+c))*ep(j)
        else
        etas        =exp(-p*log(dt+c)-q*log(ds/eg(j)+d))*ep(j)
        end if 
        end if
        inc        =inc+etas
        end do
        l(i)        =inc*k
        end do 
      deallocate(ep,eg)
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     etasfull8tintegratednew computes only the triggered intensity of the ETAS model
!     integrated on the time axis
!     on a set of points different from those observed.
!
!     input magnitude vector must be  m(j)-m0
!     27-7-2017: with covariates!!!!!
!
!     called by plot.etasclass.R
!
      subroutine etasfull8tintegratednew(n,mu,k,c,p,g,d,q,x,y,t,m,predictor,l,ngridtot, xgrid, ygrid,tmax)
      INTEGER n,ngridtot
      double precision mu,k,c,p,g,d,q,x(n),y(n),xgrid(ngridtot),ygrid(ngridtot),t(n),m(n),predictor(n),l(ngridtot)
      double precision dx,dy,ds,dt,xi,yi,etas,inc,tmax,integrt,eps,dum
      double precision, allocatable :: wt(:), eg(:)
      eps=1D-10
      dum=mu
      allocate(wt(n),eg(n))
      do j=1,n
        dt        =tmax-t(j)
        wt(j)     =0
        if (dt > 0)  then
!                check for p=1
          if (abs(p-1)<eps) then
                integrt =log(c+dt)-log(c)
          else
                integrt =((c+dt)**(1-p)- c**(1-p))/(1-p)
          end if
          wt(j)   =integrt*exp(predictor(j))
        end if
        eg(j)     =exp(g*m(j))
      end do

      do i=1,ngridtot
        inc        =0
        xi        =xgrid(i)
        yi        =ygrid(i)        
      do j=1,n
      dt        =tmax-t(j)
        etas        =0
        if (dt > 0)  then
        dx        =xi-x(j)
        dy        =yi-y(j)
        ds        =dx*dx+dy*dy
        etas        =wt(j)*exp(-q*log(ds/eg(j)+d))
        end if
        inc        =inc+etas
        end do
        l(i)        =inc*k
        end do
      deallocate(wt,eg)
       end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     etasfull8tfixednew computes only the triggered intensity of the ETAS model
!     for a given day tfixed on a set of space x-y points different from those observed.
!
!     input magnitude vector must be  m(j)-m0
!     27-7-2017 version with predictor
!
!     called by plot.etasclass.R
!
!
      subroutine etasfull8tfixednew(n,mu,k,c,p,g,d,q,x,y,t,m,predictor,l,ngridtot, xgrid, ygrid,tfixed)
      INTEGER n,ngridtot
      double precision mu,k,c,p,g,d,q,x(n),y(n),xgrid(ngridtot),ygrid(ngridtot),t(n),m(n),predictor(n),l(ngridtot)
      double precision dx,dy,ds,dt,xi,yi,etas,inc,tfixed,eps,dum
      double precision, allocatable :: wt(:), eg(:)
      dum=mu
      eps=1D-10
      allocate(wt(n),eg(n))
      do j=1,n
        dt        =tfixed-t(j)
        wt(j)     =0
        if (dt > 0)  wt(j)=((dt+c)**(-p))*exp(predictor(j))
        eg(j)     =exp(g*m(j))
      end do
      do i=1,ngridtot
        inc        =0
        xi        =xgrid(i)
        yi        =ygrid(i)        
      do j=1,n
      dt        =tfixed-t(j)
        etas        =0
        if (dt > 0)  then
        dx        =xi-x(j)
        dy        =yi-y(j)
        ds        =dx*dx+dy*dy
        etas        =wt(j)*exp(-q*log(ds/eg(j)+d))
        end if
        inc        =inc+etas

        end do
        l(i)=inc*k
        end do
      deallocate(wt,eg)
       end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        integrkdweighted
!       internal called by deltafl1kspacevar       
!
        subroutine integrkdweighted(rangex,xkern,w,nkern,k,h,kintegral)
        INTEGER nkern,k
!         input variables declaration:
      double precision xkern(nkern,k),w(nkern),rangex(k,2),h(k)
!         output variables declaration:
      double precision kintegral
!         work variables declaration:
      double precision stdx1,stdx2,ix1,ix2,wtot, hx(k),parz
        hx=h
        kintegral=0
        wtot=sum(w)
        do i=1,nkern
        parz=1
        do  j=1,k
        stdx2=((rangex(j,2)-xkern(i,j)))/hx(j)
        stdx1=((rangex(j,1)-xkern(i,j)))/hx(j)
        call probnorm(stdx2,ix2)
        call probnorm(stdx1,ix1)
        parz=parz*(ix2-ix1)
        end do
        kintegral=kintegral+parz*w(i)
        end do
        kintegral=kintegral/wtot
        return
        end

        subroutine intensitykweighted(x,m,k,xkern,w,nkern,h,dens)
!         k dimensional isotropic normal kernel
        double precision x(m,k),xkern(nkern,k),w(nkern),h(k),dens(m)
        double precision xj(k),wtot
        double precision stdj(k),fac,factot
        DATA               fac/0.3989422804014327D0/
        INTEGER nkern,m,k
!         call density3(x,y,z,n,xkern,ykern,zkern,n,hdef,dens)
        factot=(fac**k)/(PRODUCT(h))
        wtot=sum(w)
        do j=1,m
            dens(j)=0.
            xj        =x(j,1:k)
            do i=1,nkern
                stdj        =(xj-xkern(i,1:k))/h
                dens(j)=dens(j)+w(i)*exp(-0.5*dot_product(stdj,stdj))
            end do
            dens(j)=dens(j)*factot/wtot
        end do
        return
        end
        
        
        
!    deltafl1kspacevar
! #            deltafl1kspace computes flp increments with respect to 
! #                space only
! #                weighted version
! #                only isotropic use only with indanis=0
! #                2 stage variable kernel with indanis=1
! #             
! #
!        BEWARE:
!       submatrices (rowwise) MUST be passed dinamically; 
!       assign them to an allocatable array!!!!! 
!
!        k= number of variables (tipically two)
!       nh = number of smoothing parameters (at least k)
!! #######################################################################################
!
! called by flpkspace.R (kern.var=0 or 1)
!
      subroutine deltafl1kspacevar(x,t,w,n,k,m1,m2,nh,rangex,h,hdef,dens,integr,delta,expweight,indweight,allocationerr)
! input output  variables declaration:
        INTEGER n,m1,m2,indweight,k,nh,allocationerr
      double precision dens(n),delta(n),integr(n),hdef(nh)
      double precision t(n),x(n,k),w(n ),rangex(k,2),h(nh),expweight
! work variables declaration:  
      INTEGER one,two,four,m,i,err1
      double precision lambda(1),mean(k),sd(k), kintegral,deltat
! allocatable arrays       
      double precision :: x9,w9,x1,xmat9,xmat10,xkern,ykern
      ALLOCATABLE :: x9(:,:), w9(:), x1(:,:),xmat9(:,:),xmat10(:,:), xkern(:),ykern(:)
        one=1
        two=2
        four=4
 !        write(*,*) "indanis = ",indanis. 
        delta(1:n)=0
        do m=m1,m2
!        write(*,*)"start deltafl1kspacevar step m,m1,m2",m,m1,m2
        deltat=t(m+1)-t(m)
        ALLOCATE(x1(1,k), x9(m,k), w9(m),xmat9(m,four), xmat10(m,four),xkern(m),ykern(m),STAT=err1) 
!        WRITE(*,*) "x9,xmat9,err1",SHAPE(x9),SHAPE(xmat9),err1
!        write(*,*) " ----  "
        if (err1/=0)then
              allocationerr=1
              return
            else
              allocationerr=0
        x9=x(1:m,1:k)
        w9  =w(1:m)
        x1  =x((m+1):(m+1),1:k)
        if (indweight==1) then
        do i=1,k
        call univariatew(x(1:m,i),w(1:m),m,mean(i),sd(i))
        hdef(i)=h(i)*sd(i)*(m*1D0)**expweight
        end do
        else
        hdef=h
        end if
        
        call intensitykweighted(x1,one,k,x9,w(1:m),m,hdef,lambda)
        call integrkdweighted(rangex,x9,w(1:m),m,k,hdef,kintegral)
        dens(m)        =lambda(1)
        integr(m)=kintegral*deltat
        delta(m)=LOG(dens(m))-integr(m)
        end if
!       if block for allocate
        DEALLOCATE(x1, x9, w9,xmat9,xmat10, STAT=err1) 
        DEALLOCATE(xkern,ykern, STAT=err1) 

        end do
        return
        end

        subroutine univariatew(x,w,n,mean,sd)
!         input variables declaration:
        INTEGER n
        double precision x(n),w(n)
        double precision mean,sd
        double precision ww
!         computes weighted mean and standard deviation
        ww=sum(w)
        mean=sum(x*w)/ww
        sd=sqrt(sum((x-mean)*(x-mean)*w)/ww)
        return
        end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     probnorm normal probability integral
! -----------------------------------------------------------------------
! 
      subroutine probnorm(Y,P)
      double precision   P,Y
      double precision   ERFC,sqrthalf
      DATA               sqrthalf/.7071067811865475D0/
      P = -Y * sqrthalf
      IF (DABS(P) .LE. 13.2) GO TO 1
      P = 0.0
      IF (Y .LT. 0.0) RETURN
      P = 1.0
      RETURN
    1 P = 0.5 * ERFC(P)
      RETURN
      END
!
