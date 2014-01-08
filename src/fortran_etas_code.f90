!	bivariate isotropic density estimation routines
!
!	xkern,ykern, vectors of lenght nkern (centers for  kernels)
! 	x,y vectors of lenght  m (where we want to compute the density)
! 	dens vector of lenght  m with the densities in x,y
! 	h vector of 2 bandwidths
! 	w	vector of weights

      subroutine density2(x,y,m,xkern,ykern,nkern,h,w,dens)
	integer*4 nkern,m

!	input variables declaration:
 
      double precision xkern(nkern),ykern(nkern),x(m),y(m),w(nkern),ww,h(2)

!	output variables declaration:
      
      double precision dens(m)

!	work variables declaration:

      double precision stdx,stdy,fac
      DATA               fac/0.3989422804014327D0/

!	fac is 1/sqrt(2*PI)
	ww=sum(w)
	hx=h(1)
	hy=h(2)
      do 120 j=1,m
	dens(j)=0
 
      do 130 i=1,nkern
	stdx	=((x(j)-xkern(i)))/hx	
	stdy	=((y(j)-ykern(i)))/hy	
        dens(j)=dens(j)+w(i)*exp(-(stdx*stdx+stdy*stdy)/2)
 130    continue
	dens(j)=dens(j)*fac*fac/(ww*hx*hy)
 120    continue
       end

!
!	etasfull8 computes only the triggered intensity of the ETAS model
!       with 8 parametrs.
!       called by etas.mod2.R
!
      subroutine etasfull8(tflag,n,mu,k,c,p,a,g,d,q,x,y,t,m,l)
	integer*4 n,tflag
      double precision mu,k,c,p,a,g,d,q,x(n),y(n),t(n),m(n),l(n)
      double precision dx,dy,ds,dt,xi,yi,ti,etas,inc

      do 100 i=2,n
	inc	=0
	ti	=t(i)
	xi	=x(i)
	yi	=y(i)	
      do 110 j=1,i-1
      dt	=ti-t(j)
	etas	=0
	if (dt > 0)  then
	dx	=xi-x(j)
	dy	=yi-y(j)
	ds	=dx*dx+dy*dy
	if (tflag>0) then
	etas	=((dt+c)**(-p))*exp(a*m(j))
	else
	etas	=((dt+c)**(-p))*exp((a-g)*m(j))*(ds/exp(g*m(j))+d )**(-q)
	end if 
	end if
	inc	=inc+etas

 110    continue
	l(i)	=inc*k

 100    continue
       end

!
!     etasfull8tintegrated computes only the triggered intensity of the ETAS model
!     integrated on the time axis
!     on a set of points different from those observed.
! 
!

      subroutine etasfull8tintegrated(n,mu,k,c,p,a,g,d,q,x,y,t,m,l,ngridtot, xgrid, ygrid,tmax)
      integer*4 n,ngridtot
      double precision mu,k,c,p,a,g,d,q,x(n),y(n),xgrid(ngridtot),ygrid(ngridtot),t(n),m(n),l(ngridtot)
      double precision dx,dy,ds,dt,xi,yi,ti,etas,inc,tmax,integrt,eps
      eps=1D-10
      do 1100 i=1,ngridtot
	inc	=0
	xi	=xgrid(i)
	yi	=ygrid(i)	
      do 1110 j=1,n
      dt	=tmax-t(j)
	etas	=0
	if (dt > 0)  then
	dx	=xi-x(j)
	dy	=yi-y(j)
	ds	=dx*dx+dy*dy
! t- integration	
!	check for p=1
	if (abs(p-1)<eps) then
		integrt =log(c+dt)-log(c)
		 else
		integrt =(p-1)*(c**(1-p)-(c+dt)**(1-p))
	end if
	etas	=integrt*exp((a-g)*m(j))*(ds/exp(g*m(j))+d )**(-q)
	end if
	inc	=inc+etas

 1110    continue
	l(i)	=inc*k

 1100   continue
       end



!    deltafl1kspace
! #            deltafl1kspace computes flp increments with respect to 
! #		space only
! #		weighted version
! #		only isotropic use only with indanis=0
! #             
! #
! #######################################################################################
!	BEWARE:
!       submatrices (rowwise) MUST be passed dinamically; 
!       assign them to an allocatable array!!!!! 
!
!
      subroutine deltafl1kspace(x,t,w,n,k,m1,m2,rangex,h,hdef,dens,integr,delta,indanis,expweight,indweight,sigma)
! 	input variables declaration:
	integer*4 n,m1,m2,indweight,indanis,k,one,err,m,i
      double precision t(n),x(n,k),w(n),rangex(k,2),h(k),sigma(k,k),expweight
      double precision deltat
      double precision dens(n),delta(n),integr(n),hdef(k)

      double precision lambda(1),mean(k),sd(k), ranget(2),kintegral
      double precision :: x9,w9,x1
      ALLOCATABLE :: x9(:,:), w9(:),x1(:,:)
	one=1
	delta(1:n)=0
! 	hdef=h*sd(tprec)*(i^(-1/5))	
	do 200 m=m1,m2
	if (indweight==1) then
	do 195 i=1,k
	call univariatew(x(1:m,i),w(1:m),m,mean(i),sd(i))
195	hdef(i)=h(i)*sd(i)*(m*1D0)**expweight
	else
	hdef=h
	end if
	deltat=t(m+1)-t(m)
	ALLOCATE(x1(1,k),x9(m,k),w9(m),STAT=err) 
	x9=x(1:m,1:k)
        w9  =w(1:m)
	x1  =x((m+1):(m+1),1:k)
	call intensitykweighted(x1,one,k,x9,w(1:m),m,hdef,lambda)
	call integrkdweighted(rangex,x9,w(1:m),m,k,hdef,kintegral)
	dens(m)	=lambda(1)
	integr(m)=kintegral*deltat
	delta(m)=LOG(dens(m))-integr(m)

200	continue 
	return
	end

!
!	integrkdweighted
!       internal called by deltafl1kspace       
!
	subroutine integrkdweighted(rangex,xkern,w,nkern,k,h,kintegral)
	integer*4 nkern,k

! 	input variables declaration:

      double precision xkern(nkern,k),w(nkern)
      double precision rangex(k,2),h(k)

! 	output variables declaration:

      double precision kintegral

! 	work variables declaration:

      double precision stdx1
      double precision stdx2
      double precision ix1,ix2,wtot
      double precision hx(k),parz

	hx=h
	kintegral=0
	wtot=sum(w)
      do 160 i=1,nkern
	parz=1
	do 155 j=1,k
	stdx2=((rangex(j,2)-xkern(i,j)))/hx(j)
	stdx1=((rangex(j,1)-xkern(i,j)))/hx(j)
	call probnorm(stdx2,ix2)
	call probnorm(stdx1,ix1)
155	parz=parz*(ix2-ix1)
      kintegral=kintegral+parz*w(i)
 160    continue
	kintegral=kintegral/wtot
	return
	end
!
!	intensitykweighted
!       internal called by deltafl1kspace       
!

        subroutine intensitykweighted(x,m,k,xkern,w,nkern,h,dens)
! 	k dimensional isotropic normal kernel
	double precision x(m,k),xkern(nkern,k),w(nkern),h(k),dens(m)
	double precision ris(k,k),xi(k),xj(k),alphafac(k),wtot
	double precision stdi(k),stdj(k),wout(k,k),fac,factot
        DATA               fac/0.3989422804014327D0/
	integer*4 nkern,m,k

! 	call density3(x,y,z,n,xkern,ykern,zkern,n,hdef,dens)
	factot=(fac**k)/(PRODUCT(h))
	wtot=sum(w)
	do 100 j=1,m
	dens(j)=0.
	xj	=x(j,1:k)
	do 200 i=1,nkern
	stdj	=(xj-xkern(i,1:k))/h
	dens(j)=dens(j)+w(i)*exp(-0.5*dot_product(stdj,stdj))
200	continue
	dens(j)=dens(j)*factot/wtot
100	continue
	return
	end

	
!	univariatew
!	computes univariate weighted mean and standard deviation
!
	
	
	
	subroutine univariatew(x,w,n,mean,sd)
! 	input variables declaration:
	integer*4 n
        double precision x(n),w(n)
	double precision mean,sd
	double precision dev,ww
! 	computes weighted mean and standard deviation
	ww=sum(w)
	mean=sum(x*w)/ww
	sd=sqrt(sum((x-mean)*(x-mean)*w)/ww)
	return
	end
	
	
!     probnorm normal probability integral
! -----------------------------------------------------------------------
! 
      subroutine probnorm(Y,P)
      double precision   P,Y
      double precision   DERFC,sqrthalf
      DATA               sqrthalf/.7071067811865475D0/
      P = -Y * sqrthalf
      IF (DABS(P) .LE. 13.2) GO TO 1
      P = 0.0
      IF (Y .LT. 0.0) RETURN
      P = 1.0
      RETURN
    1 P = 0.5 * DERFC(P)
      RETURN
      END

