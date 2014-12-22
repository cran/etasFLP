!	bivariate adaptive (2-stages)  density estimation routine
!
!	xkern,ykern, vectors of lenght nkern (centers for  kernels)
! 	x,y vectors of lenght  m (where we want to compute the density)
! 	w	vector of input  weights
! 	h vector of 2 (or 3) bandwidths for the first stage kernel dens
!	h(1) is on x-scale, h(2) on y scale, h(3) is a correlation,
! 	
!       densdef vector of lenght  m with the second-stage densities in x,y

! 	wmat	matrix (nkern*4) of output  weights

!	bivariate adaptive (2-stages)  density estimation routine
!       divided  in two subroutines:
!       density2adaptfirst for the first stage computation of weights (wmat)
!       
! 	wmat	matrix (nkern*4) of output  weights
!       given only !	xkern,ykern, vectors of lenght nkern (centers for  kernels)
!       and density2adaptsecond for the computation of the final densities densdef (given wmat)
!
      subroutine density2adapt(x,y,m,xkern,ykern,nkern,h,w,densdef,wmat,wcomb)
	integer*4 nkern,m
!	input variables declaration:
      double precision xkern(nkern),ykern(nkern),x(m),y(m),w(nkern),h(4)
!	output variables declaration:     
      double precision wn(nkern),wx(nkern),wy(nkern),wxy(nkern),wtot
      double precision dens(nkern),densdef(m),wmat(nkern,4),wcomb(nkern,4),ww(nkern),dd
!	work variables declaration:
      double precision fac,rxy,fac2,hx,hy,hxy,stdx(nkern),stdy(nkern),detr,dbiv,alpha(2)
      logical mask(nkern),maski(nkern)
      DATA               fac/0.3989422804014327D0/       !	fac is 1/sqrt(2*PI)
!	ww=sum(w)
!      write(*,*)" in density2adapt h per wmat",h 

      call density2adaptfirst(xkern,ykern,nkern,h,wmat)
      alpha=h(3:4)

      call density2adaptsecondextended(x,y,m,xkern,ykern,nkern,h,w,densdef,wmat,wcomb,alpha)
       end
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine density2adaptfirst(xkern,ykern,nkern,h,wmat)
!input variables declaration:
      integer*4 nkern
      double precision xkern(nkern),ykern(nkern),h(4)
!output variables declaration:
       double precision wmat(nkern,4)
!work variables declaration:
      double precision wn(nkern),wx(nkern),wy(nkern),wxy(nkern),wtot,htrans(2),w(nkern)
      double precision l(nkern),dens(nkern),ww(nkern),dd,mlog,mgeom,determ(nkern),cost(nkern)
      double precision fac,rxy,fac2,hx,hy,stdx(nkern),stdy(nkern),detr,dbiv
      logical maski(nkern)
      DATA               fac/0.3989422804014327D0/     !	fac is 1/sqrt(2*PI)
	hx=h(1)
	hy=h(2)
	htrans=h(1:2)
	w=1D0/nkern
	fac2=fac*fac/(hx*hy)
        
      do  i=1,nkern
	stdx	=((xkern-xkern(i)))/hx	
	stdy	=((ykern-ykern(i)))/hy	
        dens	=fac2*exp(-(stdx*stdx+stdy*stdy)/2)
        wx(i) =SUM(dens*(xkern-xkern(i))**2)
        wy(i) =SUM(dens*(ykern-ykern(i))**2)
        wxy(i)=SUM(dens*(xkern-xkern(i))*(ykern-ykern(i)))
	wmat(i,1)=sum(dens)
         
       	end do

 !       write(*,*)"wx ----  ",wx(1:6)
 !       write(*,*)"wy ----  ",wy(1:6)
 !       write(*,*)"wxy ----  ",wxy(1:6)
       	ww=wmat(:,1)
       	wtot=sum(ww)
       	wx = wx/wtot
       	wy = wy/wtot
       	wxy=wxy/wtot     	
! computation of constants so that sqrt(volume[i])=GeometricMean(f)/f[i]
        mgeom =EXP(SUM(LOG(ww))/nkern)
        determ=wx*wy-wxy*wxy

        l  	=mgeom/(ww)
        cost	=l/sqrt(determ)
        
  !      write(*,*)"determ ----  ",determ(1:6)
   !     write(*,*)"ww ----  ",ww(1:6)
   !     write(*,*)"wx ----  ",wx(1:6)
   !     write(*,*)"wy ----  ",wy(1:6)
   !     write(*,*)"wxy ----  ",wxy(1:6)
  
        wmat(:,2)=sqrt(cost*wx)
       	wmat(:,3)=sqrt(cost*wy)
       	wmat(:,4)=wxy/sqrt(wx*wy)

        end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     subroutine density2adaptsecond(x,y,m,xkern,ykern,nkern,w,densdef,wmat)
	integer*4 nkern,m
!	input variables declaration:
      double precision xkern(nkern),ykern(nkern),x(m),y(m),w(nkern),wmat(nkern,4),h(4)
!	output variables declaration:
      double precision densdef(m)
!	work variables declaration:
      double precision wtot,dd
      double precision fac,rxy,fac2,detr,dbiv
      logical mask(nkern),maski(nkern)
      DATA               fac/0.3989422804014327D0/!	fac is 1/sqrt(2*PI)        
	wtot=sum(w)
! begin computation of kernel densities with variable mmetrics already defined in wmat
        do 180 i=1,m
        dd=0
        do 170 j=1,nkern
        call densbivnorm(x(i),y(i),xkern(j),ykern(j),wmat(j,2),wmat(j,3),wmat(j,4),dbiv) 
170        dd=dd+dbiv*w(j)
         densdef(i)=dd/wtot        
180       continue   
       return
       end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the same as density2adaptsecond but the variable kernel matrix phi(i) is 
! obtained combining individual matrices and a global average matrix
! phi(i)=(alpha*sigma(i)+(1-alpha)*sigma)/2
! alpha=0: constan(2)t kernel
! alpha=1: fully variable kernel
       subroutine density2adaptsecondextended(x,y,m,xkern,ykern,nkern,h,w,densdef,wmat,wcomb,alpha)
	integer*4 nkern,m
!	input variables declaration:
      double precision xkern(nkern),ykern(nkern),x(m),y(m),w(nkern),wmat(nkern,4),wcomb(nkern,4),h(4),alpha(2)
!	output variables declaration:    
      double precision densdef(m)
!	work variables declaration:
      double precision wtot,dd,fac,rxy,fac2,detr,dbiv
      logical mask(nkern),maski(nkern)
      DATA               fac/0.3989422804014327D0/!	fac is 1/sqrt(2*PI)
      wcomb=wmat

      do j=2,3
      wcomb(:,j)=sqrt(wmat(:,j)**2+h(j+1)**2)
      end do
      call density2adaptsecond(x,y,m,xkern,ykern,nkern,w,densdef,wcomb)
      return
      end 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!       	the same as density2adapt 
!       but input are  bivariate matrices and not 2 vectors.
!	xkernmat (nkern,2) 
! 	xmat      (n,2)
! 	w	vector of input  weights
! 	h vector of 2 (or 3) bandwidths for the first stage kernel dens
!	h(1) is on x-scale, h(2) on y scale, h(3) is a correlation
! 	
!       densdef vector of lenght  m with the second-stage densities in x,y

!       densdef vector of lenght  nkern with the second-stage densities in x,y
! 	wmat	matrix (nkern*4) of output  weights
      subroutine density2adaptmat(xmat,m,xkernmat,nkern,h,w,densdef,wmat,wcomb)
	integer*4 nkern,m
!	input variables declaration:
      double precision xkernmat(nkern,2),xmat(m,2),w(nkern),h(4)
      double precision xkern(nkern),ykern(nkern),x(m),y(m)
!	output variables declaration:    
      double precision wn(nkern),wx(nkern),wy(nkern),wxy(nkern),wtot
      double precision dens(nkern),densdef(m),wmat(nkern,4),wcomb(nkern,4),ww(nkern),dd
!	work variables declaration:
      double precision fac,rxy,fac2,hx,hy,hxy,stdx(nkern),stdy(nkern),detr,dbiv
      DATA               fac/0.3989422804014327D0/
	xkern=xkernmat(:,1)
	ykern=xkernmat(:,2)
	x=xmat(:,1)
	y=xmat(:,2)
       call density2adapt(x,y,m,xkern,ykern,nkern,h,w,densdef,wmat,wcomb)	
       end   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
       subroutine densbivnorm(x,y,mx,my,sx,sy,r,dens)
       double precision x,y,mx,my,sx,sy,r,dens,fac,detr
       data fac/4.442882938158366/ !      fac is sqrt(2)*pi       
	stdx	=(x-mx)/sx	
	stdy	=(y-my)/sy
	detr    =2*(1-r*r)    
        dens	=exp(-(stdx*stdx+stdy*stdy-2*r*stdx*stdy)/detr)/(sx*sy*sqrt(detr)*fac)
	return
	end
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine densbivnormstd(x,y,r,dens)
       double precision x,y,r,dens,fac,detr
       data fac/4.442882938158366/ !      fac is sqrt(2)*pi       
	detr    =2*(1-r*r)       
        dens	=exp(-(x*x+y*y-2*r*x*y)/detr)/(sqrt(detr)*fac)
	return
	end
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      DATA               fac/0.3989422804014327D0/ !	fac is 1/sqrt(2*PI)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	etasfull8 computes only the triggered intensity of the ETAS model
!       with 8 parametrs.
!       called by etas.mod2.R
!       input magnitude vector must be  m(j)-m0
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
       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	etasfull8 computes only the triggering intensity of the ETAS model with 8 parametrs.of each  mainshock 
!       sum the contribution of each event on the following ones
!       called by plot.etasclass
!       input magnitude vector must be  m(j)-m0
!
      subroutine etasfull8reversed(tflag,n,mu,k,c,p,a,g,d,q,x,y,t,m,ltot,l)
	integer*4 n,tflag
      double precision mu,k,c,p,a,g,d,q,x(n),y(n),t(n),m(n),ltot(n),l(n)
      double precision dx,dy,ds,dt,xi,yi,ti,etas,inc
      do 100 i=1,n-1
	inc	=0
	ti	=t(i)
	xi	=x(i)
	yi	=y(i)	
      do 110 j=i+1,n
      dt	=t(j)-ti
	etas	=0
	if (dt > 0)  then
	dx	=xi-x(j)
	dy	=yi-y(j)
	ds	=dx*dx+dy*dy
	if (tflag>0) then
	etas	=((dt+c)**(-p))*exp(a*m(i))
	else
	etas	=((dt+c)**(-p))*exp((a-g)*m(i))*(ds/exp(g*m(i))+d )**(-q)
	end if 
	end if
	inc	=inc+etas/ltot(j)
 110    continue
	l(i)	=inc*k
 100    continue
       end
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     etasfull8tintegrated computes only the triggered intensity of the ETAS model
!     integrated on the time axis
!     on a set of points different from those observed.
!
!     input magnitude vector must be  m(j)-m0
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
!		check for p=1
	if (abs(p-1)<eps) then
		integrt =log(c+dt)-log(c)
		 else
		integrt =((c+dt)**(1-p)- c**(1-p))/(1-p)
	end if
	etas	=integrt*exp((a-g)*m(j))*(ds/exp(g*m(j))+d )**(-q)
	end if
	inc	=inc+etas
 1110    continue
	l(i)	=inc*k
 1100   continue
       end
       
       
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     etasfull8tfixed computes only the triggered intensity of the ETAS model
!     for a given day tfixed on a set of space x-y points different from those observed.
!
!     input magnitude vector must be  m(j)-m0
!
!
 
      subroutine etasfull8tfixed(n,mu,k,c,p,a,g,d,q,x,y,t,m,l,ngridtot, xgrid, ygrid,tfixed)
      integer*4 n,ngridtot
      double precision mu,k,c,p,a,g,d,q,x(n),y(n),xgrid(ngridtot),ygrid(ngridtot),t(n),m(n),l(ngridtot)
      double precision dx,dy,ds,dt,xi,yi,ti,etas,inc,tfixed,eps
      eps=1D-10
      do 1100 i=1,ngridtot
	inc	=0
	xi	=xgrid(i)
	yi	=ygrid(i)	
      do 1110 j=1,n
      dt	=tfixed-t(j)
	etas	=0
	if (dt > 0)  then
	dx	=xi-x(j)
	dy	=yi-y(j)
	ds	=dx*dx+dy*dy
	etas	=((dt+c)**(-p))*exp((a-g)*m(j))*(ds/exp(g*m(j))+d )**(-q)
	end if
	inc	=inc+etas
!	if (dt < 2) write(*,*) i,j,xi,x(j),yi,y(j)

 1110    continue
	l(i)=inc*k
 1100   continue
       end

       
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	integrkdweighted
!       internal called by deltafl1kspacevar       
!
	subroutine integrkdweighted(rangex,xkern,w,nkern,k,h,kintegral)
	integer*4 nkern,k
! 	input variables declaration:
      double precision xkern(nkern,k),w(nkern),rangex(k,2),h(k)
! 	output variables declaration:
      double precision kintegral
! 	work variables declaration:
      double precision stdx1,stdx2,ix1,ix2,wtot, hx(k),parz
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
!    deltafl1kspacevar
! #            deltafl1kspace computes flp increments with respect to 
! #		space only
! #		weighted version
! #		only isotropic use only with indanis=0
! #		2 stage variable kernel with indanis=1
! #             
! #
!	BEWARE:
!       submatrices (rowwise) MUST be passed dinamically; 
!       assign them to an allocatable array!!!!! 
!
!	k= number of variables (tipically two)
!       nh = number of smoothing parameters (at least k)
!! #######################################################################################

      subroutine deltafl1kspacevar(x,t,w,n,k,m1,m2,nh,rangex,h,hdef,dens,integr,delta,indanis,expweight,indweight,allocationerr)
! input output  variables declaration:
	integer*4 n,m1,m2,indweight,indanis,k,nh,allocationerr
      double precision dens(n),delta(n),integr(n),hdef(nh)
      double precision t(n),x(n,k),w(n ),rangex(k,2),h(nh),expweight
! work variables declaration:  
      integer*4 one,two,four,err,m,i,err1,err2
      double precision deltat, wmat(n,4),xx1(1),yy1(1),alpha(2)
      double precision lambda(1),mean(k),sd(k), ranget(2),kintegral
! allocatable arrays       
      double precision :: x9,w9,x1,xmat9,xmat10,xkern,ykern
      ALLOCATABLE :: x9(:,:), w9(:), x1(:,:),xmat9(:,:),xmat10(:,:), xkern(:),ykern(:)
	one=1
	two=2
	four=4
 !	write(*,*) "indanis = ",indanis
	delta(1:n)=0
 	if (indanis==1) then
 	hdef(3:4)=h(3:4)
        alpha=h(3:4)
	ALLOCATE(xkern(n),ykern(n), STAT=err1) 
 	xkern=x(:,1)
 	ykern=x(:,2)
	if (indweight==1) then	
	call univariatew(xkern,w,n,mean(1),sd(1))
	hdef(1)=h(1)*sd(1)*(n*1D0)**expweight
	call univariatew(ykern,w,n,mean(2),sd(2))
	hdef(2)=h(2)*sd(2)*(n*1D0)**expweight
	end if
!	write(*,*)" dentro flp; kernvar; h per wmat: ",hdef 
        call density2adaptfirst(xkern,ykern,n,hdef,wmat)
! wmat contains the local weight matrices fixed for the whole process of computation of the flp contributions
!        write(*,*)"wmat "
	DEALLOCATE(xkern,ykern, STAT=err1) 	
	end if 
	idebug=0
	if(idebug==1)then
	  return
	end if	
	do 200 m=m1,m2
!	write(*,*)"start deltafl1kspacevar step m,m1,m2",m,m1,m2
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
	if (indanis==0) then
	
	call intensitykweighted(x1,one,k,x9,w(1:m),m,hdef,lambda)
	call integrkdweighted(rangex,x9,w(1:m),m,k,hdef,kintegral)
	else
        xx1(1)=x1(1,1)
        yy1(1)=x1(1,2)
       	xmat9=wmat(1:m,1:4)
       	xmat10=xmat9
	xkern=x9(:,1)
 	ykern=x9(:,2)
 	do i=1,2
	hdef(i+2)=alpha(i)*sd(i)*(m*1D0)**expweight
	end do
!! 	alpha=0.5
!       introduce weighting breaking the call to density2adaptsecondextended
        xmat10(:,2)=(alpha*xmat9(:,2)+hdef(1)*(1-alpha))
        xmat10(:,3)=(alpha*xmat9(:,3)+hdef(2)*(1-alpha))
        xmat10(:,4)=alpha*xmat9(:,4)
        
        call density2adaptsecond(xx1,yy1,one,xkern,ykern,m,w(1:m),lambda,xmat10)
!       call density2adaptmat(x1,one,x9,m,hdef,w(1:m),lambda,xmat9)
!	write(*,*)"deltafl1kspacevar after density2adaptmat and before integr call"
	call integrkdweightedvar(rangex,x9,w(1:m),m,two,xmat10,kintegral)
	endif
	dens(m)	=lambda(1)
	integr(m)=kintegral*deltat
	delta(m)=LOG(dens(m))-integr(m)
	end if
!       if block for allocate
	DEALLOCATE(x1, x9, w9,xmat9,xmat10, STAT=err1) 
	DEALLOCATE(xkern,ykern, STAT=err1) 

200	continue 
	return
	end
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine integrkdweightedvar(rangex,xkern,w,nkern,k,wmat,kintegral)
	integer*4 nkern,k
! 	input variables declaration:
      double precision xkern(nkern,k),w(nkern),rangex(k,2),wmat(nkern,4)
! 	output variables declaration:
      double precision kintegral
! 	work variables declaration:
      double precision stdx1,stdx2,stdy1,stdy2,rho,ix1,ix2,ix3,ix4,wtot, parz
  	kintegral=0
	wtot=sum(w)
      do 160 i=1,nkern
!        write(*,*)"integration step", i
	rho	=wmat(i,4)  
	stdx2	=((rangex(1,2)-xkern(i,1)))/wmat(i,2)
	stdx1	=((rangex(1,1)-xkern(i,1)))/wmat(i,2)
	stdy2	=((rangex(2,2)-xkern(i,2)))/wmat(i,3)
	stdy1	=((rangex(2,1)-xkern(i,2)))/wmat(i,3)
!	call MDBNOR(stdx2,stdy2,rho,ix4,IER)
!	call MDBNOR(stdx2,stdy1,rho,ix3,IER)
!	call MDBNOR(stdx1,stdy2,rho,ix2,IER)
!	call MDBNOR(stdx1,stdy1,rho,ix1,IER)
	parz=ix4-ix3-ix2+ix1
      kintegral=kintegral+parz*w(i)
 160    continue
	kintegral=kintegral/wtot
	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!