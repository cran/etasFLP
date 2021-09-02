#' Title daily.etasclass
#' @description A daily  estimation on a space grid is made 
#' @param x an etasclass object
#' @param ngrid  subdivisions of x and y axis for grid computation of intensities
#' @param nclass number of class for horizontal and vertical axes of the output grid
#' @param tfixed day of computation
#' @param flag.log if log intensity must be used
#' @param ... other optional parameters
#'
#' @return a grid with daily theoretical intensities
#' @export
#' 
#'
#' 
daily.etasclass <-
function(x, ngrid=201,nclass=20,tfixed=0,flag.log=FALSE,...){
	if(class(x)!="etasclass")stop("argument must be an etasclass object")
	nsimps=1+(ngrid-1)/nclass
        simpson=(nsimps==trunc(nsimps))
	if(!simpson)stop("argument nclass must divide exactly ngrid-1 in order to use Simpson's rule to compute theoretical frequencies")

 scaleres=TRUE
 ellipse=FALSE
# kern.var=x$kern.var
 kern.var=FALSE
#tfixed=0 add in the rd file
cat("Computation of the ETAS space intensity integrated intensity on a grid","\n")
cat("for large catalogs computation can take some minutes. Please wait ...","\n")
#####################################################
# plot time intensity # with points
typegraph=1

#else
#{
#dev.new()
#}
#####################################################
### computation of back intensity and triggered intensity on a given x-y grid
### back.grid computed with the same kernel used for back.dens on observed  points
## back.grid computed on km coordinates

# chamge the is.backconstant switch
magnitude=x$cat$magnitude

#if(!x$is.backconstant){


xcat.km	=x$cat$xcat.work
ycat.km	=x$cat$ycat.work
w	=x$rho.weights
hdef	=x$hdef
n	=length(xcat.km)
tmax    =x$tmax
params	=x$params
        mu= params[1]
        k0= params[2]
        c= params[3]
        p= params[4]
#        a= params[5]
        gamma= params[5]
        d= params[6]
        q= params[7]
	#
	    eps=1/n

	    rangex		=range(xcat.km)
	    rangey		=range(ycat.km)
	    eps.x   =diff(rangex)*eps/2
            eps.y   =diff(rangey)*eps/2
            rangex.km  =c(min(rangex)-eps.x,max(rangex)+eps.x)
            rangey.km  =c(min(rangey)-eps.y,max(rangey)+eps.y)

ranget=diff(range(x$cat$time.work))

alpha	=0
if (kern.var) alpha=hdef[3]

space.grid	=xy.grid(rangex.km,rangey.km,ngrid)
ngridtot	=length(space.grid[,1])
#kde=kde2dnew.fortran(xcat.km,ycat.km,space.grid[,1],space.grid[,2],w=w,factor.xy=1,h=hdef,kern.var=kern.var,alpha=alpha)
if(!x$is.backconstant){
kde=kde2dnew.fortran(xcat.km,ycat.km,space.grid[,1],space.grid[,2],w=w,factor.xy=1,h=hdef,hvarx=x$hvarx,hvary=x$hvary)
#wmat1		=matrix(kde$wmat,n,4)
back.grid	=ranget*mu*kde$z
}
else
{
kde=kde2dnew.fortran(xcat.km,ycat.km,space.grid[,1],space.grid[,2],w=w,factor.xy=1,h=c(1,1))
#wmat1		=matrix(kde$wmat,n,4)
back.grid	=ranget*mu*kde$z^0
}

# maps of triggered intensity for a single day
totfixed.grid=back.grid*0.
if(tfixed>0){
ind=x$cat$time.work<tfixed
ris=.Fortran("etasfull8tfixednew" ,NAOK=TRUE,
			n=as.integer(sum(ind)),
			mu=as.double(mu),k=as.double(k0),
			c=as.double(c),p=as.double(p),
			g=as.double(gamma),
			d=as.double(d),q=as.double(q),
			x=as.double(xcat.km[ind]),y=as.double(ycat.km[ind]), t=as.double(x$cat$time.work[ind]),m=as.double(magnitude[ind]),
			predictor=as.double(x$predictor),
                    
			l=back.grid,
			ngridtot=as.integer(ngridtot), xgrid=as.double(space.grid[,1]), 
			ygrid=as.double(space.grid[,2]),tfixed=as.double(tfixed))
		        totfixed.grid	=ris$l+back.grid/ranget

if(flag.log)totfixed.grid=log(totfixed.grid)
}


teo1=matrix(0,nclass,nclass)

coeff.integr		=simpson.kD(nsimps,2)

dx      =	diff(rangex.km)/(ngrid-1)
dy      =	diff(rangey.km)/(ngrid-1)
coeff.integr	=	as.vector(coeff.integr*dx*dy)

x0=min(rangex.km)
y0=min(rangey.km)


for (i in 1:nclass){
  for (j in 1:nclass){
    
    ivert1=(nsimps-1)*(i-1)+1
    ivert2=(nsimps-1)*i+1
    jvert1=(nsimps-1)*(j-1)+1
    jvert2=(nsimps-1)*j+1
    xymat=expand.grid(ivert1:ivert2,jvert1:jvert2)
    xyvec=ngrid*(xymat[,2]-1)+xymat[,1] 
    
    
    teovec1=totfixed.grid[xyvec]
    teo1[i,j]=sum(coeff.integr*teovec1)

  }	
}


return(list(
tfixed=tfixed,
totfixed.grid=totfixed.grid,
teo1=teo1
))
}


