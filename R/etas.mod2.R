etas.mod2 <-
function(params=c(1,1,1,1,1,1,1,1),
				params.fix	=array(1,8),
				params.ind	=array(1,8),
				params.lim	=c(0,0,0,0,0,0,0,0),
				onlytime	=FALSE,
				back.dens,
				back.integral,
				cat,
				rho.s2,
				ntheta=100,
				magn.threshold,
				trace=TRUE) {
	iprint=FALSE
	params.e	=params.fix
	params.e[params.ind==1]=exp(params)+params.lim[params.ind==1]
#	print(c(length(params.e),length(params)))
        lambda	= params.e[1]
        k0  	= params.e[2]
        c   	= params.e[3]
        p   	= params.e[4]
        a   	= params.e[5]
        gamma   = params.e[6]
        d   	= params.e[7]
        q   	= params.e[8]
	x		=cat$long
	y		=cat$lat
	magnitudes	=cat$magn1-magn.threshold
	times		=cat$time
	tmax		=max(times)
	range.t		=diff(range(times))
	time.init<-	Sys.time()
    	n   	=length(times)
	etas.comp	=as.double(array(0,n))

##################### begin of the computation Fortran etasfull8 +R  ############################

ris=.Fortran("etasfull8" ,NAOK=TRUE,
			tflag=as.integer(onlytime),
			n=as.integer(n),
			mu=as.double(lambda),k=as.double(k0),
			c=as.double(c),p=as.double(p),
			a=as.double(a),g=as.double(gamma),
			d=as.double(d),q=as.double(q),
			x=as.double(x),y=as.double(y), t=as.double(times),m=as.double(magnitudes),l=etas.comp)

	timenow		<-	Sys.time()
	timeelapsed	<-	difftime(timenow,time.init,units="secs")
	etas.comp	=ris$l
	       
	       if(iprint) print(timeelapsed)
		

   	ci	=lambda*back.dens + etas.comp
	logL.l	=sum(log(ci))

##############################################################################
# starting time integration#
##############################################################################


if(p==1){
			 it=log(c+tmax-times)-log(c)
			}
			else
			{
			 it=((c+tmax-times)^(1-p)-c^(1-p))/(1-p)
			}

if(onlytime){
		integral=k0*sum(exp(a*magnitudes)*it)
	}
else

##############################################################################
# starting space integration (polar transformation)
##############################################################################
{

m1	=as.vector(exp((a-gamma)*magnitudes))
m2	=as.vector(exp(gamma*magnitudes))
### approximate polar integration to whole space
	time.init		<-	Sys.time()

### approximate polar integration to whole space by division in ntheta triangles centered in xi,yi
etas	=rowSums(m1*m2*((rho.s2*rho.s2/m2+d)^(1-q)-d^(1-q)))

space	=(pi/((1-q)*ntheta))*etas
integral=sum(k0*it*space)

if(iprint){
print("integral polar")
print(integral)
}

}

##############################################################################
#
# end space integration
#
##############################################################################
	timenow		<-	Sys.time()
	timeelapsed	<-	difftime(timenow,time.init,units="secs")

	integraltot	=integral+lambda*range.t*back.integral
	logL		= -logL.l+integraltot

if(iprint){	print("whole Integral computation")
		print(timeelapsed)
			print(params)
			print(exp(params))
			}

attr(logL,"etas.vec")=etas.comp
attr(logL,"lambda.vec")=ci
attr(logL,"integraltot")=integraltot

	if(iprint){
	cat("ML step: likelihoods and integral",(c(logL,logL.l,integraltot)),"\n")
      }
if(trace) cat("*")

	     	return(logL)
}
