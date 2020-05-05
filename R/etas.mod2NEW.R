etas.mod2NEW <-
function(params=c(1,1,1,1,1,1,1),etas.obj,
                                params.lim	=c(0,0,0,0,0,0,0),
				trace=TRUE,
				iprint=FALSE) 
				{
 		params.ind   =etas.obj$params.ind
 		params.fix   =etas.obj$params.fix
                tmax            =etas.obj$tmax
                cat         =etas.obj$cat
                magn.threshold=etas.obj$magn.threshold
		back.dens     =etas.obj$back.dens
		back.integral =etas.obj$back.integral
		onlytime      =etas.obj$onlytime
		rho.s2        =etas.obj$rho.s2
		nparams.etas  =etas.obj$nparams.etas
		nparams       =etas.obj$nparams
                ntheta          =etas.obj$ntheta		
			
				
				
				
                params.etas    =params[1:nparams.etas]
                betacov        =params[(nparams.etas+1):nparams]
	params.e	=params.fix
	params.e[params.ind==1]=exp(params.etas)+params.lim[params.ind==1]
        lambda	= params.e[1]
        k0  	= params.e[2]
        c   	= params.e[3]
        p   	= params.e[4]
#        a   	= params.e[5]
        gamma   = params.e[5]
        d   	= params.e[6]
        q   	= params.e[7]
	x		=cat$xcat.work
	y		=cat$ycat.work
	magnitudes	=cat$magn1.work

#if(onlytime){
#predictor       =as.vector((etas.obj$cov.matrix)*(a))}
#else
#{predictor       =as.vector((etas.obj$cov.matrix)*(a-gamma))
#}
predictor       =as.matrix(etas.obj$cov.matrix)%*%as.vector(betacov)+as.vector(etas.obj$offset)



times		=cat$time.work
	range.t		=diff(range(times))
	time.init<-	Sys.time()
    	n   	=length(times)
	etas.comp	=as.double(array(0,n))

	
	
##################### begin of the computation Fortran etasfull8 +R  ############################

ris=.Fortran("etasfull8newserial" ,NAOK=TRUE,
			tflag=as.integer(onlytime),
			n=as.integer(n),
			mu=as.double(lambda),k=as.double(k0),
			c=as.double(c),p=as.double(p),
			g=as.double(gamma),
			d=as.double(d),q=as.double(q),
			x=as.double(x),y=as.double(y), t=as.double(times),m=as.double(magnitudes),
			predictor=as.double(predictor),
			l=etas.comp)



			
			
			
	timenow		<-	Sys.time()
	timeelapsed	<-	difftime(timenow,time.init,units="secs")
	etas.comp	=ris$l
	       
	       if(iprint) cat(timeelapsed,"\n")
		

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
##		integral=k0*sum(exp(a*magnitudes)*it)
	
		integralNEW=k0*sum(exp(predictor)*it)
	
	}
else

##############################################################################
# starting space integration (polar transformation)
##############################################################################
{

m1NEW	=as.vector(exp(predictor))
##m1	=as.vector(exp((a-gamma)*magnitudes))
m2	=as.vector(exp(gamma*magnitudes))

time.init		<-	Sys.time()
    ci	=lambda*back.dens + etas.comp

### approximate polar integration to whole space by division in ntheta triangles centered in xi,yi




etasNEW	=rowSums(m1NEW*m2*((rho.s2*rho.s2*etas.obj$dettrasf/m2+d)^(1-q)-d^(1-q)))
spaceNEW	=(pi/((1-q)*ntheta))*etasNEW
integralNEW=sum(k0*it*spaceNEW)

if(iprint){
cat("integral polar","\n")
cat(integral,"\n")
}


}

#print(summary(ris$l))	
#print(summary(risNEW$l))	
#print(integral)
#print(integralNEW)

integral=integralNEW
#stop()	

##############################################################################
#
# end space integration
#
##############################################################################


timenow		<-	Sys.time()
	timeelapsed	<-	difftime(timenow,time.init,units="secs")

	integraltot	=integral+lambda*range.t*back.integral
	logL		= -logL.l+integraltot

if(iprint){	cat("whole Integral computation","\n")
		cat(timeelapsed,"\n")
			cat(params,"\n")
			cat(exp(params),"\n")
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
