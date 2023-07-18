#####################################################################################
####################################################################################
#
#	profile likelihood for etasclass objects
#
#
#####################################################################################
#####################################################################################
profile.etasclass<-function(fitted,iprofile		=1,
				 nprofile	=7,
				 kprofile	=3,
				 profile.approx	=FALSE,...){
maxprofile=7 # in this version profile only for etas parameters, not for covariates parameters 27-7-2017				 
				 
if (!inherits(fitted,"etasclass")) stop("object is not of the required class etasclass")
iprofile=trunc(iprofile)
if(iprofile>maxprofile | iprofile <1 ) stop("'iprofile' must be an integer between 1 and maxprofile")

nprofile=trunc(nprofile)
if(nprofile<1) stop("'nprofile' must be at least 1")
fitted$params.ind=as.numeric(fitted$params.ind)
if(fitted$params.ind[iprofile]==0) stop("cannot compute profile for a parameter not estimated")


		logl.vec	=0
		param.vec	=0
		
		params.indprof	=fitted$params.ind
		params.fix	=fitted$params.fix
		n.params	=sum(params.indprof)
		iterlim		=100
		n		=length(fitted$cat$time)
	trace		=TRUE # controls the level of 	intermediate printing can be deleted in future versions
#####################################################################################
if(profile.approx){
#	approximation of second order for initial values for non-profile estimators			

		ind.prof	=sum(params.indprof[1:iprofile])
		ind.noprof	=setdiff(1:n.params,ind.prof)
		delta.psi	=solve(fitted$risult$hessian[ind.noprof,ind.noprof])%*%fitted$risult$hessian[ind.noprof,ind.prof]
}
#####################################################################################

		params.indprof[iprofile]	=0
		sq			=fitted$sqm[iprofile]*kprofile
		
		param.vec		=seq(fitted$params[iprofile]-sq,fitted$params[iprofile]+sq,length.out=nprofile)
		logl.vec		=array(0,nprofile)
		params.MLtot		=fitted$params
		
		for (j in 1:nprofile)	{

	cat("start profile-L computation. j=  ")
  cat(j,"\n")

		params.fix[iprofile]	=param.vec[j]
		params		=log(params.MLtot[params.indprof==1])
	
		if(profile.approx) params=params-delta.psi*(params.fix[iprofile]-params.MLtot[iprofile])
cat(params,"\n")
                fitted.prof     =fitted
                fitted.prof$params.fix=params.fix
                fitted.prof$params=params
                fitted.prof$params.ind=params.indprof
                fitted.prof$nparams.etas=fitted$nparams.etas-1
                fitted.prof$nparams=fitted$nparams-1


risult.profile= generaloptimizationNEW(fitted.prof,
		hessian	=TRUE,
                iterlim		=iterlim,
		iprint=FALSE,
		trace=trace)


logl.vec[j]=risult.profile$l.optim
cat(" profile likelihood found","\n")
cat(c(j,param.vec[j],logl.vec[j]),"\n")
	}
			general.optimum =c(fitted$params[iprofile],fitted$logl)
			names(general.optimum)[2]="-logL"
			ris=list(
			iprofile	=iprofile,
			logl.vec	=logl.vec,
			param.vec	=param.vec,
			general.optimum =general.optimum
			)

			
class(ris)=c("profile.etasclass",class(ris))

return(ris)


#############################################################################################
# end of profile computation
#############################################################################################

#############################################################################################
}
