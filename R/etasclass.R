etasclass <-
function(cat.orig,
				magn.threshold	=2.5,
				magn.threshold.back=magn.threshold+2,				
##### starting values for parameters
				mu		=1,
				k0		=1,
				c		=0.5,
				p		=1.01,
				a		=1.2,
				gamma		=.5,
				d		=1.,
				q		=1.5,
### indicators
				params.ind=replicate(8,TRUE),
				declustering	=TRUE,
				thinning	=FALSE,
				flp		=TRUE,
				ndeclust        =5,
				onlytime	=FALSE,
				is.backconstant	=FALSE,
##### end of  main input arguments. 
##### Control and secondary arguments:
				description	="",
				cat.back   	=NULL,
				back.smooth	=1.0,
				sectoday	=TRUE,
				usenlm		=TRUE,
				method		="BFGS",
				compsqm 	=TRUE,
				epsmax		= 0.0001,
				iterlim		=100,
				ntheta		=100							)   {
	this.call=match.call()
	params.ind=as.numeric(params.ind)
	flag		=eqcat(cat.orig)
	if (!flag$ok){
 			 print("WRONG EARTHQUAKE CATALOG DEFINITION")
			  return(FALSE)
		      }
	if(sum(abs(params.ind-0.5)==0.5)!=8){print("WRONG PARAMS.IND DEFINITION: ONLY 0/1 ALLOWED SEE HELP")
			  return(FALSE)}

	iter  		=0
  AIC.iter	=numeric(0)
	params.iter	=numeric(0)
	sqm.iter 	=numeric(0)
	rho.weights.iter	=numeric(0)
	hdef		=c(0,0)
	hdef.iter	=numeric(0)
	fl		=0
	fl.iter	=numeric(0)

	AIC.flag	=FALSE	
	time.start	<-	Sys.time()
	trace		=TRUE # controls the level of intermediate printing can be deleted in future versions
	region		=embedding.rect.cat.eps(cat.orig)
	nparams		=sum(params.ind)
        cat 	=   	subset(cat.orig,cat.orig$magn1>=magn.threshold)

 	if(sectoday) 	cat$time	=	cat$time/86400
#### STARTING VALUES
 	if(missing(cat.back))   cat.back=subset(cat,cat$magn1>=magn.threshold.back)
	if (missing(mu)||is.null(mu)) mu= 0.5*nrow(cat)/diff(range(cat$time))
	if (missing(k0))	k0	=mu
	if(!usenlm) compsqm=FALSE
        ord	<-	order(cat$time)
        cat 	<-	cat[ord,]
	n	=	nrow(cat)
        eps	=2*epsmax
if (onlytime || is.backconstant) {
				  declustering	=FALSE
				  flp		=FALSE
				  }

if ((declustering==FALSE)||onlytime || is.backconstant) ndeclust=1
# x.km , y.km are kilometers coordinates of back events only
# xcat.km , ycat.km are kilometers coordinates of all used events 
###############################
# In this version
# NO CORRECTION FOR SPHERICITY #######################################
#
###############################
cat.longlat=cat
rho.s2	=matrix(0,ntheta,n)

	if(onlytime){
	is.backconstant=TRUE
	}
	else
	{
		radius=6371.3
		y.km      =   radius*cat.back$lat*pi/180
		x.km      =   radius*cat.back$long*pi/180
		y.lat	=   cat.back$lat
		x.long	=cat.back$long

		ycat.km      =   radius*cat$lat*pi/180
		xcat.km      =   radius*cat$long*pi/180
cat$long=xcat.km
cat$lat=ycat.km
region=embedding.rect.cat.eps(cat)

		for(i in (1:n)){
		trasf			=region.P(region,c(xcat.km[i],ycat.km[i]),k=ntheta)
		rho.s2[,i]		=trasf$rho
							}
	rho.s2	=t(rho.s2)
	} 
	
###########################    begin clustering #################################################

while ((iter<ndeclust)&((eps>epsmax)||(eps.par>epsmax)||(AIC.flag==FALSE))){
	
	if(is.backconstant==FALSE){

if(iter>0){
 		xcat.km=cat$long
		ycat.km=cat$lat
		
    ### 	beginning of flp step for the estimation of the optimal bandwidths

		if (flp){
		etas.l		=attr(l,"etas.vec")
		etas.integral	=etasflp.integral(params=params.optim,
					params.ind=params.ind,
					params.fix=params.fix,
					cat=cat,magn.threshold=magn.threshold,
					back.dens=back.dens,
					back.integral=back.integral,
					onlytime=onlytime,	

					rho.s2=rho.s2,ntheta=ntheta,
					trace=trace)

# compute etas intensity and integral for each point
# call routine for optimization
# flp MUST be weighted 
		ris.flp=flp1.etas.nlm(cat,
			etas.params=params.MLtot,
			etas.l=etas.l,
			etas.integral=etas.integral,
			w=rho.weights,
			m1=as.integer(nrow(cat)/2),
			m2=as.integer(nrow(cat)-1),
			mh=1
			)
		hdef	=ris.flp$hdef
		fl	=ris.flp$fl
		
# h must be passed to the following 
		}
    ### 	end of flp step 

		if (thinning){
		back.ind	=runif(n)<rho.weights
		x.km	=xcat.km[back.ind]
		y.km	=ycat.km[back.ind]
		x.long	=cat.longlat$long[back.ind]
		y.lat 	=cat.longlat$lat[back.ind]
		w	=replicate(length(x.km),1)
		}
		else
		{
#		back.ind	=runif(n)<rho.weights
		x.km	=xcat.km
		y.km	=ycat.km
		x.long	=cat.longlat$long
		y.lat 	=cat.longlat$lat
		w	=rho.weights
		}
		
		
    }
    else
    {w	=replicate(length(x.km),1)

    }
		if(flp&(iter>0)){
		back.tot	=kde2dnew.fortran(x.km,y.km,xcat.km,ycat.km,factor.xy=1,h=hdef,w=w)
		back.dens	=back.tot$z
		back.integral	=kde2d.integral(x.km,y.km,xcat.km,ycat.km,factor.xy=1,w=w,h=hdef,eps=1/n)
		}
		else
		{
		back.tot	=kde2dnew.fortran(x.km,y.km,xcat.km,ycat.km,factor.xy=back.smooth,w=w)
		back.dens	=back.tot$z
		back.integral	=kde2d.integral(x.km,y.km,xcat.km,ycat.km,factor.xy=back.smooth,w=w,eps=1/n)
		hdef		=back.tot$h
		}
		
		if((iter==0)&(missing(epsmax))) epsmax=quantile(back.dens,0.01)
		}
		else
		{
#			is.backconstant==TRUE
		if(!onlytime) back.dens	=array(1,n)/(diff(range(xcat.km))*diff(range(ycat.km)))
		back.integral	=1
		}
	
			params.fix	=c(mu,k0,c,p,a,gamma,d,q)
			if (iter>0) params.fix=etas.ris$params
      namespar=c("mu","k0","c","p","a","gamma","d","q")
			params.lim=c(0,0,0,0,0,0,0,0)
	if(onlytime)	{
			back.dens=1
			params.ind[6:8]	=c(0,0,0)
			params.fix[6:8]  	=c(0,0,0)}


			params		=log((params.fix-params.lim)[params.ind==1])
cat("Start ML step; Declustering and weighting iteration number: ")
print(iter)
if (usenlm){
	risult =nlm(etas.mod2,params,
		hessian	=TRUE,
		typsize=abs(params),
		print.level=0,
		iterlim		=iterlim,
		params.ind=params.ind,
 		params.fix=params.fix,
		cat=cat,
		magn.threshold=magn.threshold,
		back.dens=back.dens,
		back.integral=back.integral,
		onlytime=onlytime,
		rho.s2=rho.s2,ntheta=ntheta,
		trace=trace)

		params.optim=risult$estimate
		l.optim	 =risult$minimum
	}
else {
	   risult =optim(params,etas.mod2,
		method	=method,
		hessian	=TRUE,
		control=list(trace=2,maxit=iterlim,fnscale=n/diff(range(cat$time)),parscale=sqrt(exp(params))),
		params.ind=params.ind,
 		params.fix=params.fix,
		cat=cat,magn.threshold=magn.threshold,
		back.dens=back.dens,back.integral=back.integral,
		onlytime=onlytime,
		rho.s2=rho.s2,
		trace=trace)
   params.optim=risult$par
   l.optim	 =risult$value

	}

l=etas.mod2(params=params.optim,
		params.ind=params.ind,
 		params.fix=params.fix,
		cat=cat,magn.threshold=magn.threshold,
		back.dens=back.dens,
		back.integral=back.integral,
		onlytime=onlytime,	

		rho.s2=rho.s2,ntheta=ntheta,
		trace=trace)

		
####################("found optimum ML step") #################################################################
det.check=det(risult$hessian)
if (abs(det.check)<1e-20) compsqm=FALSE
sqm=0

if (compsqm)	sqm	=sqrt(diag(solve(risult$hessian)))*exp(params.optim)
###### the optimization is made with respect to the log of the parameters ()
###### so the for the asymptotic standard error we use the approximation Var(exp(y))=[exp(y)]^2 var(y)
	params			=params.fix
	params[params.ind==1]	=exp(params.optim)+params.lim[params.ind==1]
	sqm.tot			=array(0,8)
	sqm.tot[params.ind==1]=sqm
	names(sqm.tot)=namespar
		params.MLtot	=params
#####################################################################################
print("found optimum; end  ML step")
print(iter)

	mu	= params[1]
        k0  	= params[2]
        c   	= params[3]
        p   	= params[4]
        a   	= params[5]
        gamma   = params[6]
        d   	= params[7]
        q   	= params[8]
#####################################################################################
#####################################################################################
#
# time.res= residuals obtained with integral tansform of time intensity function
# to be used only for time processes
#
#####################################################################################
#####################################################################################
	time.res	=array(0,n)
	if (onlytime){
	times.tot	=cat$time
	magnitudes.tot	=cat$magn1-magn.threshold
	tmin=min(times.tot)
	etas.t	=0
	for (i in 1:n) {
	tmax=times.tot[i]
	times		=subset(times.tot,times.tot<tmax)
	magnitudes	=subset(magnitudes.tot,times.tot<tmax)

	if(p==1){	
			 it=log(c+tmax-times)-log(c)
			}
			else
			{
			 it=((c+tmax-times)^(1-p)-c^(1-p))/(1-p)
			}
	
	time.res[i]=(tmax-tmin)*mu+k0*sum(exp(a*magnitudes)*it)
	}
	}

#############################################################################################
#
#		final output####   check for clustering
#
#############################################################################################

names(params.fix)=namespar
names(params.ind)=namespar
names(params.MLtot)=namespar

# rho.weights[i] is the probability that the i-th event is a background event
rho.weights=params.MLtot[1]*back.dens/attr(l,"lambda.vec")

# check convergence: back densities, AIC, parameters

	timenow		<-	Sys.time()
	time.elapsed	<-	difftime(timenow,time.start,units="secs")
	iter		=iter+1
	AIC.iter[iter]	=2*l +2*nparams
	params.iter	=rbind(params.iter,params.MLtot)
	sqm.iter	=rbind(sqm.iter,sqm.tot)
	rho.weights.iter	=rbind(rho.weights.iter,rho.weights)
	hdef.iter	=rbind(hdef.iter,hdef)
	fl.iter		=c(fl.iter,fl)
print(paste("######   ITERATION n.", iter, "  ###### AIC = ",AIC.iter[iter]))
print("Current estimates of parameters: ")
print(round(params.MLtot,5))
if (iter>1){
	eps		=max(abs(back.dens-etas.ris$back.dens))
	eps.par		=max(abs(params.iter[iter]-params.iter[iter-1]))
	AIC.flag	=(AIC.iter[iter]<=min(AIC.iter)) 
	}
rownames(params.iter)=1:nrow(params.iter)

	etas.ris=list(	this.call		=this.call,
			description		=	description,
			time.start		=	time.start,
			time.elapsed		=	time.elapsed,
			time.end		=	timenow,
			risult			=	risult,
			magn.threshold	=magn.threshold,
			magn.threshold.back=magn.threshold.back,
			params.ind	=as.logical(params.ind),
			params.fix	=params.fix,
			params		=params.MLtot,
			sqm		=sqm.tot,
			onlytime	=onlytime,
			is.backconstant	=is.backconstant,			
			time.res	=time.res,
			usenlm		=usenlm,
			AIC.iter	=AIC.iter,
			params.iter	=params.iter,
			sqm.iter	=sqm.iter,
			rho.weights.iter=rho.weights.iter,
			cat.longlat	=cat.longlat,
			cat		=cat,			
			back.integral	=back.integral,
			back.dens	=back.dens,
			declustering	=declustering,
			thinning	=thinning,
			flp 		=flp,
			hdef		=hdef,
			hdef.iter	=hdef.iter,
			fl.iter		=fl.iter,
			back.smooth	=back.smooth,
			rho.weights	=rho.weights, 
			ndeclust        =ndeclust,
			iter		=iter,
			eps 		=eps,
			logl		=l,
			l=attr(l,"lambda.vec"),
			integral	=attr(l,"integraltot"),
			rho.s2		=rho.s2,
			ntheta		=ntheta)
####   check for clustering
       }

#############################################################################################
#
#		convergence obtained. Computation of final output 
#
#############################################################################################

class(etas.ris)		="etasclass"
class(etas.ris$cat)	="eqcat"
class(etas.ris$cat.longlat)	="eqcat"

			return(etas.ris)
}
