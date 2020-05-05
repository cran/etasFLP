  etasclass <-
function(cat.orig,
         magn.threshold	=2.5,
				 magn.threshold.back=magn.threshold+2,	
         tmax		=max(cat.orig$time),

##### starting values for parameters
				mu		=1,
				k0		=1,
				c		=0.5,
				p		=1.01,
##				a		=1.2,
				gamma		=.5,
				d		=1.,
				q		=1.5,
        betacov         =0.7,
### indicators
				params.ind=replicate(7,TRUE),
### formula for covariates:
                                formula1            ="time~magn1.work-1",
                                offset              =0,
### anisotropic smoothing is not  implemented
###  ellipsoidal spatial function in triggered intensity not yet implemented
				hdef=c(1,1),
				w		=replicate(nrow(cat.orig),1),
                                hvarx            =NULL,
                                hvary            =NULL,
### flags for the kind of declustering and smoothing:
                                declustering	=TRUE,
				thinning	=FALSE,
				flp		=TRUE,
				m1		=NULL,
				ndeclust        =5,
				onlytime	=FALSE,
				is.backconstant	=FALSE,
##### end of  main input arguments. 
##### Control and secondary arguments:
				description	="",
				cat.back   	=NULL,
				back.smooth	=1.0,
				sectoday	=TRUE,
				longlat.to.km   =TRUE,
				usenlm		=TRUE,
				method		="BFGS",
				compsqm 	=TRUE,
				epsmax		= 0.0001,
				iterlim		=50,      
				ntheta		=36)   {
## unused optional arguments
                                iprint  	=FALSE
## first checks, selection of data, and initializatios	
				missingw	=missing(w)

                                this.call=match.call()
	flag		=eqcat(cat.orig)
	if (!flag$ok){
 			 cat("WRONG EARTHQUAKE CATALOG DEFINITION","\n")
			  return(FALSE)
		      }
		      
	cat.orig	=flag$cat
	
#########################################
## 
## INITIALIZATION OF SOME VARIABLES
## 
#########################################
## 
	iter  		=0
        AIC.iter	=numeric(0)
	params.iter	=numeric(0)
	sqm.iter 	=numeric(0)
	rho.weights.iter	=numeric(0)
	hdef.iter	=numeric(0)
	wmat		=numeric(0)
	fl		=0
	fl.iter		=numeric(0)

	AIC.flag	=FALSE	
	AIC.decrease	=TRUE	
	trace		=TRUE # controls the level of 	intermediate printing can be deleted in future versions

# definition of space-time region	
	
#	region		=embedding.rect.cat.eps(cat.orig)
 
	

 #       if(length(w)!=n.back)  {cat("WRONG WEIGHTS DEFINITION","\n")
#			  return(FALSE)
#			  }


        eps	=2*epsmax

        
## check of flags

		if(onlytime)	{is.backconstant=TRUE
			params.ind[5:7]	=c(0,0,0)
			gamma=0
			d=0
			q=0
			}

	if (is.backconstant) {
				  declustering	=FALSE
				  flp		=FALSE
				  }

if ((declustering==FALSE)||is.backconstant) ndeclust=1

# xback.work , yback.work are kilometers coordinates of back events only
# xcat.work , ycat.work are kilometers coordinates of all used events 
###############################
# In this version
# NO CORRECTION FOR SPHERICITY #######################################
#
###############################
## initialization of the object etasclass that will be given as output.
## added june2017
 
## fastML added on january 2018 
 
 
 ####       tmax to be corrected, after selection of catalog
         	if (missing(tmax)) is.na(tmax)=TRUE

		etas.current=list(
                        this.call		=match.call(),
			description		=	description,
			time.start		=	Sys.time(),
			magn.threshold	=magn.threshold,
			magn.threshold.back=magn.threshold.back,
			onlytime	=onlytime,
			tmax		=tmax,
			is.backconstant	=is.backconstant,			
			usenlm		=usenlm,
			method          =method,
			cat.orig        =cat.orig,
			declustering	=declustering,
			thinning	=thinning,
			flp 		=flp,
			back.smooth	=back.smooth,
			ndeclust        =ndeclust,
			eps 		=eps,
                        longlat.to.km   =longlat.to.km,
                        sectoday        =sectoday,

			ntheta		=ntheta)
class(etas.current)		="etasclass"




etas.current=cat.select(etas.current,longlat.to.km,sectoday)
etas.current
cat         =etas.current$cat # serve?
ycat.work      =  etas.current$cat$ycat.work
xcat.work      =  etas.current$cat$xcat.work
	
	n	=	nrow(cat)
	n.back	=	nrow(cat.back)

	
#### STARTING VALUES
	if(missing(m1)||is.null(m1))   m1=as.integer(nrow(cat)/2)
 	if(missing(cat.back))   cat.back=subset(cat,cat$magn1>=magn.threshold.back)
	if (missing(mu)||is.null(mu)) mu= 0.5*nrow(cat)/diff(range(cat$time.work))
	if (missing(k0))	k0	=mu
	print(missing(w))
	print(is.null(w))

## hvarx, hvary added july 2017, variable window terms correction
        if (is.null(hvarx)){ 
            hvarx=replicate(n,1)
            }
        else    
            {hvarx=cat[hvarx]
             hvarx=(prod(hvarx))^(1./n)
            }

        if (is.null(hvary)){ 
            hvary=replicate(n,1)
            }
        else    
            {hvary=cat[hvary]
             hvary=(prod(hvary))^(1./n)
            }
etas.current$hvarx=hvarx
etas.current$hvary=hvary


	if(onlytime)
                {
                rho.s2=0
                region=NULL}
                else
                {
		
		
		
		if (longlat.to.km){
		
		radius=6371.3
		
		yback.work      =   radius*cat.back$lat*pi/180
		xback.work      =   radius*cat.back$long*pi/180

		
		  }
	        else
	        {
       		yback.work      =   cat.back$lat
		xback.work      =   cat.back$long
}
	
		region    =embedding.rect.cat.epsNEW(cat)
		rho.s2	=matrix(0,ntheta,n)

               
		for(i in (1:n)){
		trasf			=region.P(region,c(xcat.work[i],ycat.work[i]),k=ntheta)
		rho.s2[,i]		=trasf$rho
							}
	rho.s2	=t(rho.s2)
              
	} 
	

##              check formula for covariates
                formula1        =as.formula(formula1)
                formula1        =update(formula1,.~.-1)
                cov.matrix      =model.matrix(formula1,data=cat)
                offset          =model.offset(model.frame(formula1,data=cat))
                if(is.null(offset)) offset=replicate(nrow(cat),0)
                if (class(cov.matrix)!="matrix"){
      	              cat("WRONG FORMULA DEFINITION","\n")
                      return(FALSE)}
       
                ncov=ncol(cov.matrix)
               if (length(betacov)!=ncov){
	cat("Wrong number of starting values for parametrs linear predictor of covariates. Correct number of zero starting values inserted","\n")
                      betacov=replicate(ncov,0)
                      
}
              if(length(params.ind)!=7){
	cat("Wrong number of elements of params.ind","\n")
        return(FALSE)}

             	params.ind=as.numeric(params.ind)
	if(sum(abs(params.ind-0.5)==0.5)!=7){
	cat("WRONG params.ind DEFINITION: ONLY FALSE/TRUE ALLOWED SEE HELP","\n")
        return(FALSE)
        }
        
                nparams.etas		=sum(params.ind)
                nparams 		=nparams.etas+ncov
 
## initialization of other element of list etas.current
## added june2017
etas.current$formula1           =formula1
etas.current$cov.matrix         =cov.matrix
etas.current$region             =region
etas.current$rho.s2		=rho.s2
etas.current$offset		=offset
	
	
#################################################	
########     begin clustering          ##########
########     and estimation            ##########
#################################################
#while ((iter<ndeclust)&((eps>epsmax)||(eps.par>epsmax)||(AIC.flag==FALSE))){
while (AIC.decrease&(iter<ndeclust)&((eps>epsmax)||(eps.par>epsmax))){
# attempt avoiding decrasing constrain:

#while ((iter<ndeclust)&((eps>epsmax)||(eps.par>epsmax))){
	
	if(is.backconstant==FALSE){

if(iter>0){
# 		xcat.work=cat$long
#		ycat.work=cat$lat
		
    ### 	beginning of flp step for the estimation of the optimal bandwidths

		if (flp){
		etas.l		=attr(l,"etas.vec")
# compute etas intensity and integral for each point and call routine for optimization
# flp MUST be weighted 
    ris.flp=flp1.etas.nlmNEW(cat,
		        h.init=hdef,
      			etas.params=params.MLtot,
      			etas.l=etas.l,
      			w=rho.weights,
      			m1=m1,
      			m2=as.integer(nrow(cat)-1),
      			mh=1
			)
		hdef	=ris.flp$hdef
		fl	=ris.flp$fl
				}
    ### 	end of flp step 

		if (thinning){
		back.ind	=runif(n)<rho.weights
		xback.work	=xcat.work[back.ind]
		yback.work	=ycat.work[back.ind]
		w	=replicate(length(xback.work),1)
		}
		else
		{
		xback.work	=xcat.work
		yback.work	=ycat.work
		w	=rho.weights
		}
				
    }
    else
    {if(missingw) w	=replicate(length(xback.work),1)
    }
# in future versions unify the two blocks    
		if(flp&(iter>0)){
#		alpha=hdef[3:4]
#		back.tot		=kde2dnew.fortran(xback.work,yback.work,xcat.work,ycat.work,factor.xy=1,eps=1/n,h=hdef,w=w,kern.var=kern.var,alpha=alpha)
print("w")
print(length(w))
print(summary(w))

		back.tot		=kde2dnew.fortran(xback.work,yback.work,xcat.work,ycat.work,factor.xy=1,eps=1/n,h=hdef,w=w,hvarx=hvarx,hvary=hvary)
                
              
		back.dens	=back.tot$z
		back.integral	=back.tot$integral
		wmat		=back.tot$wmat}
		else
		{
#		back.tot	=kde2dnew.fortran(xback.work,yback.work,xcat.work,ycat.work,factor.xy=back.smooth,eps=1/n,w=w,kern.var=kern.var,alpha=alpha)
		if(!missingw){
		xback.work	=xcat.work
		yback.work	=ycat.work
		}
print("w")
print(length(w))
print(summary(w))


		back.tot	=kde2dnew.fortran(xback.work,yback.work,xcat.work,ycat.work,factor.xy=back.smooth,eps=1/n,w=w,hvarx=hvarx,hvary=hvary)
		back.dens	=back.tot$z
		back.integral	=back.tot$integral
		wmat		=back.tot$wmat
		hdef		=back.tot$h
		}
		
		if((iter==0)&(missing(epsmax))) epsmax=quantile(back.dens,0.01)
		}
		else
		{
#			is.backconstant==TRUE
		if(!onlytime) back.dens	=array(1,n)/(diff(range(xcat.work))*diff(range(ycat.work)))
		back.integral	=1
		}
# params.fix has 7 elements (21-7-2017)	
			params.fix	=c(mu,k0,c,p,gamma,d,q)
			
			if (iter>0) params.fix=etas.ris$params
      namespar=c("mu","k0","c","p","gamma","d","q")

			params.lim=c(0,0,0,0,0,0,0)
	if(onlytime)			back.dens=1

			params		=log((params.fix-params.lim)[params.ind==1])
			
print("Start ML step; Declustering and weighting iteration number: ")

print(iter)
 
if (iter==0)rho.weights=replicate(n,0.5)
 
 
etas.current$back.integral=back.integral
etas.current$back.dens=back.dens
etas.current$params=params
etas.current$params.ind	=as.logical(params.ind)
etas.current$params.fix	=params.fix
etas.current$betacov    =betacov
etas.current$nparams    =nparams
etas.current$nparams.etas=nparams.etas
etas.current$rho.s2     =rho.s2
etas.current$rho.weights        =rho.weights


risult.opt =generaloptimizationNEW(etas.current,
		hessian	=TRUE,
                iterlim		=iterlim,
		iprint=iprint,
		trace=trace)

	
etas.current    =risult.opt$etas.obj	
params.optim    =risult.opt$params.optim[1:nparams.etas]
betacov         =risult.opt$params.optim[(nparams.etas+1):nparams]

##
##
## params.optim does not contain betacov
##
##


l.optim=risult.opt$l.optim
risult  =risult.opt$risult

l=etas.mod2NEW(params=c(params.optim,betacov),
		etas.obj=etas.current,
		trace=trace)
	
####################("found optimum ML step") #################################################################
det.check=det(risult$hessian)
if (abs(det.check)<1e-20) compsqm=FALSE
sqm=0

if (compsqm)	sqm	=sqrt(diag(solve(risult$hessian)))*c(exp(params.optim),replicate(ncov,1))
sqm.etas=sqm[1:nparams.etas]
sqm.cov=sqm[(nparams.etas+1):nparams]

###### the optimization is made with respect to the log of the parameters ()
###### so the for the asymptotic standard error we use the approximation Var(exp(y))=[exp(y)]^2 var(y)
	params			=params.fix
	params[params.ind==1]	=exp(params.optim)+params.lim[params.ind==1]
	sqm.tot			=array(0,nparams.etas)
	sqm.tot[params.ind==1]=sqm.etas
	names(sqm.tot)=namespar
	
		params.MLtot	=c(params,betacov)
		sqm.MLtot	=c(sqm.tot,sqm.cov)
		
#####################################################################################
cat("found optimum; end  ML step  ")
cat(iter,"\n")

	mu	= params[1]
        k0  	= params[2]
        c   	= params[3]
        p   	= params[4]
#        a   	= params[5]
        gamma   = params[5]
        d   	= params[6]
        q   	= params[7]
#        betacov = params[8:nparams]

#####################################################################################
# time.res= residuals obtained with integral tansform of time intensity function
# to be used only for time processes
# corrected for predictor in version 2.0.0
# maybe could be omitted later, since plot.etasclass already computes residuals
#####################################################################################
predictor       =as.matrix(etas.current$cov.matrix)%*%as.vector(betacov)+as.vector(etas.current$offset)

	time.res	=array(0,n)
	if (onlytime){
	times.tot	=cat$time
	magnitudes.tot	=cat$magn1-magn.threshold
	tmin=min(times.tot)
	etas.t	=0
	for (i in 1:n) {
	tmax.i=times.tot[i]
	times		=subset(times.tot,times.tot<tmax.i)
	magnitudes	=subset(magnitudes.tot,times.tot<tmax.i)

	if(p==1){	
			 it=log(c+tmax.i-times)-log(c)
			}
			else
			{
			 it=((c+tmax.i-times)^(1-p)-c^(1-p))/(1-p)
			}
	
#	time.res[i]=(tmax.i-tmin)*mu+k0*sum(exp(a*magnitudes)*it)
	time.res[i]=(tmax.i-tmin)*mu+k0*sum(exp(predictor)*it)
	}
	}

#############################################################################################
#		final output####   check for clustering
#############################################################################################

names(params.fix)=namespar
names(params.ind)=namespar
names(betacov)=colnames(cov.matrix)
names(params.MLtot)=c(namespar,names(betacov))
names(sqm.MLtot)=c(namespar,names(betacov))


# rho.weights[i] is the probability that the i-th event is a background event
rho.weights=params.MLtot[1]*back.dens/attr(l,"lambda.vec")

# check convergence: back densities, AIC, parameters

	timenow		<-	Sys.time()
	time.elapsed	<-	difftime(timenow,etas.current$time.start,units="secs")
	AIC.temp	=2*l +2*nparams
	AIC.decrease	=(AIC.temp<=min(AIC.iter)) 
	
#	AIC.flag	=AIC.flag1
        if (AIC.decrease){
	iter		=iter+1
	AIC.iter[iter]	=AIC.temp
	params.iter	=rbind(params.iter,params.MLtot)
	sqm.iter	=rbind(sqm.iter,sqm.tot)
	rho.weights.iter	=rbind(rho.weights.iter,rho.weights)
	hdef.iter	=rbind(hdef.iter,hdef)
	
	fl.iter		=c(fl.iter,fl)
cat(paste("######   ITERATION n.", iter, "  ###### AIC = ",AIC.iter[iter]),"\n")
cat("Current estimates of parameters: ","\n")
cat(round(params.MLtot,5))
cat("\n")
}
if (iter>1){
	eps		=max(abs(back.dens-etas.ris$back.dens))
	eps.par		=max(abs(params.iter[iter]-params.iter[iter-1]))
	AIC.flag	=(AIC.iter[iter]<=min(AIC.iter)) 
	}
rownames(params.iter)=1:nrow(params.iter)

predictor       =as.matrix(etas.current$cov.matrix)%*%as.vector(betacov)+as.vector(etas.current$offset)

names(params.MLtot)=c(namespar,names(betacov))
names(sqm.MLtot)=names(params.MLtot)

etas.current$time.elapsed   =	time.elapsed
etas.current$time.end		=	timenow
etas.current$risult		=	risult
etas.current$betacov	=betacov
etas.current$params.MLtot	=params.MLtot
etas.current$params	      =params.MLtot[1:7]
etas.current$predictor		=predictor
etas.current$sqm		=sqm.MLtot
etas.current$time.res	          =time.res
etas.current$AIC.iter	      =AIC.iter
etas.current$params.iter	=params.iter
etas.current$sqm.iter	        =sqm.iter
etas.current$rho.weights.iter=rho.weights.iter
etas.current$rho.weights        =rho.weights
etas.current$hdef		=hdef
etas.current$hdef.iter	         =hdef.iter
etas.current$wmat		=wmat
etas.current$fl.iter		=fl.iter
etas.current$iter		=iter
etas.current$logl		=l        
etas.current$l                  =attr(l,"lambda.vec")       
etas.current$integral	          =attr(l,"integraltot")        

etas.ris=etas.current
####   check for clustering
       }

#############################################################################################
#
#		convergence obtained. Computation of final output 
#
#############################################################################################
summary(etas.ris)
			return(etas.ris)
}
