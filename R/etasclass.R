etasclass <-
  function(cat.orig,
           time.update=FALSE,
           magn.threshold	=2.5,
           magn.threshold.back=magn.threshold+2,	
           tmax		=max(cat.orig$time),
           long.range=range(cat.orig$long),
           lat.range=range(cat.orig$lat),
           ##### starting values for parameters
           mu		=1,
           k0		=1,
           c		=0.5,
           p		=1.01,
           gamma	=.5,
           d		=1.,
           q		=1.5,
           betacov         =0.7,
           ### indicators: if params.ind[i] i-th parameter will be estimated
           params.ind=replicate(7,TRUE),
#           params.lim=c(0,0,0,1.0,0,0,0),
           ### formula for covariates (magnitude should always be included):
           formula1            ="time~magnitude-1",
           offset              =0,
           hdef=c(1,1),
           w		=replicate(nrow(cat.orig),1),
           hvarx  =replicate(nrow(cat.orig),1),
           hvary  =replicate(nrow(cat.orig),1),
           ### flags for the kind of declustering and smoothing:
           declustering	   =TRUE,
           thinning	       =FALSE,
           flp		         =TRUE,
           m1		           =NULL,
           ndeclust        =5,
           n.iterweight    =1,
           onlytime	=FALSE,
           is.backconstant	=FALSE,
           ##### end of  main input arguments. 
           ##### Control and secondary arguments:
           description	="",
           cat.back   	=NULL,
           back.smooth	=1.0,
           sectoday	=FALSE,
           longlat.to.km   =TRUE,
#           fastML=FALSE, #### not yet implemented
#           fast.eps=0.001, #### not yet implemented
           usenlm		=TRUE,
           method		="BFGS",
           compsqm 	=TRUE,
           epsmax		= 0.0001,
           iterlim		=50,      
           ntheta		=36)   {
    ## unused optional arguments
    iprint  	<-FALSE
    fastML    <-FALSE
    parallel   <-FALSE
    fast.eps   <-0.001
    ## first checks, selection of data, and initializatios	
          params.lim=c(0,0,0,1.0,0,0,0)
    
    this.call<-match.call()
    flag		<-eqcat(cat.orig)
    if (!flag$ok) stop("WRONG EARTHQUAKE CATALOG DEFINITION")
    
    
    cat.orig	<-flag$cat
    
    #########################################
    ## 
    ## INITIALIZATION OF SOME VARIABLES
    ## 
    #########################################
    ## 
    iter  		<-0
    AIC.iter	<-numeric(0)
    AIC.iter2	<-numeric(0)
    params.iter	<-numeric(0)
    sqm.iter 	<-numeric(0)
    rho.weights.iter	<-numeric(0)
    hdef.iter	<-numeric(0)
    wmat		<-numeric(0)
    fl		<-0
    fl.iter		<-numeric(0)
    
    AIC.decrease	<-TRUE	
    trace		<-TRUE # controls the level of 	intermediate printing can be deleted in future versions
    
    eps	<-2*epsmax
    eps.par	<-2*epsmax
    
    ## check of flags
    
    if(onlytime)	{
      is.backconstant <-TRUE
      declustering	<-FALSE
      params.ind[5:7]	<-c(0,0,0)
      gamma<-0
      d<-0
      q<-0
    }
    
    if (is.backconstant) {
      declustering	<-FALSE
    }
    
    if (!declustering) {
      ndeclust<-1
      thinning<-FALSE
      flp     <-FALSE
    }
    
    # xback.work , yback.work are kilometers coordinates of back events only
    # xcat.work , ycat.work are kilometers coordinates of all used events 
    ###############################
    # In this version
    # NO CORRECTION FOR SPHERICITY #######################################
    #
    ###############################
    ## initialization of the object etasclass that will be given as output.
    ## added june2017
    
    ####       tmax to be corrected, after selection of catalog together with spatial region
    if (missing(tmax)) is.na(tmax)<-TRUE
    if (missing(long.range)) is.na(long.range)<-TRUE
    if (missing(lat.range)) is.na(lat.range)<-TRUE
    
    
    etas.current<-list(
      parallel    =parallel,
      this.call		=match.call(),
      nstep.flp=0, #### added feb 2021
      nstep.kde=0, #### added feb 2021
      nstep.par=0, #### added feb 2021
      description		=	description,
      time.start		=	Sys.time(),
      magn.threshold	=magn.threshold,
      magn.threshold.back=magn.threshold.back,
      onlytime	=onlytime,
      tmax		=tmax,
      lat.range   =lat.range,
      long.range  =long.range,
hvarx=hvarx,
hvary=hvary,
      
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
    class(etas.current)		<-"etasclass"
    etas.current<-cat.select(etas.current,longlat.to.km,sectoday=sectoday)
    
    ################ block with hvarx,hvary,w to be checked. cat should be subsetted according to ord and ind?  
    
    #etas.current$cat<-data.frame(etas.current$cat,hvarx,hvary,w) ########   change
    etas.current$cat<-data.frame(etas.current$cat,hvarx,hvary) ########   change
    
    cat         <-etas.current$cat[etas.current$cat$ord,][etas.current$cat$ind[etas.current$cat$ord],]
    etas.current$cat<-cat
    n	<-	nrow(cat)
    missingw        <-missing(w)
    if(missingw)  w=replicate(n,1)    
    if(time.update){
    n.old<-length(w)
    w    <-c(w,replicate(n-n.old,0.5))
    hvarx    <-c(hvarx,replicate(n-n.old,1))
    hvary    <-c(hvary,replicate(n-n.old,1))
    }
    else
    {
    # cat now is ordered and selected 
    ## hvarx, hvary added july 2017, variable window terms correction
    hvarx<-cat$hvarx/(prod(cat$hvarx))^(1./n)
    hvary<-cat$hvary/(prod(cat$hvary))^(1./n)
    etas.current$w<-w
    ################
    }
    ycat.work      <-  cat$ycat.work
    xcat.work      <-  cat$xcat.work
    
    
    if(missing(cat.back)||is.null(cat.back))   ind.back<-(cat$magn1>=magn.threshold.back)  else stop("cat.back argument no more allowed, computed internally")##### back coordinates already ok
    
    xback.work      <-   cat$xcat.work[ind.back]
    yback.work      <-   cat$ycat.work[ind.back]
    if(!missingw){
      xback.work	<-xcat.work
      yback.work	<-ycat.work
    }
    if(missingw)  w=replicate(length(xback.work),1)
    
    starting<-etas.starting(cat,magn.threshold=magn.threshold,longlat.to.km = longlat.to.km,sectoday = sectoday, p.start=p,gamma.start=gamma,q.start=q,betacov.start=betacov[1],onlytime=onlytime)
    
    #### STARTING VALUES
    if(missing(m1)||is.null(m1))   m1     <-as.integer(nrow(cat)/2)
    if(missing(mu)||is.null(mu))   mu     <-starting$mu.start
    if(missing(k0)||is.null(k0))   k0	  <-starting$k0.start
    if(missing(c)||is.null(c))     c	  <-starting$c.start
    if(missing(d)||is.null(d))     d	  <-starting$d.start
    if(missing(hdef)||is.null(hdef)) hdef  <-c(bwd.nrd(xback.work,w),bwd.nrd(yback.work,w))
    print("Initial ETAS params estimates after checking: ")
    vecpar.init=c(mu,k0,c,p,gamma,d,q)
    names(vecpar.init)=c("mu","k0","c","p","gamma","d","q")
    
    print(round(vecpar.init,4))    
    
    etas.current$xback.work=xback.work
    etas.current$yback.work=yback.work
    
    etas.current$mu.start = mu
    etas.current$k0.start = k0
    etas.current$c.start = c
    etas.current$d.start = d
    etas.current$p.start = p
    etas.current$q.start = q
    etas.current$gamma.start = gamma
    etas.current$betacov.start = betacov
    etas.current$hdef.start = hdef
    etas.current$formula1           =formula1
    
    
    
    ##              check formula for covariates
    formula1        =as.formula(formula1)
    formula1        =update(formula1,.~.-1)
    cov.matrix      =model.matrix(formula1,data=cat)
    offset          =model.offset(model.frame(formula1,data=cat))
    if(is.null(offset)) offset=replicate(nrow(cat),0)
    if (!is.element("matrix",class(cov.matrix))) stop("WRONG FORMULA DEFINITION")
    
    ncov=ncol(cov.matrix)
    if (length(betacov)!=ncov){
      cat("Wrong number of starting values for parametrs linear predictor of covariates. Correct number of zero starting values inserted","\n")
      betacov=replicate(ncov,0)
    }
    if(length(params.ind)!=7) stop("Wrong number of elements of params.ind")
    params.ind  <- as.numeric(params.ind)
    if(sum(abs(params.ind-0.5)==0.5)!=7) stop("WRONG params.ind DEFINITION: ONLY FALSE/TRUE ALLOWED SEE HELP")
    
#    params.lim=c(0,0,0,0.0,0,0,0)  ###  original
#    params.lim=c(0,0,0,0.9999,0,0,0)  # ## changed 2021
    namespar=c("mu","k0","c","p","gamma","d","q")
    params.fix	<-c(mu,k0,c,p,gamma,d,q)
    nparams.etas<-sum(params.ind)
    nparams 		<-nparams.etas+ncov
    
    #       print(params.fix>=params.lim)
    ########################################### computation of the region for the integration used in the likelihood computation
    if(onlytime)
    {
      if(!(prod(params.fix[1:4]>=params.lim[1:4]))) stop("WRONG starting values for parameters: mu, k0, c must  be positive, p must be equal or greater than 1")
      rho.s2=0
      region=NULL}
    else
    {
      if(!(prod(params.fix>=params.lim))) stop("WRONG starting values for parameters: mu, k0, c, gamma, d, q must  be positive, p must be equal or greater than 1")
      region    =embedding.rect.cat.epsNEW(cat)
      rho.s2	=matrix(0,ntheta,n)
      
 #     for(i in (1:n)){
 #       trasf			=region.P(region,c(xcat.work[i],ycat.work[i]),k=ntheta)
 #       rho.s2[,i]		=trasf$rho
  #    }
      for(i in (1:n))        rho.s2[,i]			=region.P(region,c(xcat.work[i],ycat.work[i]),k=ntheta)$rho

      rho.s2	=t(rho.s2)
    } 
    ## initialization of other element of list etas.current
    ## added june2017
#    etas.current$formula1           =formula1 moved 19.2.2021
    etas.current$cov.matrix         =cov.matrix
    etas.current$region             =region
    etas.current$rho.s2		=rho.s2
    etas.current$offset		=offset
    
    
    ########### FIRST COMPUTATION OF BACKGROUND AT ITER=0 ###########
    ########### REWRITTEN JAN, 21 2021 ###########################
    if(is.backconstant)
    {
      #			is.backconstant==TRUE
      if(!onlytime) back.dens	<-array(1,n)/(diff(range(xcat.work))*diff(range(ycat.work))) else back.dens=1
      back.integral	<-1
    }
    else
    {
      #			is.backconstant==FALSE
      #           No matter if flp=TRUE or FALSE
      #           if missing(w)   xback.work and yback.work have been substituted with xcat.work e ycat.work
      back.tot	<-kde2dnew.fortran(xback.work,yback.work,xcat.work,ycat.work,h=hdef,factor.xy=back.smooth,eps=1/n,w=w,hvarx=hvarx,hvary=hvary) 
      etas.current$nstep.kde=etas.current$nstep.kde+1
      back.dens	<-back.tot$z
      back.integral	<-back.tot$integral
      
    }
    
    etas.current$back.integral <-back.integral
    etas.current$back.dens     <-back.dens
    
    
    #################################################	
    ########     begin   parameters estimation ##########
    ########        and clustering           ##########
    #################################################
    #while (AIC.decrease& ????????
    # attempt avoiding decrasing constrain:
    
    while ((iter<ndeclust)&((eps>epsmax)||(eps.par>epsmax))){
      
      
      #        while (AIC.decrease&(iter<ndeclust)&((eps>epsmax)||(eps.par>epsmax))){
      # params.fix has 7 elements (21-7-2017)	
      
      #check.if.del       if (iter>0) params.fix=etas.ris$params
      
      params		<-log((params.fix-params.lim)[params.ind==1])
      
      
      etas.current$params<-params
      etas.current$params.lim<-params.lim   ##### added 2021
      etas.current$params.ind	<-as.logical(params.ind)
      etas.current$params.fix	<-params.fix
      etas.current$betacov    <-betacov
      etas.current$nparams    <-nparams
      etas.current$nparams.etas<-nparams.etas
      #        etas.current$etas.first=etas.current
      
      cat("Start ML step number: ",iter+1,"\n")
      
      ### sperimental section for fastML added on feb. 3, 2021
      etas.current$fastML<-fastML
      etas.current$fast.eps<-fast.eps
      if(fastML) {
        
        fast<-(fastML.init(etas.current))
        etas.current$ind      <-fast$ind
        etas.current$index.tot<-fast$index.tot
 
      }
#      return(etas.current)
      risult.opt <- generaloptimizationNEW(etas.current,
                                           hessian	=TRUE,
                                           iterlim		=iterlim,
                                           iprint=iprint,
                                           trace=trace)
      
      etas.current    <-risult.opt$etas.obj	
      params.optim    <-risult.opt$params.optim[1:nparams.etas]
      betacov         <-risult.opt$params.optim[(nparams.etas+1):nparams]
      
      ## params.optim does not contain betacov
      cat("\n")
      
      for (i in 1:n.iterweight){
        
        cat("weighting step n. ",i,"\n")

        l.optim<-risult.opt$l.optim
        risult <-risult.opt$risult
        
        l<-etas.mod2NEW(params=c(params.optim,betacov),
                        etas.obj=etas.current,
                        trace=FALSE)
        
        ####################("found optimum ML step") #################################################################
        det.check<-det(risult$hessian)
        check=is.na(det.check)||(abs(det.check)<1e-20)
        if (check) compsqm<-FALSE
        sqm<-0
        
        if (compsqm)	sqm	<-sqrt(diag(solve(risult$hessian)))*c(exp(params.optim),replicate(ncov,1))
        sqm.etas<-sqm[1:nparams.etas]
        sqm.cov <-sqm[(nparams.etas+1):nparams]
        
        ###### the optimization is made with respect to the log of the parameters ()
        ###### so the for the asymptotic standard error we use the approximation Var(exp(y))=[exp(y)]^2 var(y)
        params			<-params.fix
        params[params.ind==1]	 <- exp(params.optim)+params.lim[params.ind==1]
        sqm.tot			<-array(0,nparams.etas)
        sqm.tot[params.ind==1] <-  sqm.etas
        names(sqm.tot) <- namespar
        
        params.MLtot	<-c(params,betacov)
        sqm.MLtot	    <-c(sqm.tot,sqm.cov)
        #####################################################################################
        cat("found optimum; end  ML step  ")
        cat(iter+1,"\n")
        
        mu	<- params[1]
        k0  <- params[2]
        c   <- params[3]
        p   <- params[4]
        #        a   	= params[5]
        gamma <- params[5]
        d   <- params[6]
        q   <- params[7]
        #        betacov = params[8:nparams]
        
        #####################################################################################
        # time.res= residuals obtained with integral tansform of time intensity function
        # to be used only for time processes
        # corrected for predictor in version 2.0.0
        # maybe could be omitted later, since plot.etasclass already computes residuals
        #####################################################################################
        predictor       <- as.matrix(etas.current$cov.matrix)%*%as.vector(betacov)+as.vector(etas.current$offset)
        rho.weights     <- params.MLtot[1]*back.dens/attr(l,"lambda.vec")
        w               <- rho.weights
        xback.work <- xcat.work
        yback.work <- ycat.work
        
        
        if (flp){
          etas.l		<- attr(l,"etas.vec")
          # compute etas intensity and integral for each point and call routine for optimization
          # flp MUST be weighted 
          ris.flp<-flp1.etas.nlmNEW(cat,
                                    h.init=hdef,
                                    etas.params=params.MLtot,
                                    etas.l=etas.l,
                                    w=rho.weights,
                                    m1=m1,
                                    m2=as.integer(nrow(cat)-1),
                                    mh=1
          )
          etas.current$nstep.flp=etas.current$nstep.flp+1
          
          hdef <- ris.flp$hdef
          fl	 <- ris.flp$fl
        }
        ### 	end of flp step 
        
        if (thinning){
          back.ind	 <- runif(n)<rho.weights
          xback.work <- xcat.work[back.ind]
          yback.work <- ycat.work[back.ind]
          w	<-replicate(length(xback.work),1)
        }
        
        
       if (flp) back.tot		<- kde2dnew.fortran(xback.work,yback.work,xcat.work,ycat.work,eps=1/n,h=hdef,w=w,hvarx=hvarx,hvary=hvary) else back.tot		<-kde2dnew.fortran(xback.work,yback.work,xcat.work,ycat.work,eps=1/n,w=w,hvarx=hvarx,hvary=hvary)
        etas.current$nstep.kde=etas.current$nstep.kde+1
        
        back.dens	  <- back.tot$z
        back.integral<-back.tot$integral
        wmat		<-back.tot$wmat
        
        
        ###### end of  background computation
        etas.current$rho.weights  <-rho.weights
        etas.current$back.integral  <-back.integral
        etas.current$back.dens      <-back.dens
   
        ##### 
        l<-etas.mod2NEW(params=c(params.optim,betacov),
                        etas.obj=etas.current,
                        trace=FALSE)
        AIC2   		  <- 2*l +2*nparams
        
             
      }
      
      #### time residual for time only models
      time.res	<- array(0,n)
      if (onlytime){
        times.tot	<- cat$time
        magnitudes.tot	<- cat$magn1-magn.threshold
        tmin	<- min(times.tot)
        etas.t		<- 0
        for (i in 1:n) {
          tmax.i  	<- times.tot[i]
          times			<- subset(times.tot,times.tot<tmax.i)
          magnitudes	<- subset(magnitudes.tot,times.tot<tmax.i)
          
          if(p==1){	
            it	<- log(c+tmax.i-times)-log(c)
          }
          else
          {
            it	<- ((c+tmax.i-times)^(1-p)-c^(1-p))/(1-p)
          }
          
          #	time.res[i]=(tmax.i-tmin)*mu+k0*sum(exp(a*magnitudes)*it)
          time.res[i]	<- (tmax.i-tmin)*mu+k0*sum(exp(predictor)*it)
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
      ## maybe useless dropped on 26-1-2021
      #            rho.weights=params.MLtot[1]*back.dens/attr(l,"lambda.vec")
      
      # check convergence: back densities, AIC, parameters
      
      timenow		<-	Sys.time()
      time.elapsed	<-	difftime(timenow,etas.current$time.start,units="secs")
      AIC.temp		  <- 2*l +2*nparams
#      AIC.iter[iter]		<- AIC.temp
      AIC.decrease	<- ((iter==0)||(AIC.temp<=min(AIC.iter))) 
#      cat("time.elapsed",time.elapsed, "; AIC=", AIC.temp, "; AIC.decrease=", AIC.decrease)
#      print(Sys.time())
      if (AIC.decrease){
        iter			<- iter+1
        AIC.iter[iter]		<- AIC.temp
        AIC.iter2[iter]		<- AIC2
        params.iter		<- rbind(params.iter,params.MLtot)
        sqm.iter	  	<- rbind(sqm.iter,sqm.tot)
        rho.weights.iter		<- rbind(rho.weights.iter,rho.weights)
        hdef.iter	        	<- rbind(hdef.iter,hdef)
        
        fl.iter			<- c(fl.iter,fl)
        
        cat("\n","---ITERATION n.", iter," ---; AIC = ",round(AIC.iter[iter],3),"; elapsed time: ",time.elapsed," time: ")
        print(Sys.time())
        cat("\n")
        cat("Current estimates of parameters: ","\n")
        cat(round(params.MLtot,5))
        cat("\n")
      }
      if (iter>1){
        eps			<- max(abs(back.dens-etas.ris$back.dens))
        eps.par	<- max(abs(params.iter[iter]-params.iter[iter-1]))
      }
      rownames(params.iter) <- 1:nrow(params.iter)
      
      predictor       <- as.matrix(etas.current$cov.matrix)%*%as.vector(betacov)+as.vector(etas.current$offset)
      
      names(params.MLtot)<-c(namespar,names(betacov))
      names(sqm.MLtot)   <-names(params.MLtot)
      
      etas.current$time.elapsed   <-	time.elapsed
      etas.current$time.end		    <-	timenow
      etas.current$risult		      <-	risult
      etas.current$betacov	      <-betacov
      etas.current$params.MLtot	  <-params.MLtot
      etas.current$params	        <-params.MLtot[1:7]
      etas.current$predictor	   	<-predictor
      etas.current$sqm		        <-sqm.MLtot
      etas.current$time.res	      <-time.res
      etas.current$AIC.iter	      <-AIC.iter
      etas.current$AIC.iter2	      <-AIC.iter2
      etas.current$AIC.decrease	  <-AIC.decrease
      etas.current$params.iter	  <-params.iter
      etas.current$sqm.iter	      <-sqm.iter
      etas.current$rho.weights.iter<-rho.weights.iter
      etas.current$rho.weights    <-rho.weights
      etas.current$n.iterweight    <-n.iterweight
      etas.current$hdef		        <-hdef
      etas.current$hdef.iter	    <-hdef.iter
      etas.current$back.integral  <-back.integral
      etas.current$back.dens      <-back.dens
      etas.current$wmat		        <-wmat
      etas.current$fl.iter		    <-fl.iter
      etas.current$iter		        <-iter
      etas.current$usenlm		        <-usenlm
      etas.current$method		        <-method
      etas.current$epsmax		        <-epsmax
      etas.current$iterlim		        <-iterlim
      etas.current$compsqm		        <-compsqm
      etas.current$ntheta		        <-ntheta
      etas.current$logl		        <-l        
      etas.current$l              <-attr(l,"lambda.vec")       
      etas.current$integral	      <-attr(l,"integraltot")        
      
      etas.ris                    <-etas.current
      if(!AIC.decrease){
        print("INCREASING AIC")
        return(etas.ris)
      }
    }
    # end of general iteration for declustering
    
    #############################################################################################
    #
    #		convergence obtained. Computation of final output 
    #
    #############################################################################################
    if(!time.update) {
      cat("\n")
      cat("\n")
      cat("--------- FINAL SUMMARY ---------","\n")
      
      summary(etas.ris)}
    return(etas.ris)
  }
#####   check of january 14, 2021
#####   still the last rows of rho.weights.iter,params.iter and hdef.iter do not coincide with the final output values.
#####   probably the first rho.weight is not the input one. check  




