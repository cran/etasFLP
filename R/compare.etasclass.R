## check (22-4-2020) with the new etasclass output
compare.etasclass <-
function(etas1,etas2){
# compare two etasclass objects
#check for classes
	  if (!inherits(etas1,"etasclass"))stop("first argument must be an etasclass object")
    if (!inherits(etas2,"etasclass"))stop("second argument must be an etasclass object")
	  ### check for comparability of the two input objects
	  ## same threshold
	  ## same catalog
	  ## same domain (space time)
	  ##
	  params	=0
	  npar.est	=etas1$params.ind+etas2$params.ind
	  npar.est  =c(npar.est,array(2,length(etas1$betacov)))
	  sqm		=sqrt((etas1$sqm^2+etas2$sqm^2)/npar.est)
	  params=array(0,length(etas1$params.MLtot))
	  params[sqm>0]=(etas1$params.MLtot-etas2$params.MLtot)/sqm
	  AIC=min(etas1$AIC)-min(etas2$AIC)
	  if(etas1$onlytime&etas2$onlytime){
	  weights	=0
	  weights.std	=0
	  cor.weights	=0
	  }
	  {
	  weights	=etas1$rho.weights-etas2$rho.weights
	  weights.std	=(etas1$rho.weights-etas2$rho.weights)/sqrt((etas1$rho.weights*(1-etas1$rho.weights)+etas2$rho.weights*(1-etas2$rho.weights)))
	  cor.weights	=cor(etas1$rho.weights,etas2$rho.weights)
	  }
	  cor.trig	=cor(etas1$l-etas1$params[1]*etas1$back.dens,etas2$l-etas2$params[1]*etas2$back.dens)
	  cor.back	=cor(etas1$back.dens,etas2$back.dens)
return(list(diffstd.params=params,AIC=AIC,weights=weights,weights.std=weights.std,cor.weights=cor.weights,cor.trig=cor.trig,cor.back=cor.back))
  
	  }
