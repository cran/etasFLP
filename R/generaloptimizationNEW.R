		
generaloptimizationNEW=function(etas.obj,
                hessian,
                iterlim,
  		iprint,
                trace
                ){

n       =nrow(etas.obj$cat)
params  =c(etas.obj$params,etas.obj$betacov)





print("params in general optimization")
print(params)
print(etas.obj$nparams)
print(etas.obj$nparams.etas)
                
                                  etas.obj$dettrasf=replicate(n,1)

               
 
 
 
if (etas.obj$usenlm){risult =nlm(etas.mod2NEW,params,
		typsize=abs(params),
		hessian=hessian,
		iterlim		=iterlim,
		etas.obj=etas.obj
		)
		params.optim=risult$estimate
		l.optim	 =risult$minimum}
### check optim arguments
else {	   
	
risult =optim(params,etas.mod2NEW,
		method	=etas.obj$method,
		hessian	=hessian,		control=list(trace=2,maxit=iterlim,fnscale=n/diff(range(etas.obj$cat$time.work)),parscale=sqrt(exp(params))),
 		etas.obj=etas.obj
)
		
   params.optim=risult$par

   l.optim	 =risult$value	}
   
   
return(list(params.optim=params.optim,l.optim=l.optim,risult=risult,etas.obj=etas.obj))
}
