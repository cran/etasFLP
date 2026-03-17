fastML.init <- function(etas.obj){
    params.ind   =etas.obj$params.ind
    params.fix   =etas.obj$params.fix
    nparams.etas  =etas.obj$nparams.etas
    nparams       =etas.obj$nparams
    ntheta          =etas.obj$ntheta		
    params<-etas.obj$params
    params.lim	=c(0,0,0,0.9999,0,0,0)
    fast.eps=etas.obj$fast.eps

    params.etas    =params[1:nparams.etas]
    betacov        =etas.obj$betacov
    params.e	=params.fix
    params.e[params.ind==1]=exp(params.etas)+params.lim[params.ind==1]
    lambda	= params.e[1]
    k0  	= params.e[2]
    c   	= params.e[3]
    p   	= params.e[4]
    gamma   = params.e[5]
    d   	= params.e[6]
    q   	= params.e[7]
    
      predictor       =as.matrix(etas.obj$cov.matrix)%*%as.vector(betacov)+as.vector(etas.obj$offset)

      
    eta=as.numeric(predictor)
    x		=etas.obj$cat$xcat.work
    y		=etas.obj$cat$ycat.work
    t		=etas.obj$cat$time.work
    m	=etas.obj$cat$magnitude
    n     <-length(x)
    ### without outer product
    qq=list()
   
    
    for(i in 2:n){
        dx=x[i]-x[1:(i-1)] 

        dy=y[i]-y[1:(i-1)]    
        dt1=abs(t[i]-t[1:(i-1)])    
        ds=1/((dx*dx+dy*dy)/exp(gamma*m[i])+d)^q
        dt=ds*exp(eta[i])*(dt1+c)^(-p)
        qq[[i]]=max(dt)
        print(c(i,qq[[i]]))
#        qq[[i]]=max(dt)*k0+lambda*etas.obj$back.dens[i]
    }
    
    maxq=max(unlist(qq))

    ind=array(0,n) 
    index.parz=numeric(0)
    index.tot=numeric(0)
    for(i in 2:n){
      maxq=as.numeric(qq[[i]])
        dx=x[i]-x[1:(i-1)]    
        dy=x[i]-x[1:(i-1)]    
        dt=abs(x[i]-x[1:(i-1)])    
        ds=1/((dx*dx+dy*dy)/exp(gamma*m[i])+d)^q
        dt=ds*(dt+c)^(-p)*exp(eta[i])
#        dt=dt*k0+lambda*etas.obj$back.dens[i]
        index.parz=which(dt>fast.eps*maxq)
        ind[i]=ind[i-1]+length(index.parz)
        index.tot=c(index.tot,index.parz)
        cat(c(i,length(index.parz)),";")
        
    }

    
    
    
    
    return(list(index.tot=index.tot,ind=ind))
}






