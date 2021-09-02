#' @title timeupdate.etasclass
#' @alias timeupdate.etasclass
#' Title time update for etasclass objects
#' @description New in version 2.2. A time updating of an etasclass objects: a very experimental version that can be used only on etasclass objects obtained from etasflp versions 2.2 or newer. 
#' @details It is a beta version. 
#' A new ETAS model is fitted to a previous object of class etasclass with a new catalog which must be a catalog which extends the previous one on a wider time window, that is a catalog with new observations.  
#' 
#' As a default a new quick execution is made, with one quick iteration for parameter updating and an iteration for background density estimation.
#' @param object an etasclass object  obtained from \code{etasFLP} versions 2.2 or newer that will be updated for a new time window with new events.
#' @param params.estimation logical. if \code{TRUE} parameters will be  estimated again even if quickly, with few optimizations steps. Elsewhere ETAS estimates of input object will be mantained
#' @param ... optional arguments that will override the corresponding arguments in \code{object}, possibly including a new catalog input or a new \code{tmax}
#'
#' @return a new \code{etasclass} object
#' @export
#' @seealso \code{\link{update.etasclass}}
#' @examples
timeupdate.etasclass <-
    function(object,params.estimation=FALSE, ...)   {
        if(class(object)!="etasclass")stop("argument must be an etasclass object")
        class(object$cat.orig)="data.frame"
        updated.object<-NULL
        arg.list<-"cat.orig=object$cat.orig,time.update=TRUE,magn.threshold=object$magn.threshold,magn.threshold.back=object$magn.threshold.back,long.range=object$long.range,lat.range=object$lat.range,offset=object$offset,hvarx=object$hvarx,hvary=object$hvary,declustering=object$declustering,thinning=object$thinning,flp=FALSE,ndeclust=1,n.iterweight=3,onlytime=object$onlytime,is.backconstant=object$is.backconstant,sectoday=object$sectoday,longlat.to.km=object$longlat.to.km,usenlm=object$usenlm,compsqm=object$compsqm,epsmax=object$epsmax,iterlim=8,ntheta=object$ntheta,method=object$method,mu=object$params[1],k0=object$params[2],c=object$params[3],p=object$params[4],gamma=object$params[5],d=object$params[6],q=object$params[7],hdef=object$hdef[1:2],betacov=object$betacov,w=object$rho.weights,formula1=object$formula1"
        arg.input       <-strsplit(arg.list,",")[[1]]
        call1           <-deparse(match.call(),width.cutoff=500)
        
        arg.dots        <-strsplit(call1,",",fixed=TRUE)[[1]]     
        
        name.obj        <-strsplit(arg.dots[1],split="(",fixed=TRUE )[[1]][2]
        nn              <-strsplit(name.obj,"=")[[1]]
        name.obj        <-trimws(nn)[length(nn)]
        name.obj        <-substr(name.obj,1,nchar(name.obj)-1)
        n               <-length(arg.dots)
        arg.dots[n]     <-substr(arg.dots[n],1,nchar(arg.dots[n])-1)
        arg.dots        <-trimws(arg.dots[-1])
        arg1=substr(arg.dots[1],1,10)
        if(n>1){
        if(substr(arg.dots[1],1,10)=="params.est"){
            expr.par    <-arg.dots[1]
            eval(parse(text=expr.par))        ##############
            
            arg.dots    <-arg.dots[-1]
                n       <-n-1
        }
        if(n>1){
            n            <-n-1
            name.dots    <-trimws(matrix(unlist(strsplit(arg.dots,"=")), ncol = 2, byrow = TRUE))[,1]
            name.input   <-trimws(matrix(unlist(strsplit(arg.input,"=")), ncol = 2, byrow = TRUE))[,1]
            ind          <-numeric(0)
            
            for (i in 1:n)    ind<-c(ind,which(name.input==name.dots[i]))
            arg.input[ind] <-arg.dots
        }
        }
        if(params.estimation) params.warn=" with quick reestimation of ETAS parameters " else params.warn=" without  reestimation of ETAS parameters " 
        arg.input    <-c(arg.input,paste("description= 'Time updating of ",name.obj," of date ",as.character(Sys.time()),params.warn," '"))
        e1           <-paste(arg.input,collapse=",")
        if (params.estimation){
            expr1        <-paste(e1,",params.ind=object$params.ind)",sep="")
        }
        else{
            expr1        <-paste(e1,",params.ind=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE))",sep="")
        }
        expr1        <-paste("updated.object=etasclass",expr1,sep="(")
        eval(parse(text=expr1))        ##############
        if(!params.estimation) updated.object$sqm[1:7]=object$sqm[1:7]
        summary(updated.object)
        return(updated.object)
    }
