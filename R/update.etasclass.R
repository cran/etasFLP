#' @title update.etasclass
#' @alias update.etasclass

#' Title Method update for etasclass objects
#' @description New in version 2.2. A method update for etasclass objects: a very experimental version that can be used only on etasclass objects obtained from etasflp versions 2.2 or newer. 
#' @details It is a beta version. The catalog must be the same, and options in "..." must leave unchanged the number of observations used for estimation. Arguments given in "..." will override arguments already present in \code{object}. Not all arguments are suitable for updating: among them \code{formula} and \code{params.ind} should not be included in "..." list (to update such parameters it is better to assign them to a variable and then pass the variable name) . A new etasclass execution will start, using as arguments  values of input \code{object}, eventually integrated with the list in "...". Tipically a first execution can be given with low values of iterlim, ndeclust, ntheta and high values of epsmax (e.g. iterlim=5, ndeclust=1, ntheta=24, epsmax=0.01), to obtain good starting values for parameters and for weights. Then an update can be run with better values such as iterlim=50, ndeclust=10, ntheta=60, epsmax=0.0001.
#' @param object an etasclass object  obtained from etasflp versions 2.2 or newer that will be updated.
#' @param ... optional arguments that will override the corresponding arguments in \code{object}
#'
#' @return un updated \code{etasclass} object
#' @export
#' @seealso \code{\link{timeupdate.etasclass}}
#' @examples
update.etasclass <-
function(object, ...)   {
    if(class(object)!="etasclass")stop("argument must be an etasclass object")
    class(object$cat.orig)<-"data.frame"
    arg.list<-"cat.orig=object$cat.orig,magn.threshold=object$magn.threshold,magn.threshold.back=object$magn.threshold.back,tmax=object$tmax,long.range=object$long.range,lat.range=object$lat.range,offset=object$offset,hvarx=object$hvarx,hvary=object$hvary,declustering=object$declustering,thinning=object$thinning,flp=object$flp,ndeclust=object$ndeclust,n.iterweight=object$n.iterweight,onlytime=object$onlytime,is.backconstant=object$is.backconstant,sectoday=object$sectoday,longlat.to.km=object$longlat.to.km,usenlm=object$usenlm,compsqm=object$compsqm,epsmax=object$epsmax,iterlim=object$iterlim,ntheta=object$ntheta,method=object$method,mu=object$params[1],k0=object$params[2],c=object$params[3],p=object$params[4],gamma=object$params[5],d=object$params[6],q=object$params[7],hdef=object$hdef[1:2],betacov=object$betacov,params.ind=object$params.ind,w=object$rho.weights,formula1=object$formula1"
    updated.object  <-NULL
    arg.input       <-strsplit(arg.list,",")[[1]]
    call1           <-deparse(match.call(),width.cutoff=500)
    arg.dots        <-strsplit(call1,",",fixed=TRUE)[[1]]     
    name.obj        <-strsplit(arg.dots[1],split="(",fixed=TRUE )[[1]][2]
    nn              <-strsplit(name.obj,"=")[[1]]
    name.obj        <-trimws(nn)[length(nn)]
    name.obj        <-substr(name.obj,1,length(name.obj)-1)
    n               <-length(arg.dots)
    
    if(n>1){
        arg.dots[n]  <-substr(arg.dots[n],1,nchar(arg.dots[n])-1)
        arg.dots     <-trimws(arg.dots[-1])
        n            <-n-1
        name.dots    <-trimws(matrix(unlist(strsplit(arg.dots,"=")), ncol = 2, byrow = TRUE))[,1]
        name.input   <-trimws(matrix(unlist(strsplit(arg.input,"=")), ncol = 2, byrow = TRUE))[,1]
        ind          <-numeric(0)
        for (i in 1:n)    ind<-c(ind,which(name.input==name.dots[i]))
        arg.input[ind] <-arg.dots
    }
    arg.input    <-c(arg.input,paste("description= 'updating of ",name.obj," of date ",as.character(Sys.time())," '"))
    e1           <-paste(arg.input,collapse=",")
    expr1        <-paste(e1,")",sep="")
    expr1        <-paste("updated.object=etasclass",expr1,sep="(")
    eval(parse(text=expr1))        ##############
    return(updated.object)
	
}
