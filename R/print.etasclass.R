print.etasclass<-function(x,...){
	  if(class(x)!="etasclass")stop(" argument x must be an etasclass object")

cat("Call:","\n","\n")
cat(x$this.call)
cat("\n","\n","\n")
cat(x$description,"\n")
cat("Number of observations            ",length(x$cat$time),"\n")

cat("ETAS Parameters:","\n")
ris=x$params
cat(round(ris,6),"\n")
}
