print.etasclass<-function(x,...){
cat("Call:","\n","\n")
print(x$this.call)
cat("\n","\n")
cat(x$description,"\n")
cat("Number of observations            ",length(x$cat$time),"\n")

cat("ETAS Parameters:","\n")
ris=x$params
print(round(ris,6))
}
