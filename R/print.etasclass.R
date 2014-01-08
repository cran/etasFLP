print.etasclass<-function(x,...){
print(x$description)
cat("Execution started:                 ",format(x$time.start),"\n")
cat("Elapsed time of execution (hours)  ",as.numeric(x$time.elapsed,units="hours"),"\n")
cat("Number of observations            ",length(x$cat$time),"\n")
cat("Magnitude threshold               ",x$magn.threshold,"\n")

print(summary(x$cat.longlat))

## output also declustering=TRUE|FALSE. the following is only if declustering=TRUE
declustering=TRUE
if (declustering){
cat("Number of declustering iterations  ",x$iter,"\n")
cat("Kind of declustering               ",ifelse(x$thinning,"thinning","weighting"),"\n")

print("sequence of AIC values for each iteration")
print(x$AIC.iter)
}

print("ETAS Parameters")
ris=cbind(x$params,x$sqm)
colnames(ris)=c("      Estimates","      std.err.")
print(round(ris,6))
}
