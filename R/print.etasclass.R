print.etasclass<-function(x,...){
    if(class(x)!="etasclass")stop(" argument x must be an etasclass object")
    cat(x$description,"\n")
    cat("Execution started: ",format(x$time.start))
    cat("Elapsed time of execution (hours)  ",round(as.numeric(x$time.elapsed,units="hours"),3),"\n")
    cat("Number of observations: ",length(x$cat$time))
    cat(";  Magnitude threshold ",x$magn.threshold)
    
    cat("; AIC",x$AIC.iter,"; hdef",x$hdef,"\n") ### added january 2021
    
    cat("ETAS Parameters:","\n")
    ris=x$params.MLtot
    print(round(ris,6))
}
