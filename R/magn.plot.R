magn.plot <-
function(catalog,main="Transformed plot of magnitude frequencies",...){mfreq=MLA.freq(catalog$magn1)
                        plot(mfreq$x,log(mfreq$back.cum.rel),
			    xlab="magnitude", ylab="Log-number events over a value",main=main
			    ,...)
                        lines(mfreq$x,log(mfreq$back.cum.rel) ,...)
	
			    }
