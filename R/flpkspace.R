flpkspace <-
function(theta,x,t,w,
			k	=2,
			m1=as.integer(nrow(x)/2),
			m2=as.integer(nrow(x)-1),
			etas.params,
			etas.l,
			etas.integral,
			mh=1,

			iprint=FALSE,
			indweight=TRUE)     {
# function called (for optimization) by flp.etas.nlm
# essentially is an interface to the FORTRAN  subroutine deltafl1kspace
		cat("-")	
		n	=as.integer(nrow(x))
		dens	=array(0,n)
		 s	=matrix(c(1,0,0,1),2,2)
 rangex=t(apply(x,2,range))

	ris.fl<-.Fortran("deltafl1kspace",
		x	=as.double(x),
		t	=as.double(t),
		w	=as.double(w),
		n	=as.integer(n),
		k	=as.integer(k),
		m1	=as.integer(m1),
		m2	=as.integer(m2),
		rangex	=as.double(rangex),
		h	=as.double(exp(theta[1:k])),
		hdef	=as.double(exp(theta[1:k])),
		dens	=as.double(dens),
		integr	=as.double(dens),
		delta	=as.double(dens),
		indanis	=as.integer(0),
		expweight=as.double(-0.2),
		indweight=as.integer(indweight),
		sigma	 =as.double(s),
		NAOK=TRUE
				)

		lambda	=etas.params[1]
		val	=-sum(log(lambda*ris.fl$dens[m1:m2])+log(etas.l[(m1+1):(m2+1)])
		-etas.integral[m1:m2]-lambda*ris.fl$integr[m1:m2])
                attr(val, "dens") <- ris.fl$dens
                attr(val, "delta") <- ris.fl$delta
                attr(val, "integr") <- ris.fl$integr
                attr(val, "hdef") <- ris.fl$hdef
		return(val)
 }
