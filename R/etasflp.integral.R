etasflp.integral <-
function(params=c(1,1,1,1,1,1,1,1),
				params.fix	=array(1,8),
				params.ind	=array(1,8),
				params.lim	=c(0,0,0,0,0,0,0,0),
				cat,
				rho.s2,
				ntheta=100,
				magn.threshold,
				trace=TRUE) {

	params.e	=params.fix
	params.e[params.ind==1]=exp(params)+params.lim[params.ind==1]

	lambda	= params.e[1]
        k0  	= params.e[2]
        c   	= params.e[3]
        p   	= params.e[4]
        a   	= params.e[5]
        gamma   = params.e[6]
        d   	= params.e[7]
        q   	= params.e[8]
	x.tot	=cat$long
	y.tot	=cat$lat
	magnitudes.tot	=cat$magn1-magn.threshold
	times.tot	=cat$time
	n		=length(times.tot)
	m2=n-1
	m1=n/2
intflp=0*times.tot
for (i in m1:m2){
tmax1		=times.tot[i+1]
tmax0		=times.tot[i]
ind1		=1:i
times		=times.tot[ind1]
magnitudes	=magnitudes.tot[ind1]
rho.m1		=rho.s2[ind1,]

##############################################################################
# starting time integration#
##############################################################################


if(p==1){
			 it=log(c+tmax1-times)-log(c+tmax0-times)
			}
			else
			{
			 it=((c+tmax1-times)^(1-p)-(c+tmax0-times)^(1-p))/(1-p)
			}

##############################################################################
# starting space integration (polar transformation)
##############################################################################


mag1	=as.vector(exp((a-gamma)*magnitudes))
mag2	=as.vector(exp(gamma*magnitudes))



### approximate polar integration to whole space by division in ntheta triangles centered in xi,yi

etas		=rowSums(mag1*mag2*((rho.m1*rho.m1/mag2+d)^(1-q)-d^(1-q)))
space		=(pi/((1-q)*ntheta))*etas
intflp[i]	=sum(k0*it*space)
}

return(intflp)
}
