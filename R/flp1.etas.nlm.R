flp1.etas.nlm <-
function(cat,
			logh.init=c(0,0),
			w=w,
			etas.params,
			etas.l,
			etas.integral,
			m1=as.integer(nrow(cat)/2),
			m2=as.integer(nrow(cat)-1),
			mh=1
			)
				    {
### compute the optimal bandwidth for an etas model according to the flp approach
# timescale=86400:
# for default catalogues are in seconds, our computations are in days

## in this version only isotropic x-y normal kernel
## only nlm function used to optimize flpkspace with respect to hx and hy

print("begin flp step")


time.init=Sys.time()
niter=1

k=2
x=cbind(cat$long,cat$lat)
t=cat$time


if(missing(logh.init)) logh.init=c(0,0)
else logh.init=logh.init[1:k]
    npar=k
theta.init=logh.init
print(theta.init)

gradtol=0.01
ris	=nlm(flpkspace,theta.init,x=x,t=t, w=w,m1=m1,m2=m2,mh=mh,k=k,etas.l=etas.l,
etas.params=etas.params,

etas.integral=etas.integral,hessian=TRUE
)

fl =
flpkspace( ris$estimate,x=x,t=t,
w=w,m1=m1,m2=m2,mh=mh,k=k,etas.l=etas.l,
etas.params=etas.params, etas.integral=etas.integral )

hdef=attr(fl, "hdef")
print("exit from flp step")
return(list(hdef=hdef,fl=fl))
 }



