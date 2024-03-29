flp1.etas.nlmNEW <-
function(cat,
			h.init,
			w,
			etas.params,
			etas.l,
			m1=as.integer(nrow(cat)/2),
			m2=as.integer(nrow(cat)-1),
			mh=1
			)
				    {
### compute the optimal bandwidth for an etas model according to the flp approach

## only nlm function used to optimize flpkspace with respect to hx and hy


time.init=Sys.time()
niter=1

k=2
x=cbind(cat$xcat.work,cat$ycat.work)
inew=TRUE
if(inew){
xvec=cat$xcat.work
mean1=weighted.mean(xvec,w)
s1=sqrt(weighted.mean((xvec-mean1)^2,w))
xvec=cat$ycat.work
mean1=weighted.mean(xvec,w)
s2=sqrt(weighted.mean((xvec-mean1)^2,w))
h.init[1]=h.init[1]*(m2^0.2)/s1
h.init[2]=h.init[2]*(m2^0.2)/s2

}
#h.init[3]=h.init[3]*(m2^0.2)/s1
#h.init[4]=h.init[4]*(m2^0.2)/s2

t=cat$time.work
theta.init=h.init*0
    npar=k*(1)

    theta.init[1:npar]=log(h.init[1:npar])

    cat("Start FLP step: ")


#ris	=nlm(flpkspace,theta.init,x=x,t=t, w=w,m1=m1,m2=m2,mh=mh,k=k,etas.l=etas.l,etas.params=etas.params,etas.integral=etas.integral,hessian=TRUE,kern.var=kern.var)
ris	=nlm(flpkspace,theta.init,x=x,t=t, w=w,m1=m1,m2=m2,mh=mh,k=k,etas.l=etas.l,etas.params=etas.params,hessian=TRUE)

#fl =flpkspace( ris$estimate,x=x,t=t,w=w,m1=m1,m2=m2,mh=mh,k=k,etas.l=etas.l,etas.params=etas.params, etas.integral=etas.integral,kern.var=kern.var )

fl =flpkspace(ris$estimate,x=x,t=t,w=w,m1=m1,m2=m2,mh=mh,k=k,etas.l=etas.l,etas.params=etas.params)

hdef=attr(fl, "hdef")
cat("exit from FLP step...","\n")
return(list(hdef=hdef,fl=fl))
 }



