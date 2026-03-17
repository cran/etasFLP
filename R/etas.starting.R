####### guess starting values (beta-version) 
etas.starting = function(cat.orig,
magn.threshold=2.5,
p.start=1,
gamma.start=0.5,
q.start=2,
betacov.start=.7,
longlat.to.km=TRUE,
sectoday=FALSE,
onlytime=FALSE
)
{
cat=cat.orig[cat.orig$magn1>magn.threshold,]
if(sectoday)cat$time=cat$time/86400
t=cat$time

a=cat[order(t),]
t=a$time
		if (longlat.to.km){
		
		radius=6371.3
		

		ycat.km      =   radius*cat$lat*pi/180
		xcat.km      =   radius*cat$long*pi/180


		}
else
{
    ycat.km      =   cat$lat
    xcat.km      =   cat$long
    
}

n=nrow(a)
dt=diff(t)
ds=diff(ycat.km)^2+diff(xcat.km)^2
c.start=as.numeric(quantile(dt,0.25))
d.start=as.numeric(quantile(ds,0.05))
mu.start=n*0.5/diff(range(t))
# approximate evaluation of integrals
# 
if (onlytime){
    gamma.start=0
    q.start=0
    d.start=0
}

tmax=max(t)
it=log((c.start+tmax-t)/c.start)
em=exp((betacov.start[1])*(a$magn1-magn.threshold))
is=d.start^(1-q.start)*exp(gamma.start*(a$magn1-magn.threshold))*pi/(q.start-1)

if (onlytime) k0.start=n*0.5/sum(em*it) else k0.start=n*0.5/sum(em*it*is) 

return(list(
mu.start=mu.start,
k0.start=k0.start,
c.start=c.start,
p.start=p.start,
gamma.start=gamma.start,
d.start=d.start,
q.start=q.start,
betacov.start=betacov.start,
longlat.to.km=longlat.to.km,
sectoday=sectoday
))
}


