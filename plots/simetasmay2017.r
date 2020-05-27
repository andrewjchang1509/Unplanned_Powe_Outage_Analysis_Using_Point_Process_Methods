## General ETAS simulator where the user can input densities.
    theta = list(mu=8.0,K=.3,c=.2,p=2.3,a=.1,d=.3,q=2.7,b=.7)){
    ##### If no magnitudes are desired, just let gmi = K. 
    calcbr = mean(gmi(mdensity(1000000,theta),theta)) ## calculate branching ratio. Stop if br > 1. 
    if(calcbr>1.0){
    cat("error, branching ratio = ", calcbr, " > 1.")
    return(0)
    }
    stop1 = 0
    w = y
    z1 = list()

## density = b e ^ -bm. cdf = 1 - e^-bm. m means m - m0. 
## u = unif(0,1). F(x) = u. 1 - e^-b(m-m0) = u. 
## Solve for m. 
## e^-b(m-m0) = 1-u. 
## -b(m-m0) = log(1-u). 
## m - m0 = log(1-u)/-b. 
## m = -log(1-u)/b + m0. 

expmag = function(n,theta,m0=3.5){ 
-log(1-runif(n))/theta$b + m0
}
    ## Ogata (1998). See http://wildfire.stat.ucla.edu/pdflibrary/ogata98.pdf . 
    ## ∫f(x,y)dxdy = 1 = ∫h(r)rdrdø = 2π∫h(r)rdr. 
    ## h(r) = c (r^2 + d)^(-q). 
    ## ∫ h(r)rdr = c(r^2+d)^(1-q)/(2-2q),r=0to∞. For q > 1, this is 0+cd^(1-q)/(2q-2).
    ## So c = (q-1)d^(q-1)/π. 
    v = runif(n)
    ## Here the density does not depend on magnitude of the mainshock. 
    ## To see that this is a density, 
    ## ∫f(x,y)dxdy = ∫f(r)rdrdø = 2π∫f(r)rdr 
    ## = 2alpha ∫ exp(-alpha r^2) r dr = -exp(-alpha r^2) ,r=0to∞, = 0+1, for alpha>0.  
    v = rexp(n,rate=theta$alpha)

## Examples. 
b = simhawk()
b = simhawk(T=10)
b = simhawk(T=10,theta = list(mu=8.0,K=.3,c=.2,p=2.3,a=.1,d=.3,q=2.8,b=.7))
mags = (b$m - min(b$m))/diff(range(b$m))
par(mfrow=c(1,2))
plot(b$lon,b$lat,pch=3,cex = 1+3*mags,xlab="longitude (km)",ylab="latitude (km)")
library(gplm)
library(MASS)
t2 = powergt(100000,theta=list(c=.2,p=2.3))
hist(t2[t2<1],nclass=100,prob=T,main="",xlab="time between pts")
lines(kde(t2[t2<1]),col="green")

g2 = powerxy(100000,1,theta=list(d=.3,q=2.8))
g3 = kde2d(g2[,1],g2[,2],lims=c(-1,1,-1,1))
image(g3)
contour(g3,add=T)