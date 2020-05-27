library(tidyverse)
library(spatstat)
source('/Users/andrewchang/Desktop/222/plots/simetasmay2017.r')
check <- read_csv("/Users/andrewchang/Desktop/222/data/check.csv")
### Notice that these (x, y) points were the original geometric points scaled to
### [0,1] range
points <- with(check, ppp(x,y,c(0,1),c(0,1)))

png("l.png", width = 1200, height = 1200, res = 150)
k = Kest(points, correction = "isotropic")
plot(k, sqrt(./pi)-r ~ r, ylab="L(r)-r", main="Clustering Indication Using L Function",legend=F)
dev.off()

######## parametric ETAS ########
library(ETAS)
etas_data <- check[c("StartDate","StartTime",'lon','lat',"ImpactedCustomers")]
names(etas_data) <- c('date', 'time', 'long', 'lat', 'mag')
etas_data$mag <- scaletointerval(etas_data$mag, 0, 6)
hist(etas_data$mag)
etas_data <- etas_data[etas_data$mag <= 2.5 & etas_data$mag > 0.1 ,]
hist(etas_data$mag)
power.cat <- catalog(etas_data, mag.threshold = 0.1)
plot(power.cat)
param01 <- c(0.5, 0.2, 0.05, 2.7, 1.2, 0.001, 2.3, 1)
power.fit <- etas(power.cat, param0 = param01, nthreads = 4)
plot(power.fit)
power.fit
rates(power.fit)
pr <- probs(power.fit)
summary(pr$prob)
plot(power.cat$longlat.coord[pr$target & (1 - pr$prob > 0.95), 1:2])
points(power.cat$longlat.coord[pr$target & (pr$prob > 0.95), 1:2],pch = 3, col = 2)
map("world", add = TRUE, col = "grey")
legend("bottomleft", c("background", "triggered"), pch = c(1, 3),col = 1:2)

power.res <- resid.etas(power.fit)
summary(power.res$tres)
ks.test(power.res$U, punif)

#####Psuedo likelihood model
x1 <- check$x[check$ImpactedCustomers >= 20]
y1 <- check$y[check$ImpactedCustomers >= 20]
data <- cbind(x1,y1) %>% unique()
x1 <- data[,1]
y1 <- data[,2]
d1 = as.matrix(dist(data)) ## matrix of distances between pts
n1 = length(x1)
f = function(p){  
  ## returns the negative pseudo log-likelihood
  ## p = (mu,alpha,beta,gamma,a1)
  if(p[1] < 0) return(99999)
  if(p[1] + p[2] < 0) return(99999)
  if(p[1] + p[3] < 0) return(99999)
  if(p[1] + p[2] + p[3] < 0) return(99999)
  if(p[4] < 0) return(99999)
  if(p[4] > 1) return(99999)
  if(p[5] < 0) return(99999)
  lam = p[1] + p[2] * x1 + p[3] * y1 
   for(i in 1:n1){
     for(j in c(1:n1)[-i]){
       lam[i] = lam[i] + p[4] * p[5] * exp(-p[5] * d1[i,j]) / (2*pi*d1[i,j])
     }
  }
  if (min(lam) < 0) return (99999)
  int2 = p[1] + p[2]/2 + p[3]/2 + p[4]*n1
  ## Note that this above is for a window of [0,1] x [0,1]
  cat("integral = ",int2," negative loglikelihood = ",
      int2-sum(log(lam)), "\n"," p = ",p,"\n") 
  ## integral should be roughly n when it's done
  return(int2-sum(log(lam)))
}

pstart = c(32, -20, -10, 1, 88)
fit1 = optim(pstart,f,control=list(maxit=200), hessian = T)
pend = fit1$par
b3 = sqrt(diag(solve(fit1$hess))) #standard error
#p =  31.72847 -20.5359 -10.71343 0.9711267 88.52974 
# se: 57.74315066 50.33926533 62.42760021  0.04402628  6.78584797

### Plot the Model's Background Rate
pend = c(31.72847, -20.5359, -10.71343, 0.9711267, 88.52974)
par(mfrow=c(1,2))
plot(c(0,1),c(0,1),type="n",xlab="x-coordinate",ylab="y-coordinate",
     main="background rate")
x2 = seq(0.05,0.95,length=10)
y2 = seq(0.05,0.95,length=10)
z2 = matrix(rep(0,(10*10)),ncol=10)
z3 = matrix(rep(0,(10*10)),ncol=10)
for(i in 1:10){
  for(j in 1:10){
    z2[i,j] = pend[1] + pend[2]*x2[i] + pend[3]*y2[j]
    z3[i,j] = pstart[1] + pstart[2]*x2[i] + pstart[3]*y2[j]
  }}
zmin = min(c(z2,z3))
zmax = max(c(z2,z3))
image(x2,y2,z2,col=gray((64:20)/64),zlim=c(zmin,zmax),add=T)
points(x1,y1)

### LEGEND
zrng = zmax - zmin
zmid = zmin + zrng/2
plot(c(0,10),c(zmid-2*zrng/3,zmid+2*zrng/3),type="n",axes=F,xlab="",ylab="")
zgrid = seq(zmin,zmax,length=100)
## zgrid = vector of 100 equally-spaced numbers spanning range of the values.
image(c(-1:1),zgrid,matrix(rep(zgrid,2),ncol=100,byrow=T),add=T,col=gray((64:20)/64))
text(2.5,zmin,as.character(signif(zmin,2)),cex=1)
text(2.5,zmax,as.character(signif(zmax,2)),cex=1)
text(2.5,zmid,as.character(signif(zmid,2)),cex=1)
text(4.5,zmid,"pts/unit area",srt=-90)
### PLOT LAMBDA_p on a 10 x 10 grid.
par(mfrow=c(1,2)) ## change this 3 to 2 for your projects.
plot(c(0,1),c(0,1),type="n",xlab="x-coordinate",ylab="y-coordinate",
     main="lambda_p")
x2 = seq(0.05,0.95,length=10)
y2 = seq(0.05,0.95,length=10)
zz2 = matrix(rep(0,(10*10)),ncol=10)
zz3 = matrix(rep(0,(10*10)),ncol=10) 
for(i in 1:10){
  for(j in 1:10){
    zz2[i,j] = pend[1] + pend[2] * x2[i] + pend[3] * y2[j]
    zz3[i,j] = pstart[1] + pstart[2] * x2[i] + pstart[3] * y2[j]
    for(k in c(1:n1)){
      zz2[i,j] = zz2[i,j] + pend[4] * pend[5] * exp(-pend[5] * 
                                                      sqrt((x2[i]-x1[k])^2+(y2[j]-y1[k])^2))
      zz3[i,j] = zz3[i,j] + pstart[4] * pstart[5] * exp(-pstart[5] * 
                                                          sqrt((x2[i]-x1[k])^2+(y2[j]-y1[k])^2))
    }
  }
}

zmin = min(c(zz2,zz3))
zmax = max(c(zz2,zz3))
image(x2,y2,zz2,col=gray((64:20)/64),zlim=c(zmin,zmax),add=T)
points(x1,y1,cex = 0.1)
##### Cox
# library(rgdal)
# library(raster)
# library(sf)
# library(maptools)
# library(lgcp)
# ca_df <- st_read("/Users/andrewchang/Desktop/222/tmp/CA_Counties/CA_Counties_TIGER2016.shp")
# W <- as.owin(ca_df)
# cox_data <- check[c('xp','yp', "StartDateTime")]
# names(cox_data) <- c('x','y','t')
# baseTime <- as.POSIXct("2020-01-01 00:00:00")
# cox_data$t<- difftime(cox_data$t,baseTime, units = "days") %>% as.numeric()
# tlim <- c(0, max(cox_data$t))
# xyt <- cox_data %>% as.matrix()
# xyt <- stppp(list(data = xyt, tlim = tlim, window = W))
# xyt <- integerise(xyt)
# xyt
# OW <- selectObsWindow(xyt, cellwidth = 2)
# den <- lambdaEst(xyt, axes = TRUE)
# plot(den)
# sar <- spatialAtRisk(den)
# mut1 <- muEst(xyt)
# mut <- constantInTime(xyt)
# gin <- ginhomAverage(xyt, spatial.intensity = sar, temporal.intensity = mut)
# sigmaphi1 <- spatialparsEst(gin, sigma.range = c(0, 10),phi.range = c(0, 10), spatial.covmodel = "exponential")
# kin <- KinhomAverage(xyt, spatial.intensity = sar, temporal.intensity = mut)
# attr(kin, "fname") <- "K[inhom]"
# sigmaphi2 <- spatialparsEst(kin, sigma.range = c(0, 10),phi.range = c(0, 10), spatial.covmodel = "exponential")
#baseTime <- as.POSIXct("1900-01-01 00:00:00")
#etas_data$time <- difftime(etas_data$StartDateTime,baseTime, units = "secs") %>% as.numeric()
#etas_data$time <- etas_data$time *1e-9
#etas_data$m <- scaletointerval(etas_data$m, 3, 6)
#etas_data <- data.frame(etas_data)
library(nphawkes)
hawke_data <- check[,c("StartDateTime",'x','y', "ImpactedCustomers")]
colnames(hawke_data) <- c('time','x','y','m')
hawke_data <- hawke_data[hawke_data$m >= 20,]
hawke_data$m <- hawke_data$m/sum(hawke_data$m)
baseTime <- as.POSIXct("2020-03-01 00:00:00")
hawke_data$time <- difftime(hawke_data$time,baseTime, units = "secs") %>% as.numeric()
hawke_data$time <- hawke_data$time *1e-9
hawke_data <- data.frame(hawke_data)
mydata = nphData(data = hawke_data, time_var = "time")
output = nphawkesTNS(data = mydata, nbins_t = 15, nbins_mu = 50, bw=100, eps = 1e-5)
plot(output,type="time")
plot(output, type = 'background', smooth = TRUE)

mydata = nphData(data = hawke_data, time_var = "time", x_var = "x", y_var = "y", mag = "m")
output =  nphawkesMSTNH(data = mydata)
plot(output, type = 'time')
plot(output, type = 'space')
plot(output, type = 'magnitude')
plot(output, type = 'background')



### Hawkes no magnitude
m3 <- function(x) format(round(x, 2), nsmall = 2)
theta0 = list(mu=.08,K=.75,alpha=12.5,beta=13.5,b=1) 
T <- 10^4
X1 <- 1
Y1 <- 1
z <- list('lon' = hawke_data$y, 'lat' = hawke_data$y,
          't' = hawke_data$time, 'n' = nrow(hawke_data))
loglk = function(theta,draw=0){
  # mu = theta[1]; K = theta[2]; alpha = theta[3]; beta = theta[4] 
  mu = theta0$mu
  K = theta[1]
  alpha = theta0$alpha
  beta = theta0$beta
  cat("\n mu = ",m3(mu),", K = ",m3(K),", alpha = ",m3(alpha),", beta = ",m3(beta),".\n") 
  if(min(mu,K,alpha,beta)<0.000000001) return(99999) 
  if(K>.99999) return(99999)
  if(draw){
    j1 <<- j1+1
    points(j1,theta[1],pch=3,col="green")
    #r = seq(0,3,length=100)
    #t = alpha/pi * exp(-alpha * r^2)
    #lines(r,t,col="orange",lty=2) 
  }
  sumlog = log(mu/X1/Y1) 
  intlam = mu*T + K*z$n
  const = K*alpha/pi*beta 
  for(j in 2:(z$n)){
    gij = 0
    for(i in 1:(j-1)){
      r2 = (z$lon[j]-z$lon[i])^2+(z$lat[j]-z$lat[i])^2
      gij = gij + exp(-beta*(z$t[j]-z$t[i])-alpha*r2)
    }
    lamj = mu / X1 / Y1 + const*gij
    if(lamj < 0){
      cat("lambda ",j," is less than 0.")
      return(99999)
    }
    sumlog = sumlog + log(lamj)
  }
  loglik = sumlog - intlam
  cat("loglike is ", loglik, ". sumlog = ", sumlog,". integral = ", intlam,".\n")
  # if(draw) lines(r,t,col="white",lty=2) 
  return(-1.0*loglik)
}

b1 = optim(theta0,loglk,method = "L-BFGS-B")
b2 = optim(b1$par,loglk,hessian=T,method ="L-BFGS-B")
theta2 = b2$par
sqrt(diag(solve(b2$hess))) ## for SEs 
