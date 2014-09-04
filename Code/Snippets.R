
# Problem 1.43a
nObs <- 100000
x <- rexp(n = nObs, rate = 1)
y <- rexp(n = nObs, rate = 1)
u <- (x-y)/(x+y)
min(u)
max(u)
hist(u, freq = FALSE, breaks = 50, ylim = c(0,1), main = 'Distr. of (X-Y)/(X+Y) if X and Y are indep Exp(1)')
gridVals <- seq(-1,1,length=1000)
lines(gridVals, 1/2*rep(1,length(gridVals)), col = "red")
legend(x = 0.25, y = 0.8, legend = c("Simulated", "True"),  col = c(1, 2), lty =c(1,1))


# Problem 33
nSim <- 100000
a <- runif(nSim) # Simulate A ~ U(0,1)
X <- runif(n = nSim, min = -a, max = a) # X | (A = a) ~ U(-a,a)
hist(X,freq = FALSE)
xGrid <- seq(-1,1,length=1000)
lines(xGrid,-0.5*log(abs(xGrid)), col = "red", lwd = 2)
legend(x = -1, y = 1.2, lty = c(1,1), col = c("black","red"), lwd = c(1,2), legend = c('Simulated', 'True pdf'))
mean(X)
var(X)

#install.packages("VGAM")
library(VGAM)
# Extremes of sample from laplace distribution
n <- 10
nSim <- 10000
b <- 2
xVal <- 1 # Computes Pr(maximum < xVal)
extremes <- matrix(0,nSim,2)
for (i in 1:nSim){
  x<-rlaplace(n, loc = 0, scale = b)
  extremes[i,1] <- min(x)
  extremes[i,2] <- max(x)
}
sum(extremes[,2]<xVal)/nSim                                       # The simulated answer
if (xVal<0) (1/(2^n))*exp(n*xVal/b) else (1-(1/2)*exp(-xVal/b))^n # The mathematical answer (exact)

hist(extremes[,1], 50)
hist(extremes[,2], 50)


# Chapter 5 Multivariate normal
library(manipulate)
install.packages("mvtnorm")
library(mvtnorm)

############ Play with the normal distribution using manipulate ############
PlotBivariateNormalDens <- function(mu = c(0,0), sigma = c(1,1), rho = 0){
  mu <- c(1,2)
  Lambda <- matrix(c(sigma[1]^2,sigma[1]*sigma[2]*rho,sigma[1]*sigma[2]*rho,sigma[2]^2),2,2)
  xGrid <- seq(mu[1]-3*sqrt(Lambda[1,1]),mu[1]+3*sqrt(Lambda[1,1]),length = 100)
  yGrid <- seq(mu[2]-3*sqrt(Lambda[2,2]),mu[2]+3*sqrt(Lambda[2,2]),length = 100)
  dens <- matrix(0,length(xGrid),length(yGrid))
  ix <- 0
  for (x in xGrid){
    ix <- ix + 1
    iy <- 0
    for (y in yGrid){
      iy <- iy + 1
      dens[ix,iy] <- dmvnorm(c(x,y), mu, Lambda)
    }
  }
  contour(xGrid, yGrid, dens, xlim = c(min(xGrid), max(xGrid)), ylim = c(min(yGrid), max(yGrid)), xlab ="x1", ylab = "x2", main = "Multivariate Normal Contours")
}

manipulate(
  PlotBivariateNormalDens(c(m1,mu2), c(sigma1, sigma2), rho),
  mu1 = slider(-10, 10, step=0.5, initial = 0, label = "Mean of X1"),
  mu2 = slider(-10, 10, step=0.5, initial = 0, label = "Mean of X2"),
  sigma1 = slider(0.01, 10, step=0.25, initial = 1, label = "Standard deviation of X1"),
  sigma2 = slider(0.01, 10, step=0.25, initial = 1, label = "Standard deviation of X2"),
  rho = slider(-1, 1, step=0.1, initial = 0, label = "The correlation coefficient, rho")
)

########################################################################

mu <- c(1,2)
Lambda <- matrix(c(1,-0.5,-0.5,2),2,2)
xGrid <- seq(mu[1]-3*sqrt(Lambda[1,1]),mu[1]+3*sqrt(Lambda[1,1]),length = 100)
yGrid <- seq(mu[2]-3*sqrt(Lambda[2,2]),mu[2]+3*sqrt(Lambda[2,2]),length = 100)
dens <- matrix(0,length(xGrid),length(yGrid))
ix <- 0
for (x in xGrid){
  ix <- ix + 1
  iy <- 0
  for (y in yGrid){
    iy <- iy + 1
    dens[ix,iy] <- dmvnorm(c(x,y), mu, Lambda)
  }
}
persp(xGrid, yGrid, dens, xlab ="x", ylab = "y", zlab = "Density", phi = 45, theta = 45)
contour(xGrid, yGrid, dens, xlab ="x", ylab = "y", zlab = "Density", main = "Marginals are normal, joint is normal")


# Normal density - carved out
dens <- matrix(0,length(xGrid),length(yGrid))
ix <- 0
for (x in xGrid){
  ix <- ix + 1
  iy <- 0
  for (y in yGrid){
    iy <- iy + 1
    if ((x>mu[1] & y>mu[2]) | (x<mu[1] & y<mu[2])) dens[ix,iy] <- dmvnorm(c(x,y), mu, Lambda)
  }
}
persp(xGrid, yGrid, dens, xlab ="x", ylab = "y", zlab = "Density", phi = 80, theta = 70)
contour(xGrid, yGrid, dens, xlab ="x", ylab = "y", zlab = "Density", main = "Marginals are normal, joint is not normal")

# Plotting the principal components axes
mu <- c(1,2)
Lambda <- matrix(c(1,0.9*sqrt(2),0.9*sqrt(2),2),2,2)
xGrid <- seq(mu[1]-3*sqrt(Lambda[1,1]),mu[1]+3*sqrt(Lambda[1,1]),length = 100)
yGrid <- seq(mu[2]-3*sqrt(Lambda[2,2]),mu[2]+3*sqrt(Lambda[2,2]),length = 100)
dens <- matrix(0,length(xGrid),length(yGrid))
ix <- 0
for (x in xGrid){
  ix <- ix + 1
  iy <- 0
  for (y in yGrid){
    iy <- iy + 1
    dens[ix,iy] <- dmvnorm(c(x,y), mu, Lambda)
  }
}
persp(xGrid, yGrid, dens, xlab ="x", ylab = "y", zlab = "Density", phi = 45, theta = 45)
contour(xGrid, yGrid, dens, xlab ="x", ylab = "y", zlab = "Density")
lines(seq(-2,4,length=100),seq(-2,6,length=100))
lines(seq(0,3,length=100),seq(4*0.9,-1*0.9,length=100))

# Convergence - Binomial(n, lambda/n) --> Po(lambda) in distribution

PlotBinoAndPoisson <- function(n,lambda){
  sdBinom <- sqrt(lambda*(1-lambda/n))
  xVals <- seq(0,ceiling(sdBinom*6))
  par(mfrow = c(1,2))
  barplot(names.arg = xVals, height = dbinom(x = xVals,size = n, prob = lambda/n), beside = TRUE, main = "Bin(n,lambda/n)")
  barplot(names.arg = xVals, height = dpois(xVals, lambda), beside = TRUE, main = "Pois(lambda)")
  
}
library(manipulate)
manipulate(
  PlotBinoAndPoisson(n, lambda),
  n = slider(1, 100, step=1, initial = 5, label = "n"),
  lambda = slider(0, 10, step=0.1, initial = 1, label = "Lambda")
)

# Convergence - two dice example
nOmega <- 100 # number of different omega (w) that we try, in theory we should try ALL w.
n <- 10000 # number of observations (n) for a given omega
Xn <- matrix(NA, n, nOmega)
for (omega in 1:nOmega){
  y <- matrix(NA,n,1)
  for (i in 1:n){
    y[i] <- sum(sample(1:6, 2, replace = TRUE)) # Draws of Y = sum of two dice
  }
  Xn[,omega] <- cumsum(y)/(1:n) # Computes Xn, mean of y up to time n
}
plot(Xn[,1], type = 'l', ylim = c(1,12), xlab = 'n', ylab = expression(X[n]), main = expression(omega[1]))
plot(Xn[,2], type = 'l', ylim = c(1,12), xlab = 'n', ylab = expression(X[n]), main = expression(omega[2]))

epsi <- 0.1
ProbBiggerThanEps <- matrix(NA,n,1)
for (i in 1:n){
  ProbBiggerThanEps[i] <- sum(abs((Xn[i,]-7))>epsi)/nOmega
}
plot(ProbBiggerThanEps, type = "l")

# Central limit theorem
IllustrateCLTforPoisson <- function(n = 10, mu = 1, nRep = 10000){
  means <- matrix(NA, nRep, 1)
  for (i in 1:nRep){
    x <- rpois(n,mu)
    means[i] <- (mean(x)-mu)/(sqrt(mu/n)) # Note mean=var=mu in the pois(mu) distribution.
  }
  hist(means, 20)
}
library(manipulate)
manipulate(
  IllustrateCLTforPoisson(n, mu, nRep = 50000),
  n = slider(2, 200, step=1, initial = 2, label = "n"),
  mu = slider(0.01, 5, step=0.1, initial = 1, label = "poisson mean")
)

# Empirical distribution function
PlotEmpiricalGamma <- function(n, a = 1, b = 1){
  x <- sort(rgamma(n,a,b))
  xVals <- seq(0,3, length = 100) #FIXME
  Dens <- matrix(NA,length(xVals),1)
  for (i in 1:length(xVals)){
    Dens[i] <- sum(x<=xVals[i])/n
  }
  plot(xVals, Dens, main = "Whatever", type="s")
  lines(xVals,pgamma(xVals,a,b), col = "red")
}
PlotEmpiricalGamma(n=100, 2, 2)

# Problem 6.6 X1, X2, are iid Pa(1,2). Let Yn = max(X1,X2,...,Xn). Show that n(Yn-1) converges in distribution and find the limit distribution
install.packages("VGAM") # Contains rpareto, dpareto etc
library(VGAM)
ParetoMinLimit <- function(n, scalingFact, nRep = 1000){
  # We are looking at the distribution of a(Yn-1). Only a = n will converge in distribution. You need the right scaling.
  mins <- matrix(NA,nRep,1)
  if (scalingFact == "sqrt") a <- sqrt(n)
  if (scalingFact == "n") a <- n
  if (scalingFact == "n^2") a <- n^2
  
  for (i in 1:nRep){
    mins[i] <- a*(min(rpareto(n, 1, 2))-1)
  }
 
  hist(mins,30, freq = FALSE, main = "Distribution of a*(Yn-1)")
  if (a == n){
    xVals <- seq(0,max(mins), length = 100)
    lines(xVals, dexp(xVals, 2), col = "red", lwd= 2) # Note that R's definition of Exp(lambda) is not the same as in Gut's book. In R if x = Exp(a), then E(X)=1/a
  }
}

library(manipulate)
manipulate(
  ParetoMinLimit(n, a),
  n = slider(10, 10000, step=1, initial = 100, label = "n"),
  a = picker("sqrt", "n", "n^2", label = "Scaling factor (a)")
)