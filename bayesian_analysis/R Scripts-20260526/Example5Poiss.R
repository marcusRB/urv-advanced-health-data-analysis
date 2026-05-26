
########### Example 5. Comparison of two Poisson Gamma; Goals by GironaFC 

# data, season 2022-23 on April 11th
GoalsGironaHome = c(3,0,2,3,1,1,2,2,2,0,1,6,0,2)
GoalsGironaAwayHome = c(0,1,1,1,2,1,2,2,0,0,3,2,2,0)
mean(GoalsGironaHome)
mean(GoalsGironaAwayHome)

y1=GoalsGironaHome
n1=length(y1)
y2=GoalsGironaAwayHome
n2=length(y2)


# Prior distribution. E(theta)=a/b, V(theta)=a/b^2. theta is E(y)
Home.prior = c(ah=5,bh=3)
AwayHome.prior = c(av=3,bv=3)

a1=Home.prior[1]
b1=Home.prior[2]
a2=AwayHome.prior[1]
b2=AwayHome.prior[2]

# When using the same prior for the two populations
#a1=4
#a2=4
#b1=3
#b2=3



######################### Solution through MCMC simulation (JAGS)

# YOU MUST INSTALL JAGS FROM: https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/
# You must install the library R2WinBUGS, install.packages('R2jags', repos = "http://cran.fhcrc.org/")

library(R2jags)

model.comp2pois <- "
 model
 { 
  for (i in 1:n1) {
   y1[i] ~ dpois(th1)
   }
   th1 ~ dgamma(a1,b1)
  for (j in 1:n2) {
   y2[j] ~ dpois(th2)
   }
   th2 ~ dgamma(a2,b2)
 }
"

data <- list(y1=y1, n1=n1, a1=a1, b1=b1, y2=y2, n2=n2, a2=a2, b2=b2) 

parameters <- c("th1","th2") 

examplec2Poiss <-jags(data, parameters.to.save=parameters,
                  model=textConnection(model.comp2pois),n.iter=55000,n.burnin=5000,n.thin=2,n.chains=3)

print(examplec2Poiss)

#traceplot(examplec2Poiss)

attach.jags(examplec2Poiss)

######## Inference th1-th2

par(mfrow=c(3,1))

hist(th1-th2,xlim=c(-2,3),freq=FALSE, main="Posterior th1-th2")
abline(v=c(0),lty=1)
# plot prior and posterior for th1-th2
prth1=rgamma(length(th1),a1,b1)
prth2=rgamma(length(th1),a2,b2)
plot(density(th1-th2, adjust = 1.5), ylab="",xlim=c(-2,3),xlab=expression(th1-th2), main="Posterior th1-th2")
abline(v=c(0),lty=1)
#legend("topright", c("prior","posterior"),lty=c(3,1))
plot(density(prth1-prth2,adjust=1.5),ylab="",xlim=c(-2,3),lty=3,xlab=expression(th1-th2), main="Prior th1-th2")
abline(v=c(0),lty=1)

### Estimation of posterior expected value
mean(th1-th2)
### Estimation of posterior median
median(th1-th2)
### Estimation of 90% posterior credible interval
quantile(th1-th2,probs=c(.05,.95))
### Estimation of posterior probability that 
sum(0<th1-th2)/length(th1)



###########  Inference about th1/th2

hist(th1/th2, xlim=c(0,4),freq=FALSE, xlab=expression(th1/th2), main="Posterior th1/th2")
plot(density(prth1/prth2,adjust=1.5),ylab="",xlim=c(0,4),lty=3,xlab=expression(th1/th2), main="Prior th1/th2")
plot(density(th1/th2,adjust=1.5),ylab="",xlab=expression(th1/th2),xlim=c(0,4),main="Posterior th1/th2")

### Estimation of posterior expected value
mean(th1/th2)
### Estimation of posterior median
median(th1/th2)
### Estimation of 90% posterior credible interval
quantile(th1/th2,probs=c(.05,.95))
### Estimation of posterior probability 
sum(1<th1/th2)/length(th1)

detach.jags()
