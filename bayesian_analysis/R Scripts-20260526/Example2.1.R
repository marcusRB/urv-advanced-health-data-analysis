
########### Example 2. Poisson Gamma; Goals by GironaFC 

# data 2022-23 season till April 11th
GoalsGirona = c(0,3,0,1,2,1,3,1,1,2,1,1,2,2,2,2,2,0,0,1,0,6,3,2,0,2,2,0)

y=GoalsGirona
n=length(y)
extr=as.integer(mean(y)+6*sqrt(mean(y)))


# Prior distribution. E(theta)=a/b, V(theta)=a/b^2.  theta is E(y)
overall.prior = c(ao=4,bo=3)
home.prior = c(ah=5,bh=3)
visiting.prior = c(av=3,bv=3)

a=overall.prior[1]
b=overall.prior[2]


############################## Analytic solution 

library(extraDistr)

par(mfrow=c(2,1))

#Plotting prior, std likelihood and posterior
plot(function(x)dgamma(x,a+sum(y),b+n), xlim=c(0,a/b+3*sqrt((a/b^2))),xlab=expression(theta),ylab="",main="Prior, Std likelihood and Posterior")
plot(function(x)dgamma(x,a,b),xlim=c(0,a/b+3*sqrt((a/b^2))),lty=3,add=TRUE)
likelihood = function(theta){theta^sum(y)*exp(-n*theta)}
K = integrate(likelihood,lower=0,upper=a/b+5*sqrt((a/b^2)))$value
curve(likelihood(x)/K,lty=2,add=TRUE)
legend("topright", c("prior","likelihood","posterior"),lty=c(3,2,1))


#Plotting prior and posterior predictive
plot(0:extr,dgpois(0:extr,shape=a,rate=b),type="h",xlim=c(0,extr),ylim=c(0,0.4),xlab="y",ylab="",main="Prior and posterior predictive")
points(0:extr+0.1,dgpois(0:extr,shape=a+sum(y),rate=b+n),col="red",type="h")
legend("topright", c("prior predictive","posterior predictive"),lty=1,col=c("black","red"))


# Computing point and interval estimates and probabilities of hypotheses

priormean=a/b
priormean
priorvar=a/(b*b)
priorvar
postmean=(a+sum(y))/(b+n)
postmean
postvar=(a+sum(y))/((b+n)*(b+n))
postvar
priormedian=qgamma(.5,a,b)
priormedian
postmedian=qgamma(.5,a+sum(y),b+n)
postmedian
priorint90=qgamma(c(.05,.95),a,b)
priorint90
postint90=qgamma(c(.05,.95),a+sum(y),b+n)
postint90
postprobs1=1-pgamma(1,a+sum(y),b+n)
postprobs1
postprob115=pgamma(1.5,a+sum(y),b+n)-pgamma(1,a+sum(y),b+n)
postprob115

# Computing predictive mean and variance for 

priorpredmean=a/b
priorpredmean
priorpredvar=a*(b+1)/(b*b)
priorpredvar
postpredmean=(a+sum(y))/(b+n)
postpredmean
postpredvar=(a+sum(y))*(b+n+1)/((b+n)*(b+n))
postpredvar

# Ploting possible priors

par(mfrow=c(3,1))

a=3
b=3
plot(function(x)dgamma(x,a,b),xlim=c(0,4),ylim=c(0,1.25),xlab=expression(theta),ylab="",main="Priors with E(th)=1")
plot(function(x)dgamma(x,6,6),xlim=c(0,4),lty=2,add=TRUE)
plot(function(x)dgamma(x,9,9),xlim=c(0,4),lty=3,add=TRUE)
legend("topright", c("Gamma(9,9)","Gamma(6,6)","Gamma(3,3)"),lty=c(3,2,1))

a=4
b=3
plot(function(x)dgamma(x,a,b),xlim=c(0,4),ylim=c(0,1.25),xlab=expression(theta),ylab="",main="Priors with E(th)=1.33")
plot(function(x)dgamma(x,8,6),xlim=c(0,4),lty=2,add=TRUE)
plot(function(x)dgamma(x,12,9),xlim=c(0,4),lty=3,add=TRUE)
legend("topright", c("Gamma(12,9)","Gamma(8,6)","Gamma(4,3)"),lty=c(3,2,1))

a=5
b=3
plot(function(x)dgamma(x,a,b),xlim=c(0,4),ylim=c(0,1.25),xlab=expression(theta),ylab="",main="Priors with E(th)=1.66")
plot(function(x)dgamma(x,10,6),xlim=c(0,4),lty=2,add=TRUE)
plot(function(x)dgamma(x,15,9),xlim=c(0,4),lty=3,add=TRUE)
legend("topright", c("Gamma(15,9)","Gamma(10,6)","Gamma(5,3)"),lty=c(3,2,1))







######################### Solution through MCMC simulation (JAGS)

# YOU MUST INSTALL JAGS FROM: https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/
# You must install the library R2WinBUGS, install.packages('R2jags', repos = "http://cran.fhcrc.org/")

library(R2jags)

model.pois <- "
 model
 { 
  for (i in 1:n) {
   y[i] ~ dpois(theta)
   }
   theta ~ dgamma(a,b)
 }
"

data <- list(y=y, n=n, a=a, b=b) 

parameters <- c("theta") 

examplePoiss <-jags(data, parameters.to.save=parameters,
                  model=textConnection(model.pois),n.iter=55000,n.burnin=5000,n.thin=2,n.chains=3)

print(examplePoiss)

traceplot(examplePoiss)

attach.jags(examplePoiss)

par(mfrow=c(2,1))

hist(theta,xlab=expression(theta),xlim=c(0,a/b+3*sqrt((a/b^2))),freq=FALSE, main="Posterior theta")

# plot stand likelihood, prior and posterior for theta
plot(density(theta, adjust = 1.5), ylab="", xlim=c(0,a/b+3*sqrt((a/b^2))), xlab=expression(theta), main="Prior, Std likelihood and posterior")
plot(function(x)dgamma(x,a,b),xlim=c(0,a/b+3*sqrt((a/b^2))),lty=3,add=TRUE)
likelihood = function(theta){theta^sum(y)*exp(-n*theta)}
K = integrate(likelihood,lower=0,upper=a/b+5*sqrt((a/b^2)))$value
curve(likelihood(x)/K,lty=2,add=TRUE)
legend("topright", c("prior","std likelihood","posterior"),lty=c(3,2,1))

# plot prior and posterior predictive for the number of goals scored by Girona
par(mfrow=c(2,1))
predict.post = rpois(length(theta),lambda=theta)
tbpredpost = table(predict.post)
priortheta = rgamma(length(theta),a,b)
predict.prior = rpois(length(theta),lambda=priortheta)
tbpredprior = table(predict.prior)
plot(tbpredprior/length(theta), ylab="Prob", xlab="y", main="Prior predictive")
plot(tbpredpost/length(theta), ylab="Prob", xlab="y", main="Posterior predictive")


### Estimation of posterior expected value
mean(theta)
### Estimation of posterior median
median(theta)
### Estimation of 90% posterior credible interval
quantile(theta,probs=c(.05,.95))
### Estimation of posterior probability that 1<theta
sum(1<theta)/length(theta)
### Estimation of posterior probability that 1<theta<1.5
sum((1<theta)&(theta<1.5))/length(theta)


### Posterior predictive expected value
mean(predict.post)
### Posterior predictive median
median(predict.post)
### 90% posterior predictive credible interval
quantile(predict.post,probs=c(.05,.95))
### Posterior probability of scoring more than 3 goals
sum(3<predict.post)/length(predict.post)


## Plot posterior predictive for the total number of goals scored in the 38 games of one season
## and point and interval estimation and testing of the total number of goals

nobs=as.integer(length(predict.post)/38)
predict.posttot = (1:nobs)
for (i in 1:nobs){
  predict.posttot[i]=sum(predict.post[(i-1)*38+(1:38)])
}
predict.posttot
tbpredposttot = table(predict.posttot)
plot(tbpredposttot/nobs, ylab="Prob", xlab="y", main="Posterior predictive total number goals per season")

### Posterior predictive expected value total number goals per season
mean(predict.posttot)
### Posterior predictive median total number goals per season
median(predict.posttot)
### 90% posterior predictive credible interval total number goals per season
quantile(predict.posttot,probs=c(.05,.95))
### Posterior probability of scoring more than 65 goals in one season
sum(65<predict.posttot)/length(predict.posttot)



###########  Inference about probability of scoring 0 goals in a game exp(-theta)

hist(exp(-theta),xlab=expression(exp(-theta)),xlim=c(0,1),freq=FALSE, main="Posterior Pr(y=0)")
plot(density(exp(-theta),adjust=1.5),ylab="",xlab=expression(exp(-theta)),xlim=c(0,1),main="Posterior Pr(y=0)")


### Estimation of posterior expected value
mean(exp(-theta))
### Estimation of posterior median
median(exp(-theta))
### Estimation of 90% posterior credible interval
quantile(exp(-theta),probs=c(.05,.95))
### Estimation of posterior probability that psi>.3
sum(0.3<exp(-theta))/length(theta)


###########  Inference about probability of scoring 4 or more goals in a game

hist(1-exp(-theta)*(1+theta+theta^2/2+theta^3/6),xlab=expression(1-exp(-theta)*(1+theta+theta^2/2+theta^3/6)),xlim=c(0,1),freq=FALSE, main="Posterior of Pr(y>=4)")
plot(density(1-exp(-theta)*(1+theta+theta^2/2+theta^3/6),adjust=1.5),ylab="",xlab=expression(1-exp(-theta)*(1+theta+theta^2/2+theta^3/6)),xlim=c(0,1),main="Posterior of Pr(y>=4)")

### Estimation of posterior expected value
mean(1-exp(-theta)*(1+theta+theta^2/2+theta^3/6))
### Estimation of posterior median
median(1-exp(-theta)*(1+theta+theta^2/2+theta^3/6))
### Estimation of 90% posterior credible interval
quantile(1-exp(-theta)*(1+theta+theta^2/2+theta^3/6),probs=c(.05,.95))
### Estimation of posterior probability that psi>.1
sum(0.1<1-exp(-theta)*(1+theta+theta^2/2+theta^3/6))/length(theta)


detach.jags()
