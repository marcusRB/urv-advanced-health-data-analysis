
####### Example 1. Binomial Beta; Coin, pin (point up/down) and blue eyes

# data
n=30
y=9

# Prior distribution. E(theta)=a/(a+b), V(theta)=ab/(a+b+1)(a+b)^2
coin.prior = c(ac=100,bc=100)
pin.prior = c(ap=2,bp=2)
bluee.prior = c(ab=2,bb=18)

a=pin.prior[1]
b=pin.prior[2]


############################## Analytic solution 


library(VGAM)
par(mfrow=c(2,1))

#Plotting standardized likelihood, prior and posterior all together for pins
plot(function(x)dbeta(x,a+y,b+n-y), xlab=expression(theta),ylab="",main="Std likelihood, prior and posterior")
plot(function(x)dbeta(x,a,b),lty=3,add=TRUE)
likelihood = function(theta){theta^y*(1-theta)^(n-y)}
K = integrate(likelihood,lower=0,upper=1)$value
curve(likelihood(x)/K,lty=2,add=TRUE)
legend("topright", c("prior","stand likelihood","posterior"),lty=c(3,2,1))


#Plotting prior and posterior predictive
plot(0:n,dbetabinom.ab(0:n,n,a,b),type="h",xlim=c(0,n),ylim=c(0,0.20),xlab="y",ylab="Pr",main="Prior and posterior predictive")
points(0:n+0.1,dbetabinom.ab(0:n,n,a+y,b+n-y),col="red",type="h")
legend("topright", c("prior predictive","posterior predictive"),lty=1,col=c("black","red"))

###############

#Plotting prior for coin, pin and blue eye, and then posterior for them
par(mfrow=c(3,1))
plot(function(x)dbeta(x,coin.prior[1],coin.prior[2]), xlab=expression(theta),ylab="",main="Prior for coin, pin, blue eye")
plot(function(x)dbeta(x,pin.prior[1],pin.prior[2]),lty=2,add=TRUE)
plot(function(x)dbeta(x,bluee.prior[1],bluee.prior[2]),lty=3,add=TRUE)
legend("topright", c("Coin","Pin","Blue Eyes"),lty=c(1,2,3))

#Plotting standardized likelihood
K = integrate(likelihood,lower=0,upper=1)$value
curve(likelihood(x)/K,lty=1,xlab=expression(theta),ylab="",main="Standardized likelihood")

#Plotting posterior for coin, pin and blue eye, and then posterior for them
plot(function(x)dbeta(x,coin.prior[1]+y,coin.prior[2]+n-y), xlab=expression(theta),ylab="",main="Posterior for coin, pin, blue eye")
plot(function(x)dbeta(x,pin.prior[1]+y,pin.prior[2]+n-y),lty=2,add=TRUE)
plot(function(x)dbeta(x,bluee.prior[1]+y,bluee.prior[2]+n-y),lty=3,add=TRUE)
legend("topright", c("Coin","Pin","Blue Eyes"),lty=c(1,2,3))

########################3

# Computing point and interval estimates of theta and probabilities of hypotheses about theta

priormean=a/(a+b)
priormean
priorvar=a*b/((a+b)*(a+b)*(a+b+1))
priorvar
postmean=(a+y)/(a+b+n)
postmean
postvar=((a+y)*(b+n-y))/((a+b+n)*(a+b+n)*(a+b+n+1))
postvar
priormedian=qbeta(.5,a,b)
priormedian
postmedian=qbeta(.5,a+y,b+n-y)
postmedian
priorint90=qbeta(c(.05,.95),a,b)
priorint90
postint90=qbeta(c(.05,.95),a+y,b+n-y)
postint90
postprobsm4=pbeta(.4,a+y,b+n-y)
postprobsm4
postprob57=pbeta(.7,a+y,b+n-y)-pbeta(.5,a+y,b+n-y)
postprob57

# Computing predictive mean and variance for nn new Bernouilli observations

nn=30
priorpredmean=nn*a/(a+b)
priorpredmean
priorpredvar=nn*a*b*(a+b+nn)/((a+b)*(a+b)*(a+b+1))
priorpredvar
postpredmean=nn*(a+y)/(a+b+n)
postpredmean
postpredvar=nn*((a+y)*(b+n-y)*(a+b+n+nn))/((a+b+n)*(a+b+n)*(a+b+n+1))
postpredvar

# Plot possible priors

par(mfrow=c(2,1))

plot(function(x)dbeta(x,0.5,0.5),xlim=c(0,1),lty=4,ylim=c(0,6),xlab=expression(theta),ylab="",main="Priors with E(th)=0.5")
plot(function(x)dbeta(x,1,1),xlim=c(0,1),lty=3,add=TRUE)
plot(function(x)dbeta(x,5,5),xlim=c(0,1),lty=2,add=TRUE)
plot(function(x)dbeta(x,30,30),xlim=c(0,1),lty=1,add=TRUE)
legend("topright", c("Beta(30,30)","Beta(5,5)","Beta(1,1)","Beta(.5,.5)"),lty=c(1,2,3,4))

plot(function(x)dbeta(x,2,18),xlim=c(0,1),lty=1,ylim=c(0,8),xlab=expression(theta),ylab="",main="Priors with E(th)=.1,.2,.3,.5")
plot(function(x)dbeta(x,4,16),xlim=c(0,1),lty=2,add=TRUE)
plot(function(x)dbeta(x,6,14),xlim=c(0,1),lty=3,add=TRUE)
plot(function(x)dbeta(x,10,10),xlim=c(0,1),lty=4,add=TRUE)
legend("topright", c("Beta(2,18)","Beta(4,16)","Beta(6,14)","Beta(10,10)"),lty=c(1,2,3,4))






######################### Solution through MCMC simulation (JAGS)

# YOU MUST INSTALL JAGS FROM: https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/
# You must install the library R2WinBUGS, install.packages('R2jags', repos = "http://cran.fhcrc.org/")

library(R2jags)

model.bin <- "
 model
 { 
  y ~ dbin(theta, n)  
  theta ~ dbeta(a,b)
 }
"

# data and prior must be stored in a list object
data <- list(n=n, y=y, a=a, b=b) 

# parameters to be monitored
parameters <- c("theta") 

exampleBinomial <-jags(data, parameters.to.save=parameters,
                  model=textConnection(model.bin),n.iter=55000,n.burnin=5000,n.thin=2,n.chains=3)

print(exampleBinomial)

traceplot(exampleBinomial)

attach.jags(exampleBinomial)

par(mfrow=c(2,1))

hist(theta,xlim=c(0,1),xlab=expression(theta),freq=FALSE, main="Histogram posterior theta")

# plot likelihood function, prior and posterior for theta
plot(density(theta, adjust = 1.5), ylab="", xlab=expression(theta), main="Std likelihood, prior and posterior", xlim=c(0,1))
plot(function(x)dbeta(x,a,b),lty=3 ,add=TRUE)
likelihood = function(theta){theta^y*(1-theta)^(n-y)}
K = integrate(likelihood,lower=0,upper=1)$value
curve(likelihood(x)/K,lty=2,add=TRUE)
legend("topright", c("prior","stand likelihood","posterior"),lty=c(3,2,1))


# plot prior and posterior predictive for throwing npr whatever
npr=30
par(mfrow=c(2,1))

predict.post = rbinom(length(theta),size=npr,prob=theta)
tbpredpost = table(predict.post)
#
priortheta = rbeta(length(theta),a,b)
predict.prior = rbinom(length(theta),size=npr,prob=priortheta)
tbpredprior = table(predict.prior)
plot(tbpredprior/length(theta),ylab="Prob", xlab="y", xlim=c(0,npr), main="Prior predictive")
plot(tbpredpost/length(theta),ylab="Prob", xlab="y", xlim=c(0,npr), main="Posterior predictive")


### Estimation of posterior expected value
mean(theta)
### Estimation of posterior median
median(theta)
### Estimation of 90% posterior credible interval
quantile(theta,probs=c(.05,.95))
### Estimation of posterior probability that theta<0.4
sum(0.4>theta)/length(theta)
### Estimation of posterior probability that 0.5<theta<0.7
sum((0.5<theta)&(theta<0.7))/length(theta)


### Posterior predictive expected value
mean(predict.post)
### Posterior predictive median
median(predict.post)
### 90% posterior predictive credible interval
quantile(predict.post,probs=c(.05,.95))
### Posterior probability of observing more than 15 out of 30 point up
sum(15<predict.post)/length(predict.post)




################## Inference about the odds and the logodds=log(theta/(1-theta))

plot(density(theta/(1-theta), adjust = 1.5), ylab="", xlab=expression(theta/(1-theta)), main="Posterior for odds", xlim=c(0,3))
plot(density(log(theta/(1-theta)), adjust = 1.5), ylab="", xlab=expression(log(theta/(1-theta))), main="Posterior for logodds", xlim=c(-3,3))

#hist(log(theta/(1-theta)),xlab=expression(log(theta/(1-theta))),main="Histogram simulations posterior logodds")
### Estimation of posterior expected value of logodds
mean(log(theta/(1-theta)))
### Estimation of posterior median of logodds
median(log(theta/(1-theta)))
### Estimation of 90% posterior credible interval of logodds
quantile(log(theta/(1-theta)),probs=c(.05,.95))
### Estimation of posterior probability that 0<logodds
sum(0<log(theta/(1-theta)))/length(theta)
### Estimation of posterior probability that 0.6<logodds<0.7
sum((0.6<log(theta/(1-theta)))&(log(theta/(1-theta))<0.7))/length(theta)

detach.jags()


