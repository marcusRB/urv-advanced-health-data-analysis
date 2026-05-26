
########## Example 3. Exponential Gamma; Life length slugs in years

# data. True theta is .5
BlackSlug = c(1.71,0.71,2.26,2.04,2.93,15.87,3.80,2.47,0.41,0.47)
y=BlackSlug
n=length(y)


# Prior distribution, E(theta)=a/b, V(theta)=a/b^2. E(Y)=1/theta
a=3
b=9


############################## Analytic solution 

par(mfrow=c(2,1))

#Plotting prior likelihood and posterior
plot(function(x)dgamma(x,a+n,b+sum(y)), xlim=c(0,a/b+3*sqrt((a/b^2))),xlab=expression(theta),ylab="",main="Prior, Std likelihood and Posterior")
plot(function(x)dgamma(x,a,b),lty=3,add=TRUE)
likelihood = function(theta){theta^n*exp(-theta*sum(y))}
K = integrate(likelihood,lower=0,upper=a/b+5*sqrt((a/b^2)))$value
curve(likelihood(x)/K,lty=2,add=TRUE)
legend("topright", c("prior","likelihood","posterior"),lty=c(3,2,1))


dgammaexpon = function(x,a,b){
  const = b^a/gamma(a)
  rslt = gamma(a+1)/((b+x)^(a+1))
  rslt = const*rslt
}


plot(function(x)dgammaexpon(x,a+n,b+sum(y)), xlim=c(0,mean(y)+3*sd(y)),col="red",xlab="y",ylab="",main="Prior and posterior predictive")
plot(function(x)dgammaexpon(x,a,b),xlim=c(0,mean(y)+3*sd(y)),lty=3,add=TRUE)
legend("topright", c("prior predictive","posterior predictive"),lty=c(3,1),col=c("black","red"))

# Computing point and interval estimates and probabilities of hypotheses

priormean=a/b
priormean
priorvar=a/(b*b)
priorvar
postmean=(a+n)/(b+sum(y))
postmean
postvar=(a+n)/((b+sum(y))*(b+sum(y)))
postvar
priormedian=qgamma(.5,a,b)
priormedian
postmedian=qgamma(.5,a+n,b+sum(y))
postmedian
priorint90=qgamma(c(.05,.95),a,b)
priorint90
postint90=qgamma(c(.05,.95),a+n,b+sum(y))
postint90
postprobs1=1-pgamma(1,a+n,b+sum(y))
postprobs1
postprob115=pgamma(1.5,a+n,b+sum(y))-pgamma(1,a+n,b+sum(y))
postprob115

# Computing predictive mean and variance for 

priorpredmean=b/(a-1)
priorpredmean
priorpredvar=a*b*b/((a-1)*(a-1)*(a-2))
priorpredvar
postpredmean=(b+sum(y))/(a+n-1)
postpredmean
postpredvar=(a+n)*(b+sum(y))*(b+sum(y))/((a+n-1)*(a+n-1)*(a+n-2))
postpredvar







######################### Solution through MCMC simulation (JAGS)

# YOU MUST INSTALL JAGS FROM: https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/
# You must install the library R2WinBUGS, install.packages('R2jags', repos = "http://cran.fhcrc.org/")

library(R2jags)

model.exp <- "
 model
 { 
  for (i in 1:n) {
   y[i] ~ dexp(theta)
   }
   theta ~ dgamma(a,b)
 }
"

data <- list(y=y, n=n, a=a, b=b) 

parameters <- c("theta") 

exampleExpon <-jags(data, parameters.to.save=parameters,
                   model=textConnection(model.exp),n.iter=550000,n.burnin=5000,n.thin=2,n.chains=3)

print(exampleExpon)

#traceplot(exampleExpon)

attach.jags(exampleExpon)

par(mfrow=c(3,1))

hist(theta, xlab=expression(theta), freq=FALSE, main="Posterior theta")

# plot stand likelihood function and prior and posterior for theta
plot(density(theta, adjust = 1.5), ylab="", xlim=c(0,a/b+3*sqrt((a/b^2))), xlab=expression(theta), main="Prior Std likelihood and Posterior")
plot(function(x)dgamma(x,a,b),xlim=c(0,a/b+3*sqrt((a/b^2))),lty=3,add=TRUE)
likelihood = function(th){th^n*exp(-th*sum(y))}
K = integrate(likelihood,lower=0,upper=a/b+5*sqrt((a/b^2)))$value
curve(likelihood(x)/K,lty=2,add=TRUE)
legend("topright", c("prior","std likelihood","posterior"),lty=c(3,2,1))


# plot posterior predictive for the number of slugs. For alpha<=2 prior predictive does not work
predict.post = rexp(length(theta),rate=theta)
plot(density(predict.post,adjust=.2),xlim=c(0.5,mean(y)+3*sd(y)),ylab="Prob",xlab="y",main="Posterior predictive")


### Posterior predictive expected value
mean(predict.post)
### Posterior predictive median
median(predict.post)
### 90% posterior predictive credible interval
quantile(predict.post,probs=c(.05,.95))
### Posterior probability of life lenth longer than 10y
sum(10<predict.post)/length(predict.post)


############# Inference about expected life length

par(mfrow=c(2,1))

hist(1/theta,xlab=expression(1/theta),xlim=c(0,mean(y)+2*sd(y)),freq=FALSE, main="Posterior expected life length (yr)")
plot(density(1/theta,adjust=.5),ylab="",xlab=expression(1/theta),xlim=c(0,mean(y)+2*sd(y)),main="Posterior expected life length (yr)")

### Estimation of posterior expected value
mean(1/theta)
### Estimation of posterior median
median(1/theta)
### Estimation of 90% posterior credible interval
quantile(1/theta,probs=c(.05,.95))
### Estimation of posterior probability that 3 < Expected life length
sum(3<1/theta)/length(theta)


############# Inference about Pr(life length>10years)

hist(exp(-10*theta),xlab=expression(exp(-10*theta)),xlim=c(0,1), freq=FALSE, main="Posterior Prob(lifelength>10)")
plot(density(exp(-10*theta),adjust=.5),ylab="",xlab=expression(exp(-10*theta)),xlim=c(0,1),main="Posterior Prob(lifelength>10)")

### Estimation of posterior expected value
mean(exp(-10*theta))
### Estimation of posterior median
median(exp(-10*theta))
### Estimation of 90% posterior credible interval
quantile(exp(-10*theta),probs=c(.05,.95))
### Estimation of posterior probability that
sum(0.1<exp(-10*theta))/length(theta)


detach.jags()
