
########### Example 5. Comparison of two Normal Normal+IGamma; 

# data
Bw = read.csv2("BirthWeight0.csv")
head(Bw)
ydnsmk=Bw$bwt[Bw$smoke==0]
#y1
ysmk=Bw$bwt[Bw$smoke==1]
#y2
n1=length(ydnsmk)
n1
n2=length(ysmk)
n2
summary(ydnsmk)
summary(ysmk)

# Prior distribution. 
mu1.prior = c(m1=3000,t1=0.000001)
mu2.prior = c(m2=3000,t2=0.000001)
tau1.prior = c(a1=0.01,b1=0.01)
tau2.prior = c(a2=0.01,b2=0.01)

######################### Solution through MCMC simulation (JAGS)

# YOU MUST INSTALL JAGS FROM: https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/
# You must install the library R2WinBUGS, install.packages('R2jags', repos = "http://cran.fhcrc.org/")

library(R2jags)

model.comp2norm <- "
 model
 { 
  for (i in 1:n1) {
   y1[i] ~ dnorm(mu1,tau1)
  }
   mu1 ~ dnorm(mu1.prior[1],mu1.prior[2])
   tau1 ~ dgamma(tau1.prior[1],tau1.prior[2])
  for (j in 1:n2) {
   y2[j] ~ dnorm(mu2,tau2)
   }
   mu2 ~ dnorm(mu2.prior[1],mu2.prior[2])
   tau2 ~ dgamma(tau2.prior[1],tau2.prior[2])
 }
"

data <- list(y1=ydnsmk, n1=n1, mu1.prior=mu1.prior, tau1.prior=tau1.prior, y2=ysmk, n2=n2, mu2.prior=mu2.prior, tau2.prior=tau2.prior) 

parameters <- c("mu1","tau1","mu2","tau2") 

examplec2norm <-jags(data, parameters.to.save=parameters,
                  model=textConnection(model.comp2norm),n.iter=55000,n.burnin=5000,n.thin=2,n.chains=3)

print(examplec2norm)

#traceplot(examplec2norm)

attach.jags(examplec2norm)

######## Inference th1-th2

par(mfrow=c(2,1))

hist(mu1-mu2,xlim=c(-200,800),freq=FALSE, main="Posterior mu1-mu2 in gr.")
abline(v=c(0),lty=1)
# plot pdf posterior for mu1-mu2
plot(density(mu1-mu2, adjust = 1.5), ylab="",xlim=c(-200,800),xlab=expression(mu1-mu2), main="Posterior mu1-mu2 in gr.")
#legend("topright", c("prior","posterior"),lty=c(3,1))
abline(v=c(0),lty=1)

### Estimation of posterior expected value
mean(mu1-mu2)
### Estimation of posterior median
median(mu1-mu2)
### Estimation of 90% posterior credible interval
quantile(mu1-mu2,probs=c(.05,.95))
### Estimation of posterior probability that 
sum(0<mu1-mu2)/length(mu1)


detach.jags()

