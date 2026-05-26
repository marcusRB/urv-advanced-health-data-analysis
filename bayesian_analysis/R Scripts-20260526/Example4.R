
############# Example 4. Binary response die

n=30
y=22

likel.funct = vector(length=7)
stdlikel.funct = vector(length=7)
prior.distr = vector(length=7)
post.distr = vector(length=7)


#Choice of prior distribution

# Uniform prior
uniform.prior = c(1/7,1/7,1/7,1/7,1/7,1/7,1/7)
# Binomial prior
binom.prior = vector(length=7)
for (i in 1:7){
  binom.prior[i] = dbinom(i-1,6,.5)
}

prior.distr = binom.prior


############## Compute likelihood and posterior and plot them together with prior

for (i in 1:7){
  likel.funct[i] = (((i-1)/6)^y)*(((7 -i)/6)^(n-y))
}
stdlikel.funct = likel.funct/sum(likel.funct)

post.distr = (prior.distr*likel.funct)/sum(prior.distr*likel.funct)

par(mfrow=c(3,1))
plot((0:6)/6,prior.distr,type="h",xaxt="n",ylim=c(0,0.8),xlab=expression(theta),ylab="",main="Prior distribution")
axis(1,at=round(seq(0,1,by=1/6),3),las=2)
plot((0:6)/6,stdlikel.funct,type="h",xaxt="n",ylim=c(0,0.8),xlab=expression(theta),ylab="",main="Likelihood function")
axis(1,at=round(seq(0,1,by=1/6),3),las=2)
plot((0:6)/6,post.distr,type="h",xaxt="n",ylim=c(0,0.8),xlab=expression(theta),ylab="",main="Posterior distribution")
axis(1,at=round(seq(0,1,by=1/6),3),las=2)


#### Compute prior and posterior predictive for n = tilden and plot them

tilden = 30
prior.pred = vector(length=tilden+1)
post.pred = vector(length=tilden+1)
prob.funct = vector(length=7)

for (j in 1:(tilden+1)){
  for (i in (1:7)){
    prob.funct[i] = dbinom(j-1,tilden,(i-1)/6)
  }
  prior.pred[j]=sum(prob.funct*prior.distr)
  post.pred[j]=sum(prob.funct*post.distr)
  }
#prior.pred
#post.pred

par(mfrow=c(2,1))
plot(0:tilden,prior.pred,type="h",xlim=c(0,tilden),ylim=c(0,0.15),xlab="y",ylab="",main="Prior predictive")
plot(0:tilden,post.pred,type="h",xlim=c(0,tilden),ylim=c(0,0.15),xlab="y",ylab="",main="Posterior predictive")
