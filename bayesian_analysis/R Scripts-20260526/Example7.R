
######### Example 7. Binary logistic model

Ky = read.csv2("kyphosis.csv")
head(Ky)

pairs(~Kyph+Age+Number+Start,lower.panel=NULL,Ky)

mod1 = glm(Kyph ~ Age + Number + Start, family=binomial,Ky)
summary(mod1)

library(R2jags)

model.bin00 <- "
model
{
  for (i in 1:nn) {  
#    y[i] ~ dbin(pi[i], n)
    y[i] ~ dbern(pi[i])
    pi[i] <- g[i]/(1+g[i])
      g[i] <- exp(b0 + b1*x1[i] + b2*x2[i] + b3*x3[i])
#      g[i] <- exp(max(-20,min(20,b0 + b1*x1[i] + b2*x2[i] + b3*x3[i])))
  }
  
  b0 ~ dnorm(0, 0.000001)
	b1 ~ dnorm(0, 0.000001)
	b2 ~ dnorm(0, 0.000001)
	b3 ~ dnorm(0, 0.000001)
	}
"
  nn = length(Ky$Kyph)
  nn
 
  data <- list(nn=nn, y=Ky$Kyph, x1=Ky$Age, x2=Ky$Number,x3=Ky$Start) 
  
  parameters <- c("b0","b1","b2","b3")
  
  example <-jags(data, parameters.to.save=parameters,
                 model=textConnection(model.bin00), n.iter=55000, n.burnin=5000, n.thin=2,n.chains=3)
  
  print(example)
  
# traceplot(example)
  
  attach.jags(example)
  
  par(mfrow=c(2,2))
  
  #beta0
  
  plot(density(b0, adjust = 1.5), ylab="", xlab=expression(b0), main="")
  abline(v=c(0),lty=1)
  plot(density(b1, adjust = 1.5), ylab="", xlab=expression(bage), main="")
  abline(v=c(0),lty=1)
  plot(density(b2, adjust = 1.5), ylab="", xlab=expression(bnumber), main="")
  abline(v=c(0),lty=1)
  plot(density(b3, adjust = 1.5), ylab="", xlab=expression(bstart), main="")
  abline(v=c(0),lty=1)
  
  plot(b1,b2)
  abline(v=c(0),lty=1)
  abline(h=c(0),lty=1)
  plot(b1,b3)
  abline(v=c(0),lty=1)
  abline(h=c(0),lty=1)
  plot(b2,b3)
  abline(v=c(0),lty=1)
  abline(h=c(0),lty=1)
  
  mean(b1)
  quantile(b1,probs=c(.05,.95))
  sum(0<b1)/length(b1)  
  mean(b2)
  quantile(b2,probs=c(.05,.95))
  sum(0<b2)/length(b2)  
  mean(b3)
  quantile(b3,probs=c(.05,.95))
  sum(0<b3)/length(b3)  

mean(b1) 
mean(b2) 
mean(b3)
  
quantile(b1,probs=c(.05,.95))
quantile(b2,probs=c(.05,.95))
quantile(b3,probs=c(.05,.95))

sum(0<b1)/length(b1) 
sum(0<b2)/length(b2) 
sum(0<b3)/length(b3) 
 
  
  
  ###### Inference about Prob(complications) for specific patient, (5,2,3)
  
  par(mfrow=c(2,1))
  
  prb=exp(b0+b1*5+b2*2+b3*3)/(1+exp(b0+b1*5+b2*2+b3*3))
  hist(prb,xlim=c(0,1), freq=FALSE, nclass=40, main="Posterior Prob(Compl) for (5,2,3)")
  abline(v=0.5,lty=1)
  plot(density(prb,adjust=1.5),ylab="",xlab="",xlim=c(0,1),main="Posterior Prob(Compl) for (5,2,3)")
  abline(v=0.5,lty=1)
  
  ### Estimation of posterior expected value
  mean(prb)
  ### Estimation of posterior median
  median(prb)
  ### Estimation of 90% posterior credible interval
  quantile(prb,probs=c(.05,.95))
  ### Estimation of posterior probability that 
  sum(.5 < prb)/length(b1)
  
  ###### Prediction presence complications for specific patient, (5,2,3)
  
  predict.post= rbinom(length(b1),1,prob=prb)
  tbpp=table(predict.post)
  plot(tbpp/length(b1),ylab="Prob",xlab="compl",main="Post predict for (5,2,3)")
  
  detach.jags()
  