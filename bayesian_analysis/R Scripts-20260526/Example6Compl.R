

############ Example 6. Normal linear model

Bw = read.csv2("BirthWeight0.csv")
head(Bw)

Bw$r1 = ifelse(Bw$race==1,1,0)
Bw$r2 = ifelse(Bw$race==2,1,0)
head(Bw)

pairs(~bwt+age+lwt,lower.panel=NULL,col=as.integer(Bw$smoke)+1,Bw)

# Additive model for bwt with all the x's; It assumes parallel planes.

model00 = lm(bwt~age+lwt+r1+r2+smoke+ht+ui+ptl+ftv, Bw)
summary(model00)

plot(predict(model00),rstandard(model00), col=as.integer(Bw$race),xlab="Predicted",ylab="StdRes")
abline(h=c(-2,0,2),lty=2)
qqnorm(rstandard(model00))

library(R2jags)

model.lm00 <- "
model
{
  for (i in 1:nn) {  
    y[i] ~ dnorm(mu[i],tau)
    mu[i] <- b0+b1*x1[i]+b2*x2[i]+b3*x3[i]+b4*x4[i]+b5*x5[i]+b6*x6[i]+b7*x7[i]+b8*x8[i]+b9*x9[i]
  }
  
  b0 ~ dnorm(0, 0.000001)
	b1 ~ dnorm(0, 0.000001)
	b2 ~ dnorm(0, 0.000001)
	b3 ~ dnorm(0, 0.000001)
	b4 ~ dnorm(0, 0.000001)
	b5 ~ dnorm(0, 0.000001)
	b6 ~ dnorm(0, 0.000001)
	b7 ~ dnorm(0, 0.000001)
	b8 ~ dnorm(0, 0.000001)
	b9 ~ dnorm(0, 0.000001)
	tau ~ dgamma(0.01, 0.01)
	}
"
  nn = length(Bw$bwt)
  nn
  
  data <- list(nn=nn,y=Bw$bwt,x1=Bw$age,x2=Bw$lwt,x3=Bw$r1,x4=Bw$r2,x5=Bw$smoke,x6=Bw$ht,x7=Bw$ui,x8=Bw$ptl,x9=Bw$ftv) 
  
  parameters <- c("b0","b1","b2","b3","b4","b5","b6","b7","b8","b9") 

  example <-jags(data, parameters.to.save=parameters,
                 model=textConnection(model.lm00),n.iter=55000,n.burnin=10000,n.thin=2,n.chains=3)
  
  print(example)
  
  attach.jags(example)
  
  par(mfrow=c(3,3))
  
#  plot(density(b0, adjust = 1.5), ylab="", xlab=expression(b0), main="")
  plot(density(b1, adjust = 1.5), ylab="", xlab=expression(bage), main="")
  abline(v=c(0),lty=1)
  plot(density(b2, adjust = 1.5), ylab="", xlab=expression(blwt), main="")
  abline(v=c(0),lty=1)
  plot(density(b3, adjust = 1.5), ylab="", xlab=expression(br1), main="")
  abline(v=c(0),lty=1)
  plot(density(b4, adjust = 1.5), ylab="", xlab=expression(br2), main="")
  abline(v=c(0),lty=1)
  plot(density(b5, adjust = 1.5), ylab="", xlab=expression(bsmoke), main="")
  abline(v=c(0),lty=1)
  plot(density(b6, adjust = 1.5), ylab="", xlab=expression(bht), main="")
  abline(v=c(0),lty=1)
  plot(density(b7, adjust = 1.5), ylab="", xlab=expression(bui), main="")
  abline(v=c(0),lty=1)
  plot(density(b8, adjust = 1.5), ylab="", xlab=expression(bptl), main="")
  abline(v=c(0),lty=1)
  plot(density(b9, adjust = 1.5), ylab="", xlab=expression(bftv), main="")
  abline(v=c(0),lty=1)
  
  
mean(b1)
mean(b2)
mean(b3)
mean(b4)
mean(b5)
mean(b6)
mean(b7)
mean(b8)
mean(b9)
  
quantile(b1,probs=c(.05,.95))
quantile(b2,probs=c(.05,.95))
quantile(b3,probs=c(.05,.95))
quantile(b4,probs=c(.05,.95))
quantile(b5,probs=c(.05,.95))
quantile(b6,probs=c(.05,.95))
quantile(b7,probs=c(.05,.95))
quantile(b8,probs=c(.05,.95))
quantile(b9,probs=c(.05,.95))
  
  

sum(0<b1)/length(b1)  
sum(0<b2)/length(b2)  
sum(0<b3)/length(b3)  
sum(0<b4)/length(b4)  
sum(0<b5)/length(b5)  
sum(0<b6)/length(b6)  
sum(0<b7)/length(b7)  
sum(0<b8)/length(b8)  
sum(0<b9)/length(b9)  
  
  

detach.jags()
