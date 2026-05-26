
######### Example 8. Hierarchical binomial

P06 = read.csv2("Parlament2006.csv")
head(P06)

pairs(~CIU+PSC+PP+ICV+ERC+abstencio,lower.panel=NULL,P06)

P06$Sobir = P06$CIU + P06$ERC
P06$Constit = P06$PSC + P06$PP + P06$ciudadanos
P06$BlNul = P06$blancs + P06$nuls
P06$elec = as.integer(P06$elec)
head(P06)

summary(P06)

sum(P06$abstencio<P06$elec)/length(P06$elec)



library(R2jags)

model.binhier <- "
model
{
  for (i in 1:nn) { 
    y1[i] ~ dbin(pi1[i],Vot[i])
    pi1[i] ~ dbeta(alf1,bet1)
    y2[i] ~ dbin(pi2[i],Vot[i])
    pi2[i] ~ dbeta(alf2,bet2)
    y3[i] ~ dbin(pi3[i],Vot[i])
    pi3[i] ~ dbeta(alf3,bet3)
    y4[i] ~ dbin(pi4[i],Vot[i])
    pi4[i] ~ dbeta(alf4,bet4)
  }
  alf1 ~ dexp(1)
  bet1 ~ dexp(1)
  alf2 ~ dexp(1)
  bet2 ~ dexp(1)
  alf3 ~ dexp(1)
  bet3 ~ dexp(1)
  alf4 ~ dexp(1)
  bet4 ~ dexp(1)  
  
	}
"

  
  data <- list(nn=length(P06$abstencio), y1=P06$abstencio, y2=P06$Sobir, y3=P06$Constit, y4=P06$BlNul, Vot=P06$elec) 
  
  parameters <- c("pi1","alf1","bet1","alf2","bet2","alf3","bet3","alf4","bet4")
  
  example <-jags(data, parameters.to.save=parameters,
                 model=textConnection(model.binhier), n.iter=3000, n.burnin=1000, n.thin=2,n.chains=3)
  
  print(example)
  
  traceplot(example)
  
  attach.jags(example)
  
  par(mfrow=c(2,2))

  plot(alf1,bet1,main="a1, b1 Abst")
  plot(alf2,bet2,main="a2, b2 Sobir")
  plot(alf3,bet3,main="a3, b3 Constit")
  plot(alf4,bet4,main="a4, b4 BlNul")
  

  ### Estimation of 90% posterior credible interval
  quantile(alf1,probs=c(.05,.5,.95))
  
  ### Estimation of 90% posterior credible interval
  quantile(bet1,probs=c(.05,.5,.95))
  
  
  
  par(mfrow=c(2,2))
  
  plot(density(alf1, adjust = 1.5), ylab="", xlab=expression(alpha), main="alpha")
  plot(density(bet1, adjust = 1.5), ylab="", xlab=expression(beta), main="beta")
  plot(density(alf1/(alf1+bet1), adjust = 1.5), ylab="", xlab=expression(alpha/(alpha+beta)), main="mu=E(theta_i)")
  plot(density(alf1*bet1/((alf1+bet1)^2*(alf1+bet1+1)), adjust = 1.5), ylab="", xlab=expression(V(theta)), main="V(theta_i)")
  
  
  aaa=mean(alf1)
  bbb=mean(bet1)
  plot(function(x)dbeta(x,aaa,bbb), xlim=c(0,1), xlab=expression(theta_i),ylab="",main="PDistr theta_i Abs")
  prtheta = rbeta(length(alf1),alf1,bet1)
  plot(density(prtheta, adjust = 1.5), ylab="", xlim=c(0,1), xlab=expression(theta_i), main="PDistr theta_i Abs")
  
  
  
  par(mfrow=c(2,2))
  
  plot(density(alf2, adjust = 1.5), ylab="", xlab=expression(alpha|y), main="alpha2 Sob")
  plot(density(bet2, adjust = 1.5), ylab="", xlab=expression(beta|y), main="beta2 Sob")
  plot(density(alf2/(alf2+bet2), adjust = 1.5), ylab="", xlab=expression(alpha/(alpha+beta)|y), main="mu=E(theta_i|y) Sob")
  plot(density(alf2*bet2/((alf2+bet2)^2*(alf2+bet2+1)), adjust = 1.5), ylab="", xlab=expression(V(theta|y)), main="V(theta_i|y) Sob")
  
  aaa2=mean(alf2)
  bbb2=mean(bet2)
  plot(function(x)dbeta(x,aaa2,bbb2), xlim=c(0,1), xlab=expression(theta_i),ylab="",main="PDistr theta_i Sob")
  prtheta = rbeta(length(alf2),alf2,bet2)
  plot(density(prtheta, adjust = 1.5), ylab="", xlim=c(0,1), xlab=expression(theta_i), main="PDistr theta_i")
  
  
  
  
  par(mfrow=c(2,2))
  
  plot(density(alf3, adjust = 1.5), ylab="", xlab=expression(alpha|y), main="alpha3 Const")
  plot(density(bet3, adjust = 1.5), ylab="", xlab=expression(beta|y), main="beta3 Const")
  plot(density(alf3/(alf3+bet3), adjust = 1.5), ylab="", xlab=expression(alpha/(alpha+beta)|y), main="mu=E(theta_i|y) Const")
  plot(density(alf3*bet3/((alf3+bet3)^2*(alf3+bet3+1)), adjust = 1.5), ylab="", xlab=expression(V(theta|y)), main="V(theta_i|y) Const")
  
  aaa3=mean(alf3)
  bbb3=mean(bet3)
  plot(function(x)dbeta(x,aaa3,bbb3), xlim=c(0,1), xlab=expression(theta_i),ylab="",main="PDistr theta_i Const")
  prtheta = rbeta(length(alf3),alf3,bet3)
  plot(density(prtheta, adjust = 1.5), ylab="", xlim=c(0,1), xlab=expression(theta_i), main="PDistr theta_i")
  
  
  
  
  par(mfrow=c(2,2))
  
  plot(density(alf4, adjust = 1.5), ylab="", xlab=expression(alpha|y), main="alpha4 BlNl")
  plot(density(bet4, adjust = 1.5), ylab="", xlab=expression(beta|y), main="beta4 BlNl")
  plot(density(alf4/(alf4+bet4), adjust = 1.5), ylab="", xlab=expression(alpha/(alpha+beta)|y), main="mu=E(theta_i|y)")
  plot(density(alf4*bet4/((alf4+bet4)^2*(alf4+bet4+1)), adjust = 1.5), ylab="", xlab=expression(V(theta|y)), main="V(theta_i|y)")
  
  aaa4=mean(alf4)
  bbb4=mean(bet4)
  plot(function(x)dbeta(x,aaa4,bbb4), xlim=c(0,1), xlab=expression(theta_i),ylab="",main="PDistr theta_i BlNl")
  prtheta = rbeta(length(alf4),alf4,bet4)
  plot(density(prtheta, adjust = 1.5), ylab="", xlim=c(0,1), xlab=expression(theta_i), main="PDistr theta_i")
  
  
  pi11=pi1[(1:3000)]
  pi12=pi1[3000+(1:3000)]
  pi13=pi1[2*3000+(1:3000)]
  pi14=pi1[3*3000+(1:3000)]
  pi15=pi1[4*3000+(1:3000)]
  pi16=pi1[5*3000+(1:3000)]
  pi17=pi1[6*3000+(1:3000)]
  pi18=pi1[7*3000+(1:3000)]
  pi19=pi1[8*3000+(1:3000)]
  
 
  quantile(pi11,probs=c(.05,.5,.95))
  quantile(pi12,probs=c(.05,.5,.95))  
  quantile(pi13,probs=c(.05,.5,.95))
  quantile(pi14,probs=c(.05,.5,.95))   
  quantile(pi15,probs=c(.05,.5,.95))
  quantile(pi16,probs=c(.05,.5,.95))  
  quantile(pi17,probs=c(.05,.5,.95))
  quantile(pi18,probs=c(.05,.5,.95))  
  quantile(pi19,probs=c(.05,.5,.95))  
  
  # plot stand likelihood, prior and posterior for theta
  par(mfrow=c(3,3))
  plot(density(pi11, adjust = 1.5), ylab="", xlab="Abrera", main="Posterior pi(1)")
  plot(density(pi12, adjust = 1.5), ylab="", xlab="Aguilar de S.", main="Posterior pi(2)")
  plot(density(pi13, adjust = 1.5), ylab="", xlab="Alella", main="Posterior pi(3)")
  plot(density(pi14, adjust = 1.5), ylab="", xlab="Alpens", main="Posterior pi(4)")
  plot(density(pi15, adjust = 1.5), ylab="", xlab="Ametlla de V", main="Posterior pi(5)")
  plot(density(pi16, adjust = 1.5), ylab="", xlab="Arenys de Mar 1", main="Posterior pi(6)")  
  plot(density(pi17, adjust = 1.5), ylab="", xlab="Arenys de Mar 2", main="Posterior pi(7)")
  plot(density(pi18, adjust = 1.5), ylab="", xlab="Arenys de Munt 1", main="Posterior pi(8)")
  plot(density(pi19, adjust = 1.5), ylab="", xlab="Arenys de Munt 2", main="Posterior pi(9)")  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ###### Inference about Prob(complications) for specific patient, (5,2,3)
  
  par(mfrow=c(2,2))
  
  prb=exp(b0+b1*5+b2*2+b3*3)/(1+exp(b0+b1*5+b2*2+b3*3))
  hist(prb,xlim=c(0,1),main="Post Prob(Compl) for (5,2,3)")
  plot(density(prb,adjust=1.5),ylab="",xlab="",xlim=c(0,1),main="Post Prob(Compl) for (5,2,3)")
  
  ### Estimation of posterior expected value
  mean(prb)
  ### Estimation of posterior median
  median(prb)
  ### Estimation of 90% posterior credible interval
  quantile(prb,probs=c(.05,.95))
  ### Estimation of posterior probability that 
  sum(.5 < prb)/length(b1)
  