hbgrnn.test<-function(y=NULL,n=NULL,x=NULL,M=5,adapt=4000,startburnin=1000,nChains = 2,startsample=10000,thin=20,max.time=Inf,scale=TRUE)
  {
  if (scale){
    x<-as.matrix(scale(x))
  }
  xl<-x
  xnl<-x
  N<-length(y)
  dl<-dim(x)[2]
  dnl<-dim(x)[2]
  mu.beta<-rep(0,dl)
  PHI<-diag(dl)
  dataList=list(
    n=n,
    xl=xl,
    xnl=xnl,
    mu.beta=mu.beta,
    PHI=PHI,
    y=y,
    M=M,
    dl=dl,
    dnl=dnl,
    N=N
  )
  
  
    parameters = c( "CPOinva","CPOinvb")
  
  pkgs <- c('rjags',  'runjags','coda')
  lapply(pkgs, require, character.only = T)
  modelString<- "
data
  {
for( i in 1 : N ) {
  
  ya[i]<-y[i]
  yb[i]<-y[i]
}
}
  model
  {
  
  for( j in 1 : dnl ) {
  for( k in 1 : M ) {
  # center x parameter
  mx[j , k] ~ dnorm(mx.0[j],dlambda[j])
  }
  }
  
  for( k in 1 : M ) {
  for( i in 1 : N ) {
  #hidden layer/pattern layer
  hid[i , k] <- exp(( - 0.5) * sum(z[i ,  , k]))
  }
  }
  for( k in 1 : M ) {
  #center y parameter
  v[k] ~ dnorm(v0,dlamv)
  }
  for( i in 1 : N ) {
  #output layer 
  logit(Pb[i])<-thetab[i]
  omegab[i]<-beta0b+inprod(xl[i,],betab[])+(Num[i] / Denum[i] )
  thetab[i]~dnorm(omegab[i],tau.omegab)
  logit(Pa[i])<-thetaa[i]
  omegaa[i]<-beta0a+inprod(xl[i,],betaa[])
  thetaa[i]~dnorm(omegaa[i],tau.omegaa)
  
  }
  #Prior
  betaa[1:dl]~dmnorm(mu.beta[1:dl],PHI[1:dl,1:dl])
  betab[1:dl]~dmnorm(mu.beta[1:dl],PHI[1:dl,1:dl])
  beta0a~dnorm(0,0.01)
  beta0b~dnorm(0,0.01)
  #Prior  tau dan sigma
  tau.omegaa ~ dgamma(2,1)
  tau.omegab ~ dgamma(2,1)
  
  for( i in 1 : N ) {
  
  ya[i] ~ dbin(Pa[i],n[i])
  #CPOinv
  logfya[i]<-logfact(n[i])-logfact(n[i]-ya[i])-logfact(ya[i])+ya[i]*log(Pa[i])+(n[i]-ya[i])*log(1-Pa[i])
  CPOinva[i]<-exp(-logfya[i]) 
  yb[i] ~ dbin(Pb[i],n[i])
  #CPOinv
logfyb[i]<-logfact(n[i])-logfact(n[i]-yb[i])-logfact(yb[i])+yb[i]*log(Pb[i])+(n[i]-yb[i])*log(1-Pb[i])
CPOinvb[i]<-exp(-logfyb[i]) 
  }
  #  tau ~ dgamma(0.001,0.001)
  #sigma <- 1 / sqrt(tau)
  for( j in 1 : dnl ) {
  h[j] <- 1 / sqrt(lambda[j])
  }
  for( j in 1 : dnl ) {
  for( k in 1 : M ) {
  for( i in 1 : N ) {
  z[i , j , k] <- pow((xnl[i , j] - mx[j , k]) / h[j],2)
  }
  }
  }
  for( j in 1 : dnl ) {
  lambda[j] ~ dgamma(2, 5)
  dlambda[j]<-lambda[j]/delta
  }
  
  for( k in 1 : M ) {
  wv[k] <- w[k] * v[k]
  }
  for( i in 1 : N ) {
  #node Numerator 
  Num[i] <- inprod(wv[],hid[i , ])
  }
  
  for (k in 1:M) 
  {
  w[ k] <- eta[k] / sum(eta[])
  eta[ k] ~ dgamma(alpha[k], 1)
  alpha[k] ~ dgamma(5,2)
  }
  
  for( i in 1 : N ) {
  #Node denumerator      
  Denum[i] <- inprod(w[],hid[i , ])
  }
  for( j in 1 : dnl ) {
  #mean center x parameter 
  mx.0[j] ~ dnorm( 0,1)
  }
  
  
  #mean center x parameter
  v0 ~ dnorm( 0,1)
  
  lamv ~ dgamma(1, 0.1)
  dlamv<-lamv/delta
  delta~dgamma(2,1)
}
"

model <- read.winbugs(modelString)
  
    obj<- autorun.jags( method=c("rjags","parallel")[2],
                                model=model , 
                                monitor=parameters , 
                                data=dataList ,  
                                n.chains=nChains ,
                                adapt=adapt ,
                           startburnin=startburnin , 
                           startsample=startsample ,
                           max.time=max.time,
                           thin=thin ,
                                summarise=TRUE )
  coda <- as.mcmc.list(obj)
  mcmc = as.matrix(coda,chains=TRUE)
  summaryout<-summary(obj)
  
  
  varnameCPOinva=NULL
  varnameCPOinvb=NULL
  
  for (i in 1:N){varnameCPOinva[i]<-paste("CPOinva","[",i,"]",sep="")}
  for (i in 1:N){varnameCPOinvb[i]<-paste("CPOinvb","[",i,"]",sep="")}
  
  summary<-NULL
  
  SCPOinva<-summaryout[varnameCPOinva,]
  SCPOinvb<-summaryout[varnameCPOinvb,]
  CPOa<-1/SCPOinva[,"Mean"]
  CPOb<-1/SCPOinvb[,"Mean"]
  summary$LPMLa <- sum(log(CPOa))
  summary$LPMLb <- sum(log(CPOb))
  summary$BF<-exp(summary$LPMLb-summary$LPMLa)
  summary$convergence$gbr<-gelman.diag(coda, confidence = 0.95)
  summary$convergence$heidel<-heidel.diag(coda, eps=0.1, pvalue=0.05)
  summary$convergence$geweke<-geweke.diag(coda, frac1=0.1, frac2=0.5)
  summary$convergence$raftery<-raftery.diag(coda, q=0.025, r=0.005, s=0.95, converge.eps=0.001)
  return(list(summary=summary,coda=coda)) 
}