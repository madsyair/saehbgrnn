auto.saehb<-function(y=y,n=n,x=x,adapt=4000,startburnin=1000,nChains = 2,startsample=10000,thinSteps=20,DIC=FALSE)
{
  N<-length(y)
  d<-dim(x)[2]
  mu.beta<-rep(0,d)
  PHI<-diag(d)
  
  dataList=list(
    N=N,
    x=x,
    y=y,
  
    d=d,
    n=n,
    mu.beta=mu.beta,
    PHI=PHI
  )
  
  
  if(DIC){
    parameters = c(  "CPOinv","P"  ,"beta0","beta","sigma","deviance", "pd", "popt", "dic", "ped","full.pd")
  }
    else{
  parameters=c( "CPOinv","P","beta0","beta","sigma")}
  
  pkgs <- c('rjags',  'runjags','coda')
  lapply(pkgs, require, character.only = T)
  modelString<- "
  model
{
  for(i in 1:N) {
  
  y[i]~dbin(P[i],n[i])  #Hirarki 1: sampling model (likelihood)
  #v[i]~dnorm(0,tau)  #Komponen Random Effect y
  
  #Hirarki 2: Linking Model
  logit(P[i])<-B[i]
  mu[i]<-beta0+inprod(x[i,],beta[])
  B[i]~dnorm(mu[i],tau)
#CPOinv
 logfy[i]<-logfact(n[i])-logfact(n[i]-y[i])-logfact(y[i])+y[i]*log(P[i])+(n[i]-y[i])*log(1-P[i])
  CPOinv[i]<-exp(-logfy[i]) 
  }
  
  #Prior
  beta[1:d]~dmnorm(mu.beta[1:d],PHI[1:d,1:d])
  beta0~dnorm(0,0.01)
  #Gamma sebagai hyperprior varians Random Effect
  tau~dgamma(0.01,0.01)
  sigma<-1/sqrt(tau)
} 
  
  "
  
  model <- read.winbugs(modelString)
  if(DIC){
    obj<- autorun.jags( method=c("rjags","parallel")[1] ,
                                model=model , 
                                monitor=parameters , 
                                data=dataList ,  
                                n.chains=nChains ,
                                adapt=adapt ,
                                startburnin=startburnin , 
                                startsample=startsample ,
                                thin=thinSteps ,
                                summarise=TRUE ,
                                plots=TRUE,monitor.deviance =TRUE,
                                monitor.pd = TRUE)
  }
  else{
    obj<- autorun.jags( method=c("rjags","parallel")[2] ,
                                model=model , 
                                monitor=parameters , 
                                data=dataList ,  
                                n.chains=nChains ,
                                adapt=adapt ,
                           startburnin=startburnin , 
                           startsample=startsample ,
                           thin=thinSteps ,
                                summarise=TRUE ,
                                plots=TRUE,monitor.deviance =FALSE,
                                monitor.pd = FALSE)}
  coda <- as.mcmc.list(obj)
  mcmc = as.matrix(coda,chains=TRUE)
  summaryout<-summary(obj)
  varnamepar<-varnames(coda)
  varnamep=NULL
  varnameCPOinv=NULL
  for (i in 1:N){varnamep[i]<-paste("P","[",i,"]",sep="")}
  for (i in 1:N){varnameCPOinv[i]<-paste("CPOinv","[",i,"]",sep="")}
  varnamepar<-varnamepar[ -match(varnamep,varnamepar)]
  summary<-NULL
  summary$P<-summaryout[varnamep,]
  SCPOinv<-summaryout[varnameCPOinv,]
  CPO<-1/SCPOinv[,"Mean"]
  summary$LPML <- sum(log(CPO))
  CV<-summary$P[,"SD"]/summary$P[,"Mean"]
  summary$P<-cbind(summary$P,CV)
  summary$modelpar<-summaryout[varnamepar,]
  summary$convergence$gbr<-gelman.diag(coda, confidence = 0.95)
  summary$convergence$heidel<-heidel.diag(coda, eps=0.1, pvalue=0.05)
  summary$convergence$geweke<-geweke.diag(coda, frac1=0.1, frac2=0.5)
  summary$convergence$raftery<-raftery.diag(coda, q=0.025, r=0.005, s=0.95, converge.eps=0.001)
  return(list(obj=obj,summary=summary,coda=coda)) 
}