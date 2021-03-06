saehbgrnn<-function(y=NULL,n=NULL,xl=NULL,xnl=NULL,M=5,adapt=4000,burnin=1000,nChains = 2,sample=10000,thin=1,DIC=FALSE,scale=TRUE)
{
   
  pkgs <- c('rjags',  'runjags','coda','ggmcmc')
  
  lapply(pkgs, require, character.only = T)
  nulxl=is.null(xl)
  nulxnl=is.null(xnl)
  if (nulxl==TRUE&&nulxnl==FALSE){
    if(scale) {
      
      xnl<- as.matrix(scale(xnl))
    }
    
  }else if(nulxl==FALSE&&nulxnl==TRUE){
    if(scale) {
      xl <- as.matrix(scale(xl))
      
    }
  }else if(nulxl==FALSE&&nulxnl==FALSE){
    if(scale) {
      xl <- as.matrix(scale(xl))
      xnl<- as.matrix(scale(xnl))
    }
  }
  
  if (nulxl==TRUE&&nulxnl==FALSE){
    #result<-sae.hbgrnn1(y=y,n=n,x=xnl,M=M,adapt=adapt,burnin=burnin,nChains=nChains,sample=sample,thin=thin,DIC=DIC)
    type="hbgrnn1"
    N<-length(y)
    d<-dim(xnl)[2]
    dataList=list(
      n=n,
      x=xnl,
      y=y,
      M=M,
      d=d,
      N=N
    )
        if(DIC){
      parameters = c( "CPOinv", "w","mx","P"  , "v","h","deviance", "pd", "popt", "dic", "ped","full.pd")
    }
    else{
      parameters=c(  "CPOinv","w","mx","P","v","h")}
    
    
    modelString<- "
    model
    {
    
    for( j in 1 : d ) {
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
    logit(P[i])<-theta[i]
    omega[i]<-(Num[i] / Denum[i] )
    theta[i]~dnorm(omega[i],tau.omega)
    
    
    }
    
    #Prior  tau dan sigma
    tau.omega ~ dgamma(2,1)
    
    for( i in 1 : N ) {
    
    y[i] ~ dbin(P[i],n[i])
    #CPOinv
    logfy[i]<-logfact(n[i])-logfact(n[i]-y[i])-logfact(y[i])+y[i]*log(P[i])+(n[i]-y[i])*log(1-P[i])
    CPOinv[i]<-exp(-logfy[i]) 
    }
    #  tau ~ dgamma(0.001,0.001)
    #sigma <- 1 / sqrt(tau)
    for( j in 1 : d ) {
    h[j] <- 1 / sqrt(lambda[j])
    }
    for( j in 1 : d ) {
    for( k in 1 : M ) {
    for( i in 1 : N ) {
    z[i , j , k] <- pow((x[i , j] - mx[j , k]) / h[j],2)
    }
    }
    }
    for( j in 1 : d ) {
    lambda[j] ~ dgamma(2, 1)
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
    alpha[k] ~ dgamma(2,1)
    }
    
    for( i in 1 : N ) {
    #Node denumerator      
    Denum[i] <- inprod(w[],hid[i , ])
    }
    for( j in 1 : d ) {
    #mean center x parameter 
    mx.0[j] ~ dnorm( 0,1)
    }
    
    
    #mean center x parameter
    v0 ~ dnorm( 0,1)
    
    lamv ~ dgamma(2, 1)
    dlamv<-lamv/delta
    delta~dgamma(2,1)
    }
    "
    
    model <- read.winbugs(modelString)
    if(DIC){
      obj<- run.jags( method=c("rjags","parallel")[1] ,
                          model=model , 
                          monitor=parameters , 
                          data=dataList ,  
                          n.chains=nChains ,
                          adapt=adapt ,
                          burnin=burnin , 
                          sample=sample ,
                          thin=thin ,
                          summarise=TRUE )
    }
    else{
      obj<- run.jags( method=c("rjags","parallel")[2] ,
                          model=model , 
                          monitor=parameters , 
                          data=dataList ,  
                          n.chains=nChains ,
                          adapt=adapt ,
                          burnin=burnin , 
                          sample=sample ,
                          thin=thin ,
                          summarise=TRUE)}
    
    }else if(nulxl==FALSE&&nulxnl==TRUE){
    #result<-sae.hb(y=y,n=n,x=xl,adapt=adapt,burnin=burnin,nChains=nChains,sample=sample,thin=thin,DIC=DIC)
    type="hb"
    N<-length(y)
    d<-dim(xl)[2]
    mu.beta<-rep(0,d)
    PHI<-diag(d)
    
    dataList=list(
      N=N,
      x=xl,
      y=y,
      
      d=d,
      n=n,
      mu.beta=mu.beta,
      PHI=PHI
    )
    
    
    if(DIC){
      parameters = c( "CPOinv","P"  ,"beta0","beta","sigma","deviance", "pd", "popt", "dic", "ped","full.pd")
    }
    else{
      parameters=c("CPOinv","P","beta0","beta","sigma")}
    
    
    modelString<- "
    model
    {
    for(i in 1:N) {
    
    y[i]~dbin(P[i],n[i])  #Hirarki 1: sampling model (likelihood)
    #v[i]~dnorm(0,tau)  #Komponen Random Effect y
    
    #Hirarki 2: Linking Model
    logit(P[i])<-theta[i]
    mu[i]<-beta0+inprod(x[i,],beta[])
    theta[i]~dnorm(mu[i],tau)
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
      obj<- run.jags( method=c("rjags","parallel")[1] ,
                      model=model , 
                      monitor=parameters , 
                      data=dataList ,  
                      n.chains=nChains ,
                      adapt=adapt ,
                      burnin=burnin , 
                      sample=sample ,
                      thin=thin ,
                      summarise=TRUE 
      )
    }
    else{
      obj<- run.jags( method=c("rjags","parallel")[2] ,
                      model=model , 
                      monitor=parameters , 
                      data=dataList ,  
                      n.chains=nChains ,
                      adapt=adapt ,
                      burnin=burnin , 
                      sample=sample ,
                      thin=thin ,
                      summarise=TRUE)}
    
    }else if(nulxl==FALSE&&nulxnl==FALSE){
    #result<-sae.hbgrnn2(y=y,n=n,xl=xl,xnl=xnl,M=M,adapt=adapt,burnin=burnin,nChains=nChains,sample=sample,thin=thin,DIC=DIC)
    type="hbgrnn2"
    N<-length(y)
    dl<-dim(xl)[2]
    dnl<-dim(xnl)[2]
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
    
    
    if(DIC){
      parameters = c(  "CPOinv","w","mx","P","beta0","beta"  , "v","h","deviance", "pd", "popt", "dic", "ped","full.pd")
    }
    else{
      parameters = c( "CPOinv","w","mx","P","beta0","beta", "v","h")}
    
    
    modelString<- "
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
    logit(P[i])<-theta[i]
    omega[i]<-beta0+inprod(xl[i,],beta[])+(Num[i] / Denum[i] )
    theta[i]~dnorm(omega[i],tau.omega)
    #CPOinv
    logfy[i]<-logfact(n[i])-logfact(n[i]-y[i])-logfact(y[i])+y[i]*log(P[i])+(n[i]-y[i])*log(1-P[i])
    CPOinv[i]<-exp(-logfy[i]) 
    
    }
    #Prior
    beta[1:dl]~dmnorm(mu.beta[1:dl],PHI[1:dl,1:dl])
    beta0~dnorm(0,0.01)
    #Prior  tau dan sigma
    tau.omega ~ dgamma(2,1)
    
    for( i in 1 : N ) {
    
    y[i] ~ dbin(P[i],n[i])
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
    if(DIC){
      obj<- run.jags( method=c("rjags","parallel")[1] ,
                          model=model , 
                          monitor=parameters , 
                          data=dataList ,  
                          n.chains=nChains ,
                          adapt=adapt ,
                          burnin=burnin , 
                          sample=sample ,
                          thin=thin ,
                          summarise=TRUE )
    }
    else{
      obj<- run.jags( method=c("rjags","parallel")[2],
                          model=model , 
                          monitor=parameters , 
                          data=dataList ,  
                          n.chains=nChains ,
                          adapt=adapt ,
                          burnin=burnin , 
                          sample=sample ,
                          thin=thin ,
                          summarise=TRUE )}
        }
  else{stop("input is not complete")
  }
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
  obj$N<-N
  summary$modelpar<-summaryout[varnamepar,]
  summary$convergence$gbr<-gelman.diag(coda, confidence = 0.95)
  summary$convergence$heidel<-heidel.diag(coda, eps=0.1, pvalue=0.05)
  summary$convergence$geweke<-geweke.diag(coda, frac1=0.1, frac2=0.5)
  summary$convergence$raftery<-raftery.diag(coda, q=0.025, r=0.005, s=0.95, converge.eps=0.001)
  
  ggso<-ggs(coda)
  t<-Sys.time()
  nt <- as.POSIXlt(t)
  dt<-format(t, format="%y%m%d")
  ggmcmc(ggso,file=paste(type,dt,nt[[3]],nt[[2]],floor(nt[[1]]),".pdf",sep=""))
  return(list(obj=obj,summary=summary,coda=coda,type=type)) 
}

