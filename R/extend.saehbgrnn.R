extend.saehbgrnn<-function(resultin,adapt=4000,startsample=10000,thinSteps=20,max.time=Inf,DIC=FALSE)
{
  objin<-resultin$obj
  pkgs <- c('ggmcmc','rjags',  'runjags','coda')
  lapply(pkgs, require, character.only = T)
  N<-objin$N
  
  if (resultin$type=="hbgrnn1"){
    #result<-extend.saehbgrnn1(resultin$obj,adapt=adapt,startsample=startsample,thinSteps=thinSteps,DIC=DIC)
    type="hbgrnn1"
    if(DIC){
      parameters = c( "CPOinv","w","mx","P"  , "v","h","deviance", "pd", "popt", "dic", "ped","full.pd")
    }
    else{
      parameters=c(  "CPOinv","w","mx","P","v","h")}
    
    }else if(resultin$type=="hb"){
    #result<-extend.saehb(resultin$obj,adapt=adapt,startsample=startsample,thinSteps=thinSteps,DIC=DIC)
    type="hb"
    if(DIC){
      parameters = c( "CPOinv", "P"  ,"beta0","beta","sigma","deviance", "pd", "popt", "dic", "ped","full.pd")
    }
    else{
      parameters=c( "CPOinv","P","beta0","beta","sigma")}
    }else if(resultin$type=="hbgrnn2"){
#      result<-extend.saehbgrnn2(resultin$obj,adapt=adapt,startsample=startsample,thinSteps=thinSteps,DIC=DIC)
    type="hbgrnn2"
    if(DIC){
      parameters = c(  "CPOinv","w","mx","P","beta0","beta"  , "v","h","deviance", "pd", "popt", "dic", "ped","full.pd")
    }
    else{
      parameters = c( "CPOinv","w","mx","P","beta0","beta", "v","h")}
    
    }
  else{stop("input is not complete")
  }
  N<-objin$N
  if(DIC){
  obj<-autoextend.jags(method=c("rjags","parallel")[1] ,
                       objin, 
                       add.monitor=parameters , 
                       adapt=adapt ,
                       startsample=startsample ,
                       thin=thinSteps,
                       max.time=max.time,
                       summarise=TRUE)
  }else{
    obj<-autoextend.jags(method=c("rjags","parallel")[2],
                         objin, 
                         add.monitor=parameters , 
                         adapt=adapt ,
                         startsample=startsample ,
                         thin=thinSteps,
                         max.time=max.time,
                         summarise=TRUE)
  }
  coda <- as.mcmc.list(obj)
  ggso<-ggs(coda)
  t<-Sys.time()
  nt <- as.POSIXlt(t)
  dt<-format(t, format="%y%m%d")
  ggmcmc(ggso,file=paste("extend",type,dt,nt[[3]],nt[[2]],floor(nt[[1]]),".pdf",sep=""))
  
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
  return(list(obj=obj,summary=summary,coda=coda,type=type))
}


 
