extend.saehbgrnn<-function(resultin,adapt=4000,startsample=10000,thinSteps=20,DIC=FALSE)

{

  pkgs <- c('ggmcmc')
  lapply(pkgs, require, character.only = T)

  if (resultin$type=="hbgrnn1"){
    result<-extend.saehbgrnn1(resultin$obj,adapt=adapt,startsample=startsample,thinSteps=thinSteps,DIC=DIC)
    result$type="hbgrnn1"
    }else if(resultin$type=="hb"){
      result<-extend.saehb(resultin$obj,adapt=adapt,startsample=startsample,thinSteps=thinSteps,DIC=DIC)

    result$type="hb"
    }else if(resultin$type=="hbgrnn2"){
      result<-extend.saehbgrnn2(resultin$obj,adapt=adapt,startsample=startsample,thinSteps=thinSteps,DIC=DIC)

    result$type="hbgrnn2"
    }
  else{stop("input is not complete")
  }
  ggso<-ggs(result$coda)
  t<-strftime(Sys.time(),format="%H%M%S")
  ggmcmc(ggso,file=paste("extendsaehbgrnn",t,".pdf",sep="_"))
  return(result)

}
  extend.saehb<-function(objin,adapt=4000,startsample=10000,thinSteps=20,DIC=FALSE)
  {
    
    if(DIC){
      parameters = c( "CPOinv", "P"  ,"beta0","beta","sigma","deviance", "pd", "popt", "dic", "ped","full.pd")
    }
    else{
      parameters=c( "CPOinv","P","beta0","beta","sigma")}
    
    pkgs <- c('rjags',  'runjags','coda')
    lapply(pkgs, require, character.only = T)
    obj<-autoextend.jags(objin, 
                         add.monitor=parameters , 
                         adapt=adapt ,
                         startsample=startsample ,
                         thin=thinSteps,
                         summarise=TRUE )
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
  
  extend.saehbgrnn1<-function(objin,adapt=4000,startsample=10000,thinSteps=20,DIC=FALSE)
  {
    
    
    if(DIC){
      parameters = c( "CPOinv","w","mx","P"  , "v","h","deviance", "pd", "popt", "dic", "ped","full.pd")
    }
    else{
      parameters=c(  "CPOinv","w","mx","P","v","h")}
    
    pkgs <- c('rjags',  'runjags','coda')
    lapply(pkgs, require, character.only = T)
    
    
    obj<-autoextend.jags(objin, 
                         add.monitor=parameters , 
                         adapt=adapt ,
                         startsample=startsample ,
                         thin=thinSteps,
                         summarise=TRUE)
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
  extend.saehbgrnn2<-function(objin,adapt=4000,startsample=10000,thinSteps=20,DIC=FALSE)
  {
    
    if(DIC){
      parameters = c(  "CPOinv","w","mx","P","beta0","beta"  , "v","h","deviance", "pd", "popt", "dic", "ped","full.pd")
    }
    else{
      parameters = c( "CPOinv","w","mx","P","beta0","beta", "v","h")}
    
    pkgs <- c('rjags',  'runjags','coda')
    lapply(pkgs, require, character.only = T)
    
    obj<-autoextend.jags(objin, 
                         add.monitor=parameters , 
                         adapt=adapt ,
                         startsample=startsample ,
                         thin=thinSteps,
                         summarise=TRUE )
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


 
