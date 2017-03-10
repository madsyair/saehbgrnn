sae.hbgrnn<-function(y=NULL,n=NULL,xl=NULL,xnl=NULL,M=5,adapt=4000,startburnin=1000,nChains = 2,startsample=10000,thinSteps=20,DIC=FALSE)
{

  pkgs <- c('ggmcmc')
  lapply(pkgs, require, character.only = T)
  nulxl=is.null(xl)
  nulxnl=is.null(xnl)
  if (nulxl==TRUE&&nulxnl==FALSE){
    result<-sae.hbgrnn1(y=y,n=n,x=xnl,M=5,adapt=4000,startburnin=1000,nChains = 2,startsample=10000,thinSteps=20,DIC=FALSE)
  }else if(nulxl==FALSE&&nulxnl==TRUE){
    result<-sae.hb(y=y,n=n,x=xl,adapt=4000,startburnin=1000,nChains = 2,startsample=10000,thinSteps=20,DIC=FALSE)
  }else if(nulxl==FALSE&&nulxnl==FALSE){
    result<-sae.hbgrnn2(y=y,n=n,xl=xl,xnl=xnl,M=5,adapt=4000,startburnin=1000,nChains = 2,startsample=10000,thinSteps=20,DIC=FALSE)
  }
  else{stop("input is not complete")
  }
  ggso<-ggs(result$coda)
  t<-strftime(Sys.time(),format="%H%M%S")
  ggmcmc(ggso,file=paste("saehbgrnn",t,".pdf",sep="_"))

  for (i in 1:38){varp[i]<-paste("P","[",i,"]",sep="")}

  return(list(result=result))
}
