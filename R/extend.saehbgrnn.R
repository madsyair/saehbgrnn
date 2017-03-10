extend.saehbgrnn<-function(resultin,adapt=4000,startsample=10000,thinSteps=20,DIC=FALSE)

{

  pkgs <- c('ggmcmc')
  lapply(pkgs, require, character.only = T)

  if (resultin$type=="hbgrnn1"){
    result<-extend.saehbgrnn1(resultin$obj,adapt=4000,startsample=10000,thinSteps=20,DIC=FALSE)
    result$type="hbgrnn1"
    }else if(resultin$type=="hb"){
      result<-extend.saehb(resultin$obj,adapt=4000,startsample=10000,thinSteps=20,DIC=FALSE)

    result$type="hb"
    }else if(resultin$type=="hbgrnn2"){
      result<-extend.saehbgrnn2(resultin$obj,adapt=4000,startsample=10000,thinSteps=20,DIC=FALSE)

    result$type="hbgrnn2"
    }
  else{stop("input is not complete")
  }
  ggso<-ggs(result$coda)
  t<-strftime(Sys.time(),format="%H%M%S")
  ggmcmc(ggso,file=paste("extendsaehbgrnn",t,".pdf",sep="_"))



  return(result)
}
