#' Summarize MPLE fits
#'
#' Prints out summary tables of fitted model objects from MPLE. Also
#' returns knitr::kable() summary tables.
#'@param fit a fitted MPLE object
#'@return tables based on the model fit
#'@import stats
#'@examples
#'
#' ##########generating model fit to summarize
#' #adjacency matrix A
#' A=igraph::get.adjacency(igraph::make_lattice(c(40,40)))
#' X=cbind(rep(1,1600),matrix(rnorm(1600*4),ncol=4))
#' gamma=0.6
#' beta=matrix(rnorm(5)*0.3,ncol=1)
#' y=drawSamples(beta,gamma,X,A,burnIn=10,nSamples=1)
#' fit=MPLE(X = X,y=factor(y),A = A,ciLevel = 0.99,method = "asymptotic")
#' ##########
#' 
#' ##########summarizing model fit
#' MPLE_summary(fit)
#' 
#'@export
MPLE_summary<-function(fit){

  bhat=round(fit$betaHat,digits=3)
  ciMat=round(fit$ciMatrix,digits=3)
  ciTable=bhat

  bhat=as.matrix(bhat)
  count=1
  for (q in 1:dim(bhat)[2]){
    for (p in 1:dim(bhat)[1]){

      bhat[p,q]=paste(bhat[p,q],paste("(",ciMat[count,1],",",ciMat[count,2],")",sep=""),sep=" ")
      count=count+1
    }
  }
  gammaVec=paste(round(fit$gammaHat,digits=3),paste("(",ciMat[count,1],",",ciMat[count,2],")",sep=""))
  gammaVec=c(gammaVec,rep("",dim(bhat)[1]-1))
  bhat=cbind(bhat,gammaVec)
  colnames(bhat)[dim(bhat)[2]]<-"gamma"
  print(knitr::kable(bhat,caption="Summary with confidence intervals",row.names=TRUE))
  ciTable=knitr::kable(bhat,caption="Summary with confidence intervals",row.names=TRUE)


  bhat=round(fit$betaHat,digits=3)
  bhat=as.matrix(bhat)

  count=1
  for (q in 1:dim(bhat)[2]){
    for (p in 1:dim(bhat)[1]){
      bhat[p,q]=paste(bhat[p,q],paste("(",round(fit$pValues[count],digits=3),")",sep=""),sep=" ")
      count=count+1
    }
  }
  gammaVec=paste(round(fit$gammaHat,digits=3),paste("(",round(fit$pValues[count],digits=3),")",sep=""))
  gammaVec=c(gammaVec,rep("",dim(bhat)[1]-1))
  bhat=cbind(bhat,gammaVec)
  colnames(bhat)[dim(bhat)[2]]<-"gamma"
  print(knitr::kable(bhat,caption="Summary with p-values",row.names=TRUE))
  pValueTable=knitr::kable(bhat,caption = "Summary with p-values",row.names=TRUE)
  return(list(ciTable=ciTable,pValueTable=pValueTable))
}
