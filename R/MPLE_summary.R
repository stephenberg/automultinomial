#' Summarize MPLE fits
#'
#' Prints out summary tables of fitted model objects from MPLE. Also
#' returns knitr::kable() summary tables.
#'@param fit a fitted MPLE object
#'@return tables based on the model fit
#'@import stats
#'@examples
#'
#' ##########generating coefficient values and data
#' #adjacency matrix A
#' A=igraph::get.adjacency(igraph::make_lattice(c(40,40)))
#' 
#' #design matrix
#' X=cbind(rep(1,1600),matrix(rnorm(1600*4),ncol=4))
#' 
#' #correlation parameter
#' gamma=0.6
#' 
#' #2 response categories
#' beta2=rnorm(5)*0.3
#' 
#' #This example uses a short burnIn period. Use a longer burnIn in practice.
#' y2=drawSamples(beta2,gamma,X,A,burnIn=100,nSamples=1)
#' 
#' #3 response categories (not run)
#' #beta3=rnorm(5)*0.3
#' #y3=drawSamples(beta3,gamma,X,A,nSamples=1)
#' ##########
#' 
#' ##########fitting models
#' fit2=MPLE(X = X,y=factor(y2),A = A,ciLevel = 0.99,method = "asymptotic")
#' #fit3=MPLE(X = X,y=factor(y3),A = A,ciLevel = 0.99,method = "asymptotic")
#' ##########
#' 
#' ##########summary tables
#' t2=MPLE_summary(fit2)
#' #t3=MPLE_summary(fit3)
#' ##########
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
