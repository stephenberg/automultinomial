#' Maximum pseudolikelihood estimation
#'
#' Fits an autologistic model or automultinomial model. Takes as arguments
#' a design matrix X, a response vector y (in factor form), and a square
#' symmetric adjacency matrix encoding the neighborhood structure.
#' When the number of levels of the response y is >2, the function fits a
#' multicategory generalization of the autologistic model. For a full description
#' of the models the package fits and a user guide, please see the vignette.
#'
#'@param X the n-by-p design matrix
#'@param y the response vector (required to be a factor)
#'@param A the square symmetric adjacency matrix A encoding the neighborhood structure
#'@param ciLevel the confidence level to be used for inference. Defaults to 0.95 for 95 percent intervals.
#'@param method "boot" for parametric bootstrap and "asymptotic" for asymptotic confidence intervals.
#'@param burnIn the number of burnin samples to use for the Gibbs sampler when method="boot"
#'@param nBoot the number of bootstrap samples to use when method="boot"
#'@return a fitted auto- model MPLE object
#'@import stats
#'@examples 
#' ##########generating coefficient values and data
#' A=igraph::get.adjacency(igraph::make_lattice(c(40,40))) #adjacency matrix A
#' X=cbind(rep(1,1600),matrix(rnorm(1600*4),ncol=4)) #design matrix
#' gamma=0.6 #correlation parameter
#' beta=matrix(rnorm(5)*0.3,ncol=1) #covariate parameters
#' y=drawSamples(beta,gamma,X,A,burnIn=10,nSamples=1)
#' 
#' ##########fitting model
#' fit=MPLE(X = X,y=factor(y),A = A,ciLevel = 0.99,method = "asymptotic")
#' 
#'@export
MPLE<-function(X,y,A,ciLevel=0.95,method="asymptotic",burnIn=300,nBoot=500){

  if ((ciLevel<=0) |(ciLevel>=1)){
    stop(cat("Error: ciLevel must satisfy 0<ciLevel<1\n"))
  }

  if ((burnIn<1)){
    stop(cat("Error: burnIn must be a positive integer\n"))
  }
  if ((nBoot<1)){
    stop(cat("Error: nBoot must be a positive integer\n"))
  }
  if (method!="boot"){
    if (method!="asymptotic"){
      stop(cat("Error: method should be one of \"boot\" or \"asymptotic\", for bootstrap\n or asymptotic inference, respectively.\n"))
    }
  }
  ####argument type checking
  if (!is.matrix(X)){
    stop(cat("Error: X is not a matrix. Please input a matrix for X.\n"))
  }
  if (!is.factor(y)){
    stop(cat("Error: y is not a factor. Please input y as a factor vector.\n"))
  }
  if ((!is.matrix(A)) & !(class(A)[1]=="dgCMatrix")){
    stop(cat("Error: A is not a matrix. The input A should be an nxn symmetric adjacency matrix, or matrix of type dgCMatrix.\n"))
    if (nrow(A)!=ncol(A)){
      stop(cat("Error: number of rows in A not equal to number of columns in A.\n"))
    }
  }
  ####


  #problem dimensions
  n=dim(X)[1]
  p=dim(X)[2]
  K=length(levels(y))
  parameterNames=rep("",p)

  categoryNames=levels(y)
  parameterNames=colnames(X)
  if (is.null(colnames(X))){
    parameterNames=as.character(1:p)
  }

  if (K<2){
    stop(cat("Error: need more than one response category."))
  }


  #test for bad A dimensions
  if (n!=nrow(A)){
    stop(cat("Error: number of rows of A not equal to number rows in X matrix.\n"))
  }

  #test for bad response length
  if (length(y)!=n){
    stop(cat("Error: number of responses should match number of rows in the X matrix"))
  }

  #test for zero diagonals
  if (max(abs(Matrix::diag(A)))>10^{-9}){
    stop(cat("Error: diagonal entries of A should be 0.\n"))
  }

  #symmetry check
  if (!Matrix::isSymmetric(A,tol=10^{-9})){
    stop(cat("Error: A should be a symmetric matrix\n"))
  }

  #test for bad input response vector
  if (min(table(y))<1){
    zeroObservations=which(table(y)<1)
    zeroLevels=paste(levels(y)[which(table(y)<1)],collapse=", ")
    stop(cat(paste("Error: factor level",zeroLevels,"has/have 0 observations")))
  }


  #no weighting schemes
  A=1.0*(A>0)

  #create beta matrix
  beta=matrix(0,p,K-1)
  #correlation parameter
  gamma=0

  #create response indicator matrix
  z=matrix(0,n,K)
  yNumeric=as.numeric(y)
  for (i in 1:n){
    z[i,yNumeric[i]]=1
  }

  #create neighborhood response counts
  neighborResponses=as.matrix(A%*%z)

  #########parameter estimation
  cat("Starting model fitting\n")
  betaGamma=rep(0,p*(K-1)+1)
  result=optim(par=betaGamma,fn=logPseudolikelihood,gr=logPseudolikelihoodGradient,yNumeric,X,neighborResponses,method="BFGS",control=list(fnscale=-1))
  betaGamma=result$par
  if (result$convergence!=0){
    stop(cat("Optimization failure. Potential fixes include simpler models with fewer covariates and response categories.\n"))
  }
  cat("Model fitting done, starting variance estimation\n")
  #########

  #conditional probabilities
  conditionalProbabilities=condProbs(betaGamma,yNumeric,X,neighborResponses)

  #########parameter inference
  hessian=numDeriv::jacobian(logPseudolikelihoodGradient,betaGamma,method="Richardson",side=NULL,method.args=list(),yNumeric,X,neighborResponses)
  scoreMatrix=pointWiseGradient(betaGamma,yNumeric,X,neighborResponses)
  variance=varianceComputer(hessian,scoreMatrix,A)


  #formatting and returning results
  betaHat=matrix(betaGamma[1:(p*(K-1))],ncol=K-1)
  rownames(betaHat)<-parameterNames
  betaHatColnames=rep("",K-1)
  varianceNames=rep("",p*(K-1)+1)
  for (k in 2:K){
    betaHatColnames[k-1]=paste(categoryNames[k]," vs. ", categoryNames[1],sep="")
    varianceNames[((k-2)*p+1):((k-1)*p)]=paste(parameterNames,paste("_",k,"1",sep=""),sep="")
  }
  variance=as.matrix(variance)
  varianceNames[p*(K-1)+1]="gamma"
  colnames(betaHat)<-betaHatColnames

  rownames(variance)<-varianceNames
  colnames(variance)<-varianceNames

  gammaHat=betaGamma[length(betaGamma)]
  

  cat(paste("Creating",method,"confidence intervals\n",sep=" "))
  if (method=="asymptotic"){

    ciMatrix=matrix(0,length(betaGamma),2)
    zScores=betaGamma%*%solve(diag(sqrt(diag(variance))))
    pValues=2*pnorm(-abs(zScores),lower.tail=TRUE)
    names(zScores)=varianceNames
    names(pValues)=varianceNames

    ciQuantile=qnorm(1-(1-ciLevel)/2)
    ciMatrix[,1]=betaGamma-ciQuantile*sqrt(diag(variance))
    ciMatrix[,2]=betaGamma+ciQuantile*sqrt(diag(variance))
    colnames(ciMatrix)=c("lower","upper")
    rownames(ciMatrix)=varianceNames
    return(list(betaHat=betaHat,gammaHat=gammaHat,zScores=zScores,pValues=pValues,ciMatrix=ciMatrix,variance=variance,pseudolikelihood=result$value,conditionalProbabilities=conditionalProbabilities))
  }

  #if bootstrap
  #bootStrap(betaGammaVector,yNumeric,X,neighborResponses,burnIn,nBoot)
  betaGammaMatrix=bootStrap(betaGamma,yNumeric,X,neighborResponses,burnIn,nBoot,A)
  ciMatrix=matrix(0,length(betaGamma),2)
  alpha=1-ciLevel
  ciMatrix[,1]=apply(betaGammaMatrix,1,quantile,alpha/2)
  ciMatrix[,2]=apply(betaGammaMatrix,1,quantile,1-alpha/2)
  pValues=rep(0,p*(K-1)+1)
  for (i in 1:length(pValues)){
    if (betaGamma[i]>0){
      pValues[i]=1-mean(((betaGammaMatrix[i,]-betaGamma[i])<betaGamma[i])& ((betaGammaMatrix[i,]-betaGamma[i])>(-betaGamma[i])))
    }
    if (betaGamma[i]<0){
      pValues[i]=1-mean(((betaGammaMatrix[i,]-betaGamma[i])>betaGamma[i])& ((betaGammaMatrix[i,]-betaGamma[i])<(-betaGamma[i])))
    }
  }
  colnames(betaGammaMatrix)=as.character(1:dim(betaGammaMatrix)[2])
  rownames(betaGammaMatrix)=varianceNames
  names(pValues)=varianceNames
  ciQuantile=qnorm(1-(1-ciLevel)/2)
  colnames(ciMatrix)=c("lower","upper")
  rownames(ciMatrix)=varianceNames
  return(list(betaHat=betaHat,gammaHat=gammaHat,pValues=pValues,ciMatrix=ciMatrix,variance=variance,bootStrapSamples=betaGammaMatrix,pseudolikelihood=result$value,conditionalProbabilities=conditionalProbabilities))
}



varianceComputer<-function(hessian,scoreMatrix,A){
  A_Diagonal=Matrix::Diagonal(dim(A)[1])+A
  J=t(scoreMatrix)%*%A_Diagonal%*%scoreMatrix
  I=-hessian
  variance=solve(I)
  variance=variance%*%J%*%variance
  return(variance)
}

condProbs<-function(betaGammaVector,yNumeric,X,neighborResponses){
  #problem dimensions
  K=dim(neighborResponses)[2]
  n=dim(neighborResponses)[1]
  
  betaVector=betaGammaVector[1:(length(betaGammaVector)-1)]
  gamma=betaGammaVector[length(betaGammaVector)]
  beta=matrix(betaVector,ncol=K-1)
  
  condProbs=X%*%beta
  condProbs=condProbs+gamma*(neighborResponses[,2:K,drop=FALSE]-neighborResponses[,1])
  condProbs=cbind(rep(0,n),condProbs)
  #condProbs=condProbs+gamma*neighborResponses[,2:K,drop=FALSE]
  #condProbs=cbind(gamma*neighborResponses[,1],condProbs)
  for (i in 1:n){
    condProbs[i,]=condProbs[i,]-max(condProbs[i,])
    condProbs[i,]=exp(condProbs[i,])
    condProbs[i,]=condProbs[i,]/sum(condProbs[i,])
  }
  return(condProbs)
}

logPseudolikelihood<-function(betaGammaVector,yNumeric,X,neighborResponses){

  #problem dimensions
  K=dim(neighborResponses)[2]
  n=dim(neighborResponses)[1]

  betaVector=betaGammaVector[1:(length(betaGammaVector)-1)]
  gamma=betaGammaVector[length(betaGammaVector)]
  beta=matrix(betaVector,ncol=K-1)

  condProbs=X%*%beta
  condProbs=condProbs+gamma*(neighborResponses[,2:K,drop=FALSE]-neighborResponses[,1])
  condProbs=cbind(rep(0,n),condProbs)
  #condProbs=condProbs+gamma*neighborResponses[,2:K,drop=FALSE]
  #condProbs=cbind(gamma*neighborResponses[,1],condProbs)
  for (i in 1:n){
    condProbs[i,]=condProbs[i,]-max(condProbs[i,])
    condProbs[i,]=exp(condProbs[i,])
    condProbs[i,]=condProbs[i,]/sum(condProbs[i,])
  }

  logPLike=0
  for (i in 1:length(yNumeric)){
    responseCategory=yNumeric[i]
    logPLike=logPLike+log(condProbs[i,responseCategory])
  }
  return(logPLike)
}


logPseudolikelihoodGradient<-function(betaGammaVector,yNumeric,X,neighborResponses){
  #problem dimensions
  K=dim(neighborResponses)[2]
  n=dim(neighborResponses)[1]

  betaVector=betaGammaVector[1:(length(betaGammaVector)-1)]
  gamma=betaGammaVector[length(betaGammaVector)]
  beta=matrix(betaVector,ncol=K-1)

  condProbs=X%*%beta
  condProbs=condProbs+gamma*(neighborResponses[,2:K,drop=FALSE]-neighborResponses[,1])
  condProbs=cbind(rep(0,n),condProbs)
  for (i in 1:n){
    condProbs[i,]=condProbs[i,]-max(condProbs[i,])
    condProbs[i,]=exp(condProbs[i,])
    condProbs[i,]=condProbs[i,]/sum(condProbs[i,])
  }

  naturalGradient=-condProbs[,2:K,drop=FALSE]
  for (i in 1:n){
    responseCategory=yNumeric[i]
    if (responseCategory>1){
      naturalGradient[i,responseCategory-1]=naturalGradient[i,responseCategory-1]+1
    }
  }

  betaGradient=t(X)%*%naturalGradient
  neighborResponseDifference=neighborResponses[,2:K,drop=FALSE]-neighborResponses[,1]
  gammaGradient=sum(neighborResponseDifference*naturalGradient)
  fullGradient=c(betaGradient,gammaGradient)
  return(fullGradient)
}

pointWiseGradient<-function(betaGammaVector,yNumeric,X,neighborResponses){
  #problem dimensions
  K=dim(neighborResponses)[2]
  n=dim(neighborResponses)[1]
  p=dim(X)[2]

  betaVector=betaGammaVector[1:(length(betaGammaVector)-1)]
  gamma=betaGammaVector[length(betaGammaVector)]
  beta=matrix(betaVector,ncol=K-1)

  condProbs=X%*%beta
  condProbs=condProbs+gamma*(neighborResponses[,2:K,drop=FALSE]-neighborResponses[,1])
  condProbs=cbind(rep(0,n),condProbs)
  for (i in 1:n){
    condProbs[i,]=condProbs[i,]-max(condProbs[i,])
    condProbs[i,]=exp(condProbs[i,])
    condProbs[i,]=condProbs[i,]/sum(condProbs[i,])
  }

  naturalGradient=-condProbs[,2:K,drop=FALSE]
  for (i in 1:n){
    responseCategory=yNumeric[i]
    if (responseCategory>1){
      naturalGradient[i,responseCategory-1]=naturalGradient[i,responseCategory-1]+1
    }
  }

  siteGradients=matrix(0,n,(K-1)*p+1)
  for (i in 1:n){
    siteGradients[i,1:((K-1)*p)]=c(c(X[i,])%*%t(naturalGradient[i,]))
    siteGradients[i,(K-1)*p+1]=sum(naturalGradient[i,]*(neighborResponses[i,2:K]-neighborResponses[i,1]))
  }
  return(siteGradients)
}


bootStrap<-function(betaGammaVector,yNumeric,X,neighborResponses,burnIn=300,nBoot,A){
  K=dim(neighborResponses)[2]
  n=dim(neighborResponses)[1]
  p=dim(X)[2]
  betaVector=betaGammaVector[1:(length(betaGammaVector)-1)]
  gamma=betaGammaVector[length(betaGammaVector)]
  beta=matrix(betaVector,ncol=K-1)
  z=matrix(0,n,K)
  for (i in 1:n){
    z[i,yNumeric[i]]=1
  }

  #first, do Gibbs sampling for burn-in iterations, starting from initial configuration
  #yNumeric
  A=A+Matrix::diag(1,dim(A)[1],dim(A)[2])
  indices=Matrix::summary(Matrix::Matrix(A>0,sparse=TRUE))
  A=A-Matrix::diag(1,dim(A)[1],dim(A)[2])

  nNeighbors=rep(0,n)
  for (i in 1:n){
    if (!(i%in%indices[,2])){
      nNeighbors[i]=0
    }
    if ((i%in%indices[,2])){
      nNeighbors[i]=sum(indices[,2]==i)
    }
  }
  neighbors=indices[,1]
  neighborStart=c(0,cumsum(nNeighbors)[1:(n-1)])+1
  neighborEnd=cumsum(nNeighbors)
  linPred=X%*%beta
  linPred=cbind(rep(0,n),linPred)
  cat("Burn-in samples\n")
  for (i in 1:burnIn){
    for (j in 1:n){
      neighborCount_j=rep(0,K)
      if (nNeighbors[j]>0){
        j_neighbors=neighbors[neighborStart[j]:neighborEnd[j]]
        for (q in 1:length(j_neighbors)){
          if (j_neighbors[q]!=j){
            neighborCount_j=neighborCount_j+z[j_neighbors[q],]
          }
        }
      }
      cProbs_j=linPred[j,]+gamma*neighborCount_j
      cProbs_j=cProbs_j-max(cProbs_j)
      cProbs_j=exp(cProbs_j)
      cProbs_j=cProbs_j/sum(cProbs_j)
      z[j,]=rmultinom(1,1,cProbs_j)
    }
    if (i%%100==0){
      cat(paste(i," burn-in samples so far\n"))
    }
  }


  ##############
  #now do bootstrap
  #create bootstrap coefficient matrix
  betaGammaMatrix=matrix(0,p*(K-1)+1,nBoot)
  cat("Starting bootstrap sampling and fitting\n")
  for (i in 1:nBoot){
    for (j in 1:n){
      neighborCount_j=rep(0,K)
      if (nNeighbors[j]>0){
        j_neighbors=neighbors[neighborStart[j]:neighborEnd[j]]
        for (q in 1:length(j_neighbors)){
          if (j_neighbors[q]!=j){
            neighborCount_j=neighborCount_j+z[j_neighbors[q],]
          }
        }
      }
      cProbs_j=linPred[j,]+gamma*neighborCount_j
      cProbs_j=cProbs_j-max(cProbs_j)
      cProbs_j=exp(cProbs_j)
      cProbs_j=cProbs_j/sum(cProbs_j)
      z[j,]=rmultinom(1,1,cProbs_j)
    }

    #######
    #bootstrap coefficient estimation
    yNumeric_i=apply(z,1,which.max)
    neighborResponses=as.matrix(A%*%z)
    result=optim(par=betaGammaVector,fn=logPseudolikelihood,gr=logPseudolikelihoodGradient,yNumeric_i,X,neighborResponses,method="BFGS",control=list(fnscale=-1))
    betaGammaMatrix[,i]=result$par
    if (result$convergence!=0){
      cat("Warning: bootstrap optimization failure. Confidence intervals may be unreliable.\n")
    }
    #######

    ##status printout
    if (i%%100==0){
      cat(paste(i," bootstrap samples so far\n"))
    }
  }
  return(betaGammaMatrix)
}
