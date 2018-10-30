#' Simulate data from auto- models
#'
#' Generates data from the autologistic and automultinomial
#' models via Gibbs sampling. See the vignette for an example of use.
#'
#'@param beta coefficient vector (for the autologistic model) or matrix (for the automultinomial model)
#'@param gamma the value of the autocorrelation parameter
#'@param X the design matrix
#'@param A the (square symmetric) adjacency matrix encoding the neighborhood structure
#'@param burnIn the number of burnin samples to be used. Defaults to 300
#'@param nSamples the number of samples to draw
#'@param y optional starting configuration, in factor form. Defaults to NULL
#'@return simulated samples
#'@import stats
#'@examples 
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
#' #2 response categories (1 column in coefficient matrix)
#' beta2=matrix(rnorm(5)*0.3,ncol=1)
#' #This example uses a short burnIn period. Use a longer burnIn in practice.
#' y2=drawSamples(beta2,gamma,X,A,burnIn=1,nSamples=1)
#' 
#' #3 response categories (2 columns in coefficient matrix)
#' beta3=matrix(rnorm(10)*0.3,ncol=2)
#' y3=drawSamples(beta3,gamma,X,A,burnIn=1,nSamples=1)
#' ##########
#'@export
drawSamples<-function(beta,gamma,X,A,burnIn=300,nSamples,y=NULL){
  if ((!is.double(gamma))|(length(gamma)>1)){
    stop(cat("Error: correlation parameter must be a real number\n"))
  }
  if ((burnIn<1)){
    stop(cat("Error: burnIn must be a positive integer\n"))
  }
  if ((nSamples<1) | (length(nSamples)>1)){
    stop(cat("Error: nSamples must be a positive integer\n"))
  }
  if (!is.matrix(X)){
    stop(cat("Error: X is not a matrix. Please input a matrix for X.\n"))
  }
  if ((!is.matrix(A)) & !(class(A)[1]=="dgCMatrix")){
    stop(cat("Error: A is not a matrix. The input A should be an nxn symmetric adjacency matrix, or matrix of type dgCMatrix.\n"))
    if (nrow(A)!=ncol(A)){
      stop(cat("Error: number of rows in A not equal to number of columns in A.\n"))
    }
  }
  if (!is.null(y)){
    if (!is.factor(y)){
      stop(cat("Error: if supplied, y should be a factor"))
    }
    if (length(y)!=dim(A)[1]){
      stop(cat("Error: y and A dimensions disagree"))
    }
  }
  if (!is.matrix(beta)){
    stop(cat("Error: beta must be input as a matrix"))
  }
  
  if (ncol(X)!=nrow(beta)){
    stop(cat("Error: X and beta dimensions disagree"))
  }
  
  
  #test for bad A dimensions
  n=nrow(X)
  if (n!=nrow(A)){
    stop(cat("Error: number of rows of A not equal to number rows in X matrix.\n"))
  }
  #test for zero diagonals
  if (max(abs(Matrix::diag(A)))>10^{-9}){
    stop(cat("Error: diagonal entries of A should be 0.\n"))
  }
  
  #symmetry check
  if (!Matrix::isSymmetric(A,tol=10^{-9})){
    stop(cat("Error: A should be a symmetric matrix\n"))
  }
  if (!is.vector(beta)){
    if (!is.matrix(beta)){
      stop(cat("Error: beta should be a matrix or vector\n"))
    }
  }
  beta=as.matrix(beta)
  if (nrow(beta)!=ncol(X)){
    stop(cat("Error: dimensions of beta and design matrix X disagree\n"))
  }
  
  #no weighting schemes
  A=1.0*(A>0)

  
  
  
  p=dim(X)[2]
  n=dim(X)[1]
  beta=as.matrix(beta)
  K=dim(as.matrix(beta))[2]+1
  z=matrix(0,n,K)
  
  #random initialization
  zCategories=t(rmultinom(n,1,prob=rep(1/K,K)))
  zCategories=apply(zCategories,1,which.max)
  for (i in 1:n){
    z[i,zCategories[i]]=1
  }
  
  if (!is.null(y)){
    z=z*0
    y=as.numeric(y)
    for (i in 1:n){
      z[i,y[i]]=1
    }
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
  
  sampleMatrix=matrix(0,n,nSamples)
  
  ##############
  #now draw samples
  cat("Drawing samples\n")
  for (i in 1:nSamples){
    for (j in 1:n){
      neighborCount_j=rep(0,K)
      if (nNeighbors[j]>0){
        j_neighbors=neighbors[neighborStart[j]:neighborEnd[j]]
        for (q in 1:length(j_neighbors)){
          if (j_neighbors[q]!=j){
            neighborCount_j=neighborCount_j+z[j_neighbors[q],]#
          }
        }
      }
      cProbs_j=linPred[j,]+gamma*neighborCount_j
      cProbs_j=cProbs_j-max(cProbs_j)
      cProbs_j=exp(cProbs_j)
      cProbs_j=cProbs_j/sum(cProbs_j)
      z[j,]=rmultinom(1,1,cProbs_j)
    }
    sampleMatrix[,i]=apply(z,1,which.max)
    ##status printout
    if (i%%100==0){
      cat(paste(i," samples so far\n"))
    }
  }
  return(sampleMatrix)
}
