# library(igraph)
# p=10
# K=3
# n=40
# adjMat=make_lattice(c(n,n))
# adjMat=get.adjacency(adjMat)
# beta=matrix(rnorm(p*(K-1)),ncol=K-1)*0.3
# X=matrix(rnorm(n^2*10),ncol=10)
# probs=cbind(rep(1,n^2),exp(X%*%beta))
# probs=probs/apply(probs,1,sum)
#
# z=matrix(0,n^2,K)
# for (i in 1:(n^2)){
#   z[i,]=rmultinom(1,1,probs[i,])
# }
#
#
#
# library(numDeriv)
# y=drawSamples(beta=beta,gamma=0.5,X=X,A=A,burnIn = 500,nSamples = 1)
# y=factor(c(y))
#
# fit1=MPLE(X = X,y =y,A = adjMat,ciLevel = 0.95,method = "asymptotic")
# fit2=MPLE(X = X,y =y,A = adjMat,ciLevel = 0.95,method = "boot",burnIn=300,nBoot=500)
#
