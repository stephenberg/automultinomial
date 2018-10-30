library(automultinomial)
set.seed(33)

#10 predictors
p=5


#n times n grid
n=40

#make grid and adjacency matrix
latticeGraph=igraph::make_lattice(c(n,n))
A=igraph::get.adjacency(latticeGraph)

#set coefficient values 
beta=matrix(rnorm(p),ncol=1)*0.3
beta

#set covariate values
X=matrix(rnorm(n^2*p),ncol=p)

#set the correlation parameter value (0.7 is a moderate amount of spatial correlation)
gamma=0.7
y=drawSamples(beta,gamma,X,A,nSamples = 1)
y2=drawSamples(beta,0,X,A,nSamples = 1)
im1=Matrix::image(Matrix::Matrix(matrix(y,ncol=n)))
im2=Matrix::image(Matrix::Matrix(matrix(y2,ncol=n)))
library(gridExtra)
png(filename = "vignettes/plots/plotk2.png",width = 500,height = 300)
grid.arrange(im1, im2, ncol = 2)
t=dev.off()

# responses must be input as a factor
y=factor(y)

fit1=MPLE(X=X,y=y,A=A,ciLevel = 0.95,method = "asymptotic")
fit2=MPLE(X=X,y=y,A=A,ciLevel = 0.95,method = "boot",nBoot = 500)


save(fit1,fit2,y,y2,file="vignettes/data/2category.RData")



set.seed(42)

#10 predictors
p=5

#n times n grid
n=40

#make grid and adjacency matrix
latticeGraph=igraph::make_lattice(c(n,n))
A=igraph::get.adjacency(latticeGraph)

#set coefficient values 
#with 3 categories in the response, beta is now a matrix
beta=matrix(rnorm(p*2),ncol=2)*0.3
beta

#set covariate values
X=matrix(rnorm(n^2*p),ncol=p)

#set the correlation parameter value (0.7 is a moderate amount of spatial correlation)
gamma=0.7

y=drawSamples(beta,gamma,X,A,nSamples = 1)
y2=drawSamples(beta,0,X,A,nSamples = 1)

im1=Matrix::image(Matrix::Matrix(matrix(y,ncol=n)))
im2=Matrix::image(Matrix::Matrix(matrix(y2,ncol=n)))
library(gridExtra)
png(filename = "vignettes/plots/plotk3.png",width=500,height=300)
grid.arrange(im1, im2, ncol = 2)
t=dev.off()

y=factor(y)
fit1=MPLE(X=X,y=y,A=A,ciLevel = 0.95,method = "asymptotic")
fit2=MPLE(X=X,y=y,A=A,ciLevel = 0.95,method = "boot",nBoot = 500)

save(fit1,fit2,y,y2,file="vignettes/data/3category.RData")


set.seed(33)
#making a square lattice graph and adjacency matrix with toroidal boundary conditions
#every site has 4 neighbors
t1=igraph::make_lattice(c(40,40),circular=TRUE)
a1=igraph::get.adjacency(t1)


#making a square lattice graph and adjacency matrix with free boundary conditions
#sites have 2 neighbors at the corner of the lattice, 3 neighbors on edges of the lattice,
#and 4 neighbors internal to the lattice
t2=igraph::make_lattice(c(40,40),circular=FALSE)
a2=igraph::get.adjacency(t2)

#making the X matrices: X does not have an intercept, but X_intercept has an intercept
X=matrix(rnorm(1600*2),ncol=2)
X_intercept=cbind(rep(1,1600),X)

beta=rnorm(3)
beta=matrix(beta,ncol=1)
gamma=0.5

y=drawSamples(beta=beta,gamma=gamma,X = X_intercept, A=a1,nSamples = 1)
save(y,file="vignettes/data/paramExample.RData")

