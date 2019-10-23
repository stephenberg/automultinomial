---
title: "Automultinomial Vignette"
author: "Stephen Berg"
date: "2019-10-23"
output: 
  html_document:
    keep_md: true
vignette: >
  %\VignetteIndexEntry{Automultinomial Vignette}
  %\VignetteEngine{knitr::rmarkdown}
  header-includes:
   - \usepackage{amsmath}
   - \usepackage[utf8]{inputenc}
---



This vignette explains the installation and use of the R package **automultinomial**. The **automultinomial** package is designed to be used for regressions similar to logistic or multinomial logit regression. However, unlike ordinary logistic and multinomial logit models, the autologistic/automultinomial model includes an autocorrelation parameter to account for spatial dependence between observations.

The organization of this document is:

1. Description of the problem **automultinomial** solves
2. cran and Github installation how-to
3. A data example with binary response (2 response categories)
4. A data example with 3 response categories
5. Some comments on autologistic model parameterization

\section{The problem **automultinomial** solves}

Consider a data problem where covariates (the independent variables) and categorical outcomes (the dependent variables) are observed on a spatial grid or lattice. If the outcomes are spatially autocorrelated, then our coefficient estimates and inference procedures ought to take this into account. Intuitively, when outcomes are correlated, the effective sample size of the dataset is smaller than if the outcomes were independent. The statistical information we can gain based on samples from multiple nearby sites may be less than if we were able to sample sites that are ``more independent''. 

A particular practical risk in neglecting spatial correlation is that confidence intervals may be too narrow. We might find ``statistically significant'' relationships that are in fact byproducts of noise or dependent sampling. The **automultinomial** package provides a tool for dealing with outcomes with spatial autocorrelation. It still might happen that nearby sites in a dataset are effectively independent. But we can only assess this using a method that accounts for the spatial arrangement of the samples. 

Having given a qualitative description of the data problem, we will now move on to a discussion of the probability model used by **automultinomial**.

### Mathematical description

To describe the setup in **automultinomial**, we will use the notations

* $i=1,...,n$ to index the $n$ observation sites on the spatial grid
* $K\geq 2$ to denote the number of possible values taken by the response
* $p$ to denote the number of covariates observed at each site
* $x_i\in\mathbb{R}^p$ to denote the values of the covariates at site $i$
* $z_i=1,...,K$ to denote the response at site $i$
  + $\mathbf{z}$ to denote the vector of categorical responses for all the sites
  + $z_{-i}$ to denote the entire response vector $\mathbf{z}$, *except* for site $i$
* $\boldsymbol{\beta}\in\mathbb{R}^{p\times (K-1)}$ to denote the coefficient matrix (when $K=2$, $\boldsymbol{\beta}$ is a vector of length $p$, as in logistic regression)
  + $\beta_k$ to denote the $k$-th column of $\boldsymbol{\beta}$.
* $\gamma$ to denote the autocorrelation parameter

We will assume that for each site $i$, a collection $N_i$ of neighboring sites is known. If the data are collected in a square lattice fashion, then a natural neighborhood setup is to set $N_i$ to be the 4 Up, Down, Left, and Right neighbors of site $i$ (sites at the boundary of a lattice may have less than 4 neighbors). We will assume $i\notin N_{i}$ for each $i$, so that site $i$ is not a neighbor of itself. It is allowable for sites to have no neighbors. In this case $N_i=\phi$, where $\phi$ denotes the empty set. We will use the notation $i'\sim i$ to indicate that site $i'$ is a neighbor of site $i$.

In the **automultinomial** package, we will require that the neighborhoods are *symmetric*: if $i'\in N_i$, then necessarily $i\in N_{i'}$ as well. 

We will also define an *energy* function $H(\mathbf{z}|\boldsymbol{\beta},\gamma)$. When $K=2$ and the two categories are $k=1$ and $k=2$, then
$\boldsymbol{\beta}$ is just a vector $\beta$ with length $p$, and \begin{eqnarray}
H(\mathbf{z}|\boldsymbol{\beta},\gamma)=\sum_{i=1}^{n}x_i^T\beta I(z_i=2)+\gamma \sum_{i=1}^{n}\sum_{\substack{i'\sim i\\i'>i}}\sum_{k=1}^{2}I(z_i=z_{i'}=k)\label{H}
\end{eqnarray}

In general, for $K\geq 2$ we define $H(\mathbf{z}|\boldsymbol{\beta},\gamma)$ by

\begin{eqnarray}
H(\mathbf{z}|\boldsymbol{\beta},\gamma)=\sum_{i=1}^{n}\sum_{k=1}^{K-1}x_i^T\beta_kI(z_i=k+1)+\gamma \sum_{i=1}^{n}\sum_{\substack{i'\sim i\\i'>i}}\sum_{k=1}^{K}I(z_i=z_{i'}=k)
\end{eqnarray}

With this preamble accomplished, we can define the probability model **automultinomial** seeks to estimate:

\begin{eqnarray}
p(\mathbf{z}|\boldsymbol{\beta},\gamma)=\frac{\exp\{H(\mathbf{z}|\boldsymbol{\beta},\gamma)\}}{\sum_{\mathbf{z'}}\exp\{H(\mathbf{z'}|\boldsymbol{\beta},\gamma)\}}\label{density}
\end{eqnarray}

The term $\sum_{\mathbf{z'}}$ indicates a sum over all possible categorical responses for the entire dataset.

Some intuition regarding \eqref{density} can be gained as follows: when $\gamma=0$, then $\gamma \sum_{i=1}^{n}\sum_{\substack{i'\sim i\\i'>i}}\sum_{k=1}^{2}I(z_i=z_{i'}=k)=0$ so that the $z_i$ are completely independent. In this case the autologistic (automultinomial) model is the same as an ordinary logistic (multinomial logit) model. On the other hand, when $\gamma>0$, then response configurations $\mathbf{z}$ where $z_i=z_{i'}$ for many neighbor pairs $i\sim i'$ become more likely. In this way, $\gamma$ incorporates positive spatial correlation into the responses. When $\gamma<0$, we expect neighboring $z_i,z_{i'}$ to *disagree* more frequently than if the $z_i$ were independent, and the $\gamma<0$ case will in practice be less common.

### Estimating $\boldsymbol{\beta}$ and $\gamma$

The **automultinomial** package follows the pseudolikelihood approach of (Besag, 1974) to estimate the parameters $\beta,\gamma$. We briefly explain the procedure and its motivation below.

When $\gamma$ is known to be $0$, then the responses $z_i$ are independent, and we can use common maximum likelihood techniques such as logistic or multinomial regression to estimate $\beta$. On the other hand, when spatial correlation is present, then we need to estimate $\boldsymbol{\beta}$ and $\gamma$ jointly. Unfortunately, maximum likelihood for the model in equation \eqref{density} is computationally infeasible due to the well known computational intractability of the denominator of $\eqref{density}$ when $\gamma\neq 0$.

Fortunately, a proposal of (Besag, 1974) provides a consistent estimation procedure, maximum pseudolikelihood, for the parameters of the model in equation \eqref{density}. In maximum pseudolikelihood, we maximize

\begin{eqnarray}
\ell_{PL}(\boldsymbol{\beta},\gamma)=\sum_{i=1}^{n}\log\{p(z_i|z_{-i},\boldsymbol{\beta},\gamma)\}\label{pseudolikelihood}
\end{eqnarray}

over $\boldsymbol{\beta}$ and $\gamma$. The expression $p(z_i|z_{-i},\boldsymbol{\beta},\gamma)$ refers to the conditional density of $z_i$ given the sites at all of the other grid locations. A short calculation based on equation \eqref{density} shows that for $K=2$,

\begin{eqnarray}
p(z_i=1|z_{-i},\boldsymbol{\beta},\gamma)&=&\frac{\exp\{\gamma\sum_{i'\sim i}I(z_{i'}=1)\}}{\exp\{\gamma\sum_{i'\sim i}I(z_{i'}=1)\}+\exp\{x_i^T\beta+\gamma\sum_{i'\sim i}I(z_{i'}=2)\}}\label{p1}\\
p(z_i=2|z_{-i},\boldsymbol{\beta},\gamma)&=&\frac{\exp\{x_i^T\beta+\gamma\sum_{i'\sim i}I(z_{i'}=2)\}}{\exp\{\gamma\sum_{i'\sim i}I(z_{i'}=1)\}+\exp\{x_i^T\beta+\gamma\sum_{i'\sim i}I(z_{i'}=2)\}}\label{p2}
\end{eqnarray}

For $K\geq 2$ response categories, we have 

\begin{eqnarray}
p(z_i=1|z_{-i},\boldsymbol{\beta},\gamma)&=&\frac{\exp\{\gamma\sum_{i'\sim i}I(z_{i'}=1)\}}{\exp\{\gamma\sum_{i'\sim i}I(z_{i'}=1)\}+\sum_{k=1}^{K-1}\exp\{x_i^T\beta_k+\gamma\sum_{i'\sim i}I(z_{i'}=k+1)\}}\label{p3}
\end{eqnarray}

and for $1\leq k<K$,
\begin{eqnarray}
p(z_i=k+1|z_{-i},\boldsymbol{\beta},\gamma)&=&\frac{\exp\{x_i^T\beta_k+\gamma\sum_{i'\sim i}I(z_{i'}=k+1)\}\}}{\exp\{\gamma\sum_{i'\sim i}I(z_{i'}=1)\}+\sum_{k=1}^{K-1}\exp\{x_i^T\beta_k+\gamma\sum_{i'\sim i}I(z_{i'}=k+1)\}}\label{p4}
\end{eqnarray}

In the **automultinomial** package, equations \eqref{p1}, \eqref{p2}, \eqref{p3}, and \eqref{p4} are plugged into equation \eqref{pseudolikelihood}, and \eqref{pseudolikelihood} is optimized over $\boldsymbol{\beta}$ and $\gamma$ using the **optim()** function.

\section{Installation}

### cran installation

Simply use 


```r
install.packages("automultinomial")
library(automultinomial)
```

and the package will be ready to run.

### Github installation

The most recent (development) version of the package can be installed from Github.

Installing from Github is as easy as installing from cran. First, make sure Hadley Wickham's package **devtools** is installed on your RStudio:


```r
install.packages("devtools")
```

Then, run the following:


```r
devtools::install_github(repo="stephenberg/automultinomial")
library(automultinomial)
```

Now the package is ready to run.

\section{Data example 1: K=2 response categories}

Here, we will demonstrate how to simulate data using **automultinomial**, how to fit data using **automultinomial**, and how to analyze the output.

### Simulating data

First we will use the **drawSamples()** function to simulate some data from the autologistic model. The **drawSamples()** function takes the following arguments:

* **beta**: the $p\times 1$ coefficient vector (for $K=2$ response categories) or the $p\times (K-1)$ coefficient matrix (for $K>2$ response categories)
* **gamma**: the value of the correlation parameter
* **X**: the $n\times p$ design matrix $X=[x_1,x_2,...,x_n]^T$
* **A**: a square symmetric adjacency matrix encoding the neighborhood structure
* **.***: the number of burnin iterations for the Gibbs sampler. Safe to leave at the default of 300.
* **nSamples**: the number of simulated samples to draw.

First, we will set up the simulation parameters:

```r
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
```

```
##             [,1]
## [1,] -0.04077736
## [2,] -0.01223909
## [3,]  0.30316170
## [4,] -0.04747873
## [5,] -0.64699125
```

```r
#set covariate values
X=matrix(rnorm(n^2*p),ncol=p)

#set the correlation parameter value (0.7 is a moderate amount of spatial correlation)
gamma=0.7
```




Then, we will use the **drawSamples()** function to generate simulated data.

```r
#use drawSamples to simulate data with parameters beta and gamma by Gibbs sampling
y=drawSamples(beta,gamma,X,A,nSamples = 1)
```


```r
y2=drawSamples(beta,0,X,A,nSamples = 1)
```

Figure \ref{fig:plotk2} shows plots of the responses on the grid. On the left plot, we can see "clumping" of the responses due to the positive autocorrelation parameter $\gamma=0.7$. The right hand plot is from the distribution with the same $\boldsymbol{\beta}$, and $\gamma=0$.

\begin{figure}
	\centering
	\includegraphics[width=1\textwidth]{plots/plotk2.png}
	\caption{On the left, a plot of the dataset we just simulated, where the truth is $\gamma=0.7$. For comparison, the plot on the right shows data taken from the same model but with $\gamma=0$.}
	\label{fig:plotk2}
\end{figure}


### Fitting an autologistic model to the data (K=2 categories)

Now we'll fit an autologistic model using the **MPLE()** function. The **MPLE()** function takes 3 primary arguments: 

* **X**: the design matrix $X=[x_1,x_2,...,x_n]^T$
* **A**: a square symmetric adjacency matrix encoding the neighborhood structure
* **y**: the response vector $y$, which is required by **MPLE()** to be a factor vector.

There are also several secondary arguments.

* **ciLevel**: for $xy$% confidence intervals, set **ciLevel** to $0.xy$ (the default is **ciLevel**=0.95, which produces 95% confidence intervals)
* **method**: By choosing **method**="asymptotic", MPLE will output confidence intervals based on the asymptotic distribution of the pseudolikelihood estimator. For bootstrap confidence intervals, use **method**="boot". The default is **method**="asymptotic".
* **burnIn**: the number of burnin iterations for the Gibbs sampler for **method**="boot". Safe to leave at the default of 300.
* **nBoot**: the number of bootstrap samples to use for making bootstrap confidence intervals.

After fitting the model, we will examine confidence intervals for the fitted parameters.  The asymptotic type confidence intervals can be computed very quickly, but may be less accurate in practice than bootstrap confidence intervals.

First, we will fit the model using the **MPLE()** function.


```r
# responses must be input as a factor
y=factor(y)
fit1=MPLE(X=X,y=y,A=A,ciLevel = 0.95,method = "asymptotic")
fit2=MPLE(X=X,y=y,A=A,ciLevel = 0.95,method = "boot",nBoot = 500)
```

Then, we can use the **MPLE\_summary()** function to view the model summary in table form.


```r
#to see the information contained in fit1 and fit2, use str() (not run to save space)
#str(fit1)
#str(fit2)

fitSummary1=MPLE_summary(fit1)
```

```
## 
## 
## Table: Summary with confidence intervals
## 
##      2 vs. 1                  gamma               
## ---  -----------------------  --------------------
## 1    0.001 (-0.139,0.141)     0.718 (0.641,0.795) 
## 2    -0.002 (-0.141,0.136)                        
## 3    0.277 (0.142,0.413)                          
## 4    -0.085 (-0.231,0.061)                        
## 5    -0.686 (-0.851,-0.521)                       
## 
## 
## Table: Summary with p-values
## 
##      2 vs. 1          gamma     
## ---  ---------------  ----------
## 1    0.001 (0.989)    0.718 (0) 
## 2    -0.002 (0.975)             
## 3    0.277 (0)                  
## 4    -0.085 (0.255)             
## 5    -0.686 (0)
```

```r
fitSummary2=MPLE_summary(fit2)
```

```
## 
## 
## Table: Summary with confidence intervals
## 
##      2 vs. 1                 gamma               
## ---  ----------------------  --------------------
## 1    0.001 (-0.136,0.15)     0.718 (0.651,0.805) 
## 2    -0.002 (-0.152,0.155)                       
## 3    0.277 (0.13,0.414)                          
## 4    -0.085 (-0.236,0.048)                       
## 5    -0.686 (-0.836,-0.55)                       
## 
## 
## Table: Summary with p-values
## 
##      2 vs. 1          gamma     
## ---  ---------------  ----------
## 1    0.001 (0.988)    0.718 (0) 
## 2    -0.002 (0.972)             
## 3    0.277 (0)                  
## 4    -0.085 (0.25)              
## 5    -0.686 (0)
```

![](vignetteHTML_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

\section{Data example 2: K=3 response categories}

Now, we demonstrate the use of the package when there are 3 response categories. Essentially, everything is still the same, and most of the previous code doesn't change at all. The only difference is that now, the $\boldsymbol{\beta}$ parameter is a matrix with $2=K-1$ columns, rather than a vector as in the $K=2$ case.

### Simulating data

Generating simulated data using the function **drawSamples()**.


```r
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
```

```
##            [,1]        [,2]
## [1,]  0.4112875 -0.03183735
## [2,] -0.1694095  0.45345660
## [3,]  0.1089385 -0.02839771
## [4,]  0.1898588  0.60552711
## [5,]  0.1212805 -0.01881423
```

```r
#set covariate values
X=matrix(rnorm(n^2*p),ncol=p)

#set the correlation parameter value (0.7 is a moderate amount of spatial correlation)
gamma=0.7
```




```r
#use drawSamples to simulate data with parameters beta and gamma by Gibbs sampling
y=drawSamples(beta,gamma,X,A,nSamples = 1)
```


```r
y2=drawSamples(beta,0,X,A,nSamples = 1)
```

Figure \ref{fig:plotk3} shows plots of the 3-category responses on the grid. On the left plot, we again see "clumping" of the responses due to the positive autocorrelation parameter. The right hand plot is from the distribution with the same $\boldsymbol{\beta}$, and $\gamma=0$.

\begin{figure}
	\centering
	\includegraphics[width=1\textwidth]{plots/plotk3.png}
	\caption{On the left, a plot of the 3-category dataset we just simulated, where the truth is $\gamma=0.7$. For comparison, the plot on the right shows data taken from the same model but with $\gamma=0$.}
	\label{fig:plotk3}
\end{figure}

### Fitting an automultinomial model to the data (K=3 categories)


```r
#responses must be input as a factor
y=factor(y)
fit1=MPLE(X=X,y=y,A=A,ciLevel = 0.95,method = "asymptotic")
fit2=MPLE(X=X,y=y,A=A,ciLevel = 0.95,method = "boot",nBoot = 500)
```


```r
#to see the information contained in fit1 and fit2, use str() (not run to save space)
#str(fit1)
#str(fit2)

fitSummary1=MPLE_summary(fit1)
```

```
## 
## 
## Table: Summary with confidence intervals
## 
##      2 vs. 1                  3 vs. 1                 gamma               
## ---  -----------------------  ----------------------  --------------------
## 1    0.366 (0.224,0.509)      -0.032 (-0.169,0.106)   0.696 (0.624,0.769) 
## 2    -0.204 (-0.343,-0.065)   0.301 (0.158,0.444)                         
## 3    0.079 (-0.064,0.221)     0.085 (-0.047,0.217)                        
## 4    0.181 (0.038,0.324)      0.722 (0.573,0.871)                         
## 5    0.153 (0.014,0.291)      -0.056 (-0.195,0.082)                       
## 
## 
## Table: Summary with p-values
## 
##      2 vs. 1          3 vs. 1          gamma     
## ---  ---------------  ---------------  ----------
## 1    0.366 (0)        -0.032 (0.65)    0.696 (0) 
## 2    -0.204 (0.004)   0.301 (0)                  
## 3    0.079 (0.278)    0.085 (0.206)              
## 4    0.181 (0.013)    0.722 (0)                  
## 5    0.153 (0.03)     -0.056 (0.424)
```

```r
fitSummary2=MPLE_summary(fit2)
```

```
## 
## 
## Table: Summary with confidence intervals
## 
##      2 vs. 1                  3 vs. 1                 gamma               
## ---  -----------------------  ----------------------  --------------------
## 1    0.366 (0.234,0.502)      -0.032 (-0.175,0.116)   0.696 (0.628,0.778) 
## 2    -0.204 (-0.352,-0.049)   0.301 (0.154,0.448)                         
## 3    0.079 (-0.061,0.212)     0.085 (-0.05,0.222)                         
## 4    0.181 (0.042,0.343)      0.722 (0.576,0.878)                         
## 5    0.153 (0.004,0.309)      -0.056 (-0.206,0.098)                       
## 
## 
## Table: Summary with p-values
## 
##      2 vs. 1          3 vs. 1          gamma     
## ---  ---------------  ---------------  ----------
## 1    0.366 (0)        -0.032 (0.694)   0.696 (0) 
## 2    -0.204 (0.004)   0.301 (0)                  
## 3    0.079 (0.3)      0.085 (0.246)              
## 4    0.181 (0.024)    0.722 (0)                  
## 5    0.153 (0.054)    -0.056 (0.434)
```


![](vignetteHTML_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

\section{A comparison of different autologistic parameterizations}

Finally, we'll conclude with a discussion of different parameterizations of the autologistic model. We will focus on the case $K=2$ here. The goal for now is to show by example that the parameterization used by **automultinomial** will be invariant to the choice of reference category (not guaranteed by all of the parameterizations) and to give some concrete data examples.

### Commonly used autologistic parameterizations

Frequently, autologistic parameterizations are given in terms of the conditional probability distributions at each site, rather than in terms of the joint distribution over all sites. This convention might perhaps lead to confusion about what model is actually being fit, but we will follow the common convention here.

We suppose that the response at each site is binary, so that there are $K=2$ response categories.

\textbf{Parameterization 1 (Asymmetric)}: In what we will call the asymmetric parameterization, the set of allowed responses at each site is usually taken to be $\{0,1\}$. The conditional probabilities for the set $\{0,1\}$ are

\begin{eqnarray}
p(z_i=1|z_{-i})=\frac{\exp\{x_i^T\beta+\gamma\sum_{i'\sim i}z_{i'}\}}{1+\exp\{x_i^T\beta+\gamma\sum_{i'\sim i}z_{i'}\}}\\
p(z_i=0|z_{-i})=\frac{1}{1+\exp\{x_i^T\beta+\gamma\sum_{i'\sim i}z_{i'}\}}
\end{eqnarray}

If we shift the values of the responses up by 1, so that we identify $0$ with $1$ and $1$ with $2$, then with a slight abuse of notation we have

\begin{eqnarray}
p(z_i=2|z_{-i})=\frac{\exp\{x_i^T\beta+\gamma\sum_{i'\sim i}I(z_{i'}=2)\}}{1+\exp\{x_i^T\beta+\gamma\sum_{i'\sim i}I(z_{i'}=2)\}}\\
p(z_i=1|z_{-i})=\frac{1}{1+\exp\{x_i^T\beta+\gamma\sum_{i'\sim i}I(z_{i'}=2)\}}\label{p12}
\end{eqnarray}


This form of the conditional distributions corresponds to the model 

\begin{eqnarray}
p_{asy}(\mathbf{z}|\boldsymbol{\beta},\gamma)=\frac{\exp\{H_{asy}(\mathbf{z}|\boldsymbol{\beta},\gamma)\}}{\sum_{\mathbf{z'}}\exp\{H_{asy}(\mathbf{z'}|\boldsymbol{\beta},\gamma)\}}\label{asym}
\end{eqnarray} where \begin{eqnarray}
H_{asy}(\mathbf{z}|\boldsymbol{\beta},\gamma)=\sum_{i=1}^{n}x_i^T\beta I(z_i=2)+\gamma \sum_{i=1}^{n}\sum_{\substack{i'\sim i\\i'>i}}I(z_i=z_{i'}=2)
\end{eqnarray}

By the Hammersley-Clifford theorem, the density in equation \eqref{asym} is the unique joint distribution corresponding to the conditional distributions in equation \eqref{p12}.

\textbf{Parameterization 2 (Symmetric)}: This parameterization is sometimes referred to as the $\pm 1$ parameterization. We will now show that for $K=2$, the symmetric parameterization is equivalent to the parameterization used by **automultinomial**. On the other hand, the formulation given for the **automultinomial** parameterization in section 1 is much more cleanly generalized to $K>2$ categories.

In the symmetric parameterization, the set of allowed outcomes is taken to be $\{-1,1\}$, and the conditional probabilities at each site are \begin{eqnarray}
p(z_i=1|z_{-i})=\frac{\exp\{x_i^T\beta+\gamma\sum_{i'\sim i}z_{i'}\}}{\exp\{-\gamma\sum_{i'\sim i}z_{i'}\}+\exp\{x_i^T\beta+\gamma\sum_{i'\sim i}z_{i'}\}}\\
p(z_i=-1|z_{-i})=\frac{\exp\{-\gamma\sum_{i'\sim i}z_{i'}\}}{\exp\{-\gamma\sum_{i'\sim i}z_{i'}\}+\exp\{x_i^T\beta+\gamma\sum_{i'\sim i}z_{i'}\}}
\end{eqnarray}

Now, as with the asymmetric model, we will map the response set $\{-1,1\}$ to $\{1,2\}$. We will identify $-1$ with $1$ and $1$ with $2$, so that the conditional probabilities become

\begin{eqnarray}
p(z_i=2|z_{-i})&=&\frac{\exp\left[x_i^T\beta+\gamma\sum_{i'\sim i}\{I(z_{i'}=2)-I(z_{i'}=1)\}\right]}{1+\exp\left[x_i^T\beta+\gamma\sum_{i'\sim i}\{I(z_{i'}=2)-I(z_{i'}=1)\}\right]}\\
&=&\frac{\exp\{x_i^T\beta +\gamma\sum_{i'\sim i}I(z_{i'}=2)\}}{\exp\{\gamma\sum_{i'\sim i}I(z_{i'}=1)\}+\exp\{x_i^T\beta+\gamma\sum_{i'\sim i}I(z_{i'}=2)\}}\\
p(z_i=1|z_{-i})&=&\frac{\exp\{\gamma\sum_{i'\sim i}I(z_{i'}=1)\}}{\exp\{\gamma\sum_{i'\sim i}I(z_{i'}=1)\}+\exp\{x_i^T\beta+\gamma\sum_{i'\sim i}I(z_{i'}=2)\}}\label{p22}
\end{eqnarray}

This form of the conditional distributions corresponds to the model 

\begin{eqnarray}
p_{sym}(\mathbf{z}|\boldsymbol{\beta},\gamma)=\frac{\exp\{H_{sym}(\mathbf{z}|\boldsymbol{\beta},\gamma)\}}{\sum_{\mathbf{z'}}\exp\{H_{sym}(\mathbf{z'}|\boldsymbol{\beta},\gamma)\}}\label{sym}
\end{eqnarray}

\begin{eqnarray}
H_{sym}(\mathbf{z}|\boldsymbol{\beta},\gamma)=\sum_{i=1}^{n}x_i^T\beta I(z_i=2)+\gamma \sum_{i=1}^{n}\sum_{\substack{i'\sim i\\i'>i}}\sum_{k=1}^{2}I(z_i=z_{i'}=k)\label{hsym}
\end{eqnarray}

By the Hammersley-Clifford theorem, the density in equation \eqref{sym} is the unique joint distribution corresponding to the conditional distributions in equation \eqref{p22}.


\textbf{Parameterization 3: the parameterization used in **automultinomial**}s: This parameterization is the parameterization defined in section 2 of this document.

\subsection{Comparison}
Having defined the 3 parameterizations, we will now compare them. Our first conclusion is that when $K=2$, the symmetric ($\pm 1$ ) parameterization is *exactly* equivalent to the parameterization used by **automultinomial**. This conclusion follows from determining that the energies in equations \eqref{H} and \eqref{hsym} are exactly the same, so that the densities in equation \eqref{sym} and \eqref{density} are exactly the same.

Therefore, from now on we need only compare the **automultinomial** parameterization and the asymmetric parameterization. The following statements can be shown:

* When each site has the same number of neighbors, and the covariate matrix $X=[x_1,x_2,...,x_n]^T$ includes an intercept column, then the asymmetric and the parameterization used in **automultinomial** are equivalent, in the sense that any joint distribution in the **automultinomial** parameterization has an equivalent representation in the asymmetric parameterization. The pseudolikelihoods obtained by fitting both parameterizations to a dataset will be identical. The coefficients from each parameterization can be obtained from the other.
* When the number of neighbors at each site is different, and/or the covariate matrix does not include an intercept column, then the asymmetric and the **automultinomial** parameterizations will in general *not* be equivalent. The **automultinomial** parameterization is invariant to the choice of reference category, whereas the asymmetric parameterization will in general not be invariant to the choice of reference category.

We will now give some examples on simulated data. 


```r
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
gamma=0.5
```


```r
#simulated responses
y=drawSamples(beta=beta,gamma=gamma,X = X_intercept, A=a1,nSamples = 1)
```

\textbf{Case 1: same number of neighbors, covariate matrix includes intercept}

Here, we see that the pseudolikelihoods match and that coefficients (excepting $\gamma$ and the intercept) also match. The intercept and $\gamma$ from the two models are related by a simple linear transformation.


```r
fitAsymmetric=glm(factor(y)~cbind(X,as.matrix(a1%*%(y-1))),family=binomial())
fitAutomultinomial=MPLE(X=X_intercept,y=factor(y),A = a1,method="asymptotic")
```

```
## Starting model fitting
## Model fitting done, starting variance estimation
## Creating asymptotic confidence intervals
```

```r
#comparing pseudolikelihoods: identical
fitAsymmetric$deviance/-2
```

```
## [1] -431.7837
```

```r
fitAutomultinomial$pseudolikelihood
```

```
## [1] -431.7837
```

```r
#comparing beta: identical
fitAsymmetric$coefficients[2:3]
```

```
## cbind(X, as.matrix(a1 %*% (y - 1)))1 cbind(X, as.matrix(a1 %*% (y - 1)))2 
##                           -0.1496811                           -0.1364893
```

```r
fitAutomultinomial$betaHat[2:3]
```

```
## [1] -0.1496946 -0.1364924
```

```r
#comparing intercept and gamma: not identical, but there is a mapping between them
gammaAsymmetric=fitAsymmetric$coefficients[4]
gammaAutomultinomial=fitAutomultinomial$gammaHat
gammaAsymmetric
```

```
## cbind(X, as.matrix(a1 %*% (y - 1)))3 
##                             1.162439
```

```r
gammaAutomultinomial*2
```

```
## [1] 1.162478
```

```r
interceptAutomultinomial=fitAutomultinomial$betaHat[1]
interceptAsymmetric=fitAsymmetric$coefficients[1]
interceptAsymmetric
```

```
## (Intercept) 
##   -1.671265
```

```r
interceptAutomultinomial-4*gammaAutomultinomial
```

```
## [1] -1.671386
```

\textbf{Case 2: different number of neighbors, covariate matrix includes intercept}

Here, we see that the pseudolikelihoods and coefficients no longer match between parameterizations. Additionally, while the **automultinomial** parameterization is invariant to changes in the reference category, the asymmetric parameterization is not.


```r
fitAsymmetric1=glm(factor(y)~cbind(X,as.matrix(a2%*%(y-1))),family=binomial())
fitAutomultinomial1=MPLE(X=X_intercept,y=factor(y),A = a2,method="asymptotic")
```

```
## Starting model fitting
## Model fitting done, starting variance estimation
## Creating asymptotic confidence intervals
```

```r
#refitting with different reference category
y2=as.factor(y)
y2=as.numeric(relevel(y2,ref = "2"))
fitAsymmetric2=glm(factor(y2)~cbind(X,as.matrix(a2%*%(y2-1))),family=binomial())
fitAutomultinomial2=MPLE(X=X_intercept,y=factor(y2),A = a2,method="asymptotic")
```

```
## Starting model fitting
## Model fitting done, starting variance estimation
## Creating asymptotic confidence intervals
```

```r
#comparing pseudolikelihoods: not identical
fitAsymmetric1$deviance/-2 
```

```
## [1] -444.9057
```

```r
fitAsymmetric2$deviance/-2
```

```
## [1] -435.4322
```

```r
fitAutomultinomial1$pseudolikelihood
```

```
## [1] -438.5099
```

```r
fitAutomultinomial2$pseudolikelihood
```

```
## [1] -438.5099
```

```r
#comparing beta: not identical
fitAsymmetric1$coefficients[2:3]
```

```
## cbind(X, as.matrix(a2 %*% (y - 1)))1 cbind(X, as.matrix(a2 %*% (y - 1)))2 
##                           -0.1425823                           -0.1352519
```

```r
fitAsymmetric2$coefficients[2:3]
```

```
## cbind(X, as.matrix(a2 %*% (y2 - 1)))1 
##                             0.1504164 
## cbind(X, as.matrix(a2 %*% (y2 - 1)))2 
##                             0.1323744
```

```r
fitAutomultinomial1$betaHat[2:3]
```

```
## [1] -0.1475824 -0.1334123
```

```r
fitAutomultinomial2$betaHat[2:3]
```

```
## [1] 0.1475824 0.1334123
```

```r
#comparing gamma: not identical. Asymmetric parameterization not invariant to
#choice of reference category. Automultinomial parameterization invariant to 
#choice of reference category
gammaAsymmetric1=fitAsymmetric1$coefficients[4]
gammaAsymmetric2=fitAsymmetric2$coefficients[4]
gammaAutomultinomial1=fitAutomultinomial1$gammaHat
gammaAutomultinomial2=fitAutomultinomial2$gammaHat
gammaAsymmetric1
```

```
## cbind(X, as.matrix(a2 %*% (y - 1)))3 
##                            0.9394255
```

```r
gammaAsymmetric2
```

```
## cbind(X, as.matrix(a2 %*% (y2 - 1)))3 
##                              1.114119
```

```r
gammaAutomultinomial1
```

```
## [1] 0.5378779
```

```r
gammaAutomultinomial2
```

```
## [1] 0.5378779
```

\textbf{Case 3: same number of neighbors, covariate matrix does not include intercept}

Here, we see that the pseudolikelihoods and coefficients no longer match between parameterizations. Additionally, while the **automultinomial** parameterization is invariant to changes in the reference category, the asymmetric parameterization is not.


```r
fitAsymmetric1=glm(factor(y)~cbind(X,as.matrix(a2%*%(y-1)))-1,family=binomial())
fitAutomultinomial1=MPLE(X=X,y=factor(y),A = a2,method="asymptotic")
```

```
## Starting model fitting
## Model fitting done, starting variance estimation
## Creating asymptotic confidence intervals
```

```r
#refitting with different reference category
y2=as.factor(y)
y2=as.numeric(relevel(y2,ref = "2"))
fitAsymmetric2=glm(factor(y2)~cbind(X,as.matrix(a2%*%(y2-1)))-1,family=binomial())
fitAutomultinomial2=MPLE(X=X,y=factor(y2),A = a2,method="asymptotic")
```

```
## Starting model fitting
## Model fitting done, starting variance estimation
## Creating asymptotic confidence intervals
```

```r
#comparing pseudolikelihoods: not identical
fitAsymmetric1$deviance/-2 
```

```
## [1] -447.2382
```

```r
fitAsymmetric2$deviance/-2
```

```
## [1] -1034.841
```

```r
fitAutomultinomial1$pseudolikelihood
```

```
## [1] -449.5132
```

```r
fitAutomultinomial2$pseudolikelihood
```

```
## [1] -449.5132
```

```r
#comparing beta: not identical
fitAsymmetric1$coefficients[1:2]
```

```
## cbind(X, as.matrix(a2 %*% (y - 1)))1 cbind(X, as.matrix(a2 %*% (y - 1)))2 
##                           -0.1414015                           -0.1388807
```

```r
fitAsymmetric2$coefficients[1:2]
```

```
## cbind(X, as.matrix(a2 %*% (y2 - 1)))1 
##                            0.02259219 
## cbind(X, as.matrix(a2 %*% (y2 - 1)))2 
##                            0.06077764
```

```r
fitAutomultinomial1$betaHat[1:2]
```

```
## [1] -0.1492992 -0.1240118
```

```r
fitAutomultinomial2$betaHat[1:2]
```

```
## [1] 0.1492992 0.1240118
```

```r
#comparing gamma: not identical. Asymmetric parameterization not invariant to
#choice of reference category. Automultinomial parameterization invariant to 
#choice of reference category
gammaAsymmetric1=fitAsymmetric1$coefficients[3]
gammaAsymmetric2=fitAsymmetric2$coefficients[3]
gammaAutomultinomial1=fitAutomultinomial1$gammaHat
gammaAutomultinomial2=fitAutomultinomial2$gammaHat
gammaAsymmetric1
```

```
## cbind(X, as.matrix(a2 %*% (y - 1)))3 
##                             0.698996
```

```r
gammaAsymmetric2
```

```
## cbind(X, as.matrix(a2 %*% (y2 - 1)))3 
##                             -1.021502
```

```r
gammaAutomultinomial1
```

```
## [1] 0.7831986
```

```r
gammaAutomultinomial2
```

```
## [1] 0.7831986
```
