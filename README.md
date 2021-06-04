# FGSPCA

This repository is the official implementation of **[Feature Grouping and Sparse Principle Component Analysis]()**.


Sparse principal component analysis (SPCA) attempts to find sparse weight vectors (loadings), i.e., a loading vector with only a few 'active' (nonzero) values. The SPCA approach provides better interpretability for the principal components in high-dimensional data settings. Because the principal components are formed as a linear combination of only a few of the original variables. This package provides a modified sparse principal component analysis, Feature Grouping and Sparse Principal Component Analysis (FGSPCA), which considers additional structure information among loadings (feature grouping) as well as the sparsity (feature selection) property among loadings. 


This package provides FGSPCA in R:
 
* Featrue Grouping and Sparse PCA: ``FGSPCA()``.


## Problem Formulation


FGSPCA considers an additional grouping structure where the loadings share similar coefficients (i.e., feature grouping), besides a special group with all coefficients being zero (i.e., feature selection). 
FGSPCA allows the loadings to belong to disjoint homogeneous groups, with sparsity as a special case.
More concreatly, given a data matrix ``X`` with shape ``(n,p)``, FGSPCA attemps to minimize the following
optimization problem 

```math
minimize S(A,B) = 1/2 \|X - X\cdot B\cdot A^T \|_F^2 + \lambda \|B\|_2^2 + \lambda_1 \sum_{j=1}^{k} p_1(\beta_l) + \lambda_2 \sum_{j=1}^{k}p_2(\beta_l), 
subject to A^T A = I .
```
<!-- $$ minimize S(A,B) = 1/2 \|X - X\cdot B\cdot A^T \|_F^2 + \lambda \|B\|_2^2 + \lambda_1 \sum_{j=1}^{k} p_1(\beta_l) + \lambda_2 \sum_{j=1}^{k}p_2(\beta_l), 
subject to A^T A = I . $$ -->
<!-- <pre xml:lang="latex"> minimize S(A,B) = 1/2 \|X - X\cdot B\cdot A^T \|_F^2 + \lambda \|B\|_2^2 + \lambda_1 \sum_{j=1}^{k} p_1(\beta_l) + \lambda_2 \sum_{j=1}^{k}p_2(\beta_l), 
subject to A^T A = I</pre> -->


```
minimize S(A, B), subject to A^TA=I.
```
![\large S(A,B)=1/2\|X-XBA^T\|_F^2+\lambda\|B\|_2^2+\lambda_1\sum_{j=1}^{k}p_1(\beta_l)+\lambda_2\sum_{j=1}^{k}p_2(\beta_l)](https://latex.codecogs.com/svg.latex?\large&space;S(A,B)=1/2\|X-XBA^T\|_F^2+\lambda\|B\|_2^2+\lambda_1\sum_{j=1}^{k}p_1(\beta_l)+\lambda_2\sum_{j=1}^{k}p_2(\beta_l)) 

<!-- ![\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}](https://latex.codecogs.com/svg.latex?\Large&space;x=\frac{-b\pm\sqrt{b^2-4ac}}{2a})  -->

subject to 
![\large A^TA=I](https://latex.codecogs.com/svg.latex?\large&space;A^TA=I) 

Here 
![\large p_1(\beta_j)=\sum_{l=1}^p\min(\frac{|\beta_{jl}|}{\tau},1)](https://latex.codecogs.com/svg.latex?\large&space;p_1(\beta_j)=\sum_{l=1}^p\min(\frac{|\beta_{jl}|}{\tau},1)) 

![\large p_2(\beta_j)=\sum_{l<l',(l,l')E}\min(\frac{|\beta_{jl}-\beta_{jl'}|}{\tau},1)
](https://latex.codecogs.com/svg.latex?\large&space;p_2(\beta_j)=\sum_{l<l',(l,l')E}\min(\frac{|\beta_{jl}-\beta_{jl'}|}{\tau},1)
) 

<!-- ![\large A^TA=I](https://latex.codecogs.com/svg.latex?\large&space;A^TA=I) 
p_1(\beta_j)=\sum_{l=1}^p\min(\frac{|\beta_{jl}|}{\tau},1)
p_2(\beta_j)=\sum_{l<l',(l,l')\in E}\min(\frac{|\beta_{jl}-\beta_{jl'}|}{\tau},1) -->


Here $p_1(\beta_j)$ and $p_2(\beta_j)$ are two non-convex regularization functions.
The matrix ``B`` is the sparse loadings matrix with grouping effect and ``A`` is an orthonormal matrix.
Here we use a combination of the $p_1(\beta_j)$ and $p_1(\beta_j)$ as non-convex sparsity-promoting and grouping-promoting regularizers.
Then, the principal components ``Z`` are then formed as

```
Z = X %*% B .
```

Specifically, the interface of the FGSPCA function is:

```R
FGSPCA(X, B_init, K, para, type="predictor", tau_S=0.05, lambda2=0.1, lambda3=0.1,use.corr=FALSE)

```


The description of the arguments is listed in the following:

* ``X`` is a real ``n`` by ``p`` data matrix (or data frame), where ``n`` denotes the number of observations and ``p`` the number of variables.

* ``B_init`` the initial value of matrix B.

* ``K`` specifies the target rank, i.e., number of components to be computed.

* ``para`` is a list of $\lambda_1$ with length of ``K``, the sparsity controlling parameters, tuning parameter corresponding to $p_1(\beta_l)$. Higher values lead to sparser components. 

* ``type``  type of ``X``, which can take values of ``predictor``, or ``Gram``. If `type="Gram"` the model should deal with the root matrix of ``X``. If `type="predictor"` the model should directly deal with the matrix ``X``.

* ``tau_S`` $\tau$ the controlling parameter corresponding to $p_1(\beta_l)$ and $p_2(\beta_l)$, which determines when the small values of $|\beta_l|$ will be penalized and when the small difference values of $|\beta_l - \beta_{l'}|$ will be penalized.

* ``lambda2`` $\lambda_2$ tuning parameter corresponding to $p_2(\beta_l)$.
* ``lambda3`` $\lambda$ tuning parameter corresponding to $\lambda \|B\|_2^2$, controls the amount of ridge shrinkage to apply in order to improve conditioning.
* ``use.corr`` logical value which indicates whether the variables should be shifted to be zero centered (FALSE by default).


A list with the following components is returned:

* ``loadings`` sparse loadings (weight) matrix.
* ``pev`` percentage of explained variance.
* ``var.all`` the total variance.
* ``Alpha`` the data matrix ``A``.


## Installation


To install R package from github, you need to clone this repository.


```setup
git clone https://github.com/ipapercodes/FGSPCA
```

Then just install the FGSPCA package from local files.

```R
install.packages("~/FGSPCA_0.1.0.tar.gz", repos = NULL, type = "source")
library(FGSPCA)
```

## Example: Sparse PCA


```R
library(elasticnet)  # to use the data "pitprop"
library(sparsepca)
data(pitprops)
# ## elasticnet::spca is the Sparse PCA of Zou (2006).
K <- 6
para <- c(0.06, 0.16, 0.1, 0.5, 0.5, 0.5)
out1 <- elasticnet::spca(
  pitprops, K=6, type="Gram", sparse="penalty", trace=TRUE,
  para=c(0.06,0.16,0.1,0.5,0.5,0.5))

K <- k <- 6
X <- rootmatrix(pitprops)
rspca.results <- sparsepca::rspca(X, k, alpha=1e-3, beta=1e-3, center=TRUE, scale=FALSE)

K <- 6
B_init <- rspca.results$loadings
tau_S=0.05; lambda1=0.01; lambda2=0.01; lambda3=0.1
para <- rep(lambda1, K)

out <- FGSPCA(
  pitprops, B_init, K, para, type="Gram",
  tau_S=tau_S, lambda2=lambda2, lambda3=lambda3)
NB <- out$loadings  # the normalized loadings of FGSPCA

(fgspca.pev <- round(out$pev *100, 3))
(fgspca.cpev <- round(cumsum(out$pev) * 100, 3))

```

## References

* [Zou, Hui, Trevor Hastie, and Robert Tibshirani. "Sparse principal component analysis." Journal of computational and graphical statistics 15.2 (2006): 265-286.](https://web.stanford.edu/~hastie/Papers/spc_jcgs.pdf)

