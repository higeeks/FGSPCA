# library(mvtnorm)
# library(Matrix)
# library(glmnet)
# library(lars)  ### the function lars::backsolvet, lars::downdateR
# library(elasticnet)
# library(sparsepca)


#' Calculate the root matrix of the given square matrix using eigendecomposition
#'
#' @param x the given square matrix with size \eqn{n\times n}
#' @return the root matrix
rootmatrix <- function(x){
  x.eigen <- eigen(x)
  d <- x.eigen$values
  d <- (d+abs(d))/2
  v <- x.eigen$vectors
  return (v%*%diag(sqrt(d))%*%t(v))
}


#' Convergence check of two successive matrices
#'
#' Check the convergence of two successive matrices
#'
#' The two matrices \eqn{B_1} and \eqn{B_2} may be different only in their column signs
#'   with respective to each column.
#' @param beta1 \eqn{B_1} a matrix
#' @param beta2 \eqn{B_2} a matrix
#'
#' @return the largest absolute value of the difference of their corresponding columns
convcheck <- function(beta1,beta2){
  a <- apply(abs(beta1+beta2),2,max)  # column max
  b <- apply(abs(beta1-beta2),2,max)
  d <- length(a)
  x <- rep(1, d)
  for (i in 1:d){
    x[i] <- min(a[i],b[i])
  }
  max(x)
}


updateRR <- function(
  xnew, R = NULL, xold, lambda, eps = .Machine$double.eps){
  xtx <- (sum(xnew^2)+lambda)/(1+lambda)
  norm.xnew <- sqrt(xtx)
  if(is.null(R)) {
    R <- matrix(norm.xnew, 1, 1)
    attr(R, "rank") <- 1
    return(R)
  }
  Xtx <- drop(t(xnew) %*% xold)/(1+lambda)
  r <- lars::backsolvet(R, Xtx)
  rpp <- norm.xnew^2 - sum(r^2)
  rank <- attr(R, "rank")
  if(rpp <= eps)
    rpp <- eps
  else {
    rpp <- sqrt(rpp)
    rank <- rank + 1
  }
  R <- cbind(rbind(R, 0), c(r, rpp))
  attr(R, "rank") <- rank
  R
}

#' The SPCA solver of beta for the elastic net problem
#'
#' The beta solver using the sparse PCA.
#'
#' The standard objective function of elastic net is
#' \deqn{1/(2n) \| Y-X\beta\|_2^2 + \lambda (\alpha \|\beta\|_1 + (1-\alpha )/2\|\beta\|_2^2) .}
#' But here we use the following objective function
#' \deqn{1/(2n) \| Y-X\beta\|_2^2 + \lambda_1 \|\beta\|_1 + \lambda_2/2 \|\beta\|_2^2 .}
#'
#' @param x the data matrix \eqn{X_{n\times p}}
#' @param y the response vector \eqn{Y_{n\times 1}}
#' @param paras the combination of parameters \eqn{(\lambda_2, \lambda_1)}\\
#'   `lambda = paras[1]` (which is lambda2) and `lambda1 = paras[2]` .
#' @param max.steps `100` (default) the maximum number of steps for the updating
#' @param sparse `sparse = c("penalty", "varnum")`
#' @param eps the tolerance of the stopping criterion for the termination
#' @return the solution \eqn{\beta} in each subproblem
#'
#' @references
#' \itemize{
#'  \item Zou, Hui, Trevor Hastie, and Robert Tibshirani. "Sparse principal component analysis."
#'    Journal of computational and graphical statistics 15.2 (2006): 265-286.
#' }
solvebeta <- function(
  x, y, paras, max.steps, sparse="penalty", eps = .Machine$double.eps) {
  nm <- dim(x)
  n <- nm[1]
  m <- nm[2]
  im <- seq(m)
  one <- rep(1, n)
  vn <- dimnames(x)[[2]]
  lambda <- paras[1]
  if(lambda>0){
    maxvars <- m
  }
  if (lambda==0) {
    maxvars <- min(m,n-1)
    if (m==n){
      maxvars <- m
    }
  }
  d1 <- sqrt(lambda)
  d2 <- 1/sqrt(1 + lambda)
  Cvec <- drop(t(y) %*% x) * d2
  ssy <- sum(y^2)
  residuals <- c(y, rep(0, m))
  if(missing(max.steps)) {max.steps <- 50 * maxvars}
  penalty <- max(abs(Cvec))
  if (sparse=="penalty" && penalty*2/d2<=paras[2])
  {
    beta <- rep(0,m)
  }
  else {
    beta <- rep(0, m)
    first.in <- integer(m)
    active <- NULL
    ignores <- NULL
    actions <- as.list(seq(max.steps))
    drops <- FALSE
    Sign <- NULL
    R <- NULL
    k <- 0
    while((k < max.steps) & (length(active) < maxvars - length(ignores)))
    {
      action <- NULL
      k <- k + 1
      inactive <- if(k == 1) im else im[ - c(active, ignores)]
      C <- Cvec[inactive]
      Cmax <- max(abs(C))
      if(!any(drops)) {
        new <- abs(C) == Cmax
        C <- C[!new]
        new <- inactive[new]
        for(inew in new) {
          R <- updateRR(x[, inew], R, x[, active], lambda)
          if(attr(R, "rank") == length(active)) {
            nR <- seq(length(active))
            R <- R[nR, nR, drop = FALSE]
            attr(R, "rank") <- length(active)
            ignores <- c(ignores, inew)
            action <- c(action,  - inew)
          }
          else {
            if(first.in[inew] == 0)
              first.in[inew] <- k
            active <- c(active, inew)
            Sign <- c(Sign, sign(Cvec[inew]))
            action <- c(action, inew)
          }
        }
      }
      else action <-  - dropid
      Gi1 <- backsolve(R, lars::backsolvet(R, Sign))
      A <- 1/sqrt(sum(Gi1 * Sign))
      w <- A * Gi1
      u1 <- drop(x[,active,drop=FALSE]%*%w*d2)
      u2 <- rep(0,m)
      u2[active] <- d1*d2*w
      u <- c(u1,u2)
      if(lambda>0){
        maxvars <- m-length(ignores)
      }
      if (lambda==0){
        maxvars <- min(m-length(ignores),n-1)
      }
      if(length(active) == maxvars - length(ignores)) {
        gamhat <- Cmax/A
      }
      else {
        a <- (drop(u1 %*% x[,  - c(active, ignores)]) + d1 *
                u2[ - c(active, ignores)]) * d2
        gam <- c((Cmax - C)/(A - a), (Cmax + C)/(A + a))
        gamhat <- min(gam[gam > eps], Cmax/A)
        Cdrop <- c(C - gamhat * a,  - C + gamhat * a) - (Cmax - gamhat * A)
      }
      dropid <- NULL
      b1 <- beta[active]
      z1 <-  - b1/w
      zmin <- min(z1[z1 > eps], gamhat)
      if(zmin < gamhat) {
        gamhat <- zmin
        drops <- z1 == zmin
      }
      else drops <- FALSE
      beta2 <- beta
      beta[active] <- beta[active] + gamhat * w
      residuals <- residuals - (gamhat*u)
      Cvec <- (drop(t(residuals[1:n]) %*% x) + d1 * residuals[ - (1:n)]) * d2
      penalty <- c(penalty,penalty[k]-abs(gamhat*A))
      if (sparse=="penalty" && rev(penalty)[1]*2/d2<=paras[2]){
        s1 <- rev(penalty)[1]*2/d2
        s2 <- rev(penalty)[2]*2/d2
        beta <- (s2-paras[2])/(s2-s1)*beta+(paras[2]-s1)/(s2-s1)*beta2
        beta <- beta*d2
        break
      }
      if(any(drops)) {
        dropid <- seq(drops)[drops]
        for(id in rev(dropid)) {
          R <- lars::downdateR(R, id)
        }
        dropid <- active[drops]
        beta[dropid] <- 0
        active <- active[!drops]
        Sign <- Sign[!drops]
      }
      if(sparse=="varnum" && length(active)>=paras[2]){
        break
      }
      if(!is.null(vn))
        names(action) <- vn[abs(action)]
      actions[[k]] <- action
    }
  }
  return (beta)
}


soft <- function(a, s){
  b <- abs(a)-s
  b <- (b+abs(b))/2
  b <- sign(a)*b
  b
}

#' The elastic net loss function
#'
#' The elastic net loss function (not the standard version)
#'
#' The standard version of the elastic net problem is defined as follows,
#' \deqn{1/(2n) \| Y-X\beta\|_2^2 + \lambda (\alpha \|\beta\|_1 + (1-\alpha )/2\|\beta\|_2^2) .}
#'
#' The objective function of the elastic net we use in this function is defined as
#' \deqn{1/(2n) \| Y-X\beta\|_2^2 + \lambda_1\|\beta\|_1 + \lambda_2\|\beta\|_2^2 .}
#' It is not the standard version of the elastic net problem.
#'
#' @param x the data matrix \eqn{X_{n\times p}}
#' @param y the response vector \eqn{Y_{n \times 1}}
#' @param beta \eqn{\beta} the estimation of \eqn{\beta}
#' @param lambda1 the \eqn{\lambda_1} in the loss function
#' @param lambda2 the \eqn{\lambda_2} in the loss function
#' @param verbose whether to return a list, `FALSE` (default)
#' @return the value of the loss function or a list `c(loss, pen, obj)`
enet_loss <- function(x, y, beta, lambda1, lambda2, verbose=FALSE) {
  s1 <- mean((y - x %*% beta)^2) / 2
  pen1 <- lambda1*sum(abs(beta))
  pen2 <- lambda2*sum((beta)^2)
  out <- s1 + pen1 + pen2
  res <- list(loss=s1, pen=pen1+pen2, obj=out)
  if (verbose) {return(res)} else {return(out)}
}



#' The coordinate descent for LASSO (elastic net) regression
#'
#' The coordinate descent to solve the LASSO (elastic net) regression problem.
#' It is not the standard version of the elastic net problem, so we use `coordinate_descent_LASSO`.
#' The standard version of the elastic net problem is termed `coordinate_descent_enet`.
#'
#' The objective function of the lasso (elastic net) is
#' \deqn{1/(2n) \| Y-X\beta\|_2^2 + \lambda_1\|\beta\|_1 + \lambda_2\|\beta\|_2^2 .}
#' While the soft thresholding method applied to the matrix version may not lead to convergence,
#'  the updating using the soft thresholding along each coordinate can guarantee the convergence.
#'
#' @param x the data matrix \eqn{X_{n\times p}}
#' @param y the response vector \eqn{Y_{n \times 1}}
#' @param paras the combination of parameters \eqn{\lambda_2, \lambda_1}
#' @param max.steps maximum steps, the maximum number of steps for the updating, `100` (default)
#' @param condition_tol the tolerance for the condition to stop, `1e-3` (default)
#' @param loss_return whether to return loss or not, `FALSE` (default)
#'
#' @return \eqn{\beta} the solution \eqn{\beta} or a list
#' @return a list of \eqn{\beta}, a sequence of the objective values during the CD process,
#'    a sequence of the loss values,  a sequence of the penalty values,
#'    `list(beta=beta, obj=obj, loss=loss, pen=pen)`.
coordinate_descent_LASSO <- function(
  x, y, paras, max.steps=100, condition_tol=1e-3, loss_return=FALSE) {
  # max.steps=100; condition_tol=1e-3
  lambda2 <- paras[1]
  lambda1 <- paras[2]
  n <- dim(x)[1]
  p <- dim(x)[2]
  # ## The iteration of the Coordinate Descent.
  count_k <- 0
  condition_k <- 1e3  # ## a very large value
  beta <- rep(0, p)   # ## start with some initial guess
  loss <- vector(mode = "numeric")
  pen <- vector(mode = "numeric")
  obj <- vector(mode = "numeric")
  meanxj <- 1/n*colSums(x^2)
  while ((count_k < max.steps) & (condition_k > condition_tol) ) {
    beta_old <- beta
    # ## The procedure of Coordinate Descent
    for (j in 1:p) {
      alpha <- meanxj[j] + lambda2
      bj <- y - x[, -j] %*% beta[-j]
      gamma1 <- mean(x[,j]*bj)
      beta[j] <- 1/alpha * soft(gamma1, lambda1)
      # if (alpha == 0) {
      #   beta[j] <- 0
      # } else {
      #   beta[j] <- 1/alpha * soft(gamma1, lambda1)
      # }
    }
    res <- enet_loss(x, y, beta, lambda1, lambda2, verbose = TRUE)
    loss <- c(loss, res$loss)
    pen <- c(pen, res$pen)
    obj <- c(obj, res$obj)
    count_k <- count_k + 1
    condition_k <- max(abs(beta - beta_old))  # max norm
  }
  if (loss_return) {
    return(list(beta=beta, obj=obj, loss=loss, pen=pen))
  } else {
    return(beta)
  }
}



#' The coordinate descent for the standard elastic net regression problem
#'
#' The objective function of the standard elastic net regression problem is
#' \deqn{1/(2n) \| Y-X\beta\|_2^2 + \lambda (\alpha \|\beta\|_1 + (1-\alpha )/2\|\beta\|_2^2) .}
#' It is the standard version of the elastic net problem.
#'
#' @inheritParams coordinate_descent_LASSO
#'
#' @param lambda the \eqn{\lambda} in the elastic net loss function
#' @param alpha the \eqn{\alpha} in the elastic net loss function
#'
#' @inherit  coordinate_descent_LASSO return params
coordinate_descent_enet <- function(
  x, y, lambda, alpha, max.steps=100, condition_tol=1e-3, loss_return=FALSE) {
  # max.steps=100; condition_tol=1e-3
  lambda2 <- lambda * alpha
  lambda1 <- lambda * (1-alpha )
  result <- coordinate_descent_LASSO(
    x, y, paras=c(lambda2, lambda1), max.steps, condition_tol, loss_return)
  return(result)
}


#' Different beta solvers
#'
#' Different beta solvers
#'
#' @inheritParams solvebeta
#' @inheritParams coordinate_descent_LASSO
#'
#' @param lambda2 \eqn{\lambda_2} in the elastic loss function
#' @param lambda1 \eqn{\lambda_1} in the elastic loss function
#' @param solver_type `solver_type=c("spcasolver", "lassocd", "glmnet")`
#'
#' @return the solution \eqn{\beta}
# If sparse = "penalty", then para would be a vector of penalty values;
# If sparse="varnum", para would be a vector of integers.
# sparse, paras, works only for setting with solver_type="spcasolver".
# solver_type="spcasolver", "lassocd", "glmnet".
different_beta_solver <- function(
  x, y, lambda2, lambda1, sparse, solver_type = "spcasolver"){
  if (solver_type == "spcasolver") {
    btemp <- solvebeta(x, y, paras=c(lambda2, lambda1), sparse=sparse)
  } else if (solver_type == "lassocd") {
    btemp <- coordinate_descent_LASSO(x, y, paras=c(lambda2, lambda1), loss_return=FALSE)
  } else if (solver_type == "glmnet") {
    n <- dim(x)[1]
    p <- dim(x)[2]
    y1 <- c(y, rep(0, p))
    x1 <- rbind(x, diag(sqrt(n*lambda2), p, p))
    # ## alpha=1, lasso, use transformed data.
    glmnet.result <- glmnet::glmnet(
      x1, y1, family="gaussian", lambda=lambda1, alpha=1, intercept = FALSE)
    btemp <- round(as.numeric(glmnet.result$beta), 4)
  } else {
    btemp <- solvebeta(x, y, paras=c(lambda2, lambda1), sparse=sparse)
  }
  return(btemp)
}



# type = c("predictor", "Gram")
# sparse = c("penalty", "varnum").
# ## If sparse = "penalty", then para would be penalty;
# ## if sparse="varnum", para would be integers.
# ## type, sparse, para, works only for setting with solver_type="spcasolver".
# solver_type="spcasolver", "lassocd", "glmnet"
spca_modified <- function(
  x, K, para, type="predictor", sparse="penalty", use.corr=FALSE,
  lambda=1e-6, max.iter=200, trace=FALSE, eps.conv=1e-3,
  solver_type="spcasolver") {
  vn <- dimnames(x)[[2]]
  x <- switch(
    type,
    predictor = {
      n <- dim(x)[1]
      p <- dim(x)[2]
      if (n/p >= 1e8) {
        cat("You may wish to restart and use a more efficient way \n")
        cat("let the argument x be the sample covariance/correlation matrix and set type=Gram \n")
      }
      x <- scale(x, center=TRUE, scale=use.corr)  # use.corr=FALSE
    },
    Gram = {x <- rootmatrix(x)}  # get the root of Gram.
  )
  # x <- rootmatrix(pitprops)
  svdobj <- svd(x)
  sv <- svdobj$v
  totalvariance <- sum((svdobj$d)^2)
  alpha <- as.matrix(sv[, 1:K, drop=FALSE])
  beta <- alpha
  for ( j in 1:K) {
    y <- drop(x%*%alpha[,j])
    lambda2 <- lambda
    lambda1 <- para[j]
    beta[,j] <- different_beta_solver(x, y, lambda2, lambda1, sparse, solver_type)
    # beta[,j] <- solvebeta(x,y,paras=c(lambda,para[j]),sparse=sparse)
  }
  xtx <- t(x)%*%x
  temp <- beta
  normtemp <- sqrt(apply(temp^2,2,sum))
  normtemp[normtemp==0] <- 1
  temp <- t(t(temp)/normtemp)  # the normalization.
  count_k <- 0
  mdiff <- 1
  while((count_k < max.iter) & (mdiff > eps.conv)){
    count_k <- count_k + 1
    alpha <- xtx%*%beta
    z <- svd(alpha)
    alpha <- (z$u)%*%t(z$v)
    for ( j in 1:K) {
      y <- drop(x%*%alpha[,j])
      lambda2 <- lambda
      lambda1 <- para[j]
      beta[,j] <- different_beta_solver(x, y, lambda2, lambda1, sparse, solver_type)
      # beta[,j] <- solvebeta(x,y,paras=c(lambda,para[j]),sparse=sparse)
    }
    normbeta <- sqrt(apply(beta^2,2,sum))
    normbeta[normbeta==0] <- 1
    beta2 <- t(t(beta)/normbeta)
    mdiff <- convcheck(beta2,temp)
    temp <- beta2
    if(trace){
      if (count_k %% 10 == 0){cat("iterations", count_k, fill=TRUE)}
    }
  }
  normbeta <- sqrt(apply(beta^2,2,sum))
  normbeta[normbeta==0] <- 1
  beta <- t(t(beta)/normbeta)
  dimnames(beta) <- list(vn,paste("PC",1:K,sep=""))
  Z <- x%*%beta
  R <- qr.R(qr(Z))
  pev <- diag(R^2)/totalvariance
  obj <- list(
    type=type, K=K, loadings=beta, alpha=alpha,
    pev=pev, # Percentage of Explained Variance.
    var.all=totalvariance,  # var.all, variance all
    vn=vn, # variable names
    para=para, lambda=lambda)
  class(obj) <- "spca-zou"
  return(obj)
}


#' The objective function using L2 loss with three penalty functions
#'
#' The objective function using L2-loss with three penalty functions. An extension of the elastic net regression.
#'
#' @details
#' The objective function using L2 loss is defined as follows,
#' \deqn{ \frac{1}{2n} \|Y- X\beta\|^2 + \lambda_{1}p_1(\beta) + \lambda_{2}p_2(\beta) + \lambda_{3}p_3(\beta) .}
#' The three penalties are as follows,
#' \deqn{p_1(\beta) = \sum_{j=1}^p \min\{\frac{|\beta_j|}{\tau_1}, 1\} , }
#' \deqn{p_2(\beta) = \sum_{j < j', (j, j') \in E} \min\{\frac{|\beta_j - \beta_{j'}|}{\tau_2}, 1\}, }
#' \deqn{p_3(\beta) = \sum_{j=1}^p (\min\{\beta_j, 0\})^2 . }
#'
#' @param x the data matrix \eqn{X_{n\times p}}, where \eqn{n} is the number of observations,
#'   \eqn{p} is the number of features.
#' @param y the response vector \eqn{Y_{n\times 1}} with length \eqn{n}
#' @param tau_S1 \eqn{\tau_1}, the controlling parameter corresponding to \eqn{p_1(\cdot)},
#'   which determines when the small values of \eqn{|\beta_j|} will be penalized.
#' @param tau_S2 \eqn{\tau_2}, the controlling parameter corresponding to \eqn{p_2(\cdot)},
#'   which determines when the small difference values of \eqn{|\beta_j - \beta_{j'}|} will be penalized.
#' @param lambda1 \eqn{\lambda_1}, the tuning parameter corresponding to \eqn{p_1(\cdot)}
#' @param lambda2 \eqn{\lambda_2}, the tuning parameter corresponding to \eqn{p_2(\cdot)}
#' @param lambda3 \eqn{\lambda_3}, the tuning parameter corresponding to \eqn{p_3(\cdot)}
#' @param beta \eqn{\beta}, the estimation of \eqn{\beta}
#' @param Bjj a matrix of \eqn{p\times p} with element
#'   \eqn{\beta_{jj'} = \beta_{j}-\beta_{j'} } .
#' @param SF the \eqn{\mathcal{F}} set, \eqn{p}-length vector of indicator 0-1.
#'   Its value is 1 if \eqn{|\beta_j| \leq \tau_1 }. Otherwise 0.
#' @param SFc the \eqn{\mathcal{F}^{c}} set, \eqn{p}-length vector of indicator 0-1.
#'   Its value is 1 if  \eqn{|\beta_j| > \tau_1 }. Otherwise 0.
#' @param SE the \eqn{\mathcal{E}} set, a matrix of \eqn{p\times p} with indicator 0-1.
#'   Its value is 1 if \eqn{|\beta_{j}-\beta_{j'}| \leq \tau_2}. Otherwise 0.
#'   Note if \eqn{\tau_2=0, j=j'}, then \eqn{0<=0 } is true, the diagonal of \eqn{\mathcal{E}} is 1.
#'   We should set `SEc <- (1-SE)`, `SE <- SE - diag(p)` after the calculation of \eqn{\mathcal{E}}.
#' @param SEc the \eqn{\mathcal{E}^{c}} set, a matrix of \eqn{p\times p} with indicator 0-1.
#' @param SN  the \eqn{\mathcal{N}} set, \eqn{p}-length vector of indicator 0-1.
#'   Its value is 1 when \eqn{\beta_j < 0}, corresponding to \eqn{\min(\beta_j, 0)}.
#' @return the value of the objective function
### default vector in R is formed as p*1, but shown as num[1:p].
ObjFunL2 <- function(
  x, y, tau_S1, tau_S2, lambda1, lambda2, lambda3, beta, Bjj, SF, SFc, SE, SEc, SN){
  ### The objective function of the original problem in Eq.(10), with given beta.
  ### Once beta is given, the SF, SFc, SE, SEc, SN can be calculated.
  ### In the complete graph (with C_p^2 edges), j < j', so only keep the upper Triangular Part of a Matrix.
  ### lower.tri(A) the diagonal is all False.
  SE[lower.tri(SE)] <- 0
  SEc[lower.tri(SEc)] <- 0
  s1 <- mean((y - x %*% beta)^2)/2
  # s1 <- mean(((y - x %*% beta)^2 + s_tol)^(1/2))
  # s1 <- mean(abs(y - x %*% beta))
  pen1 <- lambda1*sum(SF*(abs(beta/tau_S1))) + lambda1*sum(SFc)
  pen2 <- lambda2*sum(SE*(abs(Bjj/tau_S2))) + lambda2*sum(SEc)
  pen3 <- lambda3*sum(SN*(beta^2))
  output <- s1+pen1+pen2+pen3
  return(output)
}


#' The Feature Selection and Grouping Function FSGFun with L2 loss objective.
#'
#' The FSGFun with L2 loss is a subproblem of the FGSPCA with L2 loss.
#'
#' The parameters of `x`, `y`, `beta`, `lambda1`, `lambda2`, `lambda3` are in inherit from `ObjFunL2`.
#' @inheritParams ObjFunL2
#'
#' @param tau_S \eqn{\tau}, a global \eqn{\tau}, which is assigned to \eqn{\tau_1=\tau_2=\tau .}
#' @param v the initial value of the Lagrange multiplier, with default `v=1.0`.
#' @param c the acceleration constant to speed up the convergence procedure, default `c=1.02`.
#' @param iter_m_max the maximum number of the outer iterations (i.e. the m-iteration), default `iter_m_max=100`
#' @param iter_k_max the maximum number of the inner iterations (i.e. the k-iteration), default `iter_k_max=50`
#' @param condition_tol the conditional tolerance of both the outer and inner iterations, default `condition_tol=1e-5`
#' @param nnConstraint Boolean, indicating the non-negative constraint is true or false, default `FALSE`
#' @param sparseTruncated Boolean, indicating whether to use the truncated L1 penalty or not for sparsity,
#'   default `TRUE`
#' @param loss_return Boolean, indicating whether return the loss or not, default `FALSE`
#' @return the solution \eqn{\beta} or with the corresponding loss if `loss_return=TRUE`
FSGFunL2 <- function(
  x, y, beta, tau_S, lambda1, lambda2, lambda3, v=1.0, c=1.02,
  iter_m_max=100, iter_k_max=50, condition_tol=1e-5,
  nnConstraint=FALSE, sparseTruncated=TRUE, loss_return=FALSE){
  # beta=beta_init; tau_S=0.1; lambda1=0.1; lambda2=0.1; lambda3=10; v=1.0; c=1.02;
  # iter_m_max=100; iter_k_max=50; condition_tol=1e-5
  # ## iter_m_max for the outer loop (index by m), iter_k_max for the inner loop (by count_k).
  # n <- length(y)
  p <- length(beta)
  Ones <- as.matrix(rep(1, p), p, 1)  # p by 1 matrix
  # OnesN <- as.matrix(rep(1, n), n, 1)  # n by 1
  ObjS <- vector()
  tau_S1 <- tau_S
  tau_S2 <- tau_S
  for(m in 1:iter_m_max) {
    # ## SE, SEc, SF, SFc, SN are indicator matrices or vectors of 0-1.
    # ## 1 indicates 'yes', 0 means 'no'.
    SE <- matrix(0, p, p)  # the complete graph.
    SF <- SFc <- rep(0, p)  # the index of properties indicators with beta.
    SF[abs(beta) <= tau_S1] <- 1  # <= or <
    SFc[abs(beta) > tau_S1] <- 1
    Bjj <- beta %*% t(Ones) - Ones %*% t(beta)  # NICE! Clever. Bjj[i,k] = beta[i] - beta[j]
    SE[abs(Bjj) <= tau_S2] <- 1  # note if tau_S=0, 0<=0 is true, the diagonal is 1.
    SEc <- (1-SE)  # SE[abs(Bjj) > tau_S] <- 1
    SE <- SE - diag(p)  # remove the diagonal effect, to make the diagonal of SE with 0.
    ## this is the non-negative constraint.
    SN <- rep(0, p)
    SN[beta < 0] <- 1
    # # ## this is the L2 constraint.
    # SN <- rep(1, p)  # all the ones. lambda_3 * \|\beta\|^2
    # ### SF, for min(abs(beta)/tau_S, 1), which is stored in two vectors with 0 and 1.
    # ### SF, with codes, when abs(beta) < tau_S, the value is 1; otherwise 0.
    # ### SE, min(abs(Bjj)/tau_S, 1)
    # ### SE, Bjj is a p*p matrix, which specifies the Bjj[i,k] = beta[i]-beta[k].
    # # note, if tau_S=0, 0<0 is false, then the diagonal is 0 with SE[abs(Bjj) < tau_S] <- 1.
    # # note, if tau_S=0, 0<=0 is true, then the diagonal is 1 with SE[abs(Bjj) < tau_S] <- 1.
    # ### SN, min(beta, 0), min(beta[j], 0). The value is 1 when beta < 0.
    if (sparseTruncated == FALSE) {
      # ## method 1 - If we want to use the same structure.
      # ## If we set tau_S =1, then it causes no effect when lambda2=0; but
      # ## when lambda2!=0, the feature grouping will be affected.
      tau_S1 <- 1
      SF <- rep(1, p)
      SFc <- rep(0, p)
      # pen1 <- lambda1*sum(SF*(abs(beta/tau_S))) + lambda1*sum(SFc)
      # # ## method 2 - the straight forward way.
      # pen1 <- lambda1*sum(abs(beta))
    }
    if (nnConstraint == FALSE) {  # no non-negative constraint.
      # ## this is the L2 constraint.
      SN <- rep(1, p)  # all the ones. lambda_3 * \|\beta\|^2
    }
    # ## ObjS is only used for checking the convergence. ObjS should be Monotonically decreasing.
    ObjS[m] <- ObjFunL2(x, y, tau_S1, tau_S2, lambda1, lambda2, lambda3, beta, Bjj, SF, SFc, SE, SEc, SN)
    # ## For the inner loop for updating the beta^{m, t} to get  beta^{m} (beta^{m, t*}).
    tau <- matrix(0, p, p)
    # ## v <- 1
    condition_k <- 1e4  # error, initialized as a very large number.
    count_k <- 0
    beta_m <- beta
    while (condition_k > condition_tol & count_k <= iter_k_max){
      beta_old <- beta
      ## update beta, with the raw SE, not the upper triangle SE. alpha is a vector.
      alpha <- colMeans(x^2) + 2*lambda3*SN + v*(p*rowMeans(SE))  # rowSums(SE)
      tau_sum <- rowMeans(tau*SE)*p  # same as rowSums(tau * SE)!!!
      IndexT <- which(beta != 0)  # update only non-zero beta[j]. smaller set.
      # for (j in 1:p){
      for (j in IndexT){
        # beta[j]
        bj <- y - x[,-j] %*% beta[-j]  # ## SE %*% beta # for all the j=1, 2, ...p
        gamma0 <- mean(x[,j]*bj) - tau_sum[j] + v*(t(beta) + Bjj[j,]) %*% SE[j,]
        gammaA <- sign(gamma0)*max(abs(gamma0)-lambda1/tau_S, 0)
        gamma_final <- SFc[j]*gamma0 + SF[j]*gammaA
        beta[j] <- 1/alpha[j] * gamma_final
        ### SE %*% beta means for all js (1, 2, ...p), then get the [j]th.
        ### The s1, s2, s3 are all the same.
        # s1 <- SE[j, ] %*% beta
        # s2 <- t(beta) %*% SE[j,]
        # s3 <- (SE %*% beta)[j]
        # print(c(s1, s2, s3))
        ### to calculate Bjj[j,]) %*% SE[j,], the following are all the same.
        # a1 <- Bjj[j, ] %*% SE[j, ]
        # de <- Bjj * SE
        # a2 <- sum(de[j, ])
        # ce <- SE * Bjj
        # a3 <- sum(ce[j, ])
        # a4 <- SE[j, ] %*% Bjj[j, ]
        # a5 <- sum((Bjj * SE)[j, ])
        # print(c(a1, a2, a3, a4, a5))
      }
      ## update Bjj, $\beta_{jj'}
      Bjj_new <- beta %*% t(Ones) - Ones %*% t(beta)
      temp1 <- tau + v * Bjj_new
      index_m0 <- matrix(0, p, p)
      index_m0[(abs(temp1) - lambda2/tau_S) > 0] <- 1
      Bjj_1 <- v^(-1)*sign(temp1)*index_m0*(abs(temp1) - lambda2/tau_S)
      Bjj_2 <- Bjj  # keep the old version.
      Bjj <- SE*Bjj_1 + SEc*Bjj_2  # final update rule
      ## update tau and v
      tau <- tau + v*(Bjj_new - Bjj)
      v <- c * v
      count_k <- count_k + 1
      condition_k <- max(abs(beta - beta_old))  # max norm
      # condition_k <- mean((beta - beta_old)^2)  # l2 norm
    }  # end of the inner loop.
    # ## criterion 1, ObjS is only used for checking the convergence.
    if (m >= 2 && ObjS[m] >= ObjS[m-1]) break
    if (m >= 2 && ObjS[m-1] - ObjS[m] < condition_tol) break
    # ## criterion 2, use the norm to check the convergence.
    # current_tol <- max(abs(beta - beta_m))  # max norm
    # current_tol <- mean((beta - beta_m)^2)  # l2 norm
    # if (current_tol < condition_tol) break
  }  ## end of the outer loop.
  # ## plot(ObjS)
  if (loss_return) {
    return(list(beta=beta, ObjS=ObjS))
  } else {
    return(beta)
  }
  # if (loss_return) return(list(beta=beta, ObjS=ObjS))
  # if (!loss_return) return(beta)
}


#' The penalty function consisting of three parts.
#'
#' The penalty function consisting of three parts.
#'
#' The three penalties are as follows,
#' \deqn{p_1(\beta) = \sum_{j=1}^p \min\{\frac{|\beta_j|}{\tau_1}, 1\} , }
#' \deqn{p_2(\beta) = \sum_{j < j', (j, j') \in E} \min\{\frac{|\beta_j - \beta_{j'}|}{\tau_2}, 1\}, }
#' \deqn{p_3(\beta) = \sum_{j=1}^p (\min\{\beta_j, 0\})^2 . }
#'
#' @inheritParams ObjFunL2
#'
#' @param tau_S \eqn{\tau} a global \eqn{\tau}, which is assigned to \eqn{\tau_1=\tau_2=\tau .}
#' @param nnConstraint Boolean, indicating the non-negative constraint is true or false,
#'   default `FALSE`
#' @param sparseTruncated Boolean, indicating whether use the truncated L1 penalty or not for sparsity,
#'   default `TRUE`
#'
PenaltyFun <- function(
  beta, tau_S, lambda1, lambda2, lambda3,
  nnConstraint=FALSE, sparseTruncated=TRUE) {
  ### The objective function of the original problem in Eq.(3.7), with given beta.
  ### Once beta is given, the SF, SFc, SE, SEc, SN can be calculated.
  ### In the complete graph (with C_p^2 edges), j < j', so only keep the upper Triangular Part of a Matrix.
  ### lower.tri(A) the diagonal is all False.
  p <- length(beta)
  Ones <- as.matrix(rep(1, p), p, 1)  # p by 1 matrix
  Bjj <- beta %*% t(Ones) - Ones %*% t(beta)  # NICE! Clever. Bjj[i,k] = beta[i] - beta[j]
  # ## SE, SEc, SF, SFc, SN are indicator matrices or vectors of 0-1.
  # ## 1 indicates 'yes', 0 means 'no'.
  # ## this is the truncated grouping constraint.
  SE <- matrix(0, p, p)  # the complete graph.
  SE[abs(Bjj) <= tau_S] <- 1  # note if tau_S=0, 0<=0 is true, the diagonal is 1.
  SEc <- (1-SE)  # SE[abs(Bjj) > tau_S] <- 1
  SE <- SE - diag(p)  # remove the diagonal effect, to make the diagonal of SE with 0.
  # ## this is the non-negative constraint.
  SN <- rep(0, p)
  SN[beta < 0] <- 1
  # ## this is the truncated sparse constraint.
  SF <- SFc <- rep(0, p)  # the index of properties indicators with beta.
  SF[abs(beta) <= tau_S] <- 1  # <= or <
  SFc[abs(beta) > tau_S] <- 1
  ### In the complete graph (with C_p^2 edges), j < j', so only keep the upper Triangular Part of a Matrix.
  ### lower.tri(A) the diagonal is all False.
  SE[lower.tri(SE)] <- 0
  SEc[lower.tri(SEc)] <- 0
  pen1 <- lambda1*sum(SF*(abs(beta/tau_S))) + lambda1*sum(SFc)
  pen2 <- lambda2*sum(SE*(abs(Bjj/tau_S))) + lambda2*sum(SEc)
  pen3 <- lambda3*sum(SN*(beta^2))
  # ## if no truncated penalty, the direct calculation is easier.
  if (sparseTruncated == FALSE) {
    # # ## method 1 - If we want to use the same structure.
    # tau_S <- 1
    # SF <- rep(1, p)
    # SFc <- rep(0, p)
    # pen1 <- lambda1*sum(SF*(abs(beta/tau_S))) + lambda1*sum(SFc)
    # ## method 2 - the straight forward way.
    pen1 <- lambda1*sum(abs(beta))
  }
  if (nnConstraint == FALSE) {
    # # ## method 1 - If we want to use the same structure.
    # SN <- rep(1, p)  # all the ones. lambda_3 * \|\beta\|^2
    # pen3 <- lambda3*sum(SN*(beta^2))
    # ## method 2 - the straight forward way.
    # ## this is the corresponding L2 constraint.
    pen3 <- lambda3*sum(beta^2)
  }
  output <- pen1+pen2+pen3
  return(output)
}


#' Objective Function of Feature Grouping and Sparse PCA with L2 loss.
#'
#' Objective Function of Feature Grouping and Sparse PCA with L2 loss.
#'
#' The L2 loss is defined as follows
#' \deqn{ (X - X B A^T)^2 .}
#'
#' The penalty function with three parts is defined as follows,
#' \deqn{ \Psi(\beta) = \lambda_{1}p_1(\beta) + \lambda_{2}p_2(\beta) + \lambda_{3}p_3(\beta)}
#' where
#' \deqn{p_1(\beta) = \sum_{j=1}^p \min\{\frac{|\beta_j|}{\tau_1}, 1\} , }
#' \deqn{p_2(\beta) = \sum_{j < j', (j, j') \in E} \min\{\frac{|\beta_j - \beta_{j'}|}{\tau_2}, 1\}, }
#' \deqn{p_3(\beta) = \sum_{j=1}^p (\min\{\beta_j, 0\})^2 . }
#'
#' With the setting of `lambda2=0, nnConstraint=FALSE, sparseTruncated=FALSE`,
#' it corresponds to the same objective function of sparse PCA; \\
#'
#' `ObjFun_FG_S_PCA_L2(x, A, B, tau_S, lambda1, lambda2=0, lambda3=0.1, nnConstraint=FALSE, sparseTruncated=FALSE`
#'
#' With setting of `lambda2=0, nnConstraint=FALSE, sparseTruncated=TRUE`,
#' corresponds to Objective function of truncated sparse PCA. \\
#'
#' `ObjFun_FG_S_PCA_L2(x, A, B, tau_S, lambda1, lambda2=0, lambda3=0.1, nnConstraint=FALSE, sparseTruncated=TRUE`
#'
#' By tuning the different setting of `lambda2`, `nnConstraint`, `sparseTruncated`,
#' we get different combination of models.
#'
#' @inheritParams ObjFunL2
#' @inheritParams FSGFunL2
#' @inheritParams PenaltyFun
#'
#' @param A matrix \eqn{A} in the FGSPCA problem
#' @param B matrix \eqn{B} in the FGSPCA problem
#'
ObjFun_FG_S_PCA_L2 <- function(
  x, A, B, tau_S, lambda1, lambda2=0.1, lambda3=0.1,
  nnConstraint=FALSE, sparseTruncated=TRUE){
  ### The objective function of the original problem in Eq.(10), with given beta.
  k <- dim(B)[2]
  pen <- 0
  if (length(lambda1) == k) {
    for (j in 1:k) {
      beta <- B[, j]
      vpen <- PenaltyFun(
        beta, tau_S, lambda1[j], lambda2, lambda3, nnConstraint, sparseTruncated)
      pen <- pen + vpen
    }
  } else {
    for (j in 1:k) {
      beta <- B[, j]
      vpen <- PenaltyFun(
        beta, tau_S, lambda1, lambda2, lambda3, nnConstraint, sparseTruncated)
      pen <- pen + vpen
    }
  }
  s1 <- mean((x - x %*% B %*% t(A))^2)/2
  output <- s1 + pen
  return(output)
}


#' Feature Grouping and Sparse PCA with L2 loss function.
#'
#' Feature Grouping and Sparse PCA with L2 loss function.
#'
#' The main function of the FGSPCA problem with L2 loss using the alternating updating,
#'   which updates \eqn{A} while fixing \eqn{B} and updates \eqn{B} while fixing \eqn{A}.
#'
#' @inheritParams ObjFunL2
#' @inheritParams FSGFunL2
#' @inheritParams ObjFun_FG_S_PCA_L2
#'
#' @param x the data matrix \eqn{X}
#' @param B_init the initial value of matrix \eqn{B}
#' @param K the number of principal components
#' @param para a list of \eqn{\lambda_1} with its length being `K`
#' @param type  type of \eqn{X}, which can take values of `type=c("predictor", "Gram")`.
#'    If `type="Gram"` the model should deal with the root matrix of \eqn{X}.
#'    If `type="predictor"` the model should directly deal with the matrix \eqn{X}.
#' @param use.corr `FALSE` (default).
#'    When `type="predictor"` we need to scale the data, `scale=use.corr` .
#' @param max.iter `200` (default) the max iteration of the `A-B` updating procedure
#' @param trace `TRUE` (default) whether to do the print or not
#'   (print the number of the current iteration and the error at the iteration)
#' @param eps.conv \eqn{1e-3} (default) the convergence criterion for the `A-B` updating procedure
#' @importFrom utils tail
#'
#' @return pev Percentage of Explained Variance
#' @return var.all the total variance
#' @return loadings the final normalized loadings of \eqn{B}
#' @return Alpha the final eqn{A} from the last update
#' @return ab_errors the `A-B` updating errors
#' @return type (same as the input) the type of matrix \eqn{X}
#' @return K (same as the input) the number of principal components
#' @return para (same as the input) a list of \eqn{\lambda_1} with its length being `K`
#' @return lambda2 (same as the input) \eqn{\lambda_2}
#' @return lambda3 (same as the input) \eqn{\lambda_3}
#' @return vn variable names
#' @return Mloss a matrix of loss with shape of \eqn{i\times K},
#'   where \eqn{i} is the number of iterations.
#' @return Gloss global loss
#' @return Sloss sum loss
#' @return Tloss tail loss
#' @return LossList with each element being a list of the objective values of each subproblem.
#'   \eqn{K*(i-1) + j} .
#'
#' @examples
#' library(elasticnet)  # to use the data "pitprop"
#' library(sparsepca)
#' data(pitprops)
#' # ## elasticnet::spca is the Sparse PCA of Zou (2006).
#' K <- 6
#' para <- c(0.06, 0.16, 0.1, 0.5, 0.5, 0.5)
#' out1 <- elasticnet::spca(
#'   pitprops, K=6, type="Gram", sparse="penalty", trace=TRUE,
#'   para=c(0.06,0.16,0.1,0.5,0.5,0.5))
#'
#' K <- k <- 6
#' X <- rootmatrix(pitprops)
#' rspca.results <- sparsepca::rspca(X, k, alpha=1e-3, beta=1e-3, center=TRUE, scale=FALSE)
#'
#' K <- 6
#' B_init <- rspca.results$loadings
#' tau_S=0.05; lambda1=0.01; lambda2=0.01; lambda3=0.1
#' para <- rep(lambda1, K)
#'
#' out <- FGSPCA(
#'   pitprops, B_init, K, para, type="Gram",
#'   tau_S=tau_S, lambda2=lambda2, lambda3=lambda3)
#' NB <- out$loadings  # the normalized loadings of FGSPCA
#'
#' (fgspca.pev <- round(out$pev *100, 3))
#' (fgspca.cpev <- round(cumsum(out$pev) * 100, 3))
#'
#' @export
FGSPCA <- function(
  x, B_init, K, para, type="predictor", use.corr=FALSE,
  max.iter=200, trace=TRUE, eps.conv=1e-3,
  tau_S=0.05, lambda2=0.1, lambda3=0.1, v=1.0, c=1.02,
  iter_m_max=100, iter_k_max=50, condition_tol=1e-5,
  nnConstraint=FALSE, sparseTruncated=TRUE) {
  # type="Gram"; use.corr=FALSE;
  # max.iter=200; trace=TRUE; eps.conv=1e-3;
  # tau_S=0.05; lambda2=0.1; lambda3=0.1; v=1.0; c=1.02;
  # iter_m_max=100; iter_k_max=50; condition_tol=1e-5;
  # nnConstraint=FALSE; sparseTruncated=TRUE;
  vn <- dimnames(x)[[2]]
  x <- switch(
    type,
    predictor = {
      n <- dim(x)[1]
      p <- dim(x)[2]
      if (n/p >= 1e8) {
        cat("You may wish to restart and use a more efficient way \n")
        cat("let the argument x be the sample covariance/correlation matrix and set type=Gram \n")
      }
      x <- scale(x, center=TRUE, scale=use.corr)  # use.corr=FALSE
    },
    Gram = {x <- rootmatrix(x)}  # get the root of Gram.
  )
  # ## various types of losses
  Gloss <- vector(mode = "numeric")  # global loss
  Sloss <- vector(mode = "numeric")  # sum loss
  Tloss <- vector(mode = "numeric")  # tail loss
  Mloss <- vector(mode = "numeric")  # matrix loss, i*K
  LossList <- list()
  # ## various types of losses
  svdobj <- svd(x)
  sv <- svdobj$v
  totalvariance <- sum((svdobj$d)^2)
  Alpha <- as.matrix(sv[, 1:K, drop=FALSE])
  Beta <- Alpha
  tloss <- vector(mode = "numeric")
  sloss <- 0
  i <- 1
  for ( j in 1:K) {
    y <- drop(x %*% Alpha[,j])
    beta_init <- B_init[, j]
    res.mat <- FSGFunL2(
      x, y, beta=beta_init, tau_S, lambda1=para[j], lambda2, lambda3,
      v, c, iter_m_max, iter_k_max, condition_tol,
      nnConstraint, sparseTruncated, loss_return=TRUE)
    Beta[, j] <- res.mat$beta
    LossList[[K*(i-1) + j]] <- res.mat$ObjS
    tloss[j] <- tail(res.mat$ObjS, n=1)
    sloss <- sloss + sum(res.mat$ObjS)
  }
  # ### deal with the loss
  Mloss <- rbind(Mloss, tloss)
  Sloss[i] <- sloss
  Tloss[i] <- sum(tloss)
  Gloss[i] <- ObjFun_FG_S_PCA_L2(
    x, Alpha, Beta, tau_S, para, lambda2, lambda3, nnConstraint, sparseTruncated)
  # ### The main updating of A-B and
  xtx <- t(x)%*%x
  temp <- Beta    # the normalization.
  normtemp <- sqrt(apply(temp^2,2,sum))
  normtemp[normtemp==0] <- 1
  temp <- t(t(temp)/normtemp)  # the normalization.
  mdiff <- 1  # the initial value of mdiff.
  ab_errors <- vector(mode="numeric")
  for (i in 2:max.iter) {
    Alpha <- xtx %*% Beta
    sres <- svd(Alpha)
    Alpha <- (sres$u) %*% t(sres$v)
    for ( j in 1:K) {
      y <- drop(x %*% Alpha[,j])
      beta_init <- B_init[, j]
      mat <- FSGFunL2(
        x, y, beta=beta_init, tau_S, lambda1=para[j], lambda2, lambda3,
        v, c, iter_m_max, iter_k_max, condition_tol,
        nnConstraint, sparseTruncated, loss_return=TRUE)
      Beta[, j] <- mat$beta
      LossList[[K*(i-1) + j]] <- mat$ObjS
      tloss[j] <- tail(mat$ObjS, n=1)
      sloss <- sloss + sum(mat$ObjS)
    }
    # ### deal with the loss
    Mloss <- rbind(Mloss, tloss)
    Sloss[i] <- sloss
    Tloss[i] <- sum(tloss)
    Gloss[i] <- ObjFun_FG_S_PCA_L2(
      x, Alpha, Beta, tau_S, para, lambda2, lambda3, nnConstraint, sparseTruncated)
    # ### check the convergence
    normBeta <- sqrt(apply(Beta^2, 2, sum))
    normBeta[normBeta==0] <- 1
    Beta2 <- t(t(Beta)/normBeta)
    mdiff <- convcheck(Beta2, temp)  # convergence check.
    ab_errors[i] <- mdiff
    temp <- Beta2
    if (mdiff <= eps.conv) {break}
    if (trace) {
      if ( i %% 10 == 0 ){cat("iterations", i, "Error", mdiff, fill=TRUE)}
    }
  }
  # ## for the last step, normalization.
  normBeta <- sqrt(apply(Beta^2, 2, sum))
  normBeta[normBeta==0] <- 1
  Beta <- t(t(Beta)/normBeta)
  dimnames(Beta) <- list(vn, paste("PC", 1:K, sep=""))
  Z <- x %*% Beta
  R <- qr.R(qr(Z))
  pev <- diag(R^2)/totalvariance
  obj <- list(
    type=type, K=K, para=para,
    pev=pev, # Percentage of Explained Variance.
    var.all=totalvariance,  # var.all, variance all
    vn=vn, # variable names
    lambda2=lambda2, lambda3=lambda3,
    Mloss=Mloss, Sloss=Sloss, Tloss=Tloss, Gloss=Gloss,
    LossList=LossList,
    ab_errors=ab_errors,
    loadings=Beta, Alpha=Alpha
  )
  class(obj) <- "fgspca"
  return(obj)
}

#' Adjusted Total Variance by Zou. et.al.(2006)
#'
#' To calculate the Total Variance Explained by the new normalized loadings NB.
#' The loadings of the ordinary principal components are uncorrelated and orthogonal,
#' where the orthogonal property, and the uncorrelated property are satisfied together.
#' @param x the data matrix by \eqn{n\times p}. It should not be the sample covariance matrix,
#'   nor the correlation matrix.
#' @param NB normalized matrix of loading vectors
#' @return a list of `c(pev_svd, pev_tr, pev_spca, cum_svd, cum_tr, cum_spca)`.
#' @export
#' @references
#' \itemize{
#'  \item Zou, Hui, Trevor Hastie, and Robert Tibshirani. "Sparse principal component analysis."
#'    Journal of computational and graphical statistics 15.2 (2006): 265-286.
#' }
total_variance_explained <- function(x, NB){
  # ## the traditional method of percentage of explained variance.
  # ## Percentage of Explained Variance.
  # ## x must not be the squared matrix.
  sres <- svd(x)
  totalvariance <- sum((sres$d)^2)
  pev_svd <- ((sres$d)^2) / totalvariance
  cum_svd <- cumsum((sres$d)^2) / totalvariance
  # ## If Z are uncorrelated, trace(Z^TZ).
  Z <- x %*% NB
  zz <- t(Z) %*% Z
  s1 <- diag(zz)
  pev_tr <- s1 / totalvariance
  cum_tr <- cumsum(s1) / totalvariance
  # ## QR decomposition of Z
  R <- qr.R(qr(Z))
  pev_spca <- diag(R^2)/totalvariance
  cum_spca <- cumsum(diag(R^2)) /totalvariance
  result <- list(
    pev_svd=pev_svd*100, pev_tr=pev_tr*100, pev_spca=pev_spca*100,
    cum_svd=cum_svd*100, cum_tr=cum_tr*100, cum_spca=cum_spca*100
  )
  res <- lapply(result, round, 2)
  return(res)
}


# a <- NB[, 1]
numSparsity <- function(a, tol=1e-2) {sum(abs(a) < tol) } # which(abs(a) < tol)
numNonzero <- function(a, tol=1e-2) {sum(abs(a) > tol) } # which(abs(a) > tol)

# apply(NB, 2, numSparsity)
# apply(NB, 2, numNonzero)
# sum(apply(NB, 2, numNonzero))
# apply(NB, 2, numGroup)


numGroup <- function(a, tol=1e-2, verbose=FALSE) {
  # evaluating the grouping effect of a.
  a <- sort(a)
  da <- diff(a)
  ngroups <- 1
  igroup <- c()
  gap <- 0
  for (i in 1:length(da)) {
    d <- da[i]
    if ( (gap + d) > tol) {
      ngroups <- ngroups + 1
      igroup <- c(igroup, i)
      gap <- 0
    } else {
      gap <- gap + d
    }
  }
  # # ## to check if the groups include 0 group as a special case.
  # # ## maybe not exactly zero, but in the neighberhood of zero.
  zeroGroup = FALSE
  n <- length(a)
  if (ngroups == 1) {
    if (max(abs(a)) < tol) {zeroGroup=TRUE}
  }
  if (ngroups ==2) {
    if (max(abs(a[1:igroup[1]] - 0)) < tol) {zeroGroup=TRUE}
    if (max(abs(a[(igroup[1]+1):n] - 0)) < tol) {zeroGroup=TRUE}
  }
  if (ngroups > 2) {
    if (max(abs(a[1:igroup[1]] - 0)) < tol) {zeroGroup=TRUE}
    for (i in 1:(ngroups-2) ) {
      if (max(abs(a[(igroup[i]+1):igroup[i+1]])) < tol ) {zeroGroup=TRUE}
    }
    if (max(abs(a[(igroup[ngroups-1]+1):n] - 0)) < tol ) {zeroGroup=TRUE}
  }
  if (zeroGroup) {ngroups <- ngroups - 1}
  # ## For example, igroup=c(2,6,8,12), which means, 1-2 is a group;
  # ## 1-2, 3-6, 7-8, 9-12, 13-13 are all the groups.
  res <- list(ngroups = ngroups, igroup=igroup, zeroGroup=zeroGroup)
  if (verbose == TRUE) return(res) else return(ngroups)
  # numGroup(a, tol=0.05)
  # numGroup(a, tol=0.2)
}



