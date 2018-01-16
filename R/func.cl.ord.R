#' Composite Likelihood Calculation for Spatial Ordinal Data
#'
#' \code{func.cl.ord} calculates the composite log-likelihood for spatial ordered probit models.
#'
#' @param vec.yobs a vector of observed responses for all N sites.
#' @param mat.X regression (design) matrix, including intercepts.
#' @param mat.lattice a data matrix containing geographical information of sites. The \emph{i} th row constitutes a set of geographical coordinates.
#' @param radius weight radius.
#' @param n.cat number of categories, at least 2.
#' @param vec.par a vector of parameters consecutively as follows: a series of cutoffs (excluding -Inf, 0 and Inf) for latent responses, a vector of covariate parameters, a parameter 'sigmasq' modeling covariance matrix, 0<=sigmasq<=1, and a parameter 'rho' reflecting spatial correlation, abs(rho)<=1.
#' @return \code{func.cl.ord} returns a list: number of categories, sum of weights, composite log-likelihood, a vector of scores, and a matrix of first-order partial derivatives for \code{vec.par}.
#' @examples
#'
#' # True parameter
#' vec.cutoff <- 2; vec.beta <- c(1, 2, 1, 0, -1); sigmasq <- 0.8; rho <- 0.6; radius <- 5
#' vec.par <- c(vec.cutoff, vec.beta, sigmasq, rho)
#'
#' # Coordinate matrix
#' n.cat <- 3; n.lati <- 30; n.long <- 30
#' n.site <- n.lati * n.long
#' mat.lattice <- cbind(rep(1:n.lati, n.long), rep(1:n.long, each=n.lati))
#' mat.dist <- as.matrix(dist(mat.lattice, upper=TRUE, diag=TRUE))
#' mat.cov <- sigmasq * rho^mat.dist
#'
#' set.seed(1228)
#' # Generate regression (design) matrix with intercept
#' mat.X <- cbind(rep(1, n.site),scale(matrix(rnorm(n.site*(length(vec.beta)-1)),nrow=n.site)))
#' vec.Z <- t(chol(mat.cov)) %*% rnorm(n.site) + mat.X %*% vec.beta
#' vec.epsilon <- diag(sqrt(1-sigmasq), n.site) %*% rnorm(n.site)
#' vec.ylat <- as.numeric(vec.Z + vec.epsilon)
#'
#' # Convert to the vector of observation
#' vec.yobs <- func.obs.ord(vec.ylat, vec.alpha=c(-Inf,0,vec.cutoff,Inf))
#'
#' # Using func.cl.ord()
#' ls <- func.cl.ord(vec.yobs, mat.X, mat.lattice, radius, n.cat, vec.par)
#' ls$log.lkd
#'
#' @references Feng, Xiaoping, Zhu, Jun, Lin, Pei-Sheng, and Steen-Adams, Michelle M. (2014) Composite likelihood Estimation for Models of Spatial Ordinal Data and Spatial Proportional Data with Zero/One values. \emph{Environmetrics} 25(8): 571--583.

func.cl.ord <- function(vec.yobs, mat.X, mat.lattice, radius, n.cat, vec.par) {

  # input check
  stopifnot(radius>=0, n.cat>=2,
            vec.par[length(vec.par)-1]>=0, vec.par[length(vec.par)-1]<=1,
            vec.par[length(vec.par)]>=-1, vec.par[length(vec.par)]<=1)
  if (n.cat>2) stopifnot(length(vec.par)-NCOL(mat.X)-2>0)

  vec.dist <- c(dist(mat.lattice))
  # vec.dist: a vector of pairwise euclidean distances between
  #           sites with length equal to N*(N-1)/2
  mat.dist <- as.matrix(dist(mat.lattice, upper=TRUE, diag=TRUE))
  # mat.dist: an N by N matrix of pairwise euclidean distances between sites
  vec.idx.dist <- which(vec.dist<=radius)
  vec.dist.select <- vec.dist[vec.idx.dist]
  # vec.idx.dist:    a vector of selected indices of vector vec.dist with the
  #                  corresponding distances at most the radius
  # vec.dist.select: a vector of selected pairwise distances at most the radius
  weight.sum <- length(vec.idx.dist)
  # weight.sum: sum of weights (W_N as in (6))
  mat.idx.pair <- combn(1:length(vec.yobs), 2)[, vec.idx.dist]
  # mat.idx.pair: a matrix each col of which indexes Sites i and i' with distance
  #               at most the radius
  mat.X.1 <- mat.X[mat.idx.pair[1,],]; mat.X.2 <- mat.X[mat.idx.pair[2,],]
  # mat.X.1/mat.X.2: subsets of mat.X corresponding to the first/second rows
  # of mat.idx.pair
  vec.yobs.1 <- vec.yobs[mat.idx.pair[1,]]; vec.yobs.2 <- vec.yobs[mat.idx.pair[2,]]
  # vec.yobs.1/vec.yobs.2: subsets of vec.yobs corresponding to the first/second
  # rows in mat.idx.pair
  mat.basis <- diag(n.cat+1)
  # mat.basis: an identity matrix containing (n.cat+1) bases in R^(n.cat+1)
  inf <- 2e+50 # Note: inf must not be excessively large
  if (n.cat==2) {
    vec.alnbe <- c(-inf, 0, inf, vec.par[1:NCOL(mat.X)])
  } else {
    vec.alnbe <- c(-inf,0,vec.par[1:(n.cat-2)],inf,vec.par[(n.cat-1):(length(vec.par)-2)])
  }
  # vec.alnbe: a vector concatenating alphas and regression coefficients
  func.xi <- function(mu, rho.ab) (mu[,2] - rho.ab*mu[,1]) / sqrt(1 - rho.ab^2)
  # func.xi: see the note below
  func.der.alnbe <- function(mu, rho.ab, X.a, X.b)
    dnorm(mu[,1])*pnorm(func.xi(mu, rho.ab))*X.a +
    dnorm(mu[,2])*pnorm(func.xi(mu[,2:1], rho.ab))*X.b
  # func.der.alnbe: see the note below

  dbivnorm <- function(mu, rho.ab) # bivnorm: density of bivariate normal distribution
    1 / (2*pi*sqrt(1 - rho.ab^2)) *
    exp(-0.5*(mu[,1]^2 + mu[,2]^2 - 2*rho.ab*mu[,1]*mu[,2]) / (1 - rho.ab^2))

  # computation for each pair of ordered categories
  for (a in 0:(n.cat-1)) {
    for(b in 0:(n.cat-1)) {
      vec.idx.ab <- which(vec.yobs.1==a & vec.yobs.2==b)
      # vec.idx.ab: a vector of indices for a certain pair of categories
      n.ab <- length(vec.idx.ab)

      mat.X.a <- mat.X.1[vec.idx.ab,]; mat.X.b <- mat.X.2[vec.idx.ab,]

      mat.X.alo <- cbind(matrix(rep(mat.basis[a+1,], each=n.ab), nrow=n.ab), -mat.X.a)
      mat.X.aup <- cbind(matrix(rep(mat.basis[a+2,], each=n.ab), nrow=n.ab), -mat.X.a)
      mat.X.blo <- cbind(matrix(rep(mat.basis[b+1,], each=n.ab), nrow=n.ab), -mat.X.b)
      mat.X.bup <- cbind(matrix(rep(mat.basis[b+2,], each=n.ab), nrow=n.ab), -mat.X.b)

      mat.mu.uu <- cbind(mat.X.aup%*%vec.alnbe, mat.X.bup%*%vec.alnbe)
      mat.mu.ll <- cbind(mat.X.alo%*%vec.alnbe, mat.X.blo%*%vec.alnbe)
      mat.mu.ul <- cbind(mat.X.aup%*%vec.alnbe, mat.X.blo%*%vec.alnbe)
      mat.mu.lu <- cbind(mat.X.alo%*%vec.alnbe, mat.X.bup%*%vec.alnbe)

      vec.rho.ab <- vec.par[length(vec.par)-1] *
        vec.par[length(vec.par)] ^ vec.dist.select[vec.idx.ab]

      vec.pmf.ab <- pbivnorm::pbivnorm(x=mat.mu.uu, rho=vec.rho.ab) +
        pbivnorm::pbivnorm(x=mat.mu.ll, rho=vec.rho.ab) -
        pbivnorm::pbivnorm(x=mat.mu.ul, rho=vec.rho.ab) -
        pbivnorm::pbivnorm(x=mat.mu.lu, rho=vec.rho.ab)
      # vec.pmf.ab: a vector of PMFs as in (8)

      mat.deriv.alnbe <- func.der.alnbe(mat.mu.uu, vec.rho.ab, mat.X.aup, mat.X.bup) +
        func.der.alnbe(mat.mu.ll, vec.rho.ab, mat.X.alo, mat.X.blo) -
        func.der.alnbe(mat.mu.ul, vec.rho.ab, mat.X.aup, mat.X.blo) -
        func.der.alnbe(mat.mu.lu, vec.rho.ab, mat.X.alo, mat.X.bup)
      # mat.deriv.alnbe: a matrix of partial derivs wrt alphas and betas as in (12)

      vec.deriv.rho.ab <- dbivnorm(mat.mu.uu, vec.rho.ab) +
        dbivnorm(mat.mu.ll, vec.rho.ab) -
        dbivnorm(mat.mu.ul, vec.rho.ab) -
        dbivnorm(mat.mu.lu, vec.rho.ab)
      # vec.deriv.rho.ab: a vector of partial derivs wrt spatial corr as in (13)

      vec.deriv.sigmasq <- vec.deriv.rho.ab *
        vec.par[length(vec.par)] ^ vec.dist.select[vec.idx.ab]
      # vec.deriv.sigmasq: a vector of partial derivs wrt sigma square

      vec.deriv.rho <- vec.deriv.rho.ab * vec.dist.select[vec.idx.ab] *
        vec.par[length(vec.par)-1] *
        vec.par[length(vec.par)]^(vec.dist.select[vec.idx.ab]-1)
      # vec.deriv.rho: a vector of partial derivs wrt rho

      if (a==0 & b==0) { # pair(0,0) as the base for summation
        log.lkd <- sum(log(vec.pmf.ab)) / weight.sum
        # log.lkd: composite log-likelihood function
        mat.derivs <- unname(cbind(mat.deriv.alnbe,
                                   vec.deriv.sigmasq,vec.deriv.rho))/vec.pmf.ab
        # mat.derivs: a matrix of partial derivs for parameters including -Inf, 0, Inf
      } else {
        log.lkd <- log.lkd + sum(log(vec.pmf.ab)) / weight.sum
        mat.derivs <- unname(rbind(mat.derivs,
                                   cbind(mat.deriv.alnbe,
                                         vec.deriv.sigmasq,vec.deriv.rho)/vec.pmf.ab))
      }
    }
  }

  if (n.cat==2) {
    mat.score <- mat.derivs[, (n.cat+2):(n.cat+NCOL(mat.X)+3)]
    # mat.score: a matrix of partial derivs for parameters to be estimated
  } else {
    mat.score <- mat.derivs[, c(3:n.cat, (n.cat+2):(n.cat+NCOL(mat.X)+3))]
  }

  vec.score <- colSums(mat.score) / weight.sum
  # vec.score: a vector of score functions as in (10)
  return(list(n.cat=n.cat, weight.sum=weight.sum, log.lkd=log.lkd,
              vec.score=vec.score, mat.score=mat.score))
}
