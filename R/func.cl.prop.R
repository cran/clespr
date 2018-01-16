#' Composite Likelihood Calculation for Spatial Proportional Data
#'
#' \code{func.cl.prop} calculates the composite log-likelihood for spatial Tobit models.
#'
#' @param vec.yobs a vector of observed responses for all N sites.
#' @param mat.X regression (design) matrix, including intercepts.
#' @param mat.lattice a data matrix containing geographical information of sites. The \emph{i} th row constitutes a set of geographical coordinates.
#' @param radius weight radius.
#' @param vec.par a vector of parameters consecutively as follows: a cutoff point for latent responses, a vector of covariate parameters, a parameter 'sigmasq' modeling covariance matrix, 0<=sigmasq<=1, and a parameter 'rho' reflecting spatial correlation, abs(rho)<=1.
#' @return \code{func.cl.prop} returns a list of sum of weights, composite log-likelihood, a vector of scores, and a matrix of first-order partial derivatives for \code{vec.par}.
#'
#' @examples
#'
#' # True parameter
#' alpha <- 4; vec.beta <- c(1, 2, 1, 0, -1); sigmasq <- 0.8; rho <- 0.6; radius <- 5
#' vec.par <- c(alpha, vec.beta, sigmasq, rho)
#'
#' # Coordinate matrix
#' n.lati <- 30; n.long <- 30
#' n.site <- n.lati * n.long
#' mat.lattice <- cbind(rep(1:n.lati, n.long), rep(1:n.long, each=n.lati))
#' mat.dist <- as.matrix(dist(mat.lattice, upper=TRUE, diag=TRUE))
#' mat.cov <- sigmasq * rho^mat.dist
#'
#' set.seed(1228)
#'
#' # Generate regression (design) matrix with intercept
#' mat.X <- cbind(rep(1, n.site),scale(matrix(rnorm(n.site*(length(vec.beta)-1)),nrow=n.site)))
#' vec.Z <- t(chol(mat.cov)) %*% rnorm(n.site) + mat.X %*% vec.beta
#' vec.epsilon <- diag(sqrt(1-sigmasq), n.site) %*% rnorm(n.site)
#' vec.ylat <- as.numeric(vec.Z + vec.epsilon)
#'
#' # Convert to the vector of observation
#' vec.yobs <- func.obs.prop(vec.ylat, alpha=alpha)
#'
#' # Use func.cl.prop()
#' ls <- func.cl.prop(vec.yobs, mat.X, mat.lattice, radius, vec.par)
#' ls$log.lkd
#' @references Feng, Xiaoping, Zhu, Jun, Lin, Pei-Sheng, and Steen-Adams, Michelle M. (2014) Composite likelihood Estimation for Models of Spatial Ordinal Data and Spatial Proportional Data with Zero/One values. \emph{Environmetrics} 25(8): 571--583.

func.cl.prop <- function(vec.yobs, mat.X, mat.lattice, radius, vec.par) {

  # input check
  stopifnot(radius>=0,
            vec.par[length(vec.par)-1]>=0, vec.par[length(vec.par)-1]<=1,
            vec.par[length(vec.par)]>=-1, vec.par[length(vec.par)]<=1)


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
  vec.mu <- mat.X %*% vec.par[2:(length(vec.par)-2)]
  mat.mu.pair <- t(matrix(vec.mu[mat.idx.pair], nrow=2))
  # mat.mu.pair: a matrix each row of which denotes a pair of inner products of
  #              covariates and coefficients

  mat.derivs <- matrix(0,nrow=weight.sum,ncol=length(vec.par))
  # mat.derivs: a matrix of partial derivs for parameters of interest

  func.xi <- function(mu, rho.ab) (mu[,2] - rho.ab*mu[,1]) / sqrt(1 - rho.ab^2)

  func.der.beta <- function(mu, rho.ab, X.a, X.b)
    dnorm(mu[,1])*pnorm(func.xi(mu, rho.ab))*X.a +
    dnorm(mu[,2])*pnorm(func.xi(mu[,2:1], rho.ab))*X.b


  dbivnorm <- function(mu, rho.ab) # bivnorm: density of bivariate normal distribution
    1 / (2*pi*sqrt(1 - rho.ab^2)) *
    exp(-0.5*(mu[,1]^2 + mu[,2]^2 - 2*rho.ab*mu[,1]*mu[,2]) / (1 - rho.ab^2))

  # computation for each pair of categories
  for (a in 0:2) {
    for(b in 0:2) {
      if (a!=2 & b!=2) {
        vec.idx.ab <- which(vec.yobs.1==a & vec.yobs.2==b)
      } else if (a==2 & b!=2) {
        vec.idx.ab <- which(vec.yobs.1>0 & vec.yobs.1<1 & vec.yobs.2==b)
      } else if (a!=2 & b==2) {
        vec.idx.ab <- which(vec.yobs.1==a & vec.yobs.2>0 & vec.yobs.2<1)
      } else {
        vec.idx.ab <- which(vec.yobs.1>0 & vec.yobs.1<1 & vec.yobs.2>0 & vec.yobs.2<1)
      }
      # vec.idx.ab: a vector of indices for a certain pair of categories

      mat.X.a <- mat.X.1[vec.idx.ab,]; mat.X.b <- mat.X.2[vec.idx.ab,]
      vec.yobs.a <- vec.yobs.1[vec.idx.ab]; vec.yobs.b <- vec.yobs.2[vec.idx.ab]
      mat.mu.ab <- mat.mu.pair[vec.idx.ab, ]
      # mat.mu.ab: a matrix each row of which denotes a pair of inner products of
      #              covariates and coefficients falling into the category

      vec.rho.ab <- vec.par[length(vec.par)-1] *
        vec.par[length(vec.par)] ^ vec.dist.select[vec.idx.ab]

      if (a==0 & b==0) {
        vec.likelihood.ab <- pbivnorm::pbivnorm(-mat.mu.ab, rho=vec.rho.ab)
        vec.likelihood.ab[which(vec.likelihood.ab<=0)] <- 1e-200
        # to ensure positive likelihood

        vec.deriv.alpha <- rep(0, length(vec.idx.ab))
        mat.deriv.beta <- func.der.beta(-mat.mu.ab, vec.rho.ab, -mat.X.a, -mat.X.b)
        mat.deriv.beta <- matrix(mat.deriv.beta, ncol=length(vec.par)-3)
        # mat.deriv.beta: a matrix of partial derivs wrt regression coefficients
        vec.deriv.rho.ab <- dbivnorm(-mat.mu.ab, vec.rho.ab)
        # vec.deriv.rho.ab: a vector of partial derivs wrt spatial corr
        vec.deriv.sigmasq <- vec.deriv.rho.ab *
          vec.par[length(vec.par)] ^ vec.dist.select[vec.idx.ab]
        # vec.deriv.sigmasq: a vector of partial derivs wrt sigma square
        vec.deriv.rho <- vec.deriv.rho.ab * vec.dist.select[vec.idx.ab] *
          vec.par[length(vec.par)-1] *
          vec.par[length(vec.par)]^(vec.dist.select[vec.idx.ab]-1)
        # vec.deriv.rho: a vector of partial derivs wrt rho

        log.lkd <- sum(log(vec.likelihood.ab)) / weight.sum
        # log.lkd: composite log-likelihood function

        mat.derivs.tmp <- unname(cbind(vec.deriv.alpha,mat.deriv.beta,
                                       vec.deriv.sigmasq,vec.deriv.rho))/vec.likelihood.ab
        mat.derivs[vec.idx.ab,] <- mat.derivs.tmp


      } else if (a==0 & b==1) {
        vec.likelihood.ab <- pnorm(-mat.mu.ab[,1]) -
          pbivnorm::pbivnorm(-mat.mu.ab[,1], vec.par[1]-mat.mu.ab[,2], rho=vec.rho.ab)
        vec.likelihood.ab[which(vec.likelihood.ab<=0)] <- 1e-200
        # to ensure positive likelihood

        vec.deriv.alpha <- - dnorm(vec.par[1]-mat.mu.ab[,2]) *
          pnorm(func.xi(cbind(vec.par[1]-mat.mu.ab[,2], -mat.mu.ab[,1]),vec.rho.ab))
        mat.deriv.beta <- - dnorm(-mat.mu.ab[,1]) * mat.X.a -
          func.der.beta(cbind(-mat.mu.ab[,1],vec.par[1]-mat.mu.ab[,2]),
                        vec.rho.ab, -mat.X.a, -mat.X.b)
        mat.deriv.beta <- matrix(mat.deriv.beta, ncol=length(vec.par)-3)
        # mat.deriv.beta: a matrix of partial derivs wrt regression coefficients
        vec.deriv.rho.ab <-
          - dbivnorm(cbind(-mat.mu.ab[,1],vec.par[1]-mat.mu.ab[,2]), vec.rho.ab)
        # vec.deriv.rho.ab: a vector of partial derivs wrt spatial corr
        vec.deriv.sigmasq <- vec.deriv.rho.ab *
          vec.par[length(vec.par)] ^ vec.dist.select[vec.idx.ab]
        # vec.deriv.sigmasq: a vector of partial derivs wrt sigma square
        vec.deriv.rho <- vec.deriv.rho.ab * vec.dist.select[vec.idx.ab] *
          vec.par[length(vec.par)-1] *
          vec.par[length(vec.par)]^(vec.dist.select[vec.idx.ab]-1)
        # vec.deriv.rho: a vector of partial derivs wrt rho

        log.lkd <- log.lkd + sum(log(vec.likelihood.ab)) / weight.sum

        mat.derivs.tmp <- unname(cbind(vec.deriv.alpha,mat.deriv.beta,
                                       vec.deriv.sigmasq,vec.deriv.rho))/vec.likelihood.ab
        mat.derivs[vec.idx.ab,] <- mat.derivs.tmp

      } else if (a==1 & b==0) {
        vec.likelihood.ab <- pnorm(-mat.mu.ab[,2]) -
          pbivnorm::pbivnorm(vec.par[1]-mat.mu.ab[,1], -mat.mu.ab[,2], rho=vec.rho.ab)
        vec.likelihood.ab[which(vec.likelihood.ab<=0)] <- 1e-200
        # to ensure positive likelihood

        vec.deriv.alpha <- - dnorm(vec.par[1]-mat.mu.ab[,1]) *
          pnorm(func.xi(cbind(vec.par[1]-mat.mu.ab[,1],-mat.mu.ab[,2]),vec.rho.ab))
        mat.deriv.beta <- - dnorm(-mat.mu.ab[,2]) * mat.X.b -
          func.der.beta(cbind(vec.par[1]-mat.mu.ab[,1],-mat.mu.ab[,2]),
                        vec.rho.ab, -mat.X.a, -mat.X.b)
        mat.deriv.beta <- matrix(mat.deriv.beta, ncol=length(vec.par)-3)
        # mat.deriv.beta: a matrix of partial derivs wrt regression coefficients
        vec.deriv.rho.ab <-
          - dbivnorm(cbind(vec.par[1]-mat.mu.ab[,1],-mat.mu.ab[,2]), vec.rho.ab)
        # vec.deriv.rho.ab: a vector of partial derivs wrt spatial corr
        vec.deriv.sigmasq <- vec.deriv.rho.ab *
          vec.par[length(vec.par)] ^ vec.dist.select[vec.idx.ab]
        # vec.deriv.sigmasq: a vector of partial derivs wrt sigma square
        vec.deriv.rho <- vec.deriv.rho.ab * vec.dist.select[vec.idx.ab] *
          vec.par[length(vec.par)-1] *
          vec.par[length(vec.par)]^(vec.dist.select[vec.idx.ab]-1)
        # vec.deriv.rho: a vector of partial derivs wrt rho

        log.lkd <- log.lkd + sum(log(vec.likelihood.ab)) / weight.sum

        mat.derivs.tmp <- unname(cbind(vec.deriv.alpha,mat.deriv.beta,
                                       vec.deriv.sigmasq,vec.deriv.rho))/vec.likelihood.ab
        mat.derivs[vec.idx.ab,] <- mat.derivs.tmp

      } else if (a==1 & b==1) {
        vec.likelihood.ab <- pbivnorm::pbivnorm(mat.mu.ab-vec.par[1], rho=vec.rho.ab)
        vec.likelihood.ab[which(vec.likelihood.ab<=0)] <- 1e-200
        # to ensure positive likelihood

        vec.deriv.alpha <- func.der.beta(mat.mu.ab-vec.par[1], vec.rho.ab, -1, -1)
        mat.deriv.beta <- func.der.beta(mat.mu.ab-vec.par[1],vec.rho.ab,mat.X.a,mat.X.b)
        mat.deriv.beta <- matrix(mat.deriv.beta, ncol=length(vec.par)-3)
        # mat.deriv.beta: a matrix of partial derivs wrt regression coefficients
        vec.deriv.rho.ab <-
          dbivnorm(mat.mu.ab-vec.par[1], vec.rho.ab)
        # vec.deriv.rho.ab: a vector of partial derivs wrt spatial corr
        vec.deriv.sigmasq <- vec.deriv.rho.ab *
          vec.par[length(vec.par)] ^ vec.dist.select[vec.idx.ab]
        # vec.deriv.sigmasq: a vector of partial derivs wrt sigma square
        vec.deriv.rho <- vec.deriv.rho.ab * vec.dist.select[vec.idx.ab] *
          vec.par[length(vec.par)-1] *
          vec.par[length(vec.par)]^(vec.dist.select[vec.idx.ab]-1)
        # vec.deriv.rho: a vector of partial derivs wrt rho

        log.lkd <- log.lkd + sum(log(vec.likelihood.ab)) / weight.sum

        mat.derivs.tmp <- unname(cbind(vec.deriv.alpha,mat.deriv.beta,
                                       vec.deriv.sigmasq,vec.deriv.rho))/vec.likelihood.ab
        mat.derivs[vec.idx.ab,] <- mat.derivs.tmp

      } else if (a==0) {
        vec.u.b <- vec.par[1] * vec.yobs.b - mat.mu.ab[,2]
        vec.m <- (-mat.mu.ab[,1]-vec.rho.ab*vec.u.b)/sqrt(1-vec.rho.ab^2)
        vec.m[which(vec.m<(-10))] <- -10
        # to ensure that pnorm(vec.m) is always positive
        vec.ratio.norm <- dnorm(vec.m) / pnorm(vec.m)
        vec.likelihood.ab <- vec.par[1] * dnorm(vec.u.b) *
          pnorm((-mat.mu.ab[,1]-vec.rho.ab*vec.u.b)/sqrt(1-vec.rho.ab^2))
        vec.likelihood.ab[which(vec.likelihood.ab<=0)] <- 1e-200
        # to ensure positive likelihood

        vec.deriv.alpha <- 1/vec.par[1] - vec.u.b*vec.yobs.b -
          vec.rho.ab*vec.yobs.b*vec.ratio.norm/sqrt(1-vec.rho.ab^2)
        mat.deriv.beta <- vec.u.b*mat.X.b +
          (vec.rho.ab*mat.X.b-mat.X.a)*vec.ratio.norm/sqrt(1-vec.rho.ab^2)
        mat.deriv.beta <- matrix(mat.deriv.beta, ncol=length(vec.par)-3)
        # mat.deriv.beta: a matrix of partial derivs wrt regression coefficients
        vec.deriv.rho.ab <- - (vec.rho.ab*mat.mu.ab[,1]+vec.u.b) *
          vec.ratio.norm / (1-vec.rho.ab^2)^{3/2}
        # vec.deriv.rho.ab: a vector of partial derivs wrt spatial corr
        vec.deriv.sigmasq <- vec.deriv.rho.ab *
          vec.par[length(vec.par)] ^ vec.dist.select[vec.idx.ab]
        # vec.deriv.sigmasq: a vector of partial derivs wrt sigma square
        vec.deriv.rho <- vec.deriv.rho.ab * vec.dist.select[vec.idx.ab] *
          vec.par[length(vec.par)-1] *
          vec.par[length(vec.par)]^(vec.dist.select[vec.idx.ab]-1)
        # vec.deriv.rho: a vector of partial derivs wrt rho

        log.lkd <- log.lkd + sum(log(vec.likelihood.ab)) / weight.sum

        mat.derivs.tmp <- unname(cbind(vec.deriv.alpha,mat.deriv.beta,
                                       vec.deriv.sigmasq,vec.deriv.rho))
        mat.derivs[vec.idx.ab,] <- mat.derivs.tmp
        # mat.derivs: a matrix of partial derivs for parameters of interest

      } else if (a==1) {
        vec.u.b <- vec.par[1] * vec.yobs.b - mat.mu.ab[,2]
        vec.m <- (mat.mu.ab[,1]-vec.par[1]+vec.rho.ab*vec.u.b)/sqrt(1-vec.rho.ab^2)
        vec.m[which(vec.m<(-10))] <- -10
        # to ensure that pnorm(vec.m) is always positive
        vec.ratio.norm <- dnorm(vec.m) / pnorm(vec.m)
        vec.likelihood.ab <- vec.par[1] * dnorm(vec.u.b) *
          pnorm((mat.mu.ab[,1]-vec.par[1]+vec.rho.ab*vec.u.b)/sqrt(1-vec.rho.ab^2))
        vec.likelihood.ab[which(vec.likelihood.ab<=0)] <- 1e-200
        # to ensure positive likelihood

        vec.deriv.alpha <- 1/vec.par[1] - vec.u.b*vec.yobs.b +
          (vec.rho.ab*vec.yobs.b-1)*vec.ratio.norm/sqrt(1-vec.rho.ab^2)
        mat.deriv.beta <- vec.u.b*mat.X.b -
          (vec.rho.ab*mat.X.b-mat.X.a)*vec.ratio.norm/sqrt(1-vec.rho.ab^2)
        mat.deriv.beta <- matrix(mat.deriv.beta, ncol=length(vec.par)-3)
        # mat.deriv.beta: a matrix of partial derivs wrt regression coefficients
        vec.deriv.rho.ab <- (vec.rho.ab*(mat.mu.ab[,1]-vec.par[1])+vec.u.b) *
          vec.ratio.norm / (1-vec.rho.ab^2)^{3/2}
        # vec.deriv.rho.ab: a vector of partial derivs wrt spatial corr
        vec.deriv.sigmasq <- vec.deriv.rho.ab *
          vec.par[length(vec.par)] ^ vec.dist.select[vec.idx.ab]
        # vec.deriv.sigmasq: a vector of partial derivs wrt sigma square
        vec.deriv.rho <- vec.deriv.rho.ab * vec.dist.select[vec.idx.ab] *
          vec.par[length(vec.par)-1] *
          vec.par[length(vec.par)]^(vec.dist.select[vec.idx.ab]-1)
        # vec.deriv.rho: a vector of partial derivs wrt rho

        log.lkd <- log.lkd + sum(log(vec.likelihood.ab)) / weight.sum

        mat.derivs.tmp <- unname(cbind(vec.deriv.alpha,mat.deriv.beta,
                                       vec.deriv.sigmasq,vec.deriv.rho))
        mat.derivs[vec.idx.ab,] <- mat.derivs.tmp
        # mat.derivs: a matrix of partial derivs for parameters of interest


      } else if (b==0) {
        vec.u.a <- vec.par[1] * vec.yobs.a - mat.mu.ab[,1]
        vec.m <- (-mat.mu.ab[,2]-vec.rho.ab*vec.u.a)/sqrt(1-vec.rho.ab^2)
        vec.m[which(vec.m<(-10))] <- -10
        # to ensure that pnorm(vec.m) is always positive
        vec.ratio.norm <- dnorm(vec.m) / pnorm(vec.m)
        vec.likelihood.ab <- vec.par[1] * dnorm(vec.u.a) *
          pnorm((-mat.mu.ab[,2]-vec.rho.ab*vec.u.a)/sqrt(1-vec.rho.ab^2))
        vec.likelihood.ab[which(vec.likelihood.ab<=0)] <- 1e-200
        # to ensure positive likelihood

        vec.deriv.alpha <- 1/vec.par[1] - vec.u.a*vec.yobs.a -
          vec.rho.ab*vec.yobs.a*vec.ratio.norm/sqrt(1-vec.rho.ab^2)
        mat.deriv.beta <- vec.u.a*mat.X.a +
          (vec.rho.ab*mat.X.a-mat.X.b)*vec.ratio.norm/sqrt(1-vec.rho.ab^2)
        mat.deriv.beta <- matrix(mat.deriv.beta, ncol=length(vec.par)-3)
        # mat.deriv.beta: a matrix of partial derivs wrt regression coefficients
        vec.deriv.rho.ab <- - (vec.rho.ab*mat.mu.ab[,2]+vec.u.a) *
          vec.ratio.norm / (1-vec.rho.ab^2)^{3/2}
        # vec.deriv.rho.ab: a vector of partial derivs wrt spatial corr
        vec.deriv.sigmasq <- vec.deriv.rho.ab *
          vec.par[length(vec.par)] ^ vec.dist.select[vec.idx.ab]
        # vec.deriv.sigmasq: a vector of partial derivs wrt sigma square
        vec.deriv.rho <- vec.deriv.rho.ab * vec.dist.select[vec.idx.ab] *
          vec.par[length(vec.par)-1] *
          vec.par[length(vec.par)]^(vec.dist.select[vec.idx.ab]-1)
        # vec.deriv.rho: a vector of partial derivs wrt rho

        log.lkd <- log.lkd + sum(log(vec.likelihood.ab)) / weight.sum

        mat.derivs.tmp <- unname(cbind(vec.deriv.alpha,mat.deriv.beta,
                                       vec.deriv.sigmasq,vec.deriv.rho))
        mat.derivs[vec.idx.ab,] <- mat.derivs.tmp
        # mat.derivs: a matrix of partial derivs for parameters of interest

      } else if (b==1) {
        vec.u.a <- vec.par[1] * vec.yobs.a - mat.mu.ab[,1]
        vec.m <- (mat.mu.ab[,2]-vec.par[1]+vec.rho.ab*vec.u.a)/sqrt(1-vec.rho.ab^2)
        vec.m[which(vec.m<(-10))] <- -10
        # to ensure that pnorm(vec.m) is always positive
        vec.ratio.norm <- dnorm(vec.m) / pnorm(vec.m)
        vec.likelihood.ab <- vec.par[1] * dnorm(vec.u.a) *
          pnorm((mat.mu.ab[,2]-vec.par[1]+vec.rho.ab*vec.u.a)/sqrt(1-vec.rho.ab^2))
        vec.likelihood.ab[which(vec.likelihood.ab<=0)] <- 1e-200
        # to ensure positive likelihood

        vec.deriv.alpha <- 1/vec.par[1] - vec.u.a*vec.yobs.a +
          (vec.rho.ab*vec.yobs.a-1)*vec.ratio.norm/sqrt(1-vec.rho.ab^2)
        mat.deriv.beta <- vec.u.a*mat.X.a -
          (vec.rho.ab*mat.X.a-mat.X.b)*vec.ratio.norm/sqrt(1-vec.rho.ab^2)
        mat.deriv.beta <- matrix(mat.deriv.beta, ncol=length(vec.par)-3)
        # mat.deriv.beta: a matrix of partial derivs wrt regression coefficients
        vec.deriv.rho.ab <- (vec.rho.ab*(mat.mu.ab[,2]-vec.par[1])+vec.u.a) *
          vec.ratio.norm / (1-vec.rho.ab^2)^{3/2}
        # vec.deriv.rho.ab: a vector of partial derivs wrt spatial corr
        vec.deriv.sigmasq <- vec.deriv.rho.ab *
          vec.par[length(vec.par)] ^ vec.dist.select[vec.idx.ab]
        # vec.deriv.sigmasq: a vector of partial derivs wrt sigma square
        vec.deriv.rho <- vec.deriv.rho.ab * vec.dist.select[vec.idx.ab] *
          vec.par[length(vec.par)-1] *
          vec.par[length(vec.par)]^(vec.dist.select[vec.idx.ab]-1)
        # vec.deriv.rho: a vector of partial derivs wrt rho

        log.lkd <- log.lkd + sum(log(vec.likelihood.ab)) / weight.sum

        mat.derivs.tmp <- unname(cbind(vec.deriv.alpha,mat.deriv.beta,
                                       vec.deriv.sigmasq,vec.deriv.rho))
        mat.derivs[vec.idx.ab,] <- mat.derivs.tmp
        # mat.derivs: a matrix of partial derivs for parameters of interest


      } else { # a=2 and b=2
        mat.u.ab <- vec.par[1] * cbind(vec.yobs.a,vec.yobs.b) - mat.mu.ab
        vec.likelihood.ab <- vec.par[1]^2 * dbivnorm(mat.u.ab, vec.rho.ab)
        vec.likelihood.ab[which(vec.likelihood.ab<=0)] <- 1e-200
        # to ensure positive likelihood

        vec.deriv.alpha <- 2/vec.par[1] -
          (mat.u.ab[,1]*(vec.yobs.a-vec.rho.ab*vec.yobs.b) +
             mat.u.ab[,2]*(vec.yobs.b-vec.rho.ab*vec.yobs.a)) /
          (1-vec.rho.ab^2)
        mat.deriv.beta <- (mat.u.ab[,1]*(mat.X.a-vec.rho.ab*mat.X.b)
                           +mat.u.ab[,2]*(mat.X.b-vec.rho.ab*mat.X.a)) / (1-vec.rho.ab^2)
        mat.deriv.beta <- matrix(mat.deriv.beta, ncol=length(vec.par)-3)
        # mat.deriv.beta: a matrix of partial derivs wrt regression coefficients
        vec.deriv.rho.ab <- (mat.u.ab[,1]*mat.u.ab[,2]-
                               vec.rho.ab*(mat.u.ab[,1]^2+mat.u.ab[,2]^2-vec.rho.ab*mat.u.ab[,1]*mat.u.ab[,2]))/
          (1-vec.rho.ab^2)^2 + vec.rho.ab/(1-vec.rho.ab^2)
        # vec.deriv.rho.ab: a vector of partial derivs wrt spatial corr
        vec.deriv.sigmasq <- vec.deriv.rho.ab *
          vec.par[length(vec.par)] ^ vec.dist.select[vec.idx.ab]
        # vec.deriv.sigmasq: a vector of partial derivs wrt sigma square
        vec.deriv.rho <- vec.deriv.rho.ab * vec.dist.select[vec.idx.ab] *
          vec.par[length(vec.par)-1] *
          vec.par[length(vec.par)]^(vec.dist.select[vec.idx.ab]-1)
        # vec.deriv.rho: a vector of partial derivs wrt rho

        log.lkd <- log.lkd + sum(log(vec.likelihood.ab)) / weight.sum

        mat.derivs.tmp <- unname(cbind(vec.deriv.alpha,mat.deriv.beta,
                                       vec.deriv.sigmasq,vec.deriv.rho))
        mat.derivs[vec.idx.ab,] <- mat.derivs.tmp
        # mat.derivs: a matrix of partial derivs for parameters of interest

      }
    }
  }

  log.lkd <- ifelse(is.finite(log.lkd), log.lkd, -1e6)
  # to ensure log.lkd is finite

  vec.score <- colSums(mat.derivs) / weight.sum
  # vec.score: a vector of score functions as in (10)
  return(list(weight.sum=weight.sum, log.lkd=log.lkd, vec.score=vec.score,
              mat.score=mat.derivs))
}
