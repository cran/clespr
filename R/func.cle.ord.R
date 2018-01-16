#' Composite Likelihood Estimation for Spatial Ordinal Data
#'
#' \code{func.cle.ord} performs composite likelihood estimation of parameters and their standard errors in a spatial ordered probit model by maximizing its composite log-likelihood.
#'
#' @param vec.yobs a vector of observed responses for all N sites.
#' @param mat.X regression (design) matrix, including intercepts.
#' @param mat.lattice a data matrix containing geographical information of sites. The ith row constitutes a set of geographical coordinates.
#' @param radius weight radius.
#' @param n.cat number of categories.
#' @param n.sim number of simulations used for calculate the Godambe matrix (default: 100).
#' @param parallel logical flag indicates using parallel processing (default: \code{TRUE}).
#' @param n.core number of physical cores used for parallel processing (when \code{parallel} is \code{TRUE}, default value is \code{max(detectCores()/2,1)}).
#' @param output logical flag indicates whether printing out result (default: \code{TRUE}).
#'
#' @details Given the design matrix, the vector of observed responses, spatial lattice data, weight radius, number of categories, and the prespecified number of simulated vectors of responses used in estimating the Godambe information, this function assumes initial values of cutoff points and \eqn{\beta} as the estimates from the standard ordered probit regression with independent responses. After initial reparameterization, it first estimates parameters of interest by maximizing the composite log-likelihood using \code{optim}, then computes the reparameterized sample covariance matrix and the set of standard errors, and finally reverse the reparameterization to obtain estimates corresponding to the original parameterization.
#' @return \code{func.cle.ord} returns a list containing:
#' @return \code{vec.par}: a vector of estimator for \eqn{\theta}=(cutoff,\eqn{\beta,\sigma^2,\rho)};
#' @return \code{vec.se}: a vector of standard error for the estimator;
#' @return \code{mat.asyvar}: estimated asymptotic covariance matrix \eqn{H^{-1}(\theta)J(\theta)H^{-1}(\theta)} for the estimator; and
#' @return \code{vec.comp}: a vector of computational time for parameter and standard error estimation.
#'
#' @examples
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
#' # With parallel computing
#'
#' \dontrun{
#' ord.example <- func.cle.ord(vec.yobs, mat.X, mat.lattice, radius, n.cat,
#' n.sim=100, parallel = TRUE, n.core = 2)
#' round(ord.example$vec.par,4)
#' # [1]  1.8395  0.9550  1.9690  0.9565  0.0349 -1.0398  0.8200  0.5578
#' }
#'
#' # Without parallel computing
#'
#' \dontrun{
#' ord.example2 <- func.cle.ord(vec.yobs, mat.X, mat.lattice, radius,
#' n.cat, n.sim=100, parallel = FALSE)
#' round(ord.example2$vec.par,4)
#' # [1]  1.8395  0.9550  1.9690  0.9565  0.0349 -1.0398  0.8200  0.5578
#' }
#'
#' @references Feng, Xiaoping, Zhu, Jun, Lin, Pei-Sheng, and Steen-Adams, Michelle M. (2014) Composite likelihood Estimation for Models of Spatial Ordinal Data and Spatial Proportional Data with Zero/One values. \emph{Environmetrics} 25(8): 571--583.


func.cle.ord <- function(vec.yobs, mat.X, mat.lattice, radius, n.cat, n.sim = 100,  parallel = TRUE, n.core = max(detectCores()/2,1), output = TRUE) {

  # vec.yobs: a vector of observed responses for all N sites mat.X: regression (design) matrix, including
  # intercepts mat.lattice: a data matrix containing geographical information of sites.  The ith row constitutes
  # a set of geographical coordinates radius: weight radius n.cat: number of categories, at least 2 n.sim:
  # number of vectors simulated responses used to compute sample cov matrix

  # input check

  stopifnot(radius >= 0, n.cat >= 2, n.sim >= 0)

  options(warn=-1) # Remove warning message: glm.fit: fitted probabilities numerically 0 or 1 occurred

  vec.initial.par <- c(rep(0.5, n.cat - 2), 0, unname(coef(MASS::polr(factor(vec.yobs) ~ mat.X[, -1], method = "probit"))),
                       0.4, 0.2)
  # vec.initial.par: a vector of initial values with cutoffs and betas set to the estimates from standard
  # ordered probit model with independent responses

  options(warn=0) # Resume warning message

  vec.lower <- c(rep(0.05, n.cat - 2), rep(-10, NCOL(mat.X)), 0.01, 0.01)
  vec.upper <- c(rep(10, n.cat - 2), rep(10, NCOL(mat.X)), 0.99, 0.99)
  # vec.lower and vec.upper: vectors of lower and upper bounds for parameters of interest

  if (n.cat == 2 | n.cat == 3) {
    # if n.cat=2 or n.cat=3, no need to reparameterize
    vec.initial.repar <- vec.initial.par
  } else {
    vec.initial.repar <- c(vec.initial.par[1], diff(vec.initial.par[1:(n.cat - 2)]), vec.initial.par[(n.cat -
                                                                                                        1):length(vec.initial.par)])
    # reparameterization
  }

  fn <- function(vec.repar) -func.cl.ord.repar(vec.yobs, mat.X, mat.lattice, radius, n.cat, vec.repar)$log.lkd
  # fn: the opposite of composite log-likelihood function to be minimized

  gr <- function(vec.repar) -func.cl.ord.repar(vec.yobs, mat.X, mat.lattice, radius, n.cat, vec.repar)$vec.score
  # gr: a vector of opposite of score functions used as gradients

  t0 <- proc.time()
  vec.repar.opt <- optim(vec.initial.repar, fn = fn, gr = gr, method = "L-BFGS-B", lower = vec.lower, upper = vec.upper)$par
  # vec.repar.opt: estimated reparameterized unknowns

  if (n.cat == 2) {
    # if n.cat=2, no need to reverse reparameterization
    vec.par.opt <- vec.repar.opt
  } else {
    vec.par.opt <- c(cumsum(vec.repar.opt[1:(n.cat - 2)]), vec.repar.opt[(n.cat - 1):length(vec.repar.opt)])
    # vec.par.opt: a vector of estimates of parameters of interest
  }

  t1 <- proc.time() - t0

  if (output) cat("Parameter estimate is ",round(vec.par.opt,5), "\n")

  t0 <- proc.time()

  # estimation of standard errors
  ls.cl.repar.opt <- func.cl.ord.repar(vec.yobs, mat.X, mat.lattice, radius, n.cat, vec.repar.opt)

  func.outsq <- function(vec) vec %o% vec
  mat.H.inv <- solve(matrix(rowSums(apply(ls.cl.repar.opt$mat.score, 1, func.outsq)), nrow = length(vec.repar.opt),
                            byrow = TRUE)/ls.cl.repar.opt$weight.sum)

  # mat.H.inv: the inverse of H_N(theta) as in (14)

  mat.dist <- as.matrix(dist(mat.lattice, upper = TRUE, diag = TRUE))
  mat.cov <- vec.repar.opt[length(vec.repar.opt) - 1] * vec.repar.opt[length(vec.repar.opt)]^mat.dist

  mat.Z.sim <- t(chol(mat.cov)) %*% array(rnorm(length(vec.yobs) * n.sim), c(length(vec.yobs), n.sim)) + matrix(rep(mat.X %*%
                                                                                                                      vec.repar.opt[(n.cat - 1):(length(vec.repar.opt) - 2)], n.sim), ncol = n.sim)

  mat.eps.sim <- diag(sqrt(1 - vec.repar.opt[length(vec.repar.opt) - 1]), length(vec.yobs)) %*% array(rnorm(length(vec.yobs) *
                                                                                                              n.sim), c(length(vec.yobs), n.sim))

  mat.ylat.sim <- mat.Z.sim + mat.eps.sim

  func.obs.ord.sim <- function(vec.ylat) {
    if (n.cat == 2) {
      vec.yobs <- func.obs.ord(vec.ylat, c(-Inf, 0, Inf))
    } else {
      vec.yobs <- func.obs.ord(vec.ylat, c(-Inf, 0, vec.par.opt[1:(n.cat - 2)], Inf))
    }
    return(vec.yobs)
  }

  mat.yobs.sim <- apply(mat.ylat.sim, 2, func.obs.ord.sim)

  func.score.sim <- function(vec.yobs) func.cl.ord.repar(vec.yobs, mat.X = mat.X, mat.lattice = mat.lattice, radius = radius, n.cat = n.cat, vec.repar = vec.repar.opt)$vec.score

  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  {(abs(x - round(x)) < tol)&&(!(x== 0))}
  # Small function to check if positive integer

  if (parallel == TRUE){
    i <- 1
    options(warn=-1) # Remove warning message: closing unused connection

    if (!is.double(n.core)) {
      cat("Wrong input type for n.core, replaced by the default value, max(detectCores()/2,1) = ",max(detectCores()/2,1),"\n")
      n.core <- max(parallel::detectCores()/2,1)
    }else if(!is.wholenumber(n.core)) {
      cat("Input for n.core is not a positive integer, replaced by the default value, max(detectCores()/2,1) = ",max(detectCores()/2,1),"\n")
      n.core <- max(parallel::detectCores()/2,1)
    }

    cl <- makeCluster(n.core)
    registerDoParallel(cl)
    mat.J <- cov(foreach(i = 1:n.sim, .export="func.cl.ord.repar", .combine=rbind) %dopar% {
      func.score.sim(mat.yobs.sim[,i])
    })

    stopCluster(cl)

    options(warn=0) # Resume warning message

  } else{ mat.J <- cov(t(apply(mat.yobs.sim, 2, func.score.sim))) }
  # mat.J: the matrix hat J_N as in (16)

  if (n.cat == 2 | n.cat == 3) {
    mat.A <- diag(length(vec.par.opt))  # identity matrix
  } else {

    mat.lower.tri <- matrix(1, n.cat - 2, n.cat - 2)
    mat.lower.tri[row(mat.lower.tri) < col(mat.lower.tri)] <- 0

    mat.A <- magic::adiag(mat.lower.tri, diag(length(vec.par.opt) - n.cat + 2))
    # mat.A: a matrix of transformation as a block diagonal matrix consisting of a lower triangle matrix of ones
    # and an identity matrix
  }

  mat.asyvar <- mat.A %*% mat.H.inv %*% mat.J %*% mat.H.inv %*% t(mat.A)

  vec.se.opt <- sqrt(diag(mat.asyvar))
  # vec.se.opt: a vector of standard error based on (16) and Theorem 1

  if (output) cat("Estimated standard error is ",round(vec.se.opt ,5), "\n")

  t2 <- proc.time() - t0

  t.comp <- c(t1[3],t2[3]); names(t.comp) <- c("est", "SE")

  return(list(vec.par = vec.par.opt,
              vec.se = vec.se.opt, mat.asyvar=mat.asyvar,vec.comp=t.comp))
}
