#' Composite Likelihood Estimation for Spatial Proportional Data
#'
#' \code{func.cle.prop} performs composite likelihood estimation of parameters and their standard errors in a spatial Tobit model by maximizing its composite log-likelihood.
#'
#' @param vec.yobs a vector of observed responses for all N sites.
#' @param mat.X regression (design) matrix, including intercepts.
#' @param mat.lattice a data matrix containing geographical information of sites. The i-th row constitutes a set of geographical coordinates.
#' @param radius weight radius.
#' @param n.sim number of simulations used for calculate the Godambe matrix (default: 100).
#' @param parallel logical flag indicating using parallel processing (default: \code{TRUE})
#' @param n.core number of physical cores used for parallel processing (when \code{parallel} is \code{TRUE}), default value is \code{max(detectCores()/2,1)}).
#' @param output logical flag indicates whether printing out result (default: \code{TRUE}).
#'
#' @details Given the design matrix, the vector of observed responses, spatial lattice data, weight radius, and the prespecified number of simulated vectors of responses used in estimating the Godambe information matrix, this function assumes initial values of \eqn{\beta} as the estimates from the standard Type I Tobit model with independent responses. The initial value of \eqn{\alpha} and the right limit of the Tobit model are equally set to 1. Since there is only one cutoff point to be estimated, reparameterization is unnecessary. The function first estimates parameters of interest by maximizing the composite log-likelihood using \code{optim(...,method = "L-BFGS-B")}, then computes the simulated based standard error and asymptotic covariance matrix.
#' @return \code{func.cle.prop} returns a list containing:
#' @return \code{vec.par}: a vector of estimator for \eqn{\theta=(\alpha,\beta,\sigma^2,\rho)};
#' @return \code{vec.se}: a vector of standard error for the estimator;
#' @return \code{mat.asyvar}: estimated asymptotic covariance matrix \eqn{H^{-1}(\theta)J(\theta)H^{-1}(\theta)} for the estimator; and
#' @return \code{vec.comp}: a vector of computational time for parameter and standard error estimation.
#' @return \code{CLIC}: Composite likelihood information criterion proposed by Varin and Vidoni (2005), i.e. \eqn{-2*logCL(\theta) + 2*trace(H^{-1}(\theta)J(\theta))}
#'
#' @examples
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
#' # With parallel computing
#'
#' \dontrun{
#' prop.example <- func.cle.prop(vec.yobs, mat.X, mat.lattice, radius,
#' n.sim=100, parallel = TRUE, n.core = 2)
#'
#' round(prop.example$vec.par,4)
#' # alpha   beta0   beta1   beta2   beta3   beta4 sigma^2     rho
#' # 3.8259  0.9921  1.9679  0.9455  0.0148 -0.9871  0.8386  0.5761
#'
#' round(prop.example$vec.se ,4)
#' # alpha   beta0   beta1   beta2   beta3   beta4 sigma^2     rho
#' # 0.1902  0.1406  0.1103  0.0744  0.0385  0.0652  0.1527  0.1151
#' }
#'
#' # Without parallel computing
#'
#' \dontrun{
#' prop.example2 <- func.cle.prop(vec.yobs, mat.X, mat.lattice, radius, n.sim=100, parallel = FALSE)
#' }
#'
#' @references Feng, Xiaoping, Zhu, Jun, Lin, Pei-Sheng, and Steen-Adams, Michelle M. (2014) Composite likelihood Estimation for Models of Spatial Ordinal Data and Spatial Proportional Data with Zero/One values. \emph{Environmetrics} 25(8): 571--583.


func.cle.prop <- function(vec.yobs, mat.X, mat.lattice, radius, n.sim = 100, parallel = TRUE, n.core = max(detectCores()/2,1), output = TRUE) {

  vec.initial.par <- c(1, unname(coef(AER::tobit(vec.yobs ~ mat.X[, -1], right = 1))), 0.1, 0.2)
  vec.lower <- c(0.1, rep(-Inf, NCOL(mat.X)), 0.001, 0.001)
  vec.upper <- c(rep(Inf, NCOL(mat.X) + 1), 0.999, 0.999)

  fn <- function(vec.par) -func.cl.prop(vec.yobs, mat.X, mat.lattice, radius, vec.par)$log.lkd
  gr <- function(vec.par) -func.cl.prop(vec.yobs, mat.X, mat.lattice, radius, vec.par)$vec.score

  t0 <- proc.time()
  optim.obj <- optim(vec.initial.par, fn = fn, gr = gr, method = "L-BFGS-B", lower = vec.lower, upper = vec.upper)

  vec.par.opt <- optim.obj$par

  t1 <- proc.time() -t0

  if (output) cat("Parameter estimate is ",round(vec.par.opt,5), "\n")

  # estimation of standard errors
  ls.cl.par.opt <- func.cl.prop(vec.yobs, mat.X, mat.lattice, radius, vec.par.opt)
  func.outsq <- function(vec) vec %o% vec


  mat.Hessian <- matrix(rowSums(apply(ls.cl.par.opt$mat.score, 1, func.outsq)), nrow = length(vec.par.opt),
                            byrow = TRUE)/ls.cl.par.opt$weight.sum

  mat.H.inv <- solve(mat.Hessian)

  t0 <- proc.time()

  mat.dist <- as.matrix(dist(mat.lattice, upper = TRUE, diag = TRUE))
  mat.cov <- vec.par.opt[length(vec.par.opt) - 1] * vec.par.opt[length(vec.par.opt)]^mat.dist

  mat.Z.sim <- t(chol(mat.cov)) %*% array(rnorm(length(vec.yobs) * n.sim), c(length(vec.yobs), n.sim)) + matrix(rep(mat.X %*%
                                                                                                                      vec.par.opt[2:(length(vec.par.opt) - 2)], n.sim), ncol = n.sim)
  mat.eps.sim <- diag(sqrt(1 - vec.par.opt[length(vec.par.opt) - 1]), length(vec.yobs)) %*% array(rnorm(length(vec.yobs) *
                                                                                                          n.sim), c(length(vec.yobs), n.sim))
  mat.ylat.sim <- mat.Z.sim + mat.eps.sim

  func.obs.sim <- function(vec.ylat) func.obs.prop(vec.ylat, vec.par.opt[1])
  mat.yobs.sim <- apply(mat.ylat.sim, 2, func.obs.sim)

  func.score.sim <- function(vec.yobs) func.cl.prop(vec.yobs, mat.X = mat.X, mat.lattice = mat.lattice, radius = radius, vec.par = vec.par.opt)$vec.score

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
    mat.J <- cov(foreach(i = 1:n.sim, .export="func.cl.prop", .combine=rbind) %dopar% {
      func.score.sim(mat.yobs.sim[,i])
    })

    stopCluster(cl)

    options(warn=0) # Resume warning message

  } else {mat.J <- cov(t(apply(mat.yobs.sim, 2, func.score.sim)))}

  mat.asyvar <- mat.H.inv %*% mat.J %*% mat.H.inv

  vec.se.opt <- sqrt(diag(mat.asyvar))

  if (output) cat("Estimated standard error is ",round(vec.se.opt ,5), "\n")

  t2 <- proc.time() - t0

  t.comp <- c(t1[3],t2[3]); names(t.comp) <- c("est", "SE")

  CLIC <- clordr::clic(logCL = optim.obj$value ,mat.hessian = mat.Hessian, mat.J = mat.J)

  names(vec.par.opt) <- names(vec.se.opt) <- colnames(mat.asyvar) <-
    c("alpha",paste0("beta",0:(NCOL(mat.X)-1)),"sigma^2","rho")

  return(list(vec.par = vec.par.opt,CLIC=CLIC,
              vec.se = vec.se.opt, mat.asyvar=mat.asyvar,vec.comp=t.comp))
}
