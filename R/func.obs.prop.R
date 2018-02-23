#' Latent Response Transformation for Proportional Data
#'
#' \code{func.obs.prop} transforms a vector of latent responses into the corresponding observed ones under the spatial Tobit model.
#'
#' @param vec.ylat a vector of latent responses for all N sites.
#' @param alpha a cutoff point controlling the probability of latent reponse being one.
#' @return \code{func.obs.prop} returns a vector of observed responses.
#'
#'
#' @examples
#'
#' # A simple example for observation generation
#' a <- sample(c(0,1), 50, replace=TRUE)
#' b <- sample(runif(1000,0,10), 100, replace=TRUE)
#' alpha <- 4
#' vec.yobs <- func.obs.prop(vec.ylat=c(a, b), alpha=alpha)
#'
#' # A complex example
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
#' @references Feng, Xiaoping, Zhu, Jun, Lin, Pei-Sheng, and Steen-Adams, Michelle M. (2014) Composite likelihood Estimation for Models of Spatial Ordinal Data and Spatial Proportional Data with Zero/One values. \emph{Environmetrics} 25(8): 571--583.

func.obs.prop <- function(vec.ylat, alpha) {

    vec.yobs <- rep(0, length(vec.ylat))
    vec.yobs[which(vec.ylat > alpha)] <- 1
    vec.idx.btw <- which(vec.ylat > 0 & vec.ylat <= alpha)
    vec.yobs[vec.idx.btw] <- vec.ylat[vec.idx.btw]/alpha
    return(vec.yobs)
}
