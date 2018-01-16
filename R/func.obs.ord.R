#' Latent Response Transformation for Spatial Ordinal Data
#'
#' \code{func.obs.ord} transforms a vector of latent responses into the corresponding observed ones under the spatial Probit model.
#'
#' @param vec.ylat a vector of latent responses for all N sites.
#' @param vec.alpha a vector of prespecified cutoff points, ascending with length at least 3, including -Inf, 0, and Inf.
#' @return \code{func.obs.prop} returns a vector of observed responses.
#'
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
#' @references Feng, Xiaoping, Zhu, Jun, Lin, Pei-Sheng, and Steen-Adams, Michelle M. (2014) Composite likelihood Estimation for Models of Spatial Ordinal Data and Spatial Proportional Data with Zero/One values. \emph{Environmetrics} 25(8): 571--583.

func.obs.ord <- function(vec.ylat, vec.alpha) {

    vec.yobs <- rep(0, length(vec.ylat))
    for (i in 1:(length(vec.alpha) - 2)) {
        vec.yobs[which(vec.ylat > vec.alpha[i + 1] & vec.ylat <= vec.alpha[i + 2])] <- i
    }
    return(vec.yobs)
}
