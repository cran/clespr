#' Reparameterized Composite Likelihood Calculation for Spatial Ordinal Data
#'
#' \code{func.cl.ord} calculates the composite log-likelihood for reparameterized spatial ordered probit models. This function is internally called by \code{func.cle.ord}.
#'
#' @param vec.yobs a vector of observed responses for all N sites.
#' @param mat.X regression (design) matrix, including intercepts.
#' @param mat.lattice a data matrix containing geographical information of sites. The \emph{i} th row constitutes a set of geographical coordinates.
#' @param radius weight radius.
#' @param n.cat number of categories, at least 2.
#' @param vec.repar a vector of parameters consecutively as follows: a reparameterized vector (tau's) for latent responses, a vector of covariate parameters, a parameter 'sigmasq' modeling covariance matrix, 0<=sigmasq<=1, and a parameter 'rho' reflecting spatial correlation, abs(rho)<=1.
#' @return \code{func.cl.ord} returns a list: number of categories, sum of weights, composite log-likelihood, a vector of scores, and a matrix of first-order partial derivatives for \code{vec.par}.
#' @references Feng, Xiaoping, Zhu, Jun, Lin, Pei-Sheng, and Steen-Adams, Michelle M. (2014) Composite likelihood Estimation for Models of Spatial Ordinal Data and Spatial Proportional Data with Zero/One values. \emph{Environmetrics} 25(8): 571--583.

func.cl.ord.repar <- function(vec.yobs, mat.X, mat.lattice, radius, n.cat, vec.repar) {
  # Put func.cl.ord as local function for parallel computing
  func.cl.ord <- function(vec.yobs, mat.X, mat.lattice, radius, n.cat, vec.par) {

    # input check
    stopifnot(radius>=0, n.cat>=2,
              vec.par[length(vec.par)-1]>=0, vec.par[length(vec.par)-1]<=1,
              vec.par[length(vec.par)]>=-1, vec.par[length(vec.par)]<=1)
    if (n.cat>2) stopifnot(length(vec.par)-NCOL(mat.X)-2>0)

    vec.dist <- c(dist(mat.lattice))
    # vec.dist: a vector of pairwise euclidean distances between
    #           sites with length equal to N*(N-1)/2
    mat.dist <- as.matrix(dist(mat.lattice, upper=T, diag=T))
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

    # input check
    stopifnot(radius >= 0, n.cat >= 2, vec.repar[length(vec.repar) - 1] >= 0, vec.repar[length(vec.repar) - 1] <=
        1, vec.repar[length(vec.repar)] >= -1, vec.repar[length(vec.repar)] <= 1)
    if (n.cat > 2)
        stopifnot(length(vec.repar) - NCOL(mat.X) - 2 > 0)

    if (n.cat == 2) {
        # if n.cat=2, reparameterization is unnecessary
        ls.cl <- func.cl.ord(vec.yobs, mat.X, mat.lattice, radius, n.cat, vec.repar)
        weight.sum.repar <- ls.cl$weight.sum
        log.lkd.repar <- ls.cl$log.lkd
        vec.score.repar <- ls.cl$vec.score
        mat.score.repar <- ls.cl$mat.score
    } else {
        ls.cl <- func.cl.ord(vec.yobs, mat.X, mat.lattice, radius, n.cat, c(cumsum(vec.repar[1:(n.cat - 2)]),
            vec.repar[(n.cat - 1):length(vec.repar)]))
        # cumsum reverses reparameterization of parameters

        weight.sum.repar <- ls.cl$weight.sum
        log.lkd.repar <- ls.cl$log.lkd
        vec.score.par <- ls.cl$vec.score
        vec.score.repar <- c(rev(cumsum(rev(vec.score.par[1:(n.cat - 2)]))), vec.score.par[(n.cat - 1):length(vec.repar)])
        # vec.score.repar: a vector of reparameterized score functions

        mat.score.par <- ls.cl$mat.score

        if (n.cat == 3) {
            # if n.cat=3, reparameterization has no effect
            mat.score.repar <- mat.score.par
        } else {
            rcr <- function(vec) rev(cumsum(rev(vec)))
            mat.score.repar <- cbind(t(apply(mat.score.par[, 1:(n.cat - 2)], 1, rcr)), mat.score.par[, (n.cat -
                1):length(vec.repar)])
            # mat.score.repar: a matrix of reparameterized partial derivs
        }
    }
    return(list(n.cat = n.cat, weight.sum = weight.sum.repar, log.lkd = log.lkd.repar, vec.score = vec.score.repar,
        mat.score = mat.score.repar))
}
