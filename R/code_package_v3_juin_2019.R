#' Density of the Pareto 1 distribution
#'
#' @param x numeric, (positive) vector
#' @param mu numeric, lower bound, default is 1
#' @param alpha numeric, tail index
#' @param xi numeric, inverse of \code{alpha}
#' @return The density of the Pareto 1 distribution at points \code{x} \deqn{f(x)=alpha*mu^alpha/x^(alpha+1)} Expected value exists when \code{alpha} is larger than 1, and is \deqn{E[X]=alpha/(alpha-1)}
#' @seealso \code{\link{ppareto1}}, \code{\link{qpareto1}} and \code{\link{rpareto1}}
#' @examples
#' # density of a Pareto1(1,1.5) at point 2
#' dpareto1(2, 1, alpha=1.5)
#' 0.265165
#' # expected value of a Pareto1(1,1.5)
#' integrate(function(x) x*dpareto1(x, 1, alpha=1.5),-Inf,Inf)
#' # theoretical expected value
#' 1*(1.5)/(1.5-1)
#' @export
dpareto1 <- function (x, mu=1, alpha=1/xi, xi=1/alpha)  { (alpha*mu^alpha/abs(x)^(alpha+1)*(x>=mu)) }

#' Cumulative distribution function of the Pareto 1 distribution
#'
#' @param x numeric, (positive) vector
#' @param mu numeric, lower bound, default is 1
#' @param alpha numeric, tail index
#' @param xi numeric, inverse of \code{alpha}
#' @return the c.d.f. of the Pareto 1 distribution at points \code{x} \deqn{F(x)=1-(mu/x)^(alpha)} Expected value exists when \code{alpha} is larger than 1, and is \deqn{E[X]=alpha/(alpha-1)}
#' @seealso \code{\link{dpareto1}}, \code{\link{qpareto1}} and \code{\link{rpareto1}}
#' @examples
#' # probability to a Pareto1(1,1.5) to be smaller than 2
#' ppareto1(2, 1, 1.5)
#' 0.6464466
#' # expected value of a Pareto1(1,1.5)
#' integrate(function(x) (1-ppareto1(x, 1, alpha=1.5)),0,Inf)
#' # theoretical expected value
#' 1*(1.5)/(1.5-1)
#' @export
ppareto1 <- function (x, mu=1, alpha=1/xi, xi=1/alpha)  { (1 - ( abs(x)/mu )^(-alpha))*(x>=mu) }

#' Quantile function of the Pareto 1 distribution
#'
#' @param p numeric, (probabilities) vector
#' @param mu numeric, lower bound, default is 1
#' @param alpha numeric, tail index
#' @param xi numeric, inverse of \code{alpha}
#' @return the quantile function of the Pareto 1 distribution at points \code{p}
#' @seealso \code{\link{dpareto1}}, \code{\link{ppareto1}} and \code{\link{rpareto1}}
#' @examples
#' qpareto1(.5, 1, alpha=1.5)
#' 1.587401
#' # expected value of a Pareto1(1,1.5)
#' integrate(function(x) qpareto1(x, 1, alpha=1.5),0,1)
#' # theoretical expected value
#' 1*(1.5)/(1.5-1)
#' @export
qpareto1 <- function (p, mu, alpha=1/xi, xi=1/alpha)  { mu*(1-p)^(-1/alpha)*((p>=0)&(p<=1)) }

#' Random generation of the Pareto 1 distribution
#'
#' @param n integer
#' @param mu numeric, lower bound, default is 1
#' @param alpha numeric, tail index
#' @param xi numeric, inverse of \code{alpha}
#' @return generates \code{n} values of the Pareto 1 distribution
#' @seealso \code{\link{dpareto1}}, \code{\link{ppareto1}} and \code{\link{qpareto1}}
#' @examples
#' set.seed(123)
#' rpareto1(6, 1, alpha=1.5)
#' # expected value of a Pareto1(1,1.5)
#' mean(rpareto1(1e6, 1, alpha=1.5))
#' # theoretical expected value
#' 1*(1.5)/(1.5-1)
rpareto1 <- function (n, mu, alpha=1/xi, xi=1/alpha)  { mu*(1-runif(n))^(-1/alpha) }

#' Maximum Likelihood estimation of the Pareto 1 distribution, with weights
#'
#' @param data a vector of observations
#' @param weights a vector of weights (default = 1)
#' @param threhold the threshold parameter of the Pareto 1 distribution (\code{mu})
#' @return a list with the index \code{alpha} and \code{k}, the number of observations above \code{threshold}
#' @references Charpentier & Flachaire (2019) \emph{Pareto Models for Top Incomes } \href{https://hal.archives-ouvertes.fr/hal-02145024}{hal-02145024}
#' @examples
#' set.seed(123)
#' x <- rpareto1(100, 1, 1.5)
#' w <- rgamma(100, 10, 10)
#' estim <- MLE.pareto1(data=x, weights=w, threshold=1)
#' estim
MLE.pareto1 <- function(data, weights=rep(1,length(x)), threshold=min(x))
{
  foo <- cbind(data,weights)
  foo <- foo[foo[,1]>threshold,]
  xx <- foo[,1]
  ww <- foo[,2]/sum(foo[,2])
  m <- as.numeric(threshold)
  a <- 1/(sum(ww*log(xx))-log(m))
  k <- NROW(xx)
  return(list(alpha=a, xi=1/a, k=k))
}

#' Cumulative distribution of the Generalized Pareto distribution (GPD)
#'
#' @param x numeric, (positive) vector
#' @param mu numeric, lower bound, default is 0
#' @param alpha numeric, tail index
#' @param xi numeric, inverse of \code{alpha}
#' @param beta numeric, scaling paramater, default = 1
#' @return the c.d.f. of the Generalized Pareto distribution at points \code{x}, \deqn{F(x)=1-(1+xi*(x-mu)/beta)^(-1/xi)} (assuming \code{xi} non null). The expected value is then \deqn{E[X]=mu+beta/(1-xi)}
#' @seealso \code{\link{dgpd}}, and \code{\link{rgpd}}
#' @examples
#' pgpd(2, xi=1/1.5, 1, 1)
#' ## expected value of a GPD(1,1.5,0,1)
#' integrate(function(x) 1-pgpd(x,xi=1/1.5,mu=0,beta=1),0,Inf)
#' ## theoretical value
#' 0+1/(1-1/1.5)
pgpd <- function (x, xi=1/alpha, alpha=1/xi, mu = 0, beta = 1)  { (1 - (1+(xi*(x-mu))/beta)^(-1/xi))*(x>=mu) }

#' Density of the Generalized Pareto distribution (GPD)
#'
#' @param x numeric, (positive) vector
#' @param mu numeric, lower bound, default is 0
#' @param alpha numeric, tail index
#' @param xi numeric, inverse of \code{alpha}
#' @param beta numeric, scaling paramater, default = 1
#' @return the density of the Pareto 1 distribution at points \code{x}, \deqn{f(x)=(1/beta)*(1+xi*(x-mu)/beta)^(-1/xi-1)}
#' @seealso \code{\link{pgpd}}, and \code{\link{rgpd}}
#' @examples
#' dgpd(2, xi=1/1.5, 1, 1)
#' # expected value of a GPD(1,1.5,0,1)
#' integrate(function(x) x*dgpd(x,alpha=1.5,mu=0,beta=1),0,Inf)
#' integrate(function(x) x*dgpd(x,xi=1/1.5,mu=0,beta=1),0,Inf)
#' ## theoretical value
#' 0+1/(1-1/1.5)
dgpd <- function (x, xi=1/alpha, alpha=1/xi, mu = 0, beta = 1)  { (beta^(-1))*(1+(xi*(x-mu))/beta)^((-1/xi)-1)*(x>=mu) }

#' Random generation of the Generalized Pareto distribution (GPD)
#'
#' @param n integer
#' @param mu numeric, lower bound, default is 0
#' @param alpha numeric, tail index
#' @param xi numeric, inverse of \code{alpha}
#' @param beta numeric, scaling paramater, default = 1
#' @return generates \code{n} values of the Pareto 1 distribution at points \code{x}
#' @seealso \code{\link{dgpd}}, and \code{\link{pgpd}}
#' @examples
#' rgpd(10, xi=1/1.5, 1, 1)
#' # expected value of a GPD(1,1.5,0,1)
#' mean(rgpd(1e6,alpha=1.5,mu=0,beta=1))
rgpd <- function (n, xi=1/alpha, alpha=1/xi, mu = 0, beta = 1)  { mu + (beta/xi)*((1-runif(n))^(-xi)-1) }

#' Maximum Likelihood estimation of the Generalized Pareto distribution, with weights
#'
#' @param data numeric, vector of observations
#' @param weights numeric, vector of weights (default = 1)
#' @param threhol numeric, threshold parameter of the Generalized Pareto distribution (\code{mu})
#' @param nextrmes numeric, number of largest values considered (integer)
#' @param method character, method used for inference (\code{"ml"} for maximum likelihood)
#' @param information character (not used)
#' @return a list \code{xi} the tail index, \code{mu} the threshold, \code{"beta"} the scaling coefficient, and \code{k} the number of observations above \code{mu}
#' @references Charpentier & Flachaire (2019) \emph{Pareto Models for Top Incomes } \href{https://hal.archives-ouvertes.fr/hal-02145024}{hal-02145024}
#' @examples
#' set.seed(123)
#' x <- rpareto1(100, 1, 1.5)
#' w <- rgamma(100,10,10)
#' estim <- MLE.gpd(data=x, weights=w, threshold=1)
#' estim$par.ests
MLE.gpd <- function (data, weights=rep(1,length(x)), threshold = NA, nextremes = NA, method="ml", information = c("observed", "expected"), ...)
{
  n <- length(data)
  if (is.na(nextremes) && is.na(threshold))
    stop("Enter either a threshold or the number of upper extremes")
  if (!is.na(nextremes) && !is.na(threshold))
    stop("Enter EITHER a threshold or the number of upper extremes")
  if (!is.na(nextremes))

    data <- as.numeric(data)

  foo=cbind(data,weights)
  foo=foo[foo[,1]>threshold,]
  x=foo[,1]
  w=foo[,2]/sum(foo[,2])

  exceedances <- x
  excess <- exceedances - threshold
  Nu <- length(excess)
  xbar <- sum(w*excess)
  method <- "ml"
  s2 <- sum(w*(excess-xbar)^2)
  xi0 <- -0.5 * (((xbar * xbar)/s2) - 1)
  beta0 <- 0.5 * xbar * (((xbar * xbar)/s2) + 1)
  theta <- c(xi0, beta0)
  negloglik <- function(theta, tmp) {
    xi <- theta[1]
    beta <- theta[2]
    cond1 <- beta <= 0
    cond2 <- (xi <= 0) && (max(tmp) > (-beta/xi))
    if (cond1 || cond2)
      f <- 1e+06
    else {
      y <- logb(1 + (xi * tmp)/beta)
      y <- w*y/xi
      f <- logb(beta) + (1 + xi) * sum(y)
    }
    f
  }
  fit <- optim(theta, negloglik, hessian = TRUE, ..., tmp = excess)
  if (fit$convergence)
    warning("optimization may not have succeeded")
  par.ests <- fit$par
  converged <- fit$convergence
  nllh.final <- fit$value
  p.less.thresh <- 1 - Nu/n
  return(list(xi=par.ests[1], mu=threshold, beta=par.ests[2], k=Nu))
}

.EPDinput <- function(y, gamma, kappa, tau, kappaTau = TRUE) {
    #
    #   Modified EPD function from the ReIns R package to have:
    #        - weighted ML estimation
    #        - results from only one cutoff
    #        - direct ML estimation by default
    # original R code :  ReIns package version 1.0.7
    #    https://github.com/TReynkens/ReIns/blob/master/R/EPD.R
    #    https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R
    #
    if (!is.numeric(gamma)) {
        stop("gamma should be numeric.")
    }

    if (!is.numeric(kappa)) {
        stop("kappa should be numeric.")
    }

    if (!is.numeric(tau)) {
        stop("tau should be numeric.")
    }


    if (any(tau >= 0)) {
        stop("tau should be strictly negative.")
    }

    if (any(gamma <= 0)) {
        stop("gamma should be strictly positive.")
    }
    if (kappaTau) {
        if (any(kappa <= pmax(-1, 1/tau))) {
            stop("kappa should be larger than max(-1,1/tau).")
        }
    }

    ly <- length(y)
    lg <- length(gamma)
    lk <- length(kappa)
    lt <- length(tau)

    l <- c(ly, lg, lk, lt)

    ind <- which(l > 1)

    if (length(ind) > 1) {

        if (!length(unique(l[ind])) == 1) {
            stop("All input arguments should have length 1 or equal length.")
        }
    }
}

#' Density of the Extended Pareto distribution
#'
#' @param x mumeric, (positive) vector
#' @param gamma numeric, (strictly positive) number (the tail index)
#' @param kappa numeric, must be larger than max{-1,1/tau}
#' @param tau numeric, (negative) number (default is -1)
#' @param log logical, indicating if logarithm of density should be returned
#' @return the density of the Extended Pareto distribution at points \code{x}
#' @source \url{https://github.com/TReynkens/ReIns/blob/master/R/EPD.R} Tom Reynkens, ReIns package version 1.0.7
#' @examples
#' depd(2,.5,1,-1)
#' ## Expected value of a EPD(gamma=.5,kappa=1,tau=-1)
#' integrate(function(x) x*depd(x,.5,1,-1),-Inf,Inf)
depd <- function(x, gamma, kappa, tau = -1, log = FALSE) {
    #
    #   Modified EPD function from the ReIns R package to have:
    #        - weighted ML estimation
    #        - results from only one cutoff
    #        - direct ML estimation by default
    # original R code :  ReIns package version 1.0.7
    #    https://github.com/TReynkens/ReIns/blob/master/R/EPD.R
    #    https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R
    #

    .EPDinput(x, gamma, kappa, tau, kappaTau = TRUE)

    d <- 1 / (gamma*x^(1/gamma+1)) * (1+kappa*(1-x^tau))^(-1/gamma-1) *
    (1+kappa*(1-(1+tau)*x^tau))
    d[x <= 1] <- 0
    if (log) d <- log(d)
    return(d)
}

#' Cumulative Distribution Function of the Extended Pareto distribution
#'
#' @param x numeric, (positive) vector
#' @param gamma numeric, (strictly positive) number (the tail index)
#' @param kappa numeric, number - must be larger than max{-1,1/tau}
#' @param tau numeric, (negative) number (default is -1)
#' @param log logical, indicating if logarithm of density should be returned
#' @return the c.d.f. of the Extended Pareto distribution at points \code{x}
#' @source \url{https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R} Tom Reynkens, \code{ReIns} package version 1.0.7
#' @examples
#' pepd(2,.5,1,-1)
#' ## Expected value of a EPD(gamma=.5,kappa=1,tau=-1)
#' integrate(function(x) 1-pepd(x,.5,1,-1),0,Inf)
pepd <- function(x, gamma, kappa, tau = -1, lower.tail = TRUE, log.p = FALSE) {
    #
    #   Modified EPD function from the ReIns R package to have:
    #        - weighted ML estimation
    #        - results from only one cutoff
    #        - direct ML estimation by default
    # original R code :  ReIns package version 1.0.7
    #    https://github.com/TReynkens/ReIns/blob/master/R/EPD.R
    #    https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R
    #

    .EPDinput(x, gamma, kappa, tau, kappaTau = FALSE)

    p <- 1 - (x * (1+kappa*(1-x^tau)))^(-1/gamma)

    p[x <= 1] <- 0

    if (any(kappa <= pmax(-1, 1/tau))) {
        if (length(kappa) > 1 | length(tau) > 1) {
            p[kappa <= pmax(-1, 1/tau)] <- NA
        } else {
            p <- NA
        }
    }

    if (!lower.tail) p <- 1-p

    if (log.p) p <- log(p)

    return(p)
}

#' Quantile Function of the Extended Pareto distribution
#'
#' @param p numeric, vector of probabilities (in the interval [0,1])
#' @param gamma numeric, (strictly positive) number (the tail index)
#' @param kappa numeric, number - must be larger than max{-1,1/tau}
#' @param tau numeric, (negative) number (default is -1)
#' @param log logical, indicating if logarithm of density should be returned
#' @return the c.d.f. of the Extended Pareto distribution at points \code{x}
#' @source \url{https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R} Tom Reynkens, \code{ReIns} package version 1.0.7
#' @examples
#' qepd(.5,.5,1,-1)
qepd <-  function(p, gamma, kappa, tau = -1, lower.tail = TRUE, log.p = FALSE) {
    #
    #   Modified EPD function from the ReIns R package to have:
    #        - weighted ML estimation
    #        - results from only one cutoff
    #        - direct ML estimation by default
    # original R code :  ReIns package version 1.0.7
    #    https://github.com/TReynkens/ReIns/blob/master/R/EPD.R
    #    https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R
    #

.EPDinput(p, gamma, kappa, tau, kappaTau = TRUE)

    if (log.p) p <- exp(p)

    if (!lower.tail) p <- 1-p

    if (any(p < 0 | p > 1)) {
        stop("p should be between 0 and 1.")
    }

    l <- length(p)
    Q <- numeric(l)

    endpoint <- 10

    if (any(p < 1)) {

        mx <- max(p[p < 1])

        while (pepd(endpoint, gamma, kappa, tau) <= mx) {
            endpoint <- endpoint*10
        }
    }

    for (i in 1:l) {

        if (p[i] < .Machine$double.eps) {
            # p=0 case
            Q[i] <- 1

        } else if (abs(p[i]-1) > .Machine$double.eps) {
            # 0<p<1 case

            # Function to minimise
            f <- function(x) {
                ((1-p[i])^(-gamma) - x*(1+kappa*(1-x^tau)))^2
            }
            # If minimising fails return NA
            Q[i] <- tryCatch(optimise(f, lower=1, upper=endpoint)$minimum, error=function(e) NA)

        } else {
            # p=1 case
            Q[i] <- Inf
        }

    }

    return(Q)
}

#' Random Generation of the Extended Pareto distribution
#'
#' @param n integer, number of generations
#' @param gamma numeric, (strictly positive) number (the tail index)
#' @param kappa numeric, must be larger than max{-1,1/tau}
#' @param tau numeric, (negative) number (default is -1)
#' @param quasimc logical, indicating quasi-Monte Carlo techniques (Halton sequence) should be used
#' @return a vector of \code{n} values generated from an Extended Pareto distribution
#' @source \url{https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R} Tom Reynkens, \code{ReIns} package version 1.0.7
#' @examples
#' set.seed(123)
#' repd(6, gamma=.5,kappa=1,tau=-1)
#' repd(6, gamma=.5,kappa=1,tau=-1,quasimc=TRUE)
#' ## Mean of a EPD(gamma=.5,kappa=1,tau=-1) distribution
#' mean(repd(1e4, gamma=.5,kappa=1,tau=-1))
#' mean(repd(1e4, gamma=.5,kappa=1,tau=-1),quasimc=TRUE)
repd <-  function(n, gamma, kappa, tau = -1, quasimc = FALSE) {
  if(quasimc == FALSE) v=qepd(runif(n), gamma=gamma, kappa=kappa, tau=tau)
  if(quasimc == TRUE){
    if (!require("randtoolbox")) install.packages("randtoolbox")
    require(randtoolbox)
    u= halton(n)
    v=qepd(u, gamma=gamma, kappa=kappa, tau=tau)
  }
    return(v)
}

#' Hill estimator of the tail index, with weights
#'
#' @param numeric, data the vector of observations
#' @param numeric weights the vector of weights (default is uniform weights)
#' @return Hill estimator of \code{xi} (the tail index) and \code{alpha}, as well as \code{k} the number of observations
#'
#'  Note that \code{weights} do not need to sum up to 1 (function will convert them via \code{weights/sum(weights)})
#' @references Charpentier & Flachaire (2019) \emph{Pareto Models for Top Incomes } \href{https://hal.archives-ouvertes.fr/hal-02145024}{hal-02145024}
#' @examples
#' set.seed(123)
#' x <- rpareto1(100, 1, 1.5)
#' w <- rgamma(100,10,10)
#' Hill(x,w)
Hill = function(data,weights=rep(1,length(data))){
    w <- weights/sum(weights)
    n <- length(data)
    X <- as.numeric(sort(data))
    Hill <- sum(w[2:n]*(log(X[2:n])-log(X[1])))
    return(list(xi = Hill, alpha=1/Hill, k=n))
}

#' Fit the Extended Pareto distribution to a vector of observations, with weights
#'
#' @param data numeric, vector of observations
#' @param weights numeric, vector of (positive) weights
#' @param rho numeric, parameter of Fraga Alves et al. (2003) estimate
#' @param start numeric, vector of length 2 containing the starting values for the optimisation
#' @param direct logical, indicating if the parameters are obtained by directly maximising the log-likelihood function
#' @param warnings logical indicating if possible warnings from the optimisation function are shown
#' @references Charpentier & Flachaire (2019) \emph{Pareto Models for Top Incomes } \href{https://hal.archives-ouvertes.fr/hal-02145024}{hal-02145024}
#' @return a list with \code{k} the number of observations used, \code{gamma} the vector of the corresponding estimates for the tail parameter of the EPD, \code{kappa} the vector of the corresponding MLE estimates for the kappa parameter of the EPD and \code{tau} the vector of the corresponding estimates for the second order tail index parameter of the EPD using Hill estimates and values for \code{rho}
#' @source adapted from \url{https://github.com/TReynkens/ReIns/blob/master/R/EPD.R} Tom Reynkens, ReIns package version 1.0.7
#' @examples
#' set.seed(123)
#' x <- rpareto1(100,1,.5)
#' w <- rgamma(100,10,10)
#' EPD(data=x, weights=w)
EPD <- function(data, weights=rep(1,length(data)), rho = -1, start = NULL, direct = TRUE, warnings = FALSE, ...) {
    #
    #   Modified EPD function from the ReIns R package to have:
    #        - weighted ML estimation
    #        - results from only one cutoff
    #        - direct ML estimation by default
    # original R code :  ReIns package version 1.0.7
    #    https://github.com/TReynkens/ReIns/blob/master/R/EPD.R
    #    https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R
    #

    df=data.frame(data,weights)
    df=df[order(df$data),]
    w <- df$weights/sum(df$weights)
    X <- as.numeric(df$data)
    n <- length(X)
    K <- (n-1)

    if (n == 1) {
        stop("We need at least two data points.")
    }

    if (direct) {
        EPD <- .EPDdirectMLE(data=data, weights=w, rho=rho, start=start, warnings=warnings)
    } else {
        # Select parameter using approach of Beirlant, Joosens and Segers (2009).
        EPD <- .EPDredMLE(data=data, weights=w, rho=rho)
    }

    if (length(rho) == 1) {
        EPD$gamma <- as.vector(EPD$gamma)
        EPD$kappa <- as.vector(EPD$kappa)
        EPD$tau <- as.vector(EPD$tau)

    }

    if (length(rho) == 1) {
        return(list(k=K, gamma=unlist(EPD$gamma[K]), kappa=unlist(EPD$kappa[K]), tau=unlist(EPD$tau[K])))
    } else {
        return(list(k=K, gamma=unlist(EPD$gamma[K,]), kappa=unlist(EPD$kappa[K,]), tau=unlist(EPD$tau[K,])))
    }
}

.EPDredMLE <- function(data, weights=rep(1,length(data)), rho = -1) {
    #
    #   Modified EPD function from the ReIns R package to have:
    #        - weighted ML estimation
    #        - results from only one cutoff
    #        - direct ML estimation by default
    #   original R code :  ReIns package version 1.0.7
    #    https://github.com/TReynkens/ReIns/blob/master/R/EPD.R
    #    https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R
    #
    #   Fit EPD using approach of Beirlant, Joosens and Segers (2009).

    df=data.frame(data,weights)
    df=df[order(df$data),]
    w <- df$weights/sum(df$weights)
    X <- as.numeric(df$data)

    n <- length(X)
    K <- (n-1)

    if (n == 1) {
        stop("We need at least two data points.")
    }

    nrho <- length(rho)
    rho.orig <- rho

    H <- Hill(data, w)$xi

    if (all(rho > 0) & nrho == 1) {
        rho <- .rhoEst(data, alpha=1, tau=rho, weights=w)$rho
        beta <- -rho

    } else if (all(rho < 0)) {
        beta <- -rho

    } else {
        stop("rho should be a single positive number or a vector (of length >=1) of negative numbers.")
    }


    gamma <- matrix(0, n-1, nrho)
    kappa <- matrix(0, n-1, nrho)
    tau <- matrix(0, n-1, nrho)

    beta <- -rho

    for (j in 1:nrho) {

        if (nrho == 1 & all(rho.orig > 0)) {
            tau[K, 1] <- -beta[K]/H


            rhovec <- rho
        } else {

            tau[K, j] <- -beta[j]/H

            rhovec <- rep(rho[j], n-1)
        }


        E <- numeric(n-1)

        for (k in K) {
            i <- 1:k
            E[k] <- sum( w[n-k+i] * (X[n-k+i]/X[n-k])^tau[k,j] )
        }

        kappa[K,j] <- H * (1-2*rhovec[K]) * (1-rhovec[K])^3 / rhovec[K]^4 * (E[K] - 1 / (1-rhovec[K]))

        gamma[K,j] <- H - kappa[K,j] * rhovec[K] / (1 - rhovec[K])

    }

    return(list(gamma=gamma, kappa=kappa, tau=tau))
}

.EPDdirectMLE <- function(data, weights=rep(1,length(data)), rho = -1, start = NULL,  warnings = FALSE) {
    #
    #   Modified EPD function from the ReIns R package to have:
    #        - weighted ML estimation
    #        - results from only one cutoff
    #        - direct ML estimation by default
    # original R code :  ReIns package version 1.0.7
    #    https://github.com/TReynkens/ReIns/blob/master/R/EPD.R
    #    https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R
    #


    df=data.frame(data,weights)
    df=df[order(df$data),]
    w <- df$weights/sum(df$weights)
    X <- as.numeric(df$data)
    n <- length(X)

    if (n == 1) {
        stop("We need at least two data points.")
    }

    nrho <- length(rho)
    rho.orig <- rho

    H <- Hill(data, w)$xi

    if (all(rho > 0) & nrho == 1) {
        rho <- .rhoEst(data, alpha=1, tau=rho, weights=w)$rho
        beta <- -rho

    } else if (all(rho < 0)) {
        beta <- -rho

    } else {
        stop("rho should be a single positive number or a vector (of length >=1) of negative numbers.")
    }

    gamma <- matrix(0, n-1, nrho)
    kappa <- matrix(0, n-1, nrho)
    tau <- matrix(0, n-1, nrho)

    for (j in 1:nrho) {

        for (k in (n-1):(n-1)) {

            epddf <- df[df$data > X[n-k],]
            epddata <- epddf$data/X[n-k]
            epdw    <- epddf$w

            if (nrho == 1 & all(rho.orig > 0)) {
                tau[k,1] <- -beta[k]/H

            } else {

                tau[k,j] <- -beta[j]/H
            }

            if (is.null(start)) {
                start2 <- numeric(2)
                start2[1] <- H
                start2[2] <- 0
            } else if (is.matrix(start)) {

                if (nrow(start >= n-1)) {
                    start2 <- numeric(2)
                    start2[1] <- start[k,1]
                    start2[2] <- start[k,2]
                } else {
                    stop("start does not contain enough rows.")
                }

            } else {
                start2 <- start
            }

            if (tau[k,j] < 0) {
                tmp <- EPDfit(epddata, start=start2, tau=tau[k,j], weights=epdw)
                gamma[k,j] <- tmp[1]
                kappa[k,j] <- tmp[2]
            } else {
                gamma[k,j] <- kappa[k,j] <- NA
            }

        }

    }

    return(list(gamma=gamma, kappa=kappa, tau=tau))
}

#' Fit the Extended Pareto distribution to a vector of observations, with weights, using maximum likelihood estimation
#'
#' @param data numeric, vector of observations
#' @param tau numeric, value for tau in the EPD distribution
#' @param weights numeric, vector of (positive) weights
#' @param rho numeric, parameter of Fraga Alves et al. (2003) estimate
#' @param start numeric, vector of length 2 containing the starting values for the optimisation (default are \code{c(.1,1})
#' @param warnings logical indicating if possible warnings from the optimisation function are shown
#' @return a list with \code{gamma} the vector of the corresponding estimates for the tail parameter of the EPD, \code{kappa} the vector of the corresponding MLE estimates for the kappa parameter of the EPD, as well as \code{tau}
#' @references Charpentier & Flachaire (2019) \emph{Pareto Models for Top Incomes } \href{https://hal.archives-ouvertes.fr/hal-02145024}{hal-02145024}
#' @source adapted from \url{https://github.com/TReynkens/ReIns/blob/master/R/EPD.R} Tom Reynkens, \code{ReIns} package version 1.0.7
#' @examples
#' set.seed(123)
#' x <- rpareto1(100, mu=1, alpha=.5)
#' w <- rgamma(100,10,10)
#' EPDfit(data=x, tau=-3.3, weights=w)
EPDfit <- function(data, tau, start = c(0.1, 1), warnings = FALSE, weights=rep(1,length(data))) {
    #
    #   Modified EPD function from the ReIns R package to have:
    #        - weighted ML estimation
    #        - results from only one cutoff
    #        - direct ML estimation by default
    # original R code :  ReIns package version 1.0.7
    #    https://github.com/TReynkens/ReIns/blob/master/R/EPD.R
    #    https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R
    #

    w <- weights/sum(weights)
    if (is.numeric(start) & length(start) == 2) {
        gamma_start <- start[1]
        kappa_start <- start[2]
    } else {
        stop("start should be a 2-dimensional numeric vector.")
    }


    if (ifelse(length(data) > 1, var(data) == 0, 0)) {
        sg <- c(NA, NA)
    } else {

        fit <- optim(par=c(gamma_start, kappa_start), fn=.EPDneglogL, Y=data, tau=tau, weights=w)

        sg <- fit$par

        if (fit$convergence > 0 & warnings) {
            warning("Optimisation did not complete succesfully.")
            if (!is.null(fit$message)) {
                print(fit$message)
            }
        }
    }
    listsg=list(gamma=sg[1], kappa=sg[2], tau=tau)
    return(listsg)
}


.EPDneglogL <- function(theta, Y, tau, weights) {
    #
    #   Modified EPD function from the ReIns R package to have:
    #        - weighted ML estimation
    #        - results from only one cutoff
    #        - direct ML estimation by default
    # original R code :  ReIns package version 1.0.7
    #    https://github.com/TReynkens/ReIns/blob/master/R/EPD.R
    #    https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R
    #

    w <- weights/sum(weights)
    gamma <- theta[1]
    kappa <- theta[2]

    if (kappa <= max(-1, 1/tau) | gamma <= 0) {
        logL <- -10^6
    } else {
        logL <- sum( w*log(depd(Y, gamma=gamma, kappa=kappa, tau=tau)) )
    }

    return(-logL)
}

.rhoEst <- function(data, alpha = 1, theta1 = 2, theta2 = 3, tau = 1, weights=rep(1,length(data))) {
    #
    #   Modified EPD function from the ReIns R package to have:
    #        - weighted ML estimation
    #        - results from only one cutoff
    #        - direct ML estimation by default
    # original R code :  ReIns package version 1.0.7
    #    https://github.com/TReynkens/ReIns/blob/master/R/EPD.R
    #    https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R
    #


    if (alpha <= 0) {
        stop("alpha should be strictly positive.")
    }

    if (tau <= 0) {
        stop("tau should be strictly positive.")
    }

    df=data.frame(data,weights)
    df=df[order(df$data),]
    w <- df$weights/sum(df$weights)
    X <- as.numeric(df$data)

    n <- length(X)
    rho <- numeric(n)
    Tn <- numeric(n)
    K <- (n-1)

    M_alpha <- numeric(n)
    M_alpha_theta1 <- numeric(n)
    M_alpha_theta2 <- numeric(n)

    l <- log(X[n-K+1])
    for (k in K) {
        M_alpha[k] <- sum( (l[1:k]-log(X[n-k]))^alpha ) / k
        M_alpha_theta1[k] <- sum( (l[1:k]-log(X[n-k]))^(alpha*theta1) ) / k
        M_alpha_theta2[k] <- sum( (l[1:k]-log(X[n-k]))^(alpha*theta2) ) / k
    }

    Tn[K] <- ( (M_alpha[K]/gamma(alpha+1))^tau - (M_alpha_theta1[K]/gamma(alpha*theta1+1))^(tau/theta1)  ) /
    ( (M_alpha_theta1[K]/gamma(alpha*theta1+1))^(tau/theta1) - (M_alpha_theta2[K]/gamma(alpha*theta2+1))^(tau/theta2)  )

    rho[K] <- 1 - ( 2 * Tn[K] / ( 3 - Tn[K]) ) ^ (1/alpha)

    return(list(k=K, rho=rho[K], Tn=Tn[K]))
}


ProbEPD <- function(data, q, gamma, kappa, tau, ...) {
    #
    #   Modified EPD function from the ReIns R package to have:
    #        - weighted ML estimation
    #        - results from only one cutoff
    #        - direct ML estimation by default
    # original R code :  ReIns package version 1.0.7
    #    https://github.com/TReynkens/ReIns/blob/master/R/EPD.R
    #    https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R
    #


    if ( length(gamma) != length(kappa) | length(gamma) != length(tau)) {
        stop("gamma, kappa and tau should have equal length.")
    }

    X <- as.numeric(sort(data))


    n <- length(X)
    prob <- numeric(n)
    K <- (n-1)

    K2 <- K[which(gamma[K] > 0)]

    prob[K2] <- (K2+1)/(n+1) * (1 - pepd(q/X[n-K2], gamma=gamma[K2], kappa=kappa[K2], tau=tau[K2]))
    prob[prob < 0 | prob > 1] <- NA


    return(list(k=K, P=prob[K], q=q))

}

#' Large Return Period associated to the Extended Pareto distribution
#'
#' @param data numeric, a vector of observations
#' @param q numeric, large quantile to estimate 1/P[X>q]
#' @param gamma numeric, vector of \code{n-1} estimates for the EVD obtained from [EPD]
#' @param kappa numeric, vector of \code{n-1} estimates for the EVD obtained from [EPD]
#' @param tau numeric, vector of \code{n-1} estimates for the EVD obtained from [EPD]
#' @return a list with \code{k} the vector of the values of the tail parameter k, \code{R} the vector of the corresponding return period and \code{q} the used large quantile
#' @source \url{https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R} Tom Reynkens, ReIns package version 1.0.7
#' @examples
#' set.seed(123)
#' x <- rpareto1(100, mu=1, alpha=1.5)
#' fit <- EPD(x)
#' ReturnEPD(data=x, q=.01, gamma=fit$gamma, kappa=fit$kappa, tau=fit$tau)
ReturnEPD <- function(data, q, gamma, kappa, tau, ...) {
    #
    #   Modified EPD function from the ReIns R package to have:
    #        - weighted ML estimation
    #        - results from only one cutoff
    #        - direct ML estimation by default
    # original R code :  ReIns package version 1.0.7
    #    https://github.com/TReynkens/ReIns/blob/master/R/EPD.R
    #    https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R
    #

    if ( length(gamma) != length(kappa) | length(gamma) != length(tau)) {
        stop("gamma, kappa and tau should have equal length.")
    }

    X <- as.numeric(sort(data))
    n <- length(X)
    r <- numeric(n)
    K <- (n-1)

    K2 <- K[which(gamma[K] > 0)]

    r[K2] <- (n+1)/(K2+1) / (1 - pepd(q/X[n-K2], gamma=gamma[K2], kappa=kappa[K2], tau=tau[K2]))


    r[which(gamma[K] <= 0)] <- NA

    r[r <= 0] <- NA

    return(list(k=K, R=r[K], q=q))

}

#' @title Top_Share
#' @description This function estimates Top Share
#'
#' @param data numeric, a vector of observations
#' @param weight numeric, a vector of weights (default is equal weights)
#' @param p numeric, the probability level (default 0.01 for the top 1\%)
#' @param q numeric, (possibly a vector), the probability level to model a Pareto distribution (default 0.1)
#' @param method is the distribution considered (default \code{"edf"} for the empirical distribution function, but can be \code{"pareto1"}, \code{"gpd"} or  \code{"epd"})
#' @param edp.direct logical (default \code{TRUE}) for the method used for EPD fit
#' @return estimation of the share of income/wealth owned by the top 100p\% of the population, assuming that the top 100q\% of the distribution is Pareto distributed. The list contains elements for each value of \code{q} if it is a vector, with \code{index} the share, \code{alpha} the inverse of the tail index, \code{coef},  \code{share.index} the value of \code{p}, \code{share.pareto} the value of the threshold as a percentage (.1 for the 90\% quantile) and \code{threshold} for the numerical value.
#' @references Charpentier & Flachaire (2019) \emph{Pareto Models for Top Incomes } \href{https://hal.archives-ouvertes.fr/hal-02145024}{hal-02145024}
#' @examples
#' ################ Top Share on Synthetic Data
#' url_1 <- "https://github.com/freakonometrics/TopIncomes/raw/master/dataframe_yw_1.csv"
#' df <- read.table(url_1,sep=";",header=TRUE)
#' ## via the Empirical Distribution Function
#' Top_Share(data = df$y, weights = df$w,method="edf")$index
#' ## via a Pareto 1 distribution (fitted on the top 100q%)
#' Top_Share(data = df$y, weights = df$w,method="pareto1")$index
#' ## via a Generalized Pareto distribution (fitted on the top 100q%)
#' Top_Share(data = df$y, weights = df$w,method="gpd")$index
#' ## via an Extended Pareto distribution (fitted on the top 100q%)
#' Top_Share(data = df$y, weights = df$w,method="epd")$index
Top_Share <- function(data, weights=rep(1,length(data)), p=.01, q=.1, method="edf", epd.direct=TRUE) {

	if (!require("Hmisc")) install.packages("Hmisc")
    require(Hmisc)
    #
    # p : top 100p% share
    # q : top 100q% of the distribution being Pareto
    #
    # method="edf" : sample top share, that is, based on the EDF
    # method="pareto1" : EDF+Pareto1
    #
    x = data
    if (min(weights)<0) stop('Error: weights should be positive \n\n')
     weights = weights/sum(weights)
    if (p>1) stop('Error: p should be smaller than 1 \n\n')
    if (p<0) stop('Error: p should be greater than 0 \n\n')
    ## Top Share based on Pareto I and GPD Models

    if (min(q)>1) stop('Error: q should be smaller than 1 \n\n')
    if (max(q)<0) stop('Error: q should be greater than 0 \n\n')

    vindex=valpha=vcoef=vgamma=vkappa=vshare.index=vshare.pareto=vthreshold=vtau=vbeta=vxi=vk=rep(NA,length(q))

    L="error"
    ## Top Share based on the Empirical Distribution Function (EDF)

    if(method=="edf") {
        for(s in 1:length(q)){

            up=Hmisc::wtd.quantile(x, weights=weights, probs=1-p, normwt=TRUE) # weighted (1-p)-quantile
            up=as.numeric(up)

        vindex[s]=tis=sum(weights*x*(x>up))/sum(weights*x)
        }
        L=list(index=vindex,share=p,k=length(x),method="edf")
    }


    if(method=="pareto1" || method=="pareto2" || method=="gpd") {

              for(s in 1:length(q)){
        u=Hmisc::wtd.quantile(x, weights=weights, probs=1-q[s], normwt=TRUE)
        u=as.numeric(u)  # threshold = weighted (1-q)-quantile
        datax=cbind(x,weights)
        dataq=datax[x<u,]
        xq = Hmisc::wtd.mean(dataq[,1], weights=dataq[,2])
        # Estimate the Pareto distribution with weighted data
        if(method=="pareto1") {
            coef=MLE.pareto1(x,weights=weights,threshold=u)
            sigma=u
            alpha=coef$alpha
        }
        if(method=="pareto2" || method=="gpd") {
            coef=MLE.gpd(x,weights=weights,threshold=u)
            sigma=coef$beta/coef$xi
            alpha=1/coef$xi
        }
        # Top Income Shares with weighted data
        if (alpha<1) { tis=NaN
        } else if(p<=q[s]) {
            num = (alpha/(alpha-1))*sigma*(p/q[s])^(-1/alpha) + u - sigma
            den = (1-q[s])*xq + q[s]*sigma/(alpha-1) + q[s]*u
            tis = p*num/den
        } else if(p>q[s])  {
            up=Hmisc::wtd.quantile(x, weights=weights, probs=1-p, normwt=TRUE)
            up=as.numeric(up)
            xp = Hmisc::wtd.mean(x[x<up], weights=weigths[x<up])
            den = (1-q[s])*xq + q[s]*sigma/(alpha-1) + q[s]*u
            tis = 1 - (1-p)*xp/den
        }
        vindex[s] <- tis
        valpha[s] <- alpha
        if(method=="pareto2" || method=="gpd") vxi[s] <- coef$xi
        if(method=="pareto2" || method=="gpd") vbeta[s] <- coef$beta
        vshare.index[s]=p
        vshare.pareto[s]=q[s]
        vthreshold[s]=u
        vk[s]=sum(x>u)
              }
    L=list(index=vindex,alpha=valpha,xi=vxi,beta=vbeta,share.index=p,
        share.pareto=vshare.pareto,threshold=vthreshold,k=vk, method=method)}

    ## Top Share based on EPD Model
    if(method=="epd") {
    for(s in 1:length(q)){
        u=Hmisc::wtd.quantile(x, weights=weights, probs=1-q[s], normwt=TRUE)
        u=as.numeric(u)  # threshold = weighted (1-q)-quantile

        dataqq=x[x>=u]
        coef=EPD(x[x>=u], weights=weights[x>=u], direct=epd.direct)
        delta=coef$kappa
        tau=coef$tau
        alpha <- 1/coef$gamma
         xq = Hmisc::wtd.mean(x[x < u], weights = weights[x < u])
        up=Hmisc::wtd.quantile(x, weights=weights, probs=1-p, normwt=TRUE)
        up=as.numeric(up)
        xp = Hmisc::wtd.mean(x[x<up], weights=weights[x<up])

        pextpareto=function(x, u=1, delta=0, tau=-1, alpha){#Compute CDF
            d=1-((x/u)*(1+delta-delta*(x/u)^tau))^(-alpha)
            d[x<=u] <- 0
            return(d) }
        ff_bis <- function(x) (1-pextpareto(1/x, u=u, delta=delta, tau=tau, alpha=alpha))/x^2

        if (alpha<=1) {tis=NaN    # infinite mean (tail=alpha=1/xi < 1)
        } else if(delta<max(-1,1/tau)) {tis=NaN   # kappa should be largen than max(-1,1/tau)
        } else if(p<=q[s]) {
            uprim=u*qepd(1-p/q[s], gamma=coef$gamma, kappa=coef$kappa, tau=coef$tau)
            Eup = try( integrate(ff_bis, lower=0, upper=1/uprim)$value , TRUE)
            Eu = try( integrate(ff_bis, lower=0, upper=1/u)$value , TRUE)
            if(inherits(Eup, "try-error") && inherits(Eup, "try-error")) tis=NaN
            if(!(inherits(Eup, "try-error") && inherits(Eup, "try-error"))) tis=(p*uprim+q[s]*Eup)/((1-q[s])*xq+q[s]*(u+Eu))
        } else if(p>q[s])  {
            Eu = try( integrate(ff_bis, lower=0, upper=1/u)$value , TRUE)
            Ex = ((1-q[s])*xp+q[s]*(u+Eu))
            if(inherits(Eu, "try-error")) tis=NaN
            if(!inherits(Eu, "try-error")) tis = 1-(1-p)*xp/Ex
        }

    vindex[s] <- tis
    valpha[s] <- 1/coef$gamma
    vgamma[s]=coef$gamma
    vkappa[s]=coef$kappa
    vtau[s]=coef$tau
    vshare.index[s]=p
    vshare.pareto[s]=q[s]
    vthreshold[s]=u
    vk[s]=sum(x>u)
    }
L=list(index=vindex,alpha=valpha,tau=vtau,kappa=vkappa,gamma=vgamma,share.index=p,
share.pareto=vshare.pareto,threshold=vthreshold, k=vk, method=method)
}
    return(L)
}

#' Pareto diagrams - Pareto 1, GPD and EPD
#'
#' @param data numeric, a vector of observations
#' @param weight numeric, a vector of weights (default is equal weights)
#' @param p numeric, the probability level (default 0.01 for the top 1\%)
#' @param q numeric, the probability level to model a Pareto distribution (default 0.1)
#' @param viz logical, \code{TRUE} to plot the estimates
#' @return a table with estimations of top share and a graph
#' @references Charpentier & Flachaire (2019) \emph{Pareto Models for Top Incomes } \href{https://hal.archives-ouvertes.fr/hal-02145024}{hal-02145024}
#' @examples
#' ################ Pareto diagams on Synthetic Data
#' url_1 <- "https://github.com/freakonometrics/TopIncomes/raw/master/dataframe_yw_1.csv"
#' df <- read.table(url_1,sep=";",header=TRUE)
#' \dontrun{Pareto_diagram(data = df$y, weights = df$w)}
Pareto_diagram = function(data, weights=rep(1,length(data)), p=.01, q=.1, viz=TRUE){

    if (min(weights)<0) stop('Error: weights should be positive \n\n')
     w = weights/sum(weights)
     x = data
     idx = order(x)
     x=x[idx]
     w=w[idx]
     Fw=cumsum(w)/sum(w)
# graphical parameters
        ysup <- 5
        top.x <- .25
        top.y <- .4
        top.xx <- 38
        top.yy <- .11

res1=Top_Share(data=x, weights=w, p=p, q=q, method="pareto1")
res2=Top_Share(data=x, weights=w, p=p, q=q, method="gpd")
res3=Top_Share(data=x, weights=w, p=p, q=q, method="epd", epd.direct= TRUE)
idx=which(x>0)   # Keep positive data
if(viz) par(mfrow=c(1,1), mar=c(4, 4, 4, 1))  # bottom, left, top, right
if(viz) plot(log(x[idx]), log(1-Fw)[idx], xlab="log(x)", ylab="log(1-F(x))", cex=.6, col="gray")

u=seq(log(res1$threshold), log(max(x)), length.out=500)
yhat.par1=ppareto1(exp(u),mu=res1$threshold,alpha=res1$alpha)
yhat.par2=pgpd(exp(u),xi=res2$xi,mu=res2$threshold,beta=res2$beta)
yhat.epd=pepd(exp(u)/res3$threshold,gamma=res3$gamma,kappa=res3$kappa,tau=res3$tau)
if(viz){
lines(u,log(1-yhat.par1)+log(q), col="blue", lty=2, lwd=1.5)
lines(u,log(1-yhat.epd)+log(q), col="red", lty=1, lwd=1.5)
lines(u,log(1-yhat.par2)+log(q),col="green", lty=3, lwd=1.5)
legend("topright", legend=c("Pareto 1", "GPD", "EPD"), col=c("blue","green", "red"), lty=c(2,3,1))
}
res90=Top_Share(data, p=p, q=.10, method="pareto1")
if(viz) abline(v=log(res90$threshold), col="lightgrey", lty=2)  # percentile 90
if(viz) legend(log(res90$threshold)-top.x, top.y, legend=expression(italic('q')[90]), cex=.9, bty="n")
res95=Top_Share(data, p=p, q=.05, method="pareto1")
if(viz) abline(v=log(res95$threshold), col="lightgrey", lty=2)  # percentile 95
if(viz) legend(log(res95$threshold)-top.x, top.y, legend=expression(italic('q')[95]), cex=.9, bty="n")
res99=Top_Share(data, p=p, q=.01, method="pareto1")
if(viz) abline(v=log(res99$threshold), col="lightgrey", lty=2)  # percentile 99
legend(log(res99$threshold)-top.x, top.y, legend=expression(italic('q')[99]), cex=.9, bty="n")
}

#' Table of top shares (using three thresholds)
#'
#' @param data numeric, a vector of observations
#' @param weight numeric, a vector of weights (default is equal weights)
#' @param p numeric, the probability level (default 0.01 for the top 1\%)
#' @param q numeric (possibly a vector), the probability level to model a Pareto distribution (default \code{c(0.1, 0.05, 0.01)} for the top 10\%, 5\% and 1\%)
#' @param md logical \code{TRUE} to have the table in a markdown format (default is \code{FALSE})
#' @return top shares table (in a Markdown format if \code{md=TRUE})
#' @references Charpentier & Flachaire (2019) \emph{Pareto Models for Top Incomes } \href{https://hal.archives-ouvertes.fr/hal-02145024}{hal-02145024}
#' @examples
#' ################ Top Share on 200 Simulated Data
#' set.seed(123)
#' x <- rpareto1(200, 1, alpha=1.5)
#' w <- rgamma(200, 10, 10)
#' Table_Top_Share(data = x, weights = w, q = c(.1,.05))
#' ################ Top Share on Synthetic Data
#' url_1 <- "https://github.com/freakonometrics/TopIncomes/raw/master/dataframe_yw_1.csv"
#' df <- read.table(url_1,sep=";",header=TRUE)
#' Table_Top_Share(data = df$y, weights = df$w)
Table_Top_Share = function(data, weights=rep(1,length(data)), p=.01, q=c(.1,.05,.01), md=FALSE, verbose=FALSE){

    if (min(weights)<0) stop('Error: weights should be positive \n\n')
     w = weights/sum(weights)
     x = data
     idx = order(x)
     x=x[idx]
     w=w[idx]
     Fw=cumsum(w)/sum(w)

    thresholds=q
    res0=Top_Share(data=x, weights=w, p=.01,q=thresholds, method="edf")
    res1=Top_Share(data=x, weights=w, p=.01,q=thresholds, method="pareto1")
    res2=Top_Share(data=x, weights=w, p=.01,q=thresholds, method="gpd")
    res3=Top_Share(data=x, weights=w, p=.01,q=thresholds, method="epd")

    alpha=rbind(100*(1-res1$share.pareto),res1$alpha,res2$alpha,res3$alpha)
    rownames(alpha) <- c("cutoff","Pareto 1", "GPD", "EPD")
    colnames(alpha)=paste("top",round(100*(1-q)),"%",sep="")
    res=rbind(100*(1-res1$share.pareto),res0$index,res1$index,res2$index,res3$index)
    rownames(res) <- c("cutoff","edf","Pareto 1", "GPD", "EPD")
    colnames(res)=paste("top",round(100*(1-q)),"%",sep="")


if(verbose){
cat("----- top shares ------------------\n")
print(res)
cat("----- alpha (tail index) ----------\n")
print(alpha)
}

if(md){
    if (!require("knitr")) install.packages("knitr")
    require(knitr)
    dt=data.frame(res*100)[-1,]
    names(dt)=paste("top",round(100*(1-q)),"%",sep="")
    rownames(dt)=c("EDF","Pareto_1","GPD","EPD")
    cat(kable(dt,caption="Top Share (in percent)",
    bootstrap_options = c("striped", "hover"), full_width = FALSE))

    dt=data.frame(alpha)[-1,]
    names(dt)=paste("top",round(100*(1-q)),"%",sep="")
    rownames(dt)=c("Pareto_1","GPD","EPD")
    cat(kable(dt,caption="Tail Index (alpha)",
    bootstrap_options = c("striped", "hover"), full_width = FALSE))
}

    return(list(TopShare=res,TailIndex=alpha))}

#' Top Income plot
#'
#' @param data numeric, a vector of observations
#' @param weight numeric, a vector of weights (default is equal weights)
#' @param p numeric, probability level (default 0.01)
#' @param thr numeric, vector of probability levels to model a Pareto distribution (default is \code{seq(.85,.999,by=.001)} from 0.85 up to 0.999)
#' @param tail_index logical to plot the tail index (default \code{TRUE})
#' @return one or two graphs (depending on \code{tail==TRUE})
#' @references Charpentier & Flachaire (2019) \emph{Pareto Models for Top Incomes } \href{https://hal.archives-ouvertes.fr/hal-02145024}{hal-02145024}
#' @examples
#' ################ Top Incomes on Synthetic Data
#' url_1 <- "https://github.com/freakonometrics/TopIncomes/raw/master/dataframe_yw_1.csv"
#' df <- read.table(url_1,sep=";",header=TRUE)
#' \dontrun{Top_Incomes(data = df$y, weights = df$w)}
Top_Incomes = function(data, weights=rep(1,length(data)), p=.01, thr=seq(.85,.999,by=.001), tail_index = TRUE,...){

	 if (min(weights)<0) stop('Error: weights should be positive \n\n')
     w = weights/sum(weights)
     x = data
     idx = order(x)
     x=x[idx]
     w=w[idx]
     Fw=cumsum(w)/sum(w)

thr=round(thr,10)
tail=matrix(0,NROW(thr),7)
tis.index=matrix(0,NROW(thr),7)
tis.alpha=matrix(0,NROW(thr),7)
for(i in 1:NROW(thr)) {

    res1=Top_Share(data=x, weights=w, p=p, q=1-thr[i], method="pareto1")
    res2=Top_Share(data=x, weights=w, p=p, q=1-thr[i], method="gpd")
    res3=Top_Share(data=x, weights=w, p=p, q=1-thr[i], method="epd", epd.direct=TRUE)
    res4=Top_Share(data=x, weights=w, p=p, method="edf")

    tis.index[i,1]=res1$threshold     # threshold y0
    tis.index[i,2]=res1$k          # k largest observations
    tis.index[i,3]=thr[i]             # quantile threshold
    tis.index[i,4]=res1$index
    tis.index[i,5]=res2$index
    tis.index[i,6]=res3$index
    tis.index[i,7]=res4$index

    tis.alpha[i,1]=res2$threshold           # threshold y0
    tis.alpha[i,2]=res2$k          # k largest observations
    tis.alpha[i,3]=thr[i]             # quantile threshold
    tis.alpha[i,4]=res1$alpha
    tis.alpha[i,5]=res2$alpha
    tis.alpha[i,6]=res3$alpha
    tis.alpha[i,7]=0

}

# graphical parameters
ysup <- 5
top.x <- .25
top.y <- .4
top.xx <- 38
top.yy <- .11

if(tail_index){
    plot(tis.alpha[,2],tis.alpha[,4], ylim=c(0,8), type="b", cex=.75, pch=3, main="MLE estimates of the tail index", xlab="k largest values", ylab="tail index (alpha)", col="blue")
lines(tis.alpha[,2],tis.alpha[,4], col="blue", type="l", cex=.75)
lines(tis.alpha[,2],tis.alpha[,5], col="green", type="p", cex=.75, pch=2)
lines(tis.alpha[,2],tis.alpha[,5], col="green", type="l", cex=.75)
lines(tis.alpha[,2],tis.alpha[,6], col="red", type="b", cex=.75, pch=1)
lines(tis.alpha[,2],tis.alpha[,6], col="red", type="l", cex=.75)
abline(v=tis.alpha[(tis.alpha[,3]==.90),2], col="lightgray", lty=2) # 10% top obs
abline(v=tis.alpha[(tis.alpha[,3]==.95),2], col="lightgray", lty=2) #  5% top obs
abline(v=tis.alpha[(tis.alpha[,3]==.99),2], col="lightgray", lty=2) #  1% top obs
legend("topright", legend=c("Pareto 1 (Hill estimator)","GPD", "EPD"), col=c("blue", "green", "red"), pch=c(3,2,1), lty=1,bty="n")

legend(tis.alpha[(tis.alpha[,3]==.90),2]-top.xx,top.yy, legend=expression(italic('q')[90]), cex=.9, bty="n")
legend(tis.alpha[(tis.alpha[,3]==.95),2]-top.xx,top.yy, legend=expression(italic('q')[95]), cex=.9, bty="n")
legend(tis.alpha[(tis.alpha[,3]==.99),2]-top.xx,top.yy, legend=expression(italic('q')[99]), cex=.9, bty="n")
}

plot(tis.index[,2],tis.index[,4], ylim=c(0,.4), type="b", cex=.75, pch=3, main="Top 1% share", xlab="k largest values", ylab="share", col="blue")
lines(tis.index[,2],tis.index[,4], col="blue", type="l", cex=.75)
lines(tis.index[,2],tis.index[,5], col="green", type="p", cex=.75, pch=2)
lines(tis.index[,2],tis.index[,5], col="green", type="l", cex=.75)
lines(tis.index[,2],tis.index[,6], col="red", type="b", cex=.75, pch=1)
lines(tis.index[,2],tis.index[,6], col="red", type="l", cex=.75)
lines(tis.index[,2],tis.index[,7], col="gray", type="l", cex=.75)
abline(v=tis.index[(tis.index[,3]==.90),2], col="lightgray", lty=2) # 10% top obs
abline(v=tis.index[(tis.index[,3]==.95),2], col="lightgray", lty=2) #  5% top obs
abline(v=tis.index[(tis.index[,3]==.99),2], col="lightgray", lty=2) #  1% top obs

legend("topright", legend=c("Pareto 1","GPD", "EPD"), col=c("blue", "green", "red"), pch=c(3,2,1),lty=1,bty="n")
legend(tis.index[(tis.index[,3]==.90),2]-top.xx,top.yy, legend=expression(italic('q')[90]), cex=.9, bty="n")
legend(tis.index[(tis.index[,3]==.95),2]-top.xx,top.yy, legend=expression(italic('q')[95]), cex=.9, bty="n")
legend(tis.index[(tis.index[,3]==.99),2]-top.xx,top.yy, legend=expression(italic('q')[99]), cex=.9, bty="n")
}

#' import HMisc
#' import randtoolbox
#' import knitr
#' importFrom("graphics", "abline", "legend", "lines", "par", "plot")
#' importFrom("stats", "integrate", "optim", "optimise", "runif", "var")
#' importFrom("utils", "install.packages")
