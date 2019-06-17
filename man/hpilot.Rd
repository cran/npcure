\name{hpilot}

\alias{hpilot}

\title{Compute the Pilot Bandwidth for the Nonparametric Estimators of
  Cure Probability and Latency
}

\description{This function computes local pilot bandwidths for the
  nonparametric estimators of the probability of cure and the latency
  function.
}

\usage{hpilot(x, x0, nnfrac = 0.25)
}

\arguments{
  \item{x}{A numeric vector of observed covariate values.
  }
  \item{x0}{A numeric vector specifying a grid of covariate values.
  }
  \item{nnfrac}{A numeric value giving the fraction of the sample size
    that determines the order of the nearest neighbor. This is taken as
    \code{floor(length(x)*nnfrac)}. The default is 0.25.
  }
}

\details{The function computes a data-driven local pilot bandwidth,
  required for the bootstrap bandwidth selector of the nonparametric
  estimators of the cure rate and latency functions. Simulations
  in López-Cheda et al. (2017) show that the choice of pilot bandwidth
  has small effect on the bootstrap bandwidth. This pilot bandwidth only
  depends on the sample size and the distribution of the covariate x
  (see López-Cheda, 2018):

  \deqn{g(x0) = 0.5(d_k^{+}(x0) + d_k^{-}(x0))(100/n)^{1/9}}{g(x0) =
  0.5 * (d_k^{+}(x0) + d_k^{-}(x0)) * (100/n)^{1/9}}

  where \eqn{d_k^{+}(x0)} and \eqn{d_k^{-}(x0)} are the distances from
  \eqn{x0} to the \eqn{k}-th nearest neighbor on the right and the left,
  respectively, and \eqn{k} is a suitable integer depending on the
  sample size \eqn{n}. If there are not at least \eqn{k} neighbors on
  the right or on the left, we use \eqn{d_k^{+}(x0) = d_k^{-}(x0)}. The
  default value of \eqn{k} is \eqn{n/4}. The order
  \eqn{n^{-1/9}}{n^(-1/9)} satisfies the conditions in Theorem 1 of Li
  and Datta (2001) and coincides with the order obtained by Cao and
  González-Manteiga (1993) for the uncensored case.
}

\value{A numeric vector of local pilot bandwidths corresponding to each
  one of the values of the grid of covariate values given by \code{x0}.
}

\references{Cao R., González-Manteiga W. (1993). Bootstrap methods in
  regression smoothing. Journal of Nonparametric Statistics, 2:
  379-388. \url{https://doi.org/10.1080/10485259308832566}.
  
  Li, G., Datta, S. (2001). A bootstrap approach to nonparametric
  regression for right censored data. Annals of the Institute of
  Statistical Mathematics, 53(4):
  708-729. \url{https://doi.org/10.1023/A:1014644700806}.
  
  López-Cheda A. (2018). \emph{Nonparametric Inference in Mixture Cure
    Models}. PhD dissertation, Universidade da Coruña. Spain.
  
  López-Cheda, A., Cao, R., Jácome, M. A., Van Keilegom,
  I. (2017). Nonparametric incidence estimation and bootstrap bandwidth
  selection in mixture cure models. Computational Statistics & Data
  Analysis, 105:
  144–165. \url{https://doi.org/10.1016/j.csda.2016.08.002}.
  
  López-Cheda, A., Jácome, M. A., Cao, R. (2017). Nonparametric latency
  estimation for mixture cure models. TEST, 26:
  353–376. \url{https://doi.org/10.1007/s11749-016-0515-1}.
}

\author{Ignacio López-de-Ullibarri [aut, cre],
  Ana López-Cheda [aut],
  Maria Amalia Jácome [aut]
}

\seealso{\code{\link{controlpars}}, \code{\link{latencyhboot}},
  \code{\link{probcurehboot}}
}

\examples{
## Some artificial data
set.seed(123)
n <- 50
x <- runif(n, -2, 2) ## Covariate values
y <- rweibull(n, shape = .5*(x + 4)) ## True lifetimes
c <- rexp(n) ## Censoring values
p <- exp(2*x)/(1 + exp(2*x)) ## Probability of being susceptible
u <- runif(n)
t <- ifelse(u < p, pmin(y, c), c) ## Observed times
d <- ifelse(u < p, ifelse(y < c, 1, 0), 0) ## Uncensoring indicator
data <- data.frame(x = x, t = t, d = d)

## Computing pilot bandwidths for covariate values -1, -0.8, ..., 1
## by taking the 5-th nearest neighbor
hpilot(data$x, x0 = seq(-1, 1, by = .2), nnfrac = .05)
}