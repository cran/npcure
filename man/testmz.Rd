\name{testmz}

\alias{testmz}

\title{Test of Maller-Zhou}

\description{This function carries out the nonparametric test of Maller and
  Zhou (1992).
}

\usage{testmz(t, d, dataset)
}

\arguments{
  \item{t}{If \code{dataset} is missing, a numeric object giving the
    the observed times. If \code{dataset} is a data frame, it is
    interpreted as the name of the variable corresponding to the observed
    times in the data frame.
  }
  \item{d}{If \code{dataset} is missing, an integer object giving the
    values of the uncensoring indicator. Censored observations must be
    coded as 0, uncensored ones as 1. If \code{dataset} is a data frame,
    it is interpreted as the name of the variable corresponding to the 
    uncensoring indicator.
  }
  \item{dataset}{An optional data frame in which the variables named in
    \code{x}, \code{t} and \code{d} are interpreted. If it is missing,
    \code{x}, \code{t} and \code{d} must be objects of the workspace.
  }
}

\details{The function implements Maller and Zhou's (1992) method to
  test the null hypothesis \eqn{H_{0} : \tau_{F_{0}} > \tau_{G}}{H_0 :
  tau_F_0 > tau_G} vs. \eqn{H_{1} : \tau_{F_{0}} \leq \tau_{G}}{H_1 :
  tau_F_0 <= tau_G}, where \eqn{\tau_{F_{0}}}{tau_F_0} and
  \eqn{\tau_{G}}{tau_G} are the supports of, respectively, the
  distribution function of the survival time of the uncured and the
  distribution function of the censoring time.
}

\value{An object of S3 class 'npcure'. Formally, a list of components:
  \item{type}{The constant character string c("test", "Maller-Zhou").
  }
  \item{pvalue}{The p-value of the test.
  }
  \item{aux}{A list of components: \code{statistic}, the test
    statistic, \code{n} the sample size, \code{delta}, the difference
    between the largest observed time \eqn{t_{n}}{t_n} and the largest
    uncensored time \eqn{t_{n}^{\star}}{t_n^*}, and \code{interval}, a
    vector giving the range between
    \eqn{\max(0,t_{n}^{\star}-\texttt{delta})}{max(0, t_n^* - delta)} and
    \eqn{t_{n}^{\star}}{t_n^*}.
  }
}

\references{Maller R. A., Zhou S. (1992). Estimating the proportion of
  immunes in a censored sample. Biometrika, 79:
  731-739. \url{https://doi.org/10.1093/biomet/79.4.731}. 
}

\author{Ignacio López-de-Ullibarri [aut, cre],
  Ana López-Cheda [aut],
  Maria Amalia Jácome [aut]
}

\seealso{\code{\link{latency}}, \code{\link{probcure}}
}

\examples{
## Some artificial data
set.seed(123)
n <- 50
y <- qweibull(runif(n)*pweibull(2, shape = 2), shape = 2) ## True lifetimes
c <- qexp(runif(n)*pexp(2.5)) ## Censoring values
u <- runif(n)
## Probability of being susceptible is constantly equal to .5
t <- ifelse(u < .5, pmin(y, c), c) ## Observed times
d <- ifelse(u < .5, ifelse(y < c, 1, 0), 0) ## Uncensoring indicator
data <- data.frame(t = t, d = d)

## Maller-Zhou test
testmz(t, d, data)
}

\keyword{survival}

\keyword{nonparametric}