\name{controlpars}

\alias{controlpars}

\title{Control Values for the Bootstrap or Cross-validation
}

\description{This function returns a list of values for the control
  parameters of the functions of the package that use the bootstrap or
  cross-validation.
}

\usage{controlpars(B = 999L, hbound = c(0.1, 3), hl = 100L, hsave =
FALSE, nnfrac = 0.25, fpilot = NULL, qt = 0.75, hsmooth = 1L, \dots)
}

\arguments{
  \item{B}{An integer giving the number of bootstrap resamples.
  }
  \item{hbound}{A numeric vector of length 2 specifying the minimum
    (default, 0.1) and maximum (default, 3), respectively, of the
    initial grid of bandwidths as a multiple of the standardized
    interquartile range of the covariate values.
  }
  \item{hl}{A numeric value giving the length of the initial grid of
    bandwidths. The default is 100.
  }
  \item{hsave}{A logical value specifying if the grids of bandwidths
    must be saved as a component of the list returned by the
    \code{berancv}, \code{latencyboot} and \code{probcurehboot}
    functions. The default is \code{FALSE}.
  }
  \item{nnfrac}{A numeric value giving the fraction of the sample size
    that determines the order of the nearest neighbor used when choosing
    the pilot bandwidth. The default is 0.25.
  }
  \item{fpilot}{A function name or \code{NULL}. If \code{NULL}, the
    default, the \code{hpilot} function is used for computing a pilot
    bandwidth in case that one is needed. If not \code{NULL}, it must be
    the name of a user-defined function (given as a function name or as
    a character string). This function must necessarily have an argument
    \code{x0}, playing the same role than in \code{hpilot}, and must
    return a value of the same length than \code{x0}. If \code{fpilot}
    has more arguments, they are passed through the \code{\dots}
    argument (see below).
  }
  \item{qt}{In bandwidth selection for the latency estimator (see
    \code{latencyhboot}), a numeric value specifiying the order of a
    quantile of the observed times. It determines the right boundary
    of the integration interval in the computation of the ISE (the left
    boundary is 0). The default is 0.75 (third quartile).
  }
  \item{hsmooth}{An integer. Its value controls whether the bandwidths
    selected by the \code{latencyhboot} and \code{probcurehboot} and
    \code{berancv} function should be smoothed, and, if so, the degree
    of smoothing. The smoothing consists in computing a centered moving
    average of the unsmoothed vector of bandwidths returned by default
    by \code{latencyhboot} and \code{probcurehboot}. The value of
    \code{hsmooth} is the number of terms used to compute the
    average. The default is 1L, which means that no smoothing is done.
  }
  \item{\dots}{Arguments of \code{fpilot}, if \code{fpilot} is not
    \code{NULL}.
  }
}

\details{The output of \code{controlpars} is a list of control
  parameters required by the package functions which use the bootstrap
  or cross-validation. This is mainly the case of the \code{berancv}
  function, which computes a cross-validation bandwidth for Beran's
  estimator of survival, and of the \code{latencyhboot} and
  \code{probcurehboot} functions, which compute the bootstrap bandwidth
  selectors of the estimators of the latency and the probability of
  cure, respectively. Since these functions are indirectly called by,
  respectively, the \code{beran}, \code{latency} and \code{probcure}
  functions when their \code{h} argument is missing, the output of
  \code{controlpars} is also the expected (and default) way of passing
  to them the parameters for bandwidth selection.

  Additionally, \code{controlpars} is used by \code{beran},
  \code{latency} and \code{probcure} to set the number of bootstrap
  resamples and the value of \code{nnfrac} (see above) when confidence
  intervals are computed. The \code{testcov} function also uses it for
  setting the number of bootstrap resamples.
}

\value{A list whose components are the arguments of the function, their
  defaults being replaced with the values the function was called with.
}

\author{Ignacio López-de-Ullibarri [aut, cre],
  Ana López-Cheda [aut],
  Maria Amalia Jácome [aut]
}

\seealso{\code{\link{beran}},\code{\link{berancv}},
  \code{\link{hpilot}}, \code{\link{latency}},
  \code{\link{latencyhboot}}, \code{\link{probcure}},
  \code{\link{probcurehboot}},\code{\link{testcov}}
}
