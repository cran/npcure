## Function for setting control parameters
controlpars <- function(B = 999L,
                        hbound = c(0.1, 3),
                        hl = 100L,
                        hsave = FALSE,
                        nnfrac = 0.25,
                        fpilot = NULL,
                        qt = 0.75,
                        hsmooth = 1L,
                        ...) {
    if (length(hbound) != 2 || any(hbound <= 0))
        stop("Incorrect 'hbound' parameter")
    if (hl < 1)
        stop("Incorrect 'hl' parameter")
    if (nnfrac <= 0 || nnfrac >= 1)
        stop("Incorrect 'nnfrac' parameter")
    if (qt <= 0 || nnfrac >= 1)
        stop("Incorrect 'qt' parameter")
    if (hsmooth < 1)
        stop("Incorrect 'hsmooth' parameter")
    list(B = as.integer(B), hbound = as.numeric(hbound), hl = as.integer(hl), hsave = hsave, nnfrac = nnfrac, fpilot = fpilot, qt = qt, hsmooth = as.integer(hsmooth), dots = list(...))
}

## Bootstrap bandwidth selector for the estimator of cure probability
probcurehboot <- function(x,
                          t,
                          d,
                          dataset,
                          x0,
                          bootpars = npcure::controlpars()) {
    dfr <-
        if (missing(dataset))
            na.omit(data.frame(x, t, d))
        else 
            na.omit(dataset[, c(deparse(substitute(x)), deparse(substitute(t)), deparse(substitute(d)))])
    names(dfr) <- c("x", "t", "d")
    dfr$x <- as.numeric(dfr$x)
    dfr$t <- as.numeric(dfr$t)
    dfr$d <- as.integer(dfr$d)
    dfr <- dfr[order(dfr$t, 1 - dfr$d),]
    x0 <- as.numeric(sort(x0))
    lx0 <- length(x0)
    B <- bootpars$B
    hbound <- IQR(dfr$x)/1.349*bootpars$hbound
    lhgrid <- bootpars$hl
    steph <- (hbound[2]/hbound[1])^(1/lhgrid)
    hgrid <- as.numeric(hbound[1]*steph^seq(0, lhgrid, length.out = lhgrid))
    nrow <- dim(dfr)[1]
    fpilot <- bootpars$fpilot
    if (is.null(fpilot)) {
        pilot <- npcure::hpilot(dfr$x, x0, bootpars$nnfrac)
    }
    else
        pilot <- do.call(fpilot, c(list(x0 = x0), bootpars$dots))
    probcurepilot <- as.numeric(probcure(x, t, d, dfr, x0, pilot)$q)
    h <- .Call("probcurenp0hboot",
               dfr$t,
               dfr$x,
               dfr$d,
               nrow,
               x0,
               lx0,
               hgrid,
               lhgrid,
               pilot,
               probcurepilot,
               B,
               PACKAGE = "npcure")
    result  <- list(type = c("Bootstrap bandwidth", "cure"), x0 = x0, h = h)
    sm <- bootpars$hsmooth
    if (sm > 1) {
        if (sm >= lx0)
            warning("The number of covariate values is probably too small for smoothing the selected bandwidths with hsmooth=", sm)
        result$hsmooth <- zoo::rollapply(h, sm, mean, partial = TRUE, align = "center")
    }
    if (bootpars$hsave)
        result$hgrid <- hgrid
    structure(result, class = "npcure")
}

## Bootstrap bandwidth selector for the estimator of latency
latencyhboot <- function(x,
                         t,
                         d,
                         dataset,
                         x0,
                         bootpars = npcure::controlpars()) {
    dfr <-
        if (missing(dataset))
            na.omit(data.frame(x, t, d))
        else
            na.omit(dataset[, c(deparse(substitute(x)), deparse(substitute(t)), deparse(substitute(d)))])
    names(dfr) <- c("x", "t", "d")
    dfr$x <- as.numeric(dfr$x)
    dfr$t <- as.numeric(dfr$t)
    dfr$d <- as.integer(dfr$d)
    dfr <- dfr[order(dfr$t, 1 - dfr$d),]
    x0 <- as.numeric(sort(x0))
    lx0 <- length(x0)
    B <- bootpars$B
    boundh <- IQR(dfr$x)/1.349*bootpars$hbound
    lhgrid <- bootpars$hl
    tmax <- quantile(dfr$t, bootpars$qt)
    steph <- (boundh[2]/boundh[1])^(1/lhgrid)
    hgrid <- as.numeric(boundh[1]*steph^seq(0, lhgrid, length.out = lhgrid))
    nrow <- dim(dfr)[1]
    fpilot <- bootpars$fpilot
    if (is.null(fpilot)) {
        pilot <- npcure::hpilot(dfr$x, dfr$x, bootpars$nnfrac)
        pilotx0 <- npcure::hpilot(dfr$x, x0, bootpars$nnfrac)
    }
    else {
        pilot <- do.call(fpilot, c(list(x0 = dfr$x), bootpars$dots))
        pilotx0 <- do.call(fpilot, c(list(x0 = x0), bootpars$dots))
    }
    probcurepilot <- as.numeric(probcure(x, t, d, dfr, dfr$x, pilot)$q)
    latencypilot <- latency(x, t, d, dfr, x0, pilotx0)$S
    h <- .Call("latencynp0hboot",
               dfr$t,
               dfr$x,
               dfr$d,
               nrow,
               x0,
               lx0,
               hgrid,
               lhgrid,
               pilot,
               probcurepilot,
               latencypilot,
               B,
               tmax,
               PACKAGE = "npcure")
    result  <- list(type = c("Bootstrap bandwidth", "latency"), x0 = x0, h = h)
    sm <- bootpars$hsmooth
    if (sm > 1) {
        if (sm >= lx0)
            warning("The number of covariate values is probably too small for smoothing the selected bandwidths with hsmooth=", sm)
        result$hsmooth <- zoo::rollapply(h, sm, mean, partial = TRUE, align = "center")
    }
    if (bootpars$hsave)
        result$hgrid <- hgrid
    structure(result, class = "npcure")
}

## Cross-validation bandwidth selector for Beran's survival estimator
berancv <- function(x,
                    t,
                    d,
                    dataset,
                    x0,
                    cvpars = npcure::controlpars()) {
    dfr <-
        if (missing(dataset))
            na.omit(data.frame(x, t, d))
        else 
            na.omit(dataset[, c(deparse(substitute(x)), deparse(substitute(t)), deparse(substitute(d)))])
    names(dfr) <- c("x", "t", "d")
    dfr$x <- as.numeric(dfr$x)
    dfr$t <- as.numeric(dfr$t)
    dfr$d <- as.integer(dfr$d)
    dfr <- dfr[order(dfr$t, 1 - dfr$d),]
    x0 <- as.numeric(sort(x0))
    lx0 <- length(x0)
    nrow <- dim(dfr)[1]
    boundh <- IQR(dfr$x)/1.349*cvpars$hbound
    lhgrid <- cvpars$hl
    steph <- (boundh[2]/boundh[1])^(1/lhgrid)
    hgrid <- as.numeric(boundh[1]*steph^seq(0, lhgrid, length.out = lhgrid))
    h <- .Call("berannp0cv",
               dfr$t,
               dfr$x,
               dfr$d,
               nrow,
               x0,
               lx0,
               hgrid,
               lhgrid,
               PACKAGE = "npcure")
    result  <- list(type = c("Cross-validation bandwidth", "survival"), x0 = x0, h = h)
    sm <- cvpars$hsmooth
    if (sm > 1) {
        if (sm >= lx0)
            warning("The number of covariate values is probably too small for smoothing the selected bandwidths with hsmooth=", sm)
        result$hsmooth <- zoo::rollapply(h, sm, mean, partial = TRUE, align = "center")
    }
    if (cvpars$hsave)
        result$hgrid <- hgrid
    structure(result, class = "npcure")
}
