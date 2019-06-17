## NP estimator of the conditional probability of cure
probcure <- function(x,
                     t,
                     d,
                     dataset,
                     x0,
                     h,
                     local = TRUE,
                     conflevel = 0L,
                     bootpars = if (conflevel == 0 && !missing(h)) NULL else controlpars()) {
    dfr <-
        if (missing(dataset))
            na.omit(data.frame(x, t, d))
        else 
            na.omit(dataset[, c(deparse(substitute(x)), deparse(substitute(t)), deparse(substitute(d)))])
    names(dfr) <- c("x", "t", "d")
    dfr$x <- as.numeric(dfr$x)
    dfr$t <- as.numeric(dfr$t)
    dfr$d <- as.integer(dfr$d)
    nrow <- dim(dfr)[1]
    ordx0 <- order(x0)
    x0 <- as.numeric(x0[ordx0])
    lx0 <- length(x0)
    if (!local && missing(h)) warning("Option 'local = FALSE' overridden: with missing 'h' a local bootstrap bandwidth is computed")
    if (missing(h)) {
        sm <- bootpars$hsmooth
        h <-
            if (sm > 1)
                probcurehboot(x, t, d, dfr, x0, bootpars)$hsmooth
            else
                probcurehboot(x, t, d, dfr, x0, bootpars)$h
    }
    else {
        if (local) {
            if (lx0 != length(h)) stop("When 'local = TRUE', 'x0' and 'h' must have the same length")
            h <- as.numeric(h[ordx0])
        }
        else {
            h <- as.numeric(h)
        }
    }
    lh <- length(h)
    if (conflevel < 0 | conflevel > 1) stop("'conflevel' must be a number between 0 and 1")
    dfr <- dfr[order(dfr$t, 1 - dfr$d),]
    q <- .Call("probcurenp0",
              dfr$t, ## observed times
              dfr$x, ## observed covariate values
              dfr$d, ## observed censoring indicator
              nrow, ## sample size
              x0, ## covariate values
              lx0, ## number of covariate values
              h, ## bandwidths
              lh, ## number of bandwidths
              local, ## a boolean specifying whether the bandwidths are local or not
              PACKAGE = "npcure")
    if (!local) {
        names(q) <- paste("h", as.character(round(h, 8)), sep = "")
        q <- as.list(q)
    }
    if (conflevel > 0) {
        B <- bootpars$B
        fpilot <- bootpars$fpilot
        if (is.null(fpilot)) {
            pilot <- hpilot(dfr$x, x0, bootpars$nnfrac)
        }
        else
            pilot <- do.call(fpilot, c(list(x0 = x0), bootpars$dots))
        band <- .Call("probcurenp0confband",
                      dfr$t,
                      dfr$x,
                      dfr$d,
                      nrow,
                      x0,
                      lx0,
                      h,
                      lh,
                      1 - (1 - conflevel)/2,
                      B,
                      pilot,
                      q,
                      local,
                      PACKAGE = "npcure")
        if (local)
            names(band) <- c("lower", "upper")
        else {
            names(band) <- paste("h", as.character(round(h, 8)), sep = "")
            for (i in 1:lh) {
                names(band[[i]]) <- c("lower", "upper")
            }
        }
        structure(list(type = "cure", local = local, h = h, x0 = x0, q = q, conf = band, conflevel = conflevel), class = "npcure")
    }
    else
        structure(list(type = "cure", local = local, h = h, x0 = x0, q = q), class = "npcure")
}
