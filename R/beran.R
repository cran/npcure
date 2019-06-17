## Beran's estimator of conditional survival function
beran <- function(x,
                  t,
                  d,
                  dataset,
                  x0,
                  h,
                  local = TRUE,
                  testimate = NULL,
                  conflevel = 0L,
                  cvbootpars = if (conflevel == 0 && !missing(h)) NULL else controlpars()) {
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
    if (missing(h)) {
        sm <- cvbootpars$hsmooth
        h <-
            if (sm > 1)
                berancv(x, t, d, dfr, x0, cvbootpars)$hsmooth
            else
                berancv(x, t, d, dfr, x0, cvbootpars)$h
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
    dfr <- dfr[order(dfr$t, 1 - dfr$d),]
    if (conflevel < 0 | conflevel > 1) stop("'conflevel' must be a number between 0 and 1")
    ltestimate <- length(testimate)
    if (!is.null(testimate)) {
        if (conflevel == 0)
            warning("When 'testimate' is not NULL don't use the survival estimates for plotting")
        else
            stop("For plotting confidence bands 'testimate' must be NULL")
        testimate <- as.numeric(testimate)
    }
    S <- .Call("berannp0",
               dfr$t,
               dfr$x,
               dfr$d,
               nrow,
               x0,
               lx0,
               h,
               lh,
               local,
               testimate,
               ltestimate,
               PACKAGE = "npcure")
    if (local) { 
        names(S) <- paste("x", as.character(round(x0, 8)), sep = "")
    }
    else {
        names(S) <- paste("h", as.character(round(h, 8)), sep = "")  
        for (i in 1:lh) {
            if (lx0 == 1)
                S[[i]] <- list(S[[i]])
            names(S[[i]]) <- paste("x", as.character(round(x0, 8)), sep = "")
        }
    }
    if (conflevel > 0) {
        B <- cvbootpars$B
        fpilot <- cvbootpars$fpilot
        if (is.null(fpilot)) {
            pilot <- hpilot(dfr$x, dfr$x, cvbootpars$nnfrac)
        }
        else
            pilot <- do.call(fpilot, c(list(x0 = dfr$x), cvbootpars$dots))
        probcurepilot <- as.numeric(probcure(x, t, d, dfr, dfr$x, pilot)$q)
        band <- .Call("berannp0confband",
                      dfr$t,
                      dfr$x,
                      dfr$d,
                      nrow,
                      x0,
                      lx0,
                      h,
                      lh,
                      pilot,
                      probcurepilot,
                      1 - (1 - conflevel)/2,
                      B,
                      S,
                      local,
                      PACKAGE = "npcure")
        if (local) {
            names(band) <- paste("x", as.character(round(x0, 8)), sep = "")
            for (i in 1:lx0) {
                names(band[[i]]) <- c("lower", "upper")
            }
        }
        else {
            names(band) <- paste("h", as.character(round(h, 8)), sep = "")
            for (i in 1:lh) {
                names(band[[i]]) <- paste("x", as.character(round(x0, 8)), sep = "")
                for (j in 1:lx0) {
                    names(band[[i]][[j]]) <- c("lower", "upper")
                }
            }
        }
        structure(list(type = "survival", local = local, h = h, x0 = x0, testim = dfr$t, S = S, conf = band, conflevel = conflevel), class = "npcure")
    }
    else {
        structure(list(type = "survival", local = local, h = h, x0 = x0, testim = if (is.null(testimate)) dfr$t else testimate, S = S), class = "npcure")
    }
}
