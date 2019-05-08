## Test of Maller and Zhou (1992)
testmz <- function(t, d, dataset) {
    dfr <-
        if (missing(dataset))
            na.omit(data.frame(t, d))
        else 
            na.omit(dataset[, c(deparse(substitute(t)), deparse(substitute(d)))])
    names(dfr) <- c("t", "d")
    dfr$t <- as.numeric(dfr$t)
    dfr$d <- as.integer(dfr$d)
    nrow <- dim(dfr)[1]
    dfr <- dfr[order(dfr$t, 1 - dfr$d),]
    tmax <- max(dfr$t)
    tmaxuncens <- max(dfr$t[dfr$d == 1])
    delta <- tmax - tmaxuncens
    left <- max(0, tmaxuncens - delta)
    ninterval <-  sum(dfr$t >= left & dfr$t <= tmaxuncens)
    pvalue <- (1 - ninterval/nrow)^nrow
    structure(list(type = c("test", "Maller-Zhou"), pvalue = pvalue, aux = list(statistic = ninterval, n = nrow, delta = delta, interval = c(left, tmaxuncens))), class = "npcure")
}

## Test of covariate significance: auxiliary function
testcompute <- function(x,
                        eta,
                        nrow,
                        B,
                        type) { 
    cm <- ks <- numeric(length = B)
    meaneta <- mean(eta)
    switch(type,
    { ## binary      
        mod0 <- sort(unique(x))[1]
        n0 <- length(which(x == mod0))
        testsample0 <- sum((x == mod0)*eta) - n0*meaneta
        testsample1 <- sum(eta) - nrow*meaneta
        cmsample <- n0*testsample0^2 + (nrow - n0)*testsample1^2
        kssample <- max(abs(testsample0), abs(testsample1)) 
        for (b in 1:B) { 
            xboot <- sample(x, nrow, TRUE)
            etaboot <- sample(eta, nrow, TRUE)
            meanetaboot <- mean(etaboot)
            n0boot <- length(which(xboot == mod0))
            testboot0 <- sum((xboot == mod0)*etaboot) - n0boot*meanetaboot
            testboot1 <- sum(etaboot) - nrow*meanetaboot
            cm[b] <- n0boot*testboot0^2 + (nrow - n0boot)*testboot1^2
            ks[b] <- max(abs(testboot0), abs(testboot1))
        }
    },
    { ## discrete
        gridx <- sort(unique(x))
        lgridx <- length(gridx)
        sumx <- testsample <- numeric(length = lgridx)
        for (i in 1:lgridx) {
            sumx[i] <- sum(x == gridx[i])
            indx <- as.integer(x == gridx[i])
            if (i == 1)
                indxsum <- indx
            else
                indxsum <- indx + indxsumold
            indxsumold <- indxsum
            testsample[i] <- sum((eta - meaneta)*indxsum)
        }
        cmsample <- sum(sumx*testsample^2)
        kssample <- max(abs(testsample))
        for (b in 1:B) {
            xboot <- sample(x, nrow, TRUE)
            etaboot <- sample(eta, nrow, TRUE)
            meanetaboot <- mean(etaboot)
            gridx <- sort(unique(xboot))
            lgridx <- length(gridx)
            sumxboot <- testboot <- numeric(length = lgridx)
            for (i in 1:lgridx) {
                sumxboot[i] <- sum(xboot == gridx[i])
                indx <- as.integer(xboot == gridx[i])
                if (i == 1)
                    indxsum <- indx
                else
                    indxsum <- indx + indxsumold
                indxsumold <- indxsum
                testboot[i] <- sum((etaboot - meanetaboot)*indxsum)
            }
            cm[b] <- sum(sumxboot*testboot^2)
            ks[b] <- max(abs(testboot))
        }
    },
    { ## continuous
        testsample <- testboot <- numeric(length = nrow)
        for (i in 1:nrow) {
            condx <- x <= x[i]
            testsample[i] <- sum((eta - meaneta)*condx)
        }
        cmsample <- sum(testsample^2)
        kssample <- max(abs(testsample))
        for (b in 1:B) {
            xboot <- sample(x, nrow, TRUE)
            etaboot <- sample(eta, nrow, TRUE)
            meanetaboot <- mean(etaboot)
            for (i in 1:nrow) {
                condxboot <- xboot <= xboot[i]
                testboot[i] <- sum((etaboot - meanetaboot)*condxboot)
            }
            cm[b] <- sum(testboot^2)
            ks[b] <- max(abs(testboot))
        }
    },
    { ## qualitative
        mods <- unique(x)
        lmods <- length(mods)
        perm <- rbind(1:lmods, permute::allPerms(lmods))
        nperm <- dim(perm)[1]
        cmsamplemod <- kssamplemod <- numeric(length = nperm)
        cmmod <- ksmod <- matrix(0, nperm, B)
        testsample <- matrix(0, nperm, lmods)
        indxsum <- integer(length = lmods)
        for (i in 1:nperm) {
            for (j in 1:lmods) {
                indx <- as.integer(x %in% mods[perm[i, 1:j]])
                testsample[i, j] <- sum((eta - meaneta)*indx)
                indxsum[j] <- sum(x == mods[perm[i, j]])
            }
            cmsamplemod[i] <- indxsum %*% testsample[i, ]^2
            kssamplemod[i] <- max(abs(testsample[i,]))
        }
        cmsample <- max(cmsamplemod)
        kssample <- max(kssamplemod)
        for (b in 1:B){
            xboot <- sample(x, nrow, TRUE) 
            etaboot <- sample(eta, nrow, TRUE)    
            testboot <- matrix(0, nperm, lmods)
            indxsum <- integer(length = lmods)
            meanetaboot <- mean(etaboot)
            for (i in 1:nperm) {
                for (j in 1:lmods) {
                    indx <- as.integer(xboot %in% mods[perm[i, 1:j]])
                    testboot[i, j] <- sum((etaboot - meanetaboot)*indx)
                    indxsum[j] <- sum(xboot == mods[perm[i, j]])
                }
                cmmod[i, b] <- indxsum %*% testboot[i, ]^2
                ksmod[i, b] <- max(abs(testboot[i, ]))
            }
            cm[b] <- max(cmmod[, b])
            ks[b] <- max(ksmod[, b])
        }
    })
    list(cms = cmsample, cmb = cm, kss = kssample, ksb = ks)
}

## Test of covariate significance: main function
testcov <- function(x,
                    t,
                    d,
                    dataset,
                    bootpars = npcure::controlpars()) {
    dfr <-
        if (missing(dataset))
            dfr <- na.omit(data.frame(x, t, d))
        else
            dfr <- na.omit(dataset[, c(deparse(substitute(x)), deparse(substitute(t)), deparse(substitute(d)))])
    names(dfr) <- c("x", "t", "d")
    nmods <- length(unique(dfr$x))
    xnum <- is.numeric(dfr$x) | (!is.numeric(dfr$x) & nmods == 2)
    if (!is.numeric(dfr$x))
        dfr$x <- as.numeric(as.factor(dfr$x))
    dfr$t <- as.numeric(dfr$t)
    dfr$d <- as.integer(dfr$d)
    nrow <- dim(dfr)[1]
    dfr <- dfr[order(dfr$t, 1 - dfr$d),]
    B <- bootpars$B
    tau <- dfr$t[tail(which(dfr$d == 1), 1)]
    G <- 1
    i <- 1
    while ((dfr$t[i] <= tau) && (i <= nrow)) {
        G <- G*(((nrow - i)/(nrow - i + 1))^(1 - dfr$d[i]))
        i <- i + 1
    }   
    eta <- ifelse ((dfr$d == 0) & (dfr$t > tau), 1/G, 0)
    type <- if (xnum) {
                if (nmods == 2) 1 ## binary
                else {
                    if (nmods <= 50) 2 ## discrete
                    else 3 ## continuous
                }
            }    
            else 4 ## qualitative
    test <- testcompute(dfr$x, eta, nrow, B, type)
        prenamex <- as.character(substitute(x))
        namex <-
            if (missing(dataset)) {
                if (length(prenamex) == 1)
                    prenamex[1]
                else {
                    if (prenamex[1] == "$")
                        prenamex[3] 
                    else {
                        if (exists(prenamex[4]))                   
                            eval(parse(text = prenamex[4]))
                        else {
                            partcond <- eval(parse(text = prenamex[4]), envir = eval(parse(text = prenamex[2])), enclos = .GlobalEnv)
                            if (is.character(partcond))
                                partcond
                            else
                                prenamex[4]
                        }
                    }
                }
            }
            else
                prenamex
    structure(list(type = c("test", "Covariate"), x = namex, CM = list(statistic = test$cms/nrow^2, pvalue = mean(test$cms < (test$cmb + .Machine$double.eps^0.5))), KS = list(statistic = nrow^(-0.5)*test$kss, pvalue = mean(test$kss < (test$ksb + .Machine$double.eps^0.5)))), class = "npcure")
}
