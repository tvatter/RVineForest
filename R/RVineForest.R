RVineForestCopSelect <- function(data, familyset = c(1,3:6),
                                 selectioncrit = "AIC",
                                 indeptest = FALSE, level = 0.05,
                                 trunclevel = NA, progress = FALSE,
                                 weights = NA, se = FALSE, method = "itau",
                                 rotations = TRUE, nvines = 100, cores = 1) {

    ## register parallel backend
    if (cores != 1 | is.na(cores)) {
        if (is.na(cores))
            cores <- max(1, detectCores() - 1)
        if (cores > 1) {
            cl <- makeCluster(cores)
            registerDoParallel(cl)
            on.exit(try(stopCluster(), silent = TRUE))
            on.exit(try(closeAllConnections(), silent = TRUE), add = TRUE)
        }
    }

    d <- ncol(data)
    RVFstructures <- RVineMatrixSample(d, size = nvines, naturalOrder = FALSE)

    if (cores > 1) {
        i <- NA  # dummy for CRAN checks
        RVF <- foreach(i = 1:nvines,
                       .export = c("RVineCopSelect")) %dopar%
            RVineCopSelect(reSample(data), familyset, RVFstructures[[i]],
                           selectioncrit, indeptest, level,
                           trunclevel, se, rotations, method, cores = 1)
    } else {
        RVF <- lapply(1:nvines, function(i)
            RVineCopSelect(reSample(data), familyset, RVFstructures[[i]],
                           selectioncrit, indeptest, level,
                           trunclevel, se, rotations, method, cores = 1))
    }

    class(RVF) <- "RVineForest"
    return(RVF)
}

RVineForestLogLik <- function(data, RVF, separate = FALSE) {
  loglik <- sapply(RVF, function(RVM) RVinePDF(data, RVM))
  loglik <- apply(loglik, 1, function(x) log(mean(x)))
  if (separate == TRUE) {
    return(loglik)
  } else {
    return(sum(loglik))
  }
}

RVineForestPDF <- function(data, RVF) {
    exp(RVineForestLogLik(data, RVF, separate = TRUE))
}

reSample <- function(x, n = nrow(x)) {
  return(x[sample.int(nrow(x), n, replace = TRUE),])
}
