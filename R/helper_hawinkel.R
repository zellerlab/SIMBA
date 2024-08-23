#!/usr/bin/Rscript
### SIMBA
### Simulation of Microbiome data with Biological Accuracy
### Morgan Essex - MDC Berlin
### Jakob Wirbel - EMBL Heidelberg
### 2020 GNU GPL 3.0

# Helper functions to simulate data according to Hawinkel et al.
# based on the scripts found here:
# https://users.ugent.be/~shawinke/ABrokenPromise/02_dataGeneration.html
#
# with a lot of modifications (as best as i understood the code...)

#' # wrapper for the Dirichlet simulations
#' @keywords internal
simulate.dirmult <- function(feat, meta, sim.out, sim.params){

  ab.scale <- sim.params$ab.scale
  prop.markers <- sim.params$prop.markers
  class.balance <- sim.params$class.balance
  repeats <- sim.params$repeats

  no.marker.feat <- round(prop.markers * nrow(feat))
  el.feat.names <- rownames(feat)
  num.sample <- ncol(feat)

  libSizesOrig <- colSums(feat)

  # calculate rhos/phis
  message("++ Starting to estimate parameters ",
          "for the Dirichlet distributions")
  piMoMs <- piMoM4Wald(t(feat))
  thetaMoMs <- weirMoM4Wald(t(feat))
  message("++ Finished estimating parameters")

  pb <- progress_bar$new(total = length(ab.scale)*repeats)
  for (a in seq_along(ab.scale)){
    for (r in seq_len(repeats)){
      # generate label
      label <- rep(-1, ncol(feat))
      s <- sample(num.sample,
                  size=round(num.sample*class.balance))
      names(label) <- colnames(feat)
      label[s] <- +1
      # sample markers
      marker.idx <- sample(el.feat.names, size=no.marker.feat)

      sim.feat <- simulate.markers.dirmult(piMoMs, thetaMoMs, libSizesOrig,
                                          label, marker.idx, ab.scale[a])
      # save data in H5-file
      h5.subdir <- paste0('ab', a, '_rep', r)
      stopifnot(h5createGroup(sim.out, h5.subdir))
      h5write(label, sim.out, paste0(h5.subdir, '/labels'))
      h5write(sim.feat, sim.out, paste0(h5.subdir, '/features'))
      h5write(marker.idx, sim.out, paste0(h5.subdir, '/marker_idx'))
      h5write(colnames(sim.feat), sim.out,
              paste(h5.subdir, '/sample_names', sep=''))
      h5write(rownames(sim.feat), sim.out,
              paste(h5.subdir, '/feature_names', sep=''))
      pb$tick()
    }
  }
}

#' # wrapper for the beta-binomial simulations
#' @keywords internal
simulate.betabin <- function(feat, meta, sim.out, sim.params){
  ab.scale <- sim.params$ab.scale
  prop.markers <- sim.params$prop.markers
  class.balance <- sim.params$class.balance
  repeats <- sim.params$repeats

  test.package("TailRank")
  test.package("SpiecEasi")

  no.marker.feat <- round(prop.markers * nrow(feat))
  el.feat.names <- rownames(feat)
  num.sample <- ncol(feat)

  libSizesOrig <- colSums(feat)

  message("++ Calculating the underlying correlation structure\n",
          "\tPlease note that this may take a while!")
  # calculate correlation structure
  n.cores <- parallel::detectCores()
  if (n.cores < 16){
    message("++ Only ", n.cores, " cores detected. It would be better to",
            " have more cores available!")
  }
  res <- SpiecEasi::spiec.easi(t(feat), verbose=TRUE, pulsar.params=list(ncores=n.cores))
  correlation <- SpiecEasi::getOptCov(res)

  # calculate rhos/phis
  message("++ Starting to estimate parameters ",
          "for the beta binomial distributions")
  piMoMs <- piMoM4Wald(t(feat))
  thetaMoMs <- weirMoM4Wald(t(feat))
  message("++ Finished estimating parameters")

  pb <- progress_bar$new(total = length(ab.scale)*repeats)
  for (a in seq_along(ab.scale)){
    for (r in seq_len(repeats)){
      # generate label
      label <- rep(-1, ncol(feat))
      s <- sample(num.sample,
                  size=round(num.sample*class.balance))
      names(label) <- colnames(feat)
      label[s] <- +1
      # sample markers
      marker.idx <- sample(el.feat.names, size=no.marker.feat)

      sim.feat <- simulate.markers.betabin(piMoMs, thetaMoMs, libSizesOrig,
                                           label, marker.idx, ab.scale[a],
                                           as.matrix(correlation))

      # save data in H5-file
      h5.subdir <- paste0('ab', a, '_rep', r)
      stopifnot(h5createGroup(sim.out, h5.subdir))
      h5write(label, sim.out, paste0(h5.subdir, '/labels'))
      h5write(sim.feat, sim.out, paste0(h5.subdir, '/features'))
      h5write(marker.idx, sim.out, paste0(h5.subdir, '/marker_idx'))
      h5write(colnames(sim.feat), sim.out, paste(h5.subdir, '/sample_names',
                                                 sep=''))
      h5write(rownames(sim.feat), sim.out, paste(h5.subdir, '/feature_names',
                                                 sep=''))
      pb$tick()
    }
  }
}

#' # wrapper for the negative-binomial simulations
#' @keywords internal
simulate.negbin <- function(feat, meta, sim.out, sim.params){
  
  ab.scale <- sim.params$ab.scale
  prop.markers <- sim.params$prop.markers
  class.balance <- sim.params$class.balance
  repeats <- sim.params$repeats
  correlation <- sim.params$correlation
  # check for packages
  test.package("MASS")
  if (correlation){
    test.package("SpiecEasi")
    message("++ Calculating the underlying correlation structure\n",
            "\tPlease note that this may take a while!")
    # calculate correlation structure
    n.cores <- parallel::detectCores()
    if (n.cores < 16){
      message("++ Only ", n.cores, " cores detected. It would be better to",
              " have more cores available!")
    }
    res <- SpiecEasi::spiec.easi(t(feat), verbose=TRUE, pulsar.params=list(ncores=n.cores))
    correlation <- SpiecEasi::getOptCov(res)
  } else {
    correlation <- NULL
  }

  no.marker.feat <- round(prop.markers * nrow(feat))
  el.feat.names <- rownames(feat)
  num.sample <- ncol(feat)

  # sample libsizes
  libSizesOrig <- colSums(feat)
  # calculate rhos/phis
  message("++ Starting to calculate parameters ",
          "for the negative binomial distributions")
  logLibSizes <- log(colSums(feat))
  nbfitlist <- list()
  for (x in rownames(feat)){
    y <- feat[x,]
    res <- try(MASS::glm.nb(y ~ offset(logLibSizes),
                            link = "log"), silent = TRUE)
    nbfitlist[[x]] <- res
  }
  fit.success <- vapply(nbfitlist, length, FUN.VALUE = double(1))
  if (any(fit.success == 1)) {
    nbfitlist[[which(fit.success==1)]] <- NULL
    removed.sp <- names(fit.success)[which(fit.success==1)]
    el.feat.names <- setdiff(el.feat.names, removed.sp)
    if (!is.null(correlation)){
      idx <- which(names(fit.success) %in% removed.sp)
      correlation <- correlation[-idx, -idx]
    }
  }
  phiMLEs <- vapply(nbfitlist, function(y) {1/y$theta}, FUN.VALUE = double(1))
  rhoMLEs <- vapply(nbfitlist, function(y) {y$coef[1]}, FUN.VALUE = double(1))
  rhoMLEs <- rhoMLEs/sum(rhoMLEs)
  message("++ Finished calculating parameters")

  pb <- progress_bar$new(total = length(ab.scale)*repeats)
  for (a in seq_along(ab.scale)){
    for (r in seq_len(repeats)){
      # generate label
      label <- rep(-1, ncol(feat))
      s <- sample(num.sample,
                  size=round(num.sample*class.balance))
      names(label) <- colnames(feat)
      label[s] <- +1
      # sample markers
      marker.idx <- sample(el.feat.names, size=no.marker.feat)

      sim.feat <- simulate.markers.negbin(phiMLEs, rhoMLEs, libSizesOrig,
                                          label, marker.idx, ab.scale[a],
                                          correlation)
      # save data in H5-file
      h5.subdir <- paste0('ab', a, '_rep', r)
      stopifnot(h5createGroup(sim.out, h5.subdir))
      h5write(label, sim.out, paste0(h5.subdir, '/labels'))
      h5write(sim.feat, sim.out, paste0(h5.subdir, '/features'))
      h5write(marker.idx, sim.out, paste0(h5.subdir, '/marker_idx'))
      h5write(colnames(sim.feat), sim.out, paste(h5.subdir, '/sample_names',
                                                 sep=''))
      h5write(rownames(sim.feat), sim.out, paste(h5.subdir, '/feature_names',
                                                 sep=''))
      pb$tick()
    }
  }
}

#' # wrapper for the actual data generation
#' @keywords internal
simulate.markers.negbin <- function(phis, rhos, libSizesOrig, label,
                                    markers, ab.scale,
                                    correlation=NULL){

  sampleSizes <- as.vector(table(label))
  dmData <- matrix(NA_integer_, nrow = sum(sampleSizes), ncol  = length(phis))
  rownames(dmData) <- names(sort(label))
  colnames(dmData) <- names(phis)
  samSizeSeq <- c(0L, cumsum(sampleSizes))
  libSizes <- sample(libSizesOrig, size = length(label), replace = TRUE)

  # add fold change to rhos
  rhos.all <- cbind(rhos, rhos)
  rhos.all[markers,2] <- rhos.all[markers,2] * ab.scale  # Add fold change up
  rhos.all[,2] <- rhos.all[,2]/sum(rhos.all[,2])

  # generate counts
  for(nRun in seq_along(sampleSizes)){
    indSel <- (samSizeSeq[nRun] + 1L):samSizeSeq[nRun + 1L]
    if (!is.null(correlation)){
      dmData[indSel, ] <- rmvnegbin(n=sampleSizes[nRun],
                                    mu=t(tcrossprod(rhos.all[, nRun],
                                                    libSizes[indSel])),
                                    ks=1/phis,
                                    Sigma=as.matrix(correlation))
    } else {
      dmData[indSel,] <- mapply(rhos.all[,nRun], phis,
                                FUN=function(rho, phi){
                                  rnbinom(n=sampleSizes[nRun],
                                          mu=rho*libSizes[indSel],
                                          size=1/phi)
                                  })
    }
  }

  dmData <- t(dmData)
  dmData <- dmData[,names(label)]

  return(dmData)
}

#' # wrapper for the actual data generation
#' @keywords internal
simulate.markers.dirmult <- function(pis, theta,
                                     libSizesOrig, label,
                                     markers, ab.scale){

  libSizes <- sample(libSizesOrig, size = length(label), replace = TRUE)
  dirData <- matrix(NA, ncol=length(label), nrow=length(pis),
                   dimnames=list(names(pis), names(label)))
  piDir <- dirData

  # add fold change to rhos
  pis.all <- cbind(pis, pis)
  pis.all[markers,2] <- pis.all[markers,2] * ab.scale  # Add fold change up

  gammas <- pis.all * (1-theta)/theta

  # generate count
  for (a in c(1,2)){
    piDir[,which(label==unique(label)[a])] <-
      t(rDirichlet(n = sum(label==unique(label)[a]), alpha = gammas[,a]))
  }
  x <- which(colSums(piDir) < 1e-12)
  while (length(x) > 0){
    for (i in x){
      piDir[,x] <- rDirichlet(n=1,
                              alpha=gammas[,ifelse(label[i]==unique(label)[1],
                                                   1, 2)])
    }
    x <- which(colSums(piDir) < 1e-12)
  }
  for (x in seq_len(ncol(dirData))){
    dirData[,x] <- rmultinom(n = 1L, size = libSizes[x], prob = piDir[,x])
  }
  return(dirData)
}

#' # wrapper for the actual data generation
#' @keywords internal
simulate.markers.betabin <- function(pis, theta,
                                     libSizesOrig, label,
                                     markers, ab.scale,
                                     correlation){

  sampleSizes <- as.vector(table(label))

  # add fold change to rhos
  pis.all <- cbind(pis, pis)
  pis.all[markers,2] <- pis.all[markers,2] * ab.scale  # Add fold change up
  pis.all <- pis.all/sum(pis.all[,1]) #renormalize


  dmData <- matrix(NA_integer_, nrow = sum(sampleSizes), ncol  = nrow(pis.all))

  rownames(dmData) <- names(sort(label))
  colnames(dmData) <- rownames(pis.all)

  libSizesOrig <- libSizesOrig[names(sort(label))]
  samSizeSeq <- c(0L, cumsum(sampleSizes))

  for(nRun in seq_along(sampleSizes)){
    indSel <- (samSizeSeq[nRun] + 1L):samSizeSeq[nRun + 1L]
    dmData[indSel, ] <- rmvbetabin(n=sampleSizes[nRun], pi=pis.all[, nRun],
                                   libSizes=libSizesOrig[indSel],
                                   Sigma=correlation, theta = theta)
  }

  dmData <- t(dmData)
  dmData <- dmData[,names(label)]
  return(dmData)
}

### generate Dirichlet realisations,
### taken from gtools (identical in MCMCpack)
#' @keywords internal
rDirichlet <- function(n, alpha) {
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  res <- x/as.vector(sm)
  res[res <= 10 * .Machine$double.eps] <- 0
  res
}

### pi vector estimation with Method of Moments
#' @keywords internal
piMoM4Wald <- function(data){
  totalReads <- sum(data, na.rm = TRUE)
  piMom <- colSums(data, na.rm = TRUE)/totalReads
  zeroInds <- abs(piMom) < .Machine$double.eps
  r <- sum(zeroInds)
  rr <- length(piMom) - r
  piMom[!zeroInds] <- piMom[which(piMom != 0)] - r/(rr * 2 * (totalReads + 1))
  piMom[zeroInds] <- 1/(2 * (totalReads + 1))

  return(piMom)
}

###### my version of Weir computation (Method of Moments)
#' @keywords internal
weirMoM4Wald <- function (data, se = FALSE){
  if (missing(data))
    stop("data missing.")
  K <- ncol(data)
  J <- nrow(data)

  totalsample <- sum(data, na.rm = TRUE)
  MoM <- colSums(data, na.rm = TRUE)/totalsample
  Sn <- rowSums(data, na.rm = TRUE)
  auxData <- data / Sn

  MSP <- sum(rowSums(
    (auxData - matrix(rep(MoM, J), nrow = J, ncol = K, byrow = TRUE))^2,
    na.rm = TRUE) * Sn) / (J - 1)
  MSG <- sum(rowSums(auxData * (1 - auxData),
                     na.rm = TRUE) * Sn) / (totalsample - J)

  nc <- (sum(Sn) - sum(Sn^2)/sum(Sn)) / (J - 1)
  MoM.wh <- (MSP - MSG)/(MSP + (nc - 1) * MSG)

  if (se)
  {
    std.er <- sqrt(2 * (1 - MoM.wh)^2 / (J - 1) *
                     ((1 + (nc - 1) * MoM.wh)/nc)^2)
    return(list(theta = MoM.wh, se = std.er))
  } else
  {
    return(MoM.wh)
  }
}

#' helper function
#' @keywords internal
fixInf <- function(data) {
  # hacky way of replacing infinite values with the col max + 1
  if (any(is.infinite(data))) {
    data <- apply(data, 2, function(x) {
      if (any(is.infinite(x))) {
        x[ind <- which(is.infinite(x))] <- NA
        x[ind] <- max(x, na.rm = TRUE) + 1
      }
      x
    })
  }
  data
}

#' helper function to generate correlated counts
#' @keywords internal
rmvnegbin <- function(n, mu, Sigma, ks, ...) {
  Cor <- cov2cor(Sigma)
  SDs <- sqrt(diag(Sigma))
  if (missing(mu))
    stop("mu is required")
  if (dim(mu)[2] != dim(Sigma)[2])
    stop("Sigma and mu dimensions don't match")
  if (missing(ks)) {
    # ks <- unlist(lapply(1:length(SDs), function(i) .negbin_getK(mu[i], SDs[i])))
    stop("Problem with this dataset!")
  }
  d <- dim(mu)[2]
  normd <- mvtnorm::rmvnorm(n, rep(0, d), Sigma = Cor)  #The normal-to-anything framework
  unif <- pnorm(normd)
  data <- t(qnbinom(t(unif), mu = t(mu), size = ks, ...))
  data <- fixInf(data)
  return(data)
}

# # A custom quantile beta-binomial function with `na.rm=TRUE`. Still relies
# on the Tailrank package
#' @keywords internal
qbetabin = function(p, N, u, v) {
  pp <- cumsum(TailRank::dbb(0:N, N, u, v))
  sapply(p, function(x) sum(pp < x, na.rm = TRUE))
}

# # A function to generate correlated multivariate betabinomial data pi: a
# vector of proportions, summing to 1 libSizes: library sizes theta: the
# overdispersion parameter
#' @keywords internal
rmvbetabin = function(n, pi, Sigma, theta, libSizes, ...) {
  Cor <- cov2cor(Sigma)
  SDs <- sqrt(diag(Sigma))
  if (missing(pi))
    stop("pi is required")
  if (length(pi) != dim(Sigma)[1])
    stop("Sigma and pi dimensions don't match")
  if (missing(theta)) {
    stop("No overdispersion parameter supplied")
  }
  d <- length(pi)
  normd <- mvtnorm::rmvnorm(n, rep(0, d), Sigma = Cor)  #The normal-to-anything framework
  unif <- pnorm(normd)
  data <- mapply(unif, rep(pi, each = nrow(unif)), libSizes,
                 FUN = function(u, p, l) {
    alphaPar = p * (1 - theta)/theta
    betaPar = (1 - p) * (1 - theta)/theta
    qbetabin(u, N = l, u = alphaPar, v = betaPar, ...)
  })
  data <- fixInf(data)
  return(data)
}
