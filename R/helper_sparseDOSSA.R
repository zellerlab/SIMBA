#!/usr/bin/Rscript
### SIMBA
### Simulation of Microbiome data with Biological Accuracy
### Morgan Essex - MDC Berlin
### Jakob Wirbel - EMBL Heidelberg
### 2020 GNU GPL 3.0

# Helper functions to simulate data according to sparseDOSSA
# based on the R package found under
# https://github.com/biobakery/sparseDOSSA


simulate.sparseDOSSA <- function(feat, meta, sim.out, sim.params){

  prop.markers <- sim.params$prop.markers
  class.balance <- sim.params$class.balance
  ab.scale <- sim.params$ab.scale
  repeats <- sim.params$repeats
  no.marker.feat <- round(prop.markers * nrow(feat))
  num.sample <- ncol(feat)

  parameters <- fit.parameters(feat)

  # check very low sd values
  if (any(parameters$sd <= 1)){
    d <- min(parameters$sd)
    add <- (1 + 1e-10) - d
    warning("Very low standard deviation in some cases. ",
            sprintf(fmt='%e', add),
            ' will be added to the sd of all features')
    parameters$sd <- parameters$sd + add
  }
  pb <- progress_bar$new(total = length(ab.scale)*repeats)
  for (a in seq_along(ab.scale)){
    for (r in seq_len(repeats)){
      # generate label
      label <- rep(-1, ncol(feat))
      s <- sample(num.sample,
                  size=round(num.sample*class.balance))
      names(label) <- colnames(feat)
      label[s] <- +1 # switch random half of the samples to positive class

      # sample markers
      marker.idx <- sample(rownames(feat), size=no.marker.feat)

      # generate features
      read.depth <- sample(colSums(feat), size = length(label),
                           replace = TRUE)
      sim.feat.l <- sparseDOSSA:::funcMakeFeature(
        vdMu = parameters$mu,
        vdSD = parameters$sd,
        vdPercentZero = parameters$percentZero,
        mdLogCorr = diag(length(parameters$sd)),
        iNumberSamples = num.sample,
        iMinNumberCounts = 10,
        iMinNumberSamples = 5,
        vdTruncateThreshold =
          rep(Inf, length(parameters$sd)),
        fZeroInflate = TRUE)
      sim.feat <- t(sim.feat.l$Feature)
      sim.feat <- sparseDOSSA:::funcRoundMatrix(mtrxData=sim.feat)
      rownames(sim.feat) <- rownames(feat)
      colnames(sim.feat) <- colnames(feat)

      # implant markers
      for (m.ix in marker.idx){
        new.bug <- sparseDOSSA:::funcSpikeNewBug(label,
                                                 sim.feat[m.ix,],
                                                 ab.scale[a])
        sim.feat[m.ix,] <- new.bug
      }

      # save in h5 file
      h5.subdir <- paste0('ab', a, '_rep', r)
      stopifnot(h5createGroup(sim.out, h5.subdir))
      h5write(label, sim.out, paste0(h5.subdir, '/labels'))
      h5write(sim.feat, sim.out, paste0(h5.subdir, '/features'))
      h5write(marker.idx, sim.out, paste0(h5.subdir, '/marker_idx'))
      h5write(colnames(sim.feat), sim.out,
              paste0(h5.subdir, '/sample_names'))
      h5write(rownames(sim.feat), sim.out,
              paste0(h5.subdir, '/feature_names'))

      pb$tick()
    }
  }
}

#' Generate feature matrix
#' @keywords internal
funcMakeFeature = function(vdMu, vdSD, vdPercentZero, iNumberSamples,
                           iMinNumberCounts, iMinNumberSamples,
                           mdLogCorr = diag(length(vdSD)),
                           vdTruncateThreshold = rep(Inf, length(vdSD)),
                           fZeroInflate = TRUE, fVerbose = FALSE){
  # If not zero inflated
  if(!fZeroInflate){vdPercentZero = rep(0, length(vdMu))}

  # Check that vdLogMean and vdLogSD are same length
  if (length(vdMu) != length(vdSD)){
    stop("vdMu and vdSD must have equal length")
  }
  # Expectation of the features
  vdExpCal = sapply(seq_along(vdMu), function(i) sparseDOSSA:::funcGetExp(vdMu[i],vdSD[i]))

  # Generate features
  #mdFeature_base = func_zero_inflate(vdLogMean = log(vdMu), vdPercentZeroInflated=vdPercentZero, int_number_samples=iNumberSamples, vdLogSD=log(vdSD), mdLogCorr=mdLogCorr, viThreshold=vdTruncateThreshold)
  mdFeature_base = t(exp(tmvtnorm::rtmvnorm( iNumberSamples, vdMu,
                                             sigma = diag(vdSD)%*%mdLogCorr%*%diag(vdSD),
                                             upper = log(vdTruncateThreshold), algorithm = "gibbs" )))
  mdFeature_base = mdFeature_base * t(sapply( 1:length(vdPercentZero),
                                              function(x) sample( 0:1, iNumberSamples, replace = T,
                                                                  prob = c(vdPercentZero[x],1-vdPercentZero[x])
                                              )))

  # Extra useful measurements, the true and expected means
  vdMean = apply(mdFeature_base, 2, mean)
  vdMean_base = apply(mdFeature_base, 2, mean)

  return(list(Feature = mdFeature_base, Feature_base = mdFeature_base, Exp=vdMean, ExpCal = vdExpCal, Exp_base = vdMean))
}

#' Fit the parameters from the original data
#' @keywords internal
fit.parameters <- function(feat.filt){
  # Mean of the nonzero data
  vdExp <- vector(length=nrow(feat.filt))
  # Mean of the logged nonzero data
  vdMu <- vector(length=nrow(feat.filt))
  # Standard deviation of the logged nonzero data
  vdLogSD <- vector(length=nrow(feat.filt))
  # Percent zeros in data
  vdPercentZero <- vector(length=nrow(feat.filt))

  # Calculate parameters for each feature (column)
  for(i in seq_len(nrow(feat.filt))){
    # Get the percent zero before removing zeros for other measurements
    vdCur <- as.numeric(as.vector(as.matrix(feat.filt[i,])))
    vdPercentZero[i] <- mean(vdCur == 0)

    # Measure expectation of the feature with zeros
    vdExp[i] <- mean(vdCur)

    # Remove zeros from data
    vdCur <- vdCur[which(vdCur!=0)]

    #### Note
    #### rlnorm needs a mean and sd from a logged rlnorm distribution which would match this
    #### without further manipulation. The "mean" in the formula is actually not the expectation
    #### The expectation is e^mu+.5*sd^2 so this is always more than mu.

    # Log nonzero data
    vdLogCur <- log(vdCur)
    vdLogSD[i] <- sd(vdLogCur)
    vdMu[i] <- sparseDOSSA:::funcGetMu(vdExp[i], exp(sd(vdLogCur)))
  }
  return(list(exp=vdExp, mu=vdMu, sd = exp(vdLogSD),
              percentZero = vdPercentZero,
              dAverageReadDepth = mean(colSums(feat.filt)),
              iFeatureCount = nrow(feat.filt)))
}
