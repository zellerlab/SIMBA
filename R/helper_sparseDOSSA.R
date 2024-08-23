#!/usr/bin/Rscript
### SIMBA
### Simulation of Microbiome data with Biological Accuracy
### Morgan Essex - MDC Berlin
### Jakob Wirbel - EMBL Heidelberg
### 2020 GNU GPL 3.0

# Helper functions to simulate data according to sparseDOSSA
# based on the R package found under
# https://github.com/biobakery/sparseDOSSA

#' @keywords internal
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
      sim.feat.l <- funcMakeFeature(
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
      sim.feat <- funcRoundMatrix(mtrxData=sim.feat)
      rownames(sim.feat) <- rownames(feat)
      colnames(sim.feat) <- colnames(feat)

      # implant markers
      for (m.ix in marker.idx){
        new.bug <- funcSpikeNewBug(label,
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
    vdMu[i] <- funcGetMu(vdExp[i], exp(sd(vdLogCur)))
  }
  return(list(exp=vdExp, mu=vdMu, sd = exp(vdLogSD),
              percentZero = vdPercentZero,
              dAverageReadDepth = mean(colSums(feat.filt)),
              iFeatureCount = nrow(feat.filt)))
}

# ##################
# Some functions from sparseDOSSA, because they are not exported :/
# see https://github.com/biobakery/sparseDOSSA for the source code
# The code has been provided under the MIT license:
# Copyright (c) <2014> <Huttenhower Lab,Department of Biostatistics, Harvard School of Public Health>
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ##################


#' @keywords internal
funcGetMu <- function(dEx, dSD) {
  if ((log(dSD)^2 - log(dEx)) > 745) 
    warning(paste("mu = e^(-x), but x =", 
                  round(log(dSD)^2 - log(dEx), 2), "will give mu=0."))
  return(exp(-1 * (((log(dSD, exp(1))^2)/2) - log(dEx, exp(1)))))
}

#' @keywords internal
funcGetExp <- function(dMu, dSD){
  return(exp(log(dMu, exp(1)) + 0.5 * (log(dSD, exp(1))^2)))
}

#' @keywords internal
funcRoundMatrix <- function(mtrxData, fZeroInflated = FALSE){
  if (fZeroInflated) {
    vdSubCount = intersect(which(mtrxData < 1), which(mtrxData > 0))
    mtrxData[vdSubCount] = 1
  }
  return(round(mtrxData))
}

#' @keywords internal
funcMakeFeature <- function(
    vdMu, vdSD, vdPercentZero, iNumberSamples, iMinNumberCounts, 
    iMinNumberSamples, mdLogCorr = diag(length(vdSD)), vdTruncateThreshold = NA, 
    fZeroInflate = TRUE, fVerbose = FALSE){
  if (!fZeroInflate) {
    vdPercentZero = rep(0, length(vdMu))
  }
  if (length(vdMu) != length(vdSD)) {
    stop("vdMu and vdSD must have equal length")
  }
  vdExpCal = sapply(seq_along(vdMu), function(i) funcGetExp(vdMu[i], 
                                                            vdSD[i]))
  mdFeature_base = func_zero_inflate(vdLogMean = log(vdMu), 
                                     vdPercentZeroInflated = vdPercentZero, int_number_samples = iNumberSamples, 
                                     vdLogSD = log(vdSD), mdLogCorr = mdLogCorr, viThreshold = vdTruncateThreshold)
  mdFeature = matrix(NA, ncol = ncol(mdFeature_base), nrow = nrow(mdFeature_base))
  for (k in seq_len(ncol(mdFeature))) {
    mdFeature[, k] = funcUpdateDistributionToExpectation(vdFeatures = mdFeature_base[, 
                                                                                     k], dExp = vdExpCal[k])
  }
  mdFeature = apply(mdFeature, 2, funcForceMinCountsInMinSamples, 
                    iMinNumberCounts = iMinNumberCounts, iMinNumberSamples = iMinNumberSamples)
  mdFeature_base = apply(mdFeature_base, 2, funcForceMinCountsInMinSamples, 
                         iMinNumberCounts = iMinNumberCounts, iMinNumberSamples = iMinNumberSamples)
  vdMean = apply(mdFeature, 2, mean)
  vdMean_base = apply(mdFeature_base, 2, mean)
  return(list(Feature = mdFeature, Feature_base = mdFeature_base, 
              Exp = vdMean, ExpCal = vdExpCal, Exp_base = vdMean))
}

#' @keywords internal
funcSpikeNewBug <- function(vdCurMetadata, vdCurData, 
                            multiplier, fZeroInflated = TRUE){
  metadata_average = funcGetRowMetric(vdCurMetadata, mean)
  metadata_sigma = funcGetSD(vdCurMetadata)
  data_average = mean(vdCurData[which(vdCurData > 0)])
  data_sigma = sd(vdCurData[which(vdCurData > 0)])
  if (is.na(data_sigma)) {
    return(NULL)
  }
  if (!fZeroInflated) {
    data_average = mean(vdCurData)
    data_sigma = sd(vdCurData)
  }
  metadata_sigma[metadata_sigma == 0] = 1
  if (data_sigma == 0) {
    data_sigma = 1
  }
  liOccurenceFilter = vdCurData
  liOccurenceFilter[liOccurenceFilter > 0] = 1
  if (!fZeroInflated) {
    liOccurenceFilter[seq_along(liOccurenceFilter)] = 1
  }
  scaled_metadata = (vdCurMetadata - metadata_average)/metadata_sigma
  if (data_sigma != 0) {
    scaled_metadata = scaled_metadata * data_sigma
  }
  scaled_metadata = scaled_metadata + data_average
  if (!(is.null(nrow(vdCurMetadata)) || (nrow(vdCurMetadata) == 
                                         1))) {
    scaled_metadata = colSums(scaled_metadata)
  }
  vdSpikedBug = vdCurData + (multiplier * scaled_metadata * 
                               liOccurenceFilter)
  iNumberRows = nrow(vdCurMetadata)
  if (is.null(iNumberRows)) {
    iNumberRows = 1
  }
  vdSpikedBug = vdSpikedBug/((iNumberRows * multiplier) + 1)
  vdSpikedBug[which(vdSpikedBug < 0)] = 0
  return(vdSpikedBug)
}

#' @keywords internal
funcGetRowMetric <- function(lxValues, funcMethod){
  if (is.null(dim(lxValues)[1])) {
    return(funcMethod(lxValues))
  }
  return(apply(lxValues, 1, funcMethod))
}

#' @keywords internal
funcGetSD <- function(lxValues){
  if (is.null(dim(lxValues)[1])) {
    return(sd(lxValues, na.rm = TRUE))
  }
  return(apply(lxValues, 1, sd, na.rm = TRUE))
}

#' @keywords internal
func_zero_inflate <- function(vdLogMean, 
                              vdPercentZeroInflated, 
                              int_number_samples, 
                              vdLogSD, mdLogCorr = diag(length(vdLogSD)), 
                              viThreshold = NA) {
  mdFeature = funcTruncatedRLNorm(int_number_samples, vdLogMean, 
                                  vdLogSD, mdLogCorr = mdLogCorr, viThreshold = viThreshold)
  miZeroLocations = as.logical(sapply(vdPercentZeroInflated, 
                                      rbinom, n = int_number_samples, size = 1))
  mdFeature[miZeroLocations] = 0
  if (any(colSums(mdFeature) == 0)) {
    zero_cols <- which(colSums(mdFeature) == 0)
    one_rows <- sample(nrow(mdFeature), length(zero_cols), 
                       replace = TRUE)
    mdFeature[cbind(one_rows[zero_cols], zero_cols)] = 1
  }
  return(mdFeature)
}

#' @keywords internal
funcTruncatedRLNorm <- function(iNumberMeasurements, 
                                vdLogMean, vdLogSD, 
                                mdLogCorr = diag(length(vdLogSD)), 
                                viThreshold = NA) 
{
  if (length(vdLogMean) != length(vdLogSD)) {
    stop("vdLogMean and vdLogSD must have equal length")
  }
  mdLogSD = diag(x = vdLogSD, nrow = length(vdLogSD))
  mdLogVar = mdLogSD %*% mdLogCorr %*% mdLogSD
  if (length(viThreshold) == 1 && is.na(viThreshold)) {
    viThreshold <- rep(Inf, length(vdLogSD))
  }
  mdFeature <- exp(tmvtnorm::rtmvnorm(n = iNumberMeasurements, 
                                      mean = vdLogMean, sigma = mdLogVar, 
                                      upper = log(viThreshold), 
                                      algorithm = "gibbs"))
  return(mdFeature)
}

#' @keywords internal
funcUpdateDistributionToExpectation <- function(vdFeatures, dExp){
  vdUpdateIndices = which(vdFeatures > 0)
  dDifference = dExp - mean(vdFeatures)
  if (abs(dDifference) > 0) {
    dCounts = ceiling(abs(dDifference * length(vdFeatures)))
    if (dDifference > 0) {
      if (length(vdUpdateIndices) > 1) {
        vdUpdateIndices = funcSample(
          vdUpdateIndices, 
          dCounts, replace = TRUE, 
          prob = vdFeatures[vdUpdateIndices]/sum(vdFeatures[vdUpdateIndices]))
      }
      update.freq = table(vdUpdateIndices)
      vdFeatures[as.numeric(names(update.freq))] = as.numeric(update.freq) + 
        vdFeatures[as.numeric(names(update.freq))]
    }
    else if (dDifference < 0) {
      for (iIndex in 1:dCounts) {
        viGreaterThan1 = which(vdFeatures > 1)
        if (length(viGreaterThan1) > 1) {
          iUpdateIndex = funcSample(
            viGreaterThan1, 1, 
            prob = vdFeatures[viGreaterThan1]/sum(vdFeatures[viGreaterThan1]))
          vdFeatures[iUpdateIndex] = vdFeatures[iUpdateIndex] - 
            1
        }
      }
    }
  }
  return(vdFeatures)
}


#' @keywords internal
funcSample <- function(x, size, replace = FALSE, prob = NULL){
  iLength = length(x)
  if (iLength == 0) {
    stop("funcSample:: Can not sample from length of 0 vector.")
  }
  else if (iLength == 1) {
    if (size == 1) {
      return(x)
    }
    if (replace) {
      return(rep(x, size))
    }
    stop("funcSample:: Can not create a vector of size > 1 from 1 entry without replacement.")
  }
  else {
    return(sample(x = x, size = size, replace = replace, 
                  prob = prob))
  }
}

#' @keywords internal
funcForceMinCountsInMinSamples <- function(vdFeature, 
                                           iMinNumberCounts = 0, 
                                           iMinNumberSamples = 0){
  if ((iMinNumberCounts + iMinNumberSamples) == 0) {
    return(vdFeature)
  }
  iNonZeroSamplesNeeded = min(0, iMinNumberSamples - length(which(vdFeature == 
                                                                    0)))
  if (iNonZeroSamplesNeeded > 0) {
    dSignalMean = round(mean(vdFeature[which(vdFeature > 
                                               0)]))
    viUpdate = funcSample(which(vdFeature == 0), iNonZeroSamplesNeeded)
    vdFeature[viUpdate] = dSignalMean
  }
  iNeededExtraValues = iMinNumberSamples - length(which(vdFeature >= 
                                                          iMinNumberCounts))
  vdProbabilities = vdFeature/sum(vdFeature)
  viAddIndices = which(vdFeature > 0)
  while (iNeededExtraValues > 0) {
    iIndexAdd = viAddIndices
    if (length(iIndexAdd) > 1) {
      iIndexAdd = funcSample(viAddIndices, 1, prob = vdProbabilities[viAddIndices])
    }
    vdFeature[iIndexAdd] = vdFeature[iIndexAdd] + 1
    iNeededExtraValues = iMinNumberSamples - length(which(vdFeature >= 
                                                            iMinNumberCounts))
  }
  return(vdFeature)
}
