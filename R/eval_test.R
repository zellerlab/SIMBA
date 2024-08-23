#!/usr/bin/Rscript
### metaGsim package
### Morgan Essex - MDC Berlin
### Jakob Wirbel - EMBL Heidelberg
### 2020 GNU GPL 3.0

#' @title Evaluate the results of different tests
#'
#' @description This function takes the results of the \link{apply.test}
#' function and calculates different evaluation measures (see Details).
#'
#' @usage eval.test(sim.location, group, res.mat, adjust='BH', alpha = 0.05)
#'
#' @param sim.location file name for the .h5 file containing the simulations
#'
#' @param group character, name of the group within the .h5 file containing
#' the simulated features
#'
#' @param res.mat matrix, output from the \link{apply.test} function
#'
#' @param adjust character, indicate the multiple hypothesis testing to be
#' performed on the P-values, defaults to "BH"
#' 
#' @param alpha numeric, significance threshold at which the test will be
#' evaluated, defaults to \code{0.05}
#'
#' @keywords SIMBA eval.test
#'
#' @details This function will check that tests ran successfully and calculate
#' how well the supplied P-values can distinguish between the true differential
#' abudance features and background features. The measures computed are:
#' \itemize{
#' \item AUROC: measure for the separation of true and background features
#' \item TP: number of true positives detected at \code{alpha}
#' \item FP: number of false positives detected at \code{alpha}
#' \item TN: number of true negatives at \code{alpha}
#' \item FN: number of false negatives detected at \code{alpha}
#' \item PR: precision, calculated as TP/(FP+TP)
#' \item R: recall, calculated as TP/(TP+FN)
#' \item FDR: false discovery rate, calculated as FP/(FP+TP)
#' }
#'
#' @return A dataframe containing the evaluation measures at the given alpha
#' value for the different test runs.
#'
#' @export
eval.test <- function(sim.location, group, res.mat, adjust='BH', alpha = 0.05) {

  res <- check.eval.parameters(sim.location, group, res.mat, adjust, alpha)
  # browser()

  # create data.frame to store the evaluation results
  df.res <- data.frame()
  marker <- res[,'marker']
  stopifnot(length(unique(marker)) == 2)
  stopifnot(all(unique(marker) %in% c(0,1)))

  # loop over columns to evaluate tests
  for (x in setdiff(colnames(res), 'marker')){

    info <- strsplit(x, split='_')[[1]]
    if (length(setdiff('rep', info)) > 0){
      stop("Column name in test results need the info 'rep'!")
    }
    r <- info[which(info == 'rep')+1]

    # named vector of p.vals for test rep x
    p.val <- res[,x]
    p.val <- p.adjust(p.val, method=adjust)
    # auc
    auroc <- roc(predictor = -log10(p.val + 1e-50), response = marker,
                 levels = c(0, 1), direction = '<')
    # p value of alpha as cutoff
    # TP
    TP <- sum(marker[p.val < alpha] == 1)
    # TN
    TN <- sum(marker[p.val > alpha] == 0)
    # Type1
    FP <- sum(marker[p.val < alpha] == 0)
    # Type2
    FN <- sum(marker[p.val > alpha] == 1)
    # FDR
    FDR <- FP/(FP+TP)
    FDR <- ifelse(is.na(FDR), 0, FDR)
    # precision
    PR <- TP/(FP+TP)
    PR <- ifelse(is.na(PR), 0, PR)
    # recall
    R <- TP/(TP+FN)
    R <- ifelse(is.na(R), 0, R)
    df.res <- rbind(df.res, data.frame(rep=r, auroc=as.double(auroc$auc),
                                       TP=TP, FP=FP,
                                       TN=TN, FN=FN, PR=PR,
                                       R=R, FDR=FDR))
  }
  return(df.res)
}


#' @keywords internal
check.eval.parameters <- function(sim.location, group, 
                                  res.mat, adjust, alpha) {

  # check simulation file
  if (!file.exists(sim.location)) {
    stop("No such simulation exists!") }

  if (mode(adjust)!='character' | length(adjust)!=1){
    stop("Parameter 'adjust' should be a character and of length 1!")
  }
  if (!adjust %in% stats::p.adjust.methods){
    stop("Parameter 'adjust' must be one of these: c('", 
         paste(stats::p.adjust.methods, collapse = "', '"), "')")
  }
  # get marker features
  this <- h5read(file = sim.location, name = group)
  marker.feat <- this$marker_idx

  # add marker column
  marker <- ifelse((rownames(res.mat) %in% marker.feat), 1, 0)
  res <- cbind(res.mat, marker)

  # check results
  if (mode(res) == 'numeric'){
    if (!is.matrix(res)){
      stop("Parameter 'test.results' should be a matrix!")
    }
    if (is.null(colnames(res))){
      stop("Test result matrix needs column names!")
    }
    if (!('marker' %in% colnames(res))){
      stop("Column 'marker' is missing in the test result matrix!")
    }
    if (any(is.na(res))) {
      stop("NAs present in data! Check input of 'test.results'!")
    }
  } else {
    stop("Cannot interpret the input of 'test.results'!")
  }

  # check that alpha makes sense
  if (mode(alpha) != 'numeric'){
    stop("Parameter 'alpha' has be numeric!")
  }
  if (length(alpha) > 1){
    stop("Parameter 'alpha' has to be of length 1!")
  }
  if (alpha > 0.4){
    stop("Parameter 'alpha' is unreasonably large (alpha = ", alpha, ")!")
  } else if (alpha < 0.0001){
    stop("Parameter 'alpha' is unreasonably small (alpha = ", alpha, ")!")
  }
  return(res)
}


eval.real.test <- function(sim.location, group, res.mat, alpha = 0.05) {

  res <- check.eval.parameters(sim.location, group, res.mat, alpha)

  # real biological names and 'marker' species
  feat.names <- h5read(sim.location, '/original_data/filt_feature_names_real')
  # browser()
  if (nrow(res) != length(feat.names))
    stop('+ check that filtered biological names are stored correctly')
  # idx <- c(grep('bartlettii', feat.names),
  #          grep('gnavus', feat.names),
  #          grep('Escherichia', feat.names),
  #          grep('Enterococc', feat.names),
  #          grep('Streptococc', feat.names),
  #          grep('plebeius', feat.names),
  #          grep('Rothia', feat.names))
  rownames(res) <- feat.names
  # res[idx,'marker'] <- 1
  # use to calculate statistics for marker 'signatures' ?
  # remove marker col -- invalid for real data
  # res <- res[,-ncol(res)]

  df.res <- data.frame()
  # df.tax <- data.frame()
  bin.mat <- matrix(nrow = nrow(res), ncol = ncol(res)-1)
  colnames(bin.mat) <- setdiff(colnames(res), 'marker')
  rownames(bin.mat) <- feat.names

  # loop over columns (reps)
  for (x in setdiff(colnames(res), 'marker')) {

    info <- strsplit(x, split='_')[[1]]
    if (length(setdiff('rep', info)) > 0){
      stop("Column name in test results need the info 'rep'!")
    }
    r <- info[which(info == 'rep')+1]

    # named vector of p.vals for test rep x
    p.val <- res[,x]
    # df.tax <- rbind(df.tax, t(p.val[idx]))
    n.sig <- sum(p.val < alpha)
    pct.sig <- n.sig/length(p.val)*100
    df.res <- rbind(df.res, data.frame(rep = r, pct.sig = pct.sig))
    p.bin <- ifelse(p.val <= alpha, 1, 0)
    bin.mat[,x] <- p.bin
  }
  # browser()
  # check for significant features -- return subset
  # idx <- rowSums(bin.mat)[which(rowSums(bin.mat) >= ncol(bin.mat)/2)]
  # return(cbind(df.res, t(res[idx, -ncol(res)])))
  # return(cbind(df.res, df.tax))
  return(cbind(df.res, t(res[, -ncol(res)])))
}
