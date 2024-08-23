#!/usr/bin/Rscript
### SIMBA
### Simulation of Microbiome data with Biological Accuracy
### Morgan Essex - MDC Berlin
### Jakob Wirbel - EMBL Heidelberg
### 2020 GNU GPL 3.0

#' @title Check simulations for realism
#'
#' @description This function will compare simulated data with the
#' real data on which the simulations are based
#'
#' @usage reality.check(sim.location, group, ml=FALSE)
#'
#' @param sim.location file name for the .h5 file containing the simulations
#'
#' @param group character, name of the group within the .h5 file containing
#' the simulated features \strong{OR} \code{"original"} to get the same
#' measures for the original data set
#'
#' @param ml Boolean, indicating if a machine-learning model should be trained
#' to distinguish real and simulated data (see Details below).
#' Defaults to \code{FALSE}
#'
#' @export
#'
#' @keywords SIMBA reality.check
#'
#' @details This function computes a list containing different measures for how
#' realistic the simulations are, containing
#' \itemize{
#' \item \code{sparsity} measures for the sparsity of samples, computed as
#' L0-norm and as Gini-coefficient
#' \item \code{variance} measures for the variance of features, containing
#' additionally measures for differences between groups, such as generalized
#' fold change and prevalence difference
#' \item \code{AUC-BC} AUROC value for how well real and simulated data can be
#' distinguished based on the Bray-Curtis distance
#' \item \code{AUC-E} AUROC value for how well real and simulated data can be
#' distinguished based on the log-Euclidean distance
#' \item \code{AUC-ML} AUROC value for how well real and simulated data can be
#' distinguished based on a \link[SIAMCAT]{SIAMCAT}
#' machine learning model
#' }
#'
#' @return Returns a list of measures
#'
reality.check <- function(sim.location, group, ml=FALSE){
  message("+ Checking parameters")
  log.n0 <- reality.check.parameters(sim.location, group, ml)

  # finally run simulation checks :D
  # get the original dataset
  feat.filt <- h5read(file=sim.location,
                      name='/original_data/filt_features')
  colnames(feat.filt) <- h5read(file=sim.location,
                                name='/original_data/sample_names')
  colnames(feat.filt) <- paste0('real_', colnames(feat.filt))
  rownames(feat.filt) <- h5read(file=sim.location,
                                name='/original_data/filt_feature_names')
  feat.rel <- prop.table(feat.filt, 2)

  if (group=='original'){
    message("+ Calculating sparsity")
    df.sparsity <- data.frame(Gini=vapply(colnames(feat.filt),
                                          FUN=function(x){
                                            Gini(feat.filt[,x])
                                            }, FUN.VALUE = double(1)),
                              L0=vapply(colnames(feat.filt),
                                        FUN=function(x){sum(feat.filt[,x]==0)},
                                        FUN.VALUE = integer(1)),
                              group=group)
    rownames(df.sparsity) <- NULL
    message("+ Calculating feature variance")
    df.select <- data.frame(
      idx=rownames(feat.filt),
      median=rowMedians(feat.rel[rownames(feat.filt),]),
      quant75=rowQuantiles(feat.rel[rownames(feat.filt),],
                           probs=0.75))
    df.select$type <- ifelse(df.select$median==0 & df.select$quant75==0, 'low',
                             ifelse(df.select$median==0 & df.select$quant75!=0,
                                    'middle','high'))
    df.variance <- data.frame(var.rel=rowVars(feat.rel),
                              mean.rel=rowMeans(feat.rel),
                              features=rownames(feat.filt),
                              selection=FALSE,
                              feat.type=df.select[rownames(feat.filt), 'type'])
    message("+ Calculating background for correlation structure")
    cor.real <- cor(log10(t(prop.table(feat.filt, 2)) + log.n0))
    diag(cor.real) <- NA
    cor.real[lower.tri(cor.real)] <- NA
    cor.real <- as.vector(cor.real)
    cor.real <- cor.real[!is.na(cor.real)]

    pb <- progress_bar$new(total=50)
    cor.diffs <- c()
    for (i in seq_len(50)){
      cor.resample <- cor(log10(t(prop.table(
        feat.filt[,sample(colnames(feat.filt), ncol(feat.filt)/2)], 2)) +
          log.n0))
      diag(cor.resample) <- NA
      cor.resample[lower.tri(cor.resample)] <- NA
      cor.resample <- as.vector(cor.resample)
      cor.resample <- cor.resample[!is.na(cor.resample)]
      cor.diffs <- c(cor.diffs, mean(abs(cor.real-cor.resample)))
      pb$tick()
    }

    l.return <- list("sparsity"=df.sparsity,
                     "variance"=df.variance,
                     "Bray_Curtis"=NULL,
                     "Euclidean"=NULL,
                     "Correlation"=cor.diffs,
                     "ML"=NULL)
  } else {
    # simulated data
    feat.sim <- h5read(file=sim.location,
                       name=paste0('/',group, '/features'))
    colnames(feat.sim) <- h5read(file=sim.location,
                                 name=paste0('/',group, '/sample_names'))
    rownames(feat.sim) <- h5read(file=sim.location,
                                 name=paste0('/',group, '/feature_names'))
    markers <- h5read(file=sim.location,
                      name=paste0('/',group, '/marker_idx'))
    sim.label <- h5read(file=sim.location,
                        name=paste0('/', group, '/labels'))
    names(sim.label) <- colnames(feat.sim)
    if (any(colSums(feat.sim) == 0)){
        idx <- which(colSums(feat.sim)==0)
        feat.sim <- feat.sim[,-idx]
        sim.label <- sim.label[-idx]
    }
    feat.sim.rel <- prop.table(feat.sim, 2)

    # calculate different measures
    message("+ Calculating sparsity")
    # sparsity
    df.sparsity <- data.frame(
      Gini=vapply(colnames(feat.sim),FUN=function(x){
        Gini(feat.sim[,x])},FUN.VALUE = double(1)),
      L0=vapply(colnames(feat.sim),FUN=function(x){
        sum(feat.sim[,x]==0)},FUN.VALUE = integer(1)),
      group=group)
    rownames(df.sparsity) <- NULL

    # which types of features were selected for implantation?
    # (potential effects for compositionality)
    message("+ Calculating feature variance, gFC, and prevalence shift")
    df.select <- data.frame(
      idx=rownames(feat.filt),
      median=rowMedians(feat.rel[rownames(feat.filt),]),
      quant75=rowQuantiles(feat.rel[rownames(feat.filt),],
                           probs=0.75))
    df.select$type <- ifelse(df.select$median==0 & df.select$quant75==0, 'low',
                             ifelse(df.select$median==0 & df.select$quant75!=0,
                                    'middle','high'))

    # feature variance
    x <- rownames(feat.filt)
    y <- rownames(feat.sim)
    # problem for negbin/betabin simulations (some models might fail for 
    # some features)
    stopifnot(all(y %in% x))
    overlap <- intersect(x, y)
    feat.filt <- feat.filt[overlap,]
    feat.sim <- feat.sim[overlap,]

    # add also gFC and prevalence shift
    df.variance <- data.frame(
      var.rel=rowVars(feat.sim.rel),
      mean.rel=rowMeans(feat.sim.rel),
      features=rownames(feat.filt),
      selection=rownames(feat.sim) %in% markers,
      feat.type=df.select[rownames(feat.filt), 'type'])

    gfc.simulated <- vapply(rownames(feat.sim.rel), FUN = function(x){
      temp <- log10(feat.sim.rel[x,] + log.n0)
      q.ctr <- quantile(temp[which(sim.label==-1)], probs=seq(.05, .95, .05))
      q.case <- quantile(temp[which(sim.label==1)], probs=seq(.05, .95, .05))
      sum(q.case - q.ctr)/length(q.case)
    }, FUN.VALUE = double(1))

    prev.shift.simulated <- vapply(rownames(feat.sim.rel), FUN = function(x){
      temp <- feat.sim.rel[x,]==0
      p.ctr <- sum(temp[which(sim.label==-1)])/length(which(sim.label==-1))
      p.case <- sum(temp[which(sim.label==1)])/length(which(sim.label==-1))
      p.case-p.ctr
    }, FUN.VALUE = double(1))
    df.variance$fc.sim <- gfc.simulated
    df.variance$prev.shift.sim <- prev.shift.simulated

    # BC distances
    message("+ Calculating separation between real and simulated data")
    message("++ based on Bray-Curtis distance")
    tmp <- cbind(feat.filt, feat.sim)
    bc <- vegdist(t(tmp))
    permanova.bc <- vegan::adonis2(bc~c(colnames(tmp) %in% colnames(feat.filt)))
    bc <- as.matrix(bc)
    diag(bc) <- NA
    bc[lower.tri(bc)] <- NA
    # can the two groups be distinguished by the PCoA alone?
    within.class <- c(c(bc[colnames(feat.filt), colnames(feat.filt)]),
                      c(bc[colnames(feat.sim), colnames(feat.sim)]))
    between.class <- c(c(bc[colnames(feat.filt), colnames(feat.sim)]),
                       c(bc[colnames(feat.sim), colnames(feat.filt)]))
    res.bc <- roc(cases=within.class, controls=between.class, direction='>')

    # log-euclidean distances
    message("++ based on log-Euclidean distance")
    euc <- vegdist(t(log10(tmp + log.n0)), method='euclidean')
    permanova.euc <- vegan::adonis2(euc~c(colnames(tmp) %in% colnames(feat.filt)))
    euc <- as.matrix(euc)
    diag(euc) <- NA
    euc[lower.tri(euc)] <- NA
    # can the two groups be distinguished by the PCoA alone?
    within.class <- c(c(euc[colnames(feat.filt), colnames(feat.filt)]),
                      c(euc[colnames(feat.sim), colnames(feat.sim)]))
    between.class <- c(c(euc[colnames(feat.filt), colnames(feat.sim)]),
                       c(euc[colnames(feat.sim), colnames(feat.filt)]))
    res.euc <- roc(cases=within.class, controls=between.class, direction='>')


    message("++ compute the correlation-structure across taxa")
    cor.real <- cor(log10(t(prop.table(feat.filt, 2)) + log.n0))
    diag(cor.real) <- NA
    cor.real[lower.tri(cor.real)] <- NA
    cor.real <- as.vector(cor.real)
    cor.real <- cor.real[!is.na(cor.real)]

    cor.sim <- cor(log10(t(prop.table(feat.sim, 2)) + log.n0))
    diag(cor.sim) <- NA
    cor.sim[lower.tri(cor.sim)] <- NA
    cor.sim <- as.vector(cor.sim)
    cor.sim <- cor.sim[!is.na(cor.sim)]

    cor.diff <- mean(abs(cor.real - cor.sim))

    # ML
    if (ml){
      message("+ Using machine-learning to assess separation between real",
              " and simulated data")
      test.package('SIAMCAT')
      label.vec <- c(rep('original', ncol(feat.filt)),
                    rep('sim', ncol(feat.sim)))
      names(label.vec) <- c(paste0('real_', seq_len(ncol(feat.filt))),
                            paste0('sim_', seq_len(ncol(feat.sim))))
      tmp <- cbind(feat.filt, feat.sim)
      colnames(tmp) <- names(label.vec)
      # subsample to 100 at most in each group
      if (ncol(feat.filt) > 100){
        real <- paste0('real_', seq_len(ncol(feat.filt)))
        sim <- paste0('sim_', seq_len(ncol(feat.sim)))
        label.vec <- label.vec[c(sample(real, size=100), sample(sim, size=100))]
      }

      sc <- SIAMCAT::siamcat(feat=prop.table(tmp, 2),
                    label=label.vec,
                    case='sim',
                    verbose=0)
      sc <- SIAMCAT::normalize.features(sc,
                               norm.method = 'log.std',
                               norm.param = list(log.n0=log.n0, sd.min.q=0),
                               feature.type = 'original',
                               verbose=0)
      sc <- SIAMCAT::create.data.split(sc, 
                                       num.folds = 10, num.resample = 1,
                                       verbose=0)
      sc <- SIAMCAT::train.model(sc, verbose=0, method='lasso')
      sc <- SIAMCAT::make.predictions(sc, verbose=0)
      sc <- SIAMCAT::evaluate.predictions(sc, verbose=0)
      ml.res <- as.numeric(sc@eval_data$auroc)
    } else {
      ml.res <- NULL
    }
    l.return <- list("sparsity"=df.sparsity,
                     "variance"=df.variance,
                     "Bray_Curtis"=list('AUC'=res.bc$auc,
                                        'F'=permanova.bc$aov.tab$F.Model[1],
                                        'R2'=permanova.bc$aov.tab$R2[1],
                                        'P'=permanova.bc$aov.tab$`Pr(>F)`[1]),
                     "Euclidean"=list('AUC'=res.euc$auc,
                                      'F'=permanova.euc$aov.tab$F.Model[1],
                                      'R2'=permanova.euc$aov.tab$R2[1],
                                      'P'=permanova.euc$aov.tab$`Pr(>F)`[1]),
                     "Correlation"=cor.diff,
                     "ML"=ml.res)
  }
  return(l.return)
}


#' @keywords internal
reality.check.parameters <- function(sim.location, group, ml){
  # check the H5 file
  if (!file.exists(sim.location)){
    stop("Simulation file does not exist!")
  }
  # get and check simulation parameters of the h5file
  params <- h5read(file=sim.location,
                   name='simulation_parameters')
  if (is.null(params)){
    stop("Parameters are empty for this file!")
  }
  if (params$general_params$sim.type == 'time-course'){
    stop('Not implemented yet!')
  }
  if (params$general_params$sim.method == 'pass'){
    stop("Reality checks are not applicable for the 'pass' method!")
  }
  log.n0 <- as.vector(params$filt_params$log.n0)
  params <- params$sim_params
  needed <- setdiff(c('ab.scale', 'repeats'), names(params))
  if (length(needed) > 0){
    stop("Simulation layout is not suitable for checks!\n",
         "Missing the parameters:", paste(needed, collapse = ', '))
  }
  # ml
  if (!is.logical(ml)){
    stop("'ml' parameters must be logical!")
  }

  # info
  if (group=='original'){
    if (ml){
      message("Parameter 'ml' will be ignored for the original data measures!")
    }

  } else {
    all.groups <- h5ls(file=sim.location, recursive=FALSE)
    if (!group %in% all.groups$name){
      stop("Group '", group, "' not found in the h5 file!")
    }
  }
  return(log.n0)
}
