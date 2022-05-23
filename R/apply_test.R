#!/usr/bin/Rscript
### metaGsim package
### Morgan Essex - MDC Berlin
### Jakob Wirbel - EMBL Heidelberg
### 2020 GNU GPL 3.0

#' @title Apply different tests to the simulated data
#'
#' @description This function applies one of the implemented differential
#' abundance testing methods to a dataset within the simulation
#'
#' @usage apply.test(sim.location, group, type='default',
#' subset=NA_real_, norm='pass',
#' test='wilcoxon', conf=NULL, ...)
#'
#' @param sim.location file name for the .h5 file containing the simulations
#'
#' @param group character, name of the group within the .h5 file containing
#' the simulated features
#'
#' @param type character, name of the index group within the .h5 file
#' containing the testing or bootstrapping indices
#' Can be either \code{"default"}, \code{"biomarker"},  \code{"confounder"},
#' or \code{"random"}
#'
#' @param subset vector, subset size(s) to be tested
#'
#' @param norm character, normalization method to be applied before testing
#'
#' @param test character, differential abundance testing method
#'
#' @param conf character, name of the confounder metadata column
#'
#' @param ... additional parameters (for experts)
#'
#' @export
#'
#' @keywords SIMBA apply.test
#'
#' @return This function returns a list of P-values matrices
#'
#' @details todo
apply.test <- function(sim.location, group, type='default',
                       subset=NA_real_, norm='pass',
                       test='wilcoxon', conf=NULL, ...){

  additional.arguments <- list(...)
  # check testing parameters
  message("+ Checking testing parameters!")
  all <- check.testing.parameters(sim.location, group, type, subset, norm,
                                  test, conf, additional.arguments)

  res.all <- list()
  meta.all <- all[["meta"]]

  for (g in names(all)){

    # reconstruct abundance table
    feat.sim <- all[[g]]$features
    markers <- all[[g]]$markers
    test.idx <- all[[g]]$idx
    log.n0 <- all[[g]]$log.n0
    conf.l <- all[[g]]$conf

    if (!is.null(meta.all) & !is.null(conf.l)){
      conf.mat <- cbind(meta.all, conf.l)
      colnames(conf.mat)[ncol(conf.mat)] <- 'conf'
      conf.mat <- conf.mat[,conf, drop=FALSE]
      rownames(conf.mat) <- colnames(feat.sim)
    } else if (!is.null(meta.all) & is.null(conf.l)){
      conf.mat <- meta.all
      conf.mat <- conf.mat[,conf, drop=FALSE]
      rownames(conf.mat) <- colnames(feat.sim)
    } else if (is.null(meta.all) & !is.null(conf.l)){
      conf.mat <- data.frame(conf.l)
      colnames(conf.mat) <- 'conf'
      rownames(conf.mat) <- colnames(feat.sim)
    } else if (is.null(meta.all) & is.null(conf.l)){
      conf.mat <- NULL
    }

    res.list <- list()

    for (s in names(test.idx)){
      info.mat <- test.idx[[s]]
      label <- info.mat['label',]
      idx.mat <- info.mat[which(grepl('^rep[0-9]*$', rownames(info.mat))),,
                          drop=FALSE]
      p.val.mat <- matrix(1, ncol=nrow(idx.mat), nrow=nrow(feat.sim),
                          dimnames = list(rownames(feat.sim),
                                          paste0('rep_',
                                                 seq_len(nrow(idx.mat)))))
      pb <- progress::progress_bar$new(total = nrow(idx.mat))
      for (r in seq_len(nrow(idx.mat))){

        df.test <- feat.sim[,idx.mat[r,]]

        # normalize
        df.test <- norm.data(df.test, norm, log.n0)
        # browser()
        names(label) <- colnames(df.test)
        # run test
        p.val <- run.test(data=df.test,
                          label=label,
                          test=test,
                          conf=conf.mat)
        # conf=conf.mat[idx.mat[r,],,drop=FALSE])
        stopifnot(!is.null(names(p.val)))
        if (any(is.na(p.val))){
          p.val[is.na(p.val)] <- 1
        }
        p.val.mat[names(p.val),r] <- p.val
        pb$tick()
      }
      res.list[[s]] <- p.val.mat
      message("++ Finished calculations for ", s)
    }

    res.all[[g]] <- res.list
    message('+ Finished calculations for group ', g)
  }
  return(res.all)
}

#' @keywords internal
check.testing.parameters <- function(sim.location, group, type,
                                     subset, norm, test, conf,
                                     additional.arguments = list()){

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
  log.n0 <- as.vector(params$filt_params$log.n0)
  sim.params <- params$sim_params

  # check norm and test combination
  implemented.test <- c('wilcoxon',
                        'DESeq2',
                        'metagenomeSeq',
                        'metagenomeSeq2',
                        'mdc-FE', # metadeconfoundR: fixed effect
                        'mdc-RE', # metadeconfoundR: random effect
                        'edgeR',
                        'ZIBSeq',
                        'ZIBSeq-sqrt',
                        'ANCOM',
                        'ANCOM_old',
                        'ANCOMBC',
                        "lm",
                        "corncob",
                        "limma",
                        "lme",
                        'fisher',
                        'ALDEx2',
                        'KS',
                        'scde',
                        'MAST',
                        'distinct',
                        'mixMC',
                        'ZINQ',
                        'songbird',
                        'gFC')
  if (!is.character(test)){
    stop("Parameter 'test' should be a character!")
  }
  if (length(test) != 1){
    stop("Parameter 'test' should be of length 1!")
  }
  if (!test %in% implemented.test){
    stop("Test '", test, "' not implemented yet!")
  }
  # TODO test the norm parameter
  implemented.norms <- c('pass', 'TSS', 'TSS.log', 'clr', 'rclr',
                         'TSS.arcsin','rarefy','rarefy.TSS',
                         'rarefy.TSS.log',
                         'rank')
  if (!is.character(norm)){
    stop("Parameter 'norm' should be a character!")
  }
  if (length(norm) != 1){
    stop("Parameter 'norm' should be of length 1!")
  }
  if (!norm %in% implemented.norms){
    stop("Normalization method '", norm, "' not implemented yet!")
  }
  if (test == 'gFC'){
    if (norm != 'TSS.log'){
      stop("Parameter 'norm' should be 'TSS.log' for gFC reference!")
    }
  }
  # is the type allowed
  allowed.types <- c('default', 'biomarker', 'confounder', 'random')
  names(allowed.types) <- c('test_idx', 'bio_bias_test_idx',
                            'conf_bias_test_idx', 'random_test_idx')
  if (type %in% allowed.types){
    name <- names(allowed.types)[which(allowed.types==type)]
    if (type=='biomarker'){
      bias <- seq_along(params[[name]]$bias)
    } else if (type=='confounder'){
      bias <- seq_along(params[[name]]$bias)
    } else {
      bias <- NA_real_
    }
  } else if (grepl(pattern = '::[0-9]*$', type)){
    type.split <- strsplit(type, split = '::')[[1]]
    type <- type.split[1]
    name <- names(allowed.types)[which(allowed.types==type)]
    bias <- as.numeric(type.split[2])
    if (!type %in% c('biomarker', 'confounder', 'random')){
      stop("Parameter 'type' should be either biomarker or confounder, not ",
           type)
    }
    # check bias
    if (is.na(bias)){
      bias.ref <- params[[name]]$bias
      if (bias > length(bias.ref)){
        stop("Bias parameter not present in indices!")
      }
    }
  } else {
    stop("Parameter 'type' must be one of these:\n\t",
         paste(allowed.types, collapse = ', '))
  }

  # info
  all.groups <- h5ls(file=sim.location, recursive = FALSE)
  all.groups <- setdiff(all.groups$name,
                        c('original_data','simulation_parameters'))
  all.groups.red <- unique(gsub(all.groups, replacement = '',
                                pattern = '_rep[0-9]*$'))
  # distinguish real data ? no ab_prev pattern -- all groups unique
  if (all(all.groups == all.groups.red))
    group <- group
  else if (group %in% all.groups.red) {
    group <- all.groups[grep(group, all.groups)] }

  # confounder info
  if (!is.null(conf)){
    if (typeof(conf)!='character'){
      stop("Parameter 'confounder' has to be a character!")
    }
    if (length(conf)==1 & conf=='conf'){
      TRUE
    } else {
      meta <- h5read(sim.location, 'original_data/metadata')
      if (!all(conf %in% c('conf', colnames(meta)))){
        stop("All confounders have to be present in the metadata!")
      }
    }
  }

  # check subsets
  sizes <- params[[name]]$subsets
  if (!is.numeric(subset)){
    stop("Subsets need to be numeric")
  }
  if (length(subset) == 1){
    if (is.na(subset)){
      subset <- sizes
    }
  }
  ov <- intersect(subset, sizes)
  if (length(ov) == 0){
    stop("Supplied subset sizes do not match available subset sizes!")
  }
  if (length(ov) < length(subset)){
    warning("Not all supplied subset sizes are found in the data!")
    subset <- subset[subset %in% sizes]
  }

  # loop through groups
  all <- list()

  for (g in group){
    all[[g]] <- list()
    all[[g]]$log.n0 <- log.n0

    if (!g %in% all.groups){
      stop("Group '", g, "' not found in the h5 file!")
    }
    complete.group <- h5read(file=sim.location, name=g)

    if (!name %in% names(complete.group)){
      stop("Indices for ", type, " testing have not been created yet!")
    }

    # get indices and store them

    idx.all <- complete.group[[name]]
    if (!all(is.na(bias))){
      idx.all <- list()
      idx.all.temp <- complete.group[[name]]
      for (a in names(idx.all.temp)){
        for (b in names(idx.all.temp[[a]])){
          idx.all[[paste0(a, '-', b)]] <- idx.all.temp[[a]][[b]]
        }
      }
    }
    # rownames (important :D)
    for (a in seq_along(idx.all)){
      rownames(idx.all[[a]]) <- params[[name]][['names']]
    }

    # restrict to correct subsets
    idx.all <- idx.all[grep(paste(paste0('subset_', subset, '$'),
                                  collapse = '|'),
                            names(idx.all), value=TRUE)]

    # restrict to correct biases
    if (!all(is.na(bias))){
      if (length(bias)==1){
        idx.all <- idx.all[grep(paste0('bias_', bias),
                                names(idx.all), value=TRUE)]
      }
    }

    # check if rep_no in the additional arguments
    if (length(additional.arguments > 0)) {
      if ('rep_no' %in% names(additional.arguments)){
        if (length(subset) > 1){
          warning("Parameter 'rep_no' will be ignored, since more",
                  " than one subset is present")
          additional.arguments$rep_no <- NULL
        } else {
          for (a in names(idx.all)){
            idx.all[[a]] <- idx.all[[a]][[
              paste0('rep', additional.arguments$rep.no), , drop=FALSE]]
          }
        }
      }
    }
    all[[g]][['idx']] <- idx.all


    if (!is.null(conf) & any(conf=='conf')){
      all[[g]][['conf']] <- complete.group$conf_label
    }
    # pass the features/markers/labels nicely
    feat.sim <- complete.group$features
    colnames(feat.sim) <- complete.group$sample_names
    rownames(feat.sim) <- complete.group$feature_names
    all[[g]]$features <- feat.sim

    all[[g]]$markers <- complete.group$marker_idx
    if (!is.null(conf) & any(conf!='conf')){
      all[[meta]] <- h5read(sim.location, 'original_data/metadata')
    }
  }
  return(all)
}



#' @keywords internal
norm.data <- function(df, norm, log.n0 = NULL){
  if (norm=='pass'){
    return(df)
  } else if (norm=='TSS'){
    return(prop.table(df, 2))
  } else if (norm=='TSS.log'){
    if (is.null(log.n0)){
      stop("Parameter 'log.n0' needed for TSS.log normalization!") }
    x <- prop.table(df, 2)
    x <- x + log.n0
    x <- log10(x)
    return(x)
  } else if (norm=='clr'){
    # add pseudocount
    # compute geometric mean per sample
    # divide by geometric mean
    df.clr <- df + 1
    df.norm <- vapply(colnames(df.clr), FUN = function(x){
      temp <- df.clr[,x]
      gm <- exp(mean(log(temp)))
      return(log(temp/gm))
    }, FUN.VALUE = double(nrow(df.clr)))
    return(df.norm)
  } else if (norm=='rclr'){
    # compute geometric mean per sample of non-zero entries
    # add pseudocount
    # divide by geometric mean
    df.clr <- df
    df.norm.r <- vapply(colnames(df.clr), FUN = function(x){
      temp <- df.clr[,x]
      gm <- exp(mean(log(temp[temp!=0])))
      return(log((temp + 1)/gm))
    }, FUN.VALUE = double(nrow(df.clr)))
    return(df.norm.r)
  } else if (norm=='TSS.arcsin'){
    x <- prop.table(df, 2)
    y <- asin(sqrt(x))
    return(y)
  } else if (norm=='rarefy') {
    df.int <- vapply(colnames(df),
                     FUN=function(x){round(df[,x])},
                     FUN.VALUE = double(nrow(df)))
    rarefied <- t(vegan::rrarefy(t(df.int), quantile(colSums(df), 0.25)))
    return(rarefied)
  } else if (norm=='rarefy.TSS') {
    df.int <- vapply(colnames(df),
                     FUN=function(x){round(df[,x])},
                     FUN.VALUE = double(nrow(df)))
    rarefied <- t(vegan::rrarefy(t(df.int), quantile(colSums(df), 0.25)))
    return(prop.table(rarefied, 2))
  } else if (norm=='rarefy.TSS.log') {
    df.int <- vapply(colnames(df),
                     FUN=function(x){round(df[,x])},
                     FUN.VALUE = double(nrow(df)))
    rarefied <- t(vegan::rrarefy(t(df.int), quantile(colSums(df), 0.25)))
    rarefied.tss <- prop.table(rarefied, 2)
    return(log10(rarefied.tss + log.n0))
  } else if (norm=='rank') {
    ranked <- vapply(colnames(df), FUN = function(x) {
      rank(df[,x], ties.method='min') },
      FUN.VALUE = double(nrow(df)))
    return(ranked)
  } else {
    stop("Unknown normalization method!: ", norm)
  }
}
