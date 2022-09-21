#!/usr/bin/Rscript
### SIMBA
### Simulation of Microbiome data with Biological Accuracy
### Morgan Essex - MDC Berlin
### Jakob Wirbel - EMBL Heidelberg
### 2020 GNU GPL 3.0

#' @title Create biased indices for testing
#'
#' @description This function creates biased indices for testing with an
#' additional confounder
#'
#' @usage create.biased.idx(sim.location, subsets, repetitions,
#' lab.svar = NULL, meta.svar = NULL, strat.scheme = 'all', bias = 0.5)
#'
#' @param sim.location file name for the .h5 file containing the simulations
#'
#' @param subsets vector of intended subset sizes for resampled data
#'
#' @param repetitions number of unique sets of indices to generate for each
#' subset size deemed valid
#'
#' @param lab.svar for real data, name of the main stratification "label"
#' variable in the original metadata; defaults to \code{NULL} and is not
#' applicable for simulated data
#'
#' @param meta.svar for real data, a single string or vector of strings
#' corresponding to the name(s) of binary metadata variable(s) by which biased
#' resampling should be stratified; defaults to \code{NULL} and is not
#' applicable for simulated data
#'
#' @param strat.scheme character, name of the index group to be written within
#' the .h5 file containing different indices for downstream testing
#' Can be either \code{'biomarker'},\code{'confounder'}, \code{'random'},
#' or \code{'all'} (the default)
#'
#' @param bias float, strength of the confounder bias in the subsampling
#'
#' @keywords internal
#'
#' @return vector of real metadata variables insufficient for resampling, if
#' applicable, otherwise nothing returned & indices are written to
#' \code{'sim.location'}
#'
#' @details todo
create.biased.idx <- function(sim.location, subsets, repetitions,
                              # should both be NULL if simulated data
                              lab.svar = NULL, meta.svar = NULL,
                              strat.scheme = 'all', bias = 0.5) {
  # browser()
  if (strat.scheme == 'all') {
    clear.biased.idx(sim.location, 'bio_bias')
    clear.biased.idx(sim.location, 'conf_bias')
    clear.biased.idx(sim.location, 'random')
  } else {
    id <- c('confounder'='conf_bias', 'biomarker'='bio_bias',
            'random'='random')
    clear.biased.idx(sim.location, id[strat.scheme])
  }

  params <- h5read(file = sim.location, name = 'simulation_parameters')

  # real data biased resampling, no implanted confounder
  if (params$general_params$sim.method == 'pass') {
    if (is.null(meta.svar) | length(meta.svar) == 0)
      stop('Variable(s) for stratification needed! None provided') }

  skipped <- perform.idx.bias.checks(sim.location, subsets, repetitions,
                                     lab.svar, meta.svar,
                                     strat.scheme = strat.scheme,
                                     bias = bias)

  if (!is.null(skipped)) {
    message('+ Stratification variables are ok -- ',
            paste(skipped, collapse = ', '), 'will be removed') }

  sim.params <- params$sim_params
  groups <- h5ls(sim.location, recursive = FALSE)
  groups <- groups$name
  groups <- setdiff(groups, c('original_data', 'simulation_parameters'))

  message(paste0('+ Creating ', strat.scheme, ' biased test indices'))
  pb <- progress_bar$new(total=length(groups))
  for (group in groups) {
    switch(strat.scheme,
           'biomarker' = {
             temp.bio <- save.biomarker.idx(subsets, repetitions,
                                            sim.location,
                                            group, bias)},
           'confounder' = {
             temp.conf <- save.confounder.idx(subsets, repetitions,
                                              sim.location,
                                              group, bias)},
           'random' = {
             temp.random <- save.random.idx(subsets, repetitions,
                                            sim.location, group)},
           'all' = {
             temp.bio <- save.biomarker.idx(subsets, repetitions,
                                            sim.location,
                                            group, bias)

             temp.conf <- save.confounder.idx(subsets, repetitions,
                                              sim.location,
                                              group, bias)

             temp.random <- save.random.idx(subsets, repetitions,
                                            sim.location, group)
           })
    pb$tick()
  }

  # save rownames and bias info
  if (exists('temp.bio')) {
    h5write(temp.bio, file = sim.location,
            name = '/simulation_parameters/bio_bias_test_idx') }
  if (exists('temp.conf')) {
    h5write(temp.conf, file = sim.location,
            name = '/simulation_parameters/conf_bias_test_idx') }
  if (exists('temp.random')) {
    h5write(temp.random, file = sim.location,
            name = '/simulation_parameters/random_test_idx') }

  message('++ Finished saving all biased test indices to simulation file')
}

#' @keywords internal
#'
#' vector of real metadata variables insufficient for resampling, if
#' applicable, otherwise nothing returned & indices are written to
#' \code{'sim.location'}
perform.idx.bias.checks <- function(sim.location, subsets, repetitions,
                                    lab.svar, meta.svar,
                                    strat.scheme, bias) {

  # browser()
  message('+ Starting to check parameters for bootstrap idx creation')

  ## general parameter checks
  # check the H5 file
  if (!file.exists(sim.location)) {
    stop('Simulation file does not exist!') }
  params <- h5read(file=sim.location,
                   name='simulation_parameters')
  if (is.null(params)) { stop('Parameters are empty for this file!') }
  if (params$general_params$sim.type == 'time-course') {
    stop('Not implemented yet!') }

  # bias parameters
  bias <- sort(bias)
  if (!any(is.numeric(bias))) stop("The 'bias' parameter should be numeric!")

  if (any(bias > 1) | any(bias < 0.5)) {
    out.idx <- which((bias > 1 | bias < 0.5))
    message("+ The 'bias' parameter should range from 0.5 to 1! ",
            "Will remove values: ", paste(bias[out.idx], collapse = ', '))
    bias <- bias[-out.idx]
    if (length(bias) == 0) stop("Removed all values as too extreme! Exiting.")
  }

  # check repetitions
  if (!is.numeric(repetitions)) {
    stop('Parameter ', repetitions, ' must be numeric!')
  } else if (length(repetitions) > 1) {
    stop('Parameter ', repetitions, ' needs to be of length 1!') }

  # check if stratified bootstrap scheme is valid
  valid.schemes <- c('all','biomarker','confounder','random')
  if (!strat.scheme %in% valid.schemes)
    stop('Selected stratification scheme not available, has to be one of:\n\t',
         paste(valid.schemes, collapse = ', '))

  ## real data checks
  sim <- suppressWarnings(h5read(sim.location, '/original_data'))
  meta <- suppressWarnings(h5read(sim.location, '/original_data/metadata'))

  # get the number of eligible samples
  n.samples <- nrow(meta)

  # check subset size
  if (!is.numeric(subsets)){
    stop('Parameter ', subsets, ' must be numeric!') }
  subsets <- sort(unique(subsets))
  if (any(subsets < 0)){
    message('++ All entries in parameter ', subsets, ' need to be positive!')
    subsets <- subsets[subsets > 0]
    message('++ Parameter ', subset, ' reduced to: ',
            paste(subsets, collapse = ', ')) }
  if (any(subsets > n.samples)){
    message('++ Some entries in parameter subsets are too large!  c(',
            paste(subsets, collapse = ', '), ')')
    subsets <- subsets[subsets <= n.samples]
    message('++ Parameter subsets reduced to: ',
            paste(subsets, collapse = ', ')) }

  skipped <- NULL
  if (params$general_params$sim.method == 'pass') {
    # check if label already written to h5 or in metadata
    in.h5 <- 'labels' %in% names(sim)
    in.meta <- 'label' %in% tolower(colnames(meta))
    if (!any(in.h5, in.meta))
      stop('Label must be stored in /original_data/labels or in metadata table
           as Label/label/LABEL column!')
    if (in.meta & !in.h5) {
      label <- as.numeric(as.character(
        meta[,which(tolower(colnames(meta)) == 'label')]))
      low <- sort(unique(label))[1]
      hi <- sort(unique(label))[2]
      label[label == low] <- -1
      label[label == hi] <- 1
      h5write(label, sim.location, paste0('original_data', '/labels'))
      meta <- meta[, -which(tolower(colnames(meta)) == 'label')]
      message('++ Label written to /original_data/labels and removed from
            metadata') }
    if (all(in.meta, in.h5)) {
      meta <- meta[, -which(tolower(colnames(meta)) == 'label')]
      message('++ Label already written to .h5 -- removing from metadata') }

    # remove from meta.svar
    if ('label' %in% tolower(meta.svar))
      meta.svar <- meta.svar[-which(tolower(meta.svar) == 'label')]

    # check that meta.svar in metadata table
    excluded <- setdiff(meta.svar, colnames(meta))
    if (length(excluded) > 0) {
      message('++ The following metadata variables were already removed
            in preprocessing: ', paste(excluded, collapse = ', '))
      meta.svar <- meta.svar[-which(meta.svar %in% excluded)] }

    # remove extraneous variables from metadata table
    extra <- setdiff(colnames(meta), meta.svar)
    if (length(extra) > 0) {
      message('++ The following metadata variables were not specified as
            stratification variables and will be removed: ',
              paste(extra, collapse = ', '))
      meta <- meta[, -match(extra, colnames(meta))] }

    # individual variable checks & partitioning of h5 into 'groups'
    skipped <- c()
    labs <- h5read(sim.location, 'original_data/labels')
    for (v in meta.svar) {
      tryCatch(this.group <- h5read(sim.location, v),
               error = function(e) {
                 message(paste('Creating', v, 'group in .h5 file'))
                 h5createGroup(sim.location, v) })
      vec <- get(v, meta)
      n.levels <- setdiff(unique(vec), NA)
      # skip if not binary -- dummy vars later?
      if (length(n.levels) != 2) {
        skipped <- c(skipped, v)
        message(paste0('++ Metadata variables must contain exactly 2 groups! ',
                       v,' is not binary -- skipping for now.'))
        next() }

      # make variable specific metadata table for partitioning
      meta.tmp <- data.frame(label = labs, meta.label = vec)
      sample.names <- h5read(sim.location, '/original_data/sample_names')
      feat <- h5read(sim.location, '/original_data/filt_features')
      feat.names <- h5read(sim.location, '/original_data/filt_feature_names')

      # check NAs
      if (any(is.na(vec))) {
        prop.nas <- table(vec, useNA='ifany')/n.samples
        i <- which(is.na(names(prop.nas)))

        # keep variable if <10% NAs and n.samples large enough
        if (prop.nas[i] <= 0.1 & n.samples > 500) {
          message(paste0('++ Subsetting metadata and features
                         to remove NA samples for ', v))
          meta.tmp <- meta.tmp[-which(is.na(vec)),]

          # write subset features and samples to .h5
          h5write(obj = sample.names[-which(is.na(vec))],
                  file = sim.location,
                  name = paste0(v, '/sample_names'))

          feat <- feat[,-which(is.na(vec))]
          rowkeep <- which(rowSums(feat) >= as.numeric(
            params$filt_params$ab.cutoff))
          feat <- feat[rowkeep,]
          h5write(obj = feat,
                  file = sim.location,
                  name = paste0(v, '/features'))
          h5write(obj = feat.names[rowkeep],
                  file = sim.location,
                  name = paste0(v, '/feature_names'))
        } else {
          message(paste0('++ More than 10% of samples for ', v, ' are NAs,
                       sample size not large enough to subset -- removing'))
          skipped <- c(skipped, v)
          next() }
      } else {
        h5write(obj = sample.names,
                file = sim.location,
                name = paste0(v, '/sample_names'))
        h5write(obj = feat,
                file = sim.location,
                name = paste0(v, '/features'))
        h5write(obj = feat.names,
                file = sim.location,
                name = paste0(v, '/feature_names'))
      }

      # write main and meta labels
      h5write(obj = meta.tmp$label, file = sim.location,
              name = paste0(v, '/labels'))
      h5write(obj = meta.tmp$meta.label, file = sim.location,
              name = paste0(v, '/conf_label'))

      # update regular test_idx to reflect NA downsampling
      tryCatch(reg.idx <- h5read(sim.location, paste0(v,'/test_idx')),
               error = function(e) {
                 # message(paste('Saving test_idx for variable', v))
                 tmp.idx <- save.idx(NULL, NULL, subsets, repetitions, 0.5,
                                     sim.file = sim.location, group = v,
                                     replace = TRUE) })
    }
  }
  return(skipped)
}

#' @keywords internal
calculate.phi <- function(x, y) {

  tbl <- table(as.double(x), as.double(y))
  if (all(dim(tbl) != c(2,2)))
    return(NA)
  a <- tbl[1,1]
  b <- tbl[1,2]
  c <- tbl[2,1]
  d <- tbl[2,2]
  R1 <- as.double(a+b)
  R2 <- as.double(c+d)
  C1 <- as.double(a+c)
  C2 <- as.double(b+d)
  phi <- ((a*d)-(b*c))/sqrt(R1 * R2 * C1 * C2)
  return(phi)
}

#' @keywords internal
check.resampling.limits <- function(sim.file, group, bias) {

  # browser()
  sim <- h5read(file = sim.file, name=group)
  phi <- calculate.phi(as.double(sim$labels), as.double(sim$conf_label))
  legal <- bias[which(bias >= phi)]
  if (length(bias) == length(legal)) {
    message('\n+ bias parameters okay for ', group)
  } else {
    message('\n+ some bias parameters not okay for ', group, ', including: ',
            paste0(setdiff(bias, legal),
                   collapse = ','), '; removing!') }
  return(legal)

  # testing what's actually generated in idx.matrices
  # tmp <- c()
  # for (r in 1:50) {
  #   idx <- mat[r,]
  #   tmp <- c(tmp, calculate.phi(as.double(sim$label[idx]),
  #                               as.double(sim$conf_label[idx])))
  # }
}


#' @keywords internal
#' originally: randomize feature-metadata variable association, if any
#' currently: does not generate meaningful indices -- do not use !
save.biomarker.idx <- function(subs, reps, sim.file,
                               group, bias) {

  sim <- h5read(file = sim.file, name = group)
  params <- h5read(file = sim.file, name = 'simulation_parameters')
  # valid <- check.resampling.limits(sim.file, 'biomarker', group)

  all.idx.list <- list()
  for (b in seq_along(bias)) {
    all.idx.list[[paste0('bias_', b)]] <- list()
    bias.value <- bias[b]

    if (params$general_params$sim.method == 'pass') {
      message('+ Resampling biomarker indices for ', group,
              ' bias = ', bias[b]) }

    # generalize confounder label classes
    conf1 <- sort(unique(sim$conf_label))[1]
    conf2 <- sort(unique(sim$conf_label))[2]
    # mini df of label and confounder
    tmp <- cbind(sim$labels, sim$conf_label)
    rownames(tmp) <- sim$sample_names
    # initial class probabilities
    tbl <- table(sim$labels, sim$conf_label)/length(sim$labels)
    conf.probs <- colSums(tbl)
    # g1.probs <- c(sum(conf.probs)*(1-bias.value), sum(conf.probs)*bias.value)
    # g2.probs <- c(sum(conf.probs)*bias.value, sum(conf.probs)*(1-bias.value))

    tmp <- cbind(tmp, 0)
    tmp[which(tmp[,2]==conf1), 3] <- conf.probs[1]
    tmp[which(tmp[,2]==conf2), 3] <- conf.probs[2]

    # tmp[which(tmp[,1]==1 & tmp[,2]==conf1), 3] <- g1.probs[1]
    # tmp[which(tmp[,1]==1 & tmp[,2]==conf2), 3] <- g1.probs[2]
    # tmp[which(tmp[,1]==-1 & tmp[,2]==conf1), 3] <- g2.probs[1]
    # tmp[which(tmp[,1]==-1 & tmp[,2]==conf2), 3] <- g2.probs[2]

    # loop over subs & reps
    for (sub in subs) {
      mat <- matrix(NA, nrow = reps, ncol = sub)
      for (i in seq_len(reps)) {
        sub.label <- c(
        sample(which(tmp[,2] == conf1),
               prob = tmp[tmp[,2] == conf1,3],
               # maintain original class balance
               size=round(sub*table(sim$labels)[1]/length(sim$labels)),
               replace=TRUE),
        sample(which(tmp[,2] == conf2),
               prob = tmp[tmp[,2] == conf2,3],
               size=round(sub*table(sim$labels)[2]/length(sim$labels)),
               replace=TRUE))

        if (any(is.na(sub.label))) browser()
        if (length(sub.label) != sub) browser()

        mat[i,] <- sub.label
        # print(table(sim$labels[mat[i,]],
        #             sim$conf_label[mat[i,]]))
      }
      mat <- rbind(
        c(rep(-1, round(sub*table(sim$labels)[1]/length(sim$labels))),
          rep(1, round(sub*table(sim$labels)[2]/length(sim$labels)))), mat)
      rownames(mat) <- c('label', paste0('rep', seq_len(reps)))
      all.idx.list[[paste0('bias_', b)]][[paste0('subset_', sub)]] <- mat
    }
  }
  # browser()
  # anti-diagonal always 0 if bias = 1
  # table(sim$labels[mat[2,]], sim$conf_label[mat[2,]])
  h5write(obj = all.idx.list, file = sim.file,
          name = paste0(group, '/bio_bias_test_idx'))
  return(list('names'=rownames(mat), 'bias'=bias, 'subsets'=subs))
}

#' @keywords internal
#' originally: randomize feature-label association, if any
#' currently: bias parameter influences probability of a given label-metalabel
#' pair to be resampled (higher bias => higher phi coefficient)
save.confounder.idx <- function(subs, reps, sim.file,
                                group, bias) {

  # browser()
  params <- h5read(file = sim.file, name = 'simulation_parameters')
  sim <- h5read(sim.file, group)

  if (params$general_params$sim.method == 'pass') {
    bias <- check.resampling.limits(sim.file, group, bias)
    # just resample according to actual class balances..
    if (length(bias) < 1) {
      phi <- calculate.phi(as.double(sim$labels),
                           as.double(sim$conf_label))
      message('++ cannot create biased indices supplied for ', group,
              '; going ahead with natural phi coefficient: ', round(phi,
                                                                    digits=3))
      bias <- round(phi, digits=3) }
  }

  all.idx.list <- list()
  for (b in seq_along(bias)) {
    all.idx.list[[paste0('bias_', b)]] <- list()
    bias.value <- bias[b]

    # generalize confounder label classes
    conf1 <- sort(unique(sim$conf_label))[1]
    conf2 <- sort(unique(sim$conf_label))[2]
    # mini df of label and confounder
    tmp <- cbind(sim$labels, sim$conf_label)
    rownames(tmp) <- sim$sample_names
    # initial class probabilities
    tbl <- table(sim$labels, sim$conf_label)/length(sim$labels)
    # new probabilities come from weighting diagonal according to bias param
    conf.probs <- colSums(tbl)
    g1.probs <- c(sum(conf.probs)*bias.value, sum(conf.probs)*(1-bias.value))
    g2.probs <- c(sum(conf.probs)*(1-bias.value), sum(conf.probs)*bias.value)

    tmp <- cbind(tmp, 0)
    tmp[which(tmp[,1]==1 & tmp[,2]==conf1), 3] <- g1.probs[1]
    tmp[which(tmp[,1]==1 & tmp[,2]==conf2), 3] <- g1.probs[2]
    tmp[which(tmp[,1]==-1 & tmp[,2]==conf1), 3] <- g2.probs[1]
    tmp[which(tmp[,1]==-1 & tmp[,2]==conf2), 3] <- g2.probs[2]

    # loop over subs & reps
    for (sub in subs) {
      mat <- matrix(NA, nrow = reps, ncol = sub)
      for (i in seq_len(reps)) {
        sub.label <- c(
          sample(which(tmp[,1] == -1),
                 prob = tmp[tmp[,1]==-1,3],
                 # maintain original class balance wrt label
                 size=round(sub*table(sim$labels)[1]/length(sim$labels)),
                 replace=TRUE),
          sample(which(tmp[,1] == 1),
                 prob = tmp[tmp[,1]==1,3],
                 size=round(sub*table(sim$labels)[2]/length(sim$labels)),
                 replace=TRUE))
        if (any(is.na(sub.label))) browser()
        if (length(sub.label) != sub) browser()

        mat[i,] <- sub.label
        # print(table(sim$labels[mat[i,]],
        #             sim$conf_label[mat[i,]]))
      }
      # first row of matrix stores (main) label
      mat <- rbind(
        c(rep(-1, round(sub*table(sim$labels)[1]/length(sim$labels))),
          rep(1, round(sub*table(sim$labels)[2]/length(sim$labels)))), mat)
      rownames(mat) <- c('label', paste0('rep', seq_len(reps)))
      all.idx.list[[paste0('bias_', b)]][[paste0('subset_', sub)]] <- mat
    }
  }
  # browser()
  # mat[1,]
  # table(sim$labels[mat[2,]], sim$conf_label[mat[2,]])
  h5write(obj = all.idx.list, file = sim.file,
          name = paste0(group, '/conf_bias_test_idx'))
  if (params$general_params$sim.method == 'pass') {
    h5write(obj = as.array(bias), file = sim.file,
            name = paste0(group, '/conf_bias_values'))
    h5write(obj = names(all.idx.list), file = sim.file,
            name = paste0(group, '/conf_bias_names'))
  }
  return(list('names'=rownames(mat), 'bias'=bias, 'subsets'=subs))
}

#' @keywords internal
save.random.idx <- function(subs, reps, sim.file, group) {

  sim <- rhdf5::h5read(sim.file, group)
  params <- h5read(file = sim.file, name = 'simulation_parameters')

  idx.list <- list()

  # loop over subs & reps
  for (sub in subs) {
    mat <- matrix(NA, nrow = reps, ncol = sub)
    for (i in seq_len(reps)) {
      # take any samples, will hard code label in idx.mat anyway
      sub.label <- sample(sim$sample_names,
                          size = sub,
                          replace = TRUE)

      if (any(is.na(sub.label))) browser()
      if (length(sub.label) != sub) browser()

      sub.idx <- match(sub.label, sim$sample_names)
      mat[i,] <- sub.idx
    }
    # hard code label ignoring actual label
    mat <- rbind(c(rep(-1, round(sub*0.5)),
                   rep(1, round(sub*0.5))),
                 mat)
    rownames(mat) <- c('label', paste0('rep', seq_len(reps)))
    idx.list[[paste0('subset_', sub)]] <- mat
  }
  # browser()
  h5write(obj = idx.list, file = sim.file,
          name = paste0(group, '/random_test_idx'))
  return(list('names'=rownames(mat), 'subsets'=subs))

}

#' @keywords internal
clear.biased.idx <- function(sim.location, idx.type) {

  all <- h5ls(sim.location, recursive = 2)
  g.idx <- which(all$name == paste0(idx.type, '_test_idx'))

  if (length(g.idx) == 0) {
    message(paste0('++ Biased test indices do not exist yet for type ',
                   idx.type))
  } else {
    for (g in g.idx){
      group <- paste0(all$group[g], paste0('/', idx.type, '_test_idx'))
      h5delete(file = sim.location, name = group)
    }
    message(paste0('++ Finished removing existing biased test indices of type ',
                   idx.type))
  }
}
