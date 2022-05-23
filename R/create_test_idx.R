#!/usr/bin/Rscript
### SIMBA
### Simulation of Microbiome data with Biological Accuracy
### Morgan Essex - MDC Berlin
### Jakob Wirbel - EMBL Heidelberg
### 2020 GNU GPL 3.0

#' @title Create indices for testing
#'
#' @description This function creates a matrix of test indices for each group
#' in the simulation file on which differential abundance testing methods can
#' be run
#'
#' @usage create.test.idx(sim.location, subsets, repetitions)
#'
#' @param sim.location file name for the .h5 file containing the simulations
#'
#' @param subsets vector of subset sizes to be sampled
#'
#' @param repetitions integer, number of repetitions for each group/subset
#' combination
#'
#' @export
#'
#' @keywords SIMBA create.test.idx
#'
#' @return Returns nothing, but modifies the simulation file
create.test.idx <- function(sim.location, subsets, repetitions){

  message("+ Starting to check parameters for test idx creation")
  # check parameters
  # check the H5 file
  if (!file.exists(sim.location)){
    stop("Simulation file does not exist!")
  }
  params <- h5read(file = sim.location, name = 'simulation_parameters')
  if (params$general_params$sim.type == 'time-course'){
    stop('Not implemented yet!')
  }

  # check repetitions
  if (!is.numeric(repetitions)){
    stop("Parameter 'repetitions' must be numeric!")
  } else if (length(repetitions) > 1){
    stop("Parameter 'repetitions' needs to be of length 1!")
  }
  # remove idx if present already
  clear.idx(sim.location)

  # get original size of feature table
  n.samples <- nrow(suppressWarnings(h5read(sim.location,
                                            '/original_data/metadata')))
  sim.params <- params$sim_params
  class.balance <- 0.5
  if ('class.balance' %in% names(sim.params)){
    class.balance <- sim.params$class.balance
  }

  # check subset size
  if (!is.numeric(subsets)){
    stop("Parameter 'subsets' must be numeric!")
  }
  subsets <- sort(unique(subsets))
  if (any(subsets < 0)){
    message("++ All entries in parameter 'subsets' need to be positive!")
    subsets <- subsets[subsets > 0]
    message("++ Parameter 'subset' reduced to: ",
            paste(subsets, collapse = ', '))
  }
  if (any(subsets > n.samples)){
    message("++ Some entries in parameter 'subsets' are too large!")
    subsets <- subsets[subsets <= n.samples]
    message("++ Parameter 'subset' reduced to: ",
            paste(subsets, collapse = ', '))
  }

  message("+ Start test idx creation")
  if (params$general_params$sim.method == 'pass') {
    temp <- save.idx(NULL, NULL, subsets, repetitions,
                     class.balance, sim.location)
  } else {
    ## loop through subset sizes and repetitions and simulation repetitions
    if ('prev.scale' %in% names(sim.params)){
      pb <- progress_bar$new(total=length(sim.params$ab.scale) *
                               length(sim.params$prev.scale) *
                               sim.params$repeats)
      for (ab in seq_along(sim.params$ab.scale)){
        for (r in seq_len(sim.params$repeats)){
          for (p in seq_along(sim.params$prev.scale)){
            temp <- save.idx(ab, r, subsets, repetitions, class.balance,
                             sim.location, prev=p)
            pb$tick()
          }
        }
      }
    } else {
      pb <- progress_bar$new(total=length(sim.params$ab.scale) *
                               sim.params$repeats)
      for (ab in seq_along(sim.params$ab.scale)){
        for (r in seq_len(sim.params$repeats)){
          temp <- save.idx(ab, r, subsets, repetitions,
                           class.balance, sim.location)
          pb$tick()
        }
      }
    }
  }
  message("+ Finished test idx creation")

  # save rownames for idx matrices
  h5write(temp, file=sim.location,
          name='/simulation_parameters/test_idx')
}

#' @keywords internal
save.idx <- function(ab, rep, subs, reps, class.balance, sim.file,
                     prev=NULL, ...){
  additional.arguments <- list(...)
  # load the original or specific instance of simulated data
  if ('group' %in% names(additional.arguments)) {
    group <- additional.arguments$group }
    #if (group=='ace_inhib') browser() }
  else if (is.null(ab) & is.null(prev)) { group <- 'original_data'}
  else if (is.null(prev)) {
    group <- paste0('ab', ab, '_rep', rep)
  } else {
    group <- paste0('ab', ab, '_prev', prev, '_rep', rep)
  }
  sim.label <- h5read(sim.file, paste0('/', group, '/labels'))
  pos.samples <- which(sim.label == -1)
  neg.samples <- which(sim.label == 1)

  idx.list <- list()
  # create indices
  for (s in subs){
    mat <- matrix(NA, nrow=reps, ncol=s)
    for (r in seq_len(reps)){
      if ('replace' %in% names(additional.arguments)) {
        sub.label <- c(sample(pos.samples, s*(1-class.balance),
                              replace = TRUE),
                       sample(neg.samples, s*class.balance,
                              replace = TRUE)) }
      else {
        sub.label <- c(sample(pos.samples, s*(1-class.balance)),
                       sample(neg.samples, s*class.balance)) }
      stopifnot(length(sub.label) == s)
      mat[r,] <- sub.label
    }
    fix.label <- c(rep(-1, s*(1-class.balance)), rep(1, s*class.balance))
    mat <- rbind(fix.label, mat)
    rownames(mat) <- c('label', paste0('rep', seq_len(reps)))
    idx.list[[paste0('subset_', s)]] <- mat
  }

  # save the indices
  h5write(obj=idx.list, file=sim.file,
          name=paste0(group, '/test_idx'))
  return(list('names'=rownames(mat), 'subsets'=subs))
}

#' @keywords internal
clear.idx <- function(sim.location){
  message("++ Remove test idx, if present")
  # get test_idx groups
  all <- h5ls(sim.location, recursive = 2)

  g.idx <- which(all$name=='test_idx')
  if (length(g.idx)==0){
    message("++ Test idx do not exist yet")
  } else {
    for (g in g.idx){
      group <- paste0(all$group[g], '/test_idx')
      h5delete(file=sim.location, name=group)
    }
    message("++ Finished to remove test idx")
  }
}
