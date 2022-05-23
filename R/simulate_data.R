#!/usr/bin/Rscript
### SIMBA package
### Simulation of Microbiome data with Biological Accuracy
### Morgan Essex - MDC Berlin
### Jakob Wirbel - EMBL Heidelberg
### 2020 GNU GPL 3.0

#' @title Simulate metagenomic features
#'
#' @description This function uses a feature and metadata table
#' as well as specific parameters to create simulated metagenomic samples
#' and stores those in the indicated .h5 file
#'
#' @usage simulate.data(feat, meta, sim.location,
#' sim.method='resampling',
#' sim.type="cross-section",
#' filt.params=list(prev.cutoff=0.05, ab.cutoff=1e-04, log.n0=1e-05),
#' sim.params=list(class.balance=0.5, prop.markers=0.1,
#' ab.scale=c(1.0, 1.25, 1.5, 2.0, 5.0, 10.0, 20.0),
#' prev.scale=c(0.0, 0.1, 0.2, 0.3),
#' repeats=100, feature.type='all'))
#'
#' @param feat feature matrix, can be either a single matrix or a list of
#' matrices from different studies
#'
#' @param meta metadata about samples in the features matrix, should be a
#' dataframe or a list of dataframes
#'
#' @param sim.location filename for an .h5 file which will store the
#' simulations, original data, and simulation parameters
#'
#' @param sim.method methodology to create the simulations. Can either be
#' a string describing a pre-implemented method or a user-specified function
#' (see Details below). Defaults to \code{"resampling"}
#'
#' @param sim.type character, can be either \code{"cross-section"} or
#' \code{"time-course"}, defaults to \code{'cross-section'}
#'
#' @param sim.params list of parameters which are given to
#' the \code{sim.method} function, see Details
#'
#' @param filt.params list of parameters for data filtering, defaults
#' to \code{filt.params=list(prev.cutoff=0.05, ab.cutoff=1e-04,
#' log.n0=1e-05)}, see Details below
#'
#' @export
#'
#' @keywords SIMBA simulate.data
#'
#' @return Does not return anything but instead creates an h5 file
#' containing simulated samples
#'
#' @details This functions TODO
#'
#' \bold{Filtering}
#' \bold{Pre-implemented Simulations}
#' \bold{Simulations}
#' \bold{h5-File Organisation}
#'
simulate.data <- function(feat, meta, sim.location,
                          sim.method='resampling',
                          sim.type="cross-section",
                          filt.params=list(prev.cutoff=0.05,
                                          ab.cutoff=1e-04,
                                          log.n0=1e-05),
                          sim.params=list(class.balance=0.5, prop.markers=0.1,
                                          conf='None',
                                          ab.scale=c(1.0, 1.25, 1.5, 2.0,
                                                     5.0, 10.0, 20.0),
                                          prev.scale=c(0.0, 0.1, 0.2, 0.3),
                                          repeats=100,
                                          feature.type='all')){

  # check data
  message("+ Start checking data")
  res <- check.original.data(feat, meta, sim.type, sim.method)
  message("+ Finished checking data")
  message("")
  feat.original <- res$feat
  meta.original <- res$meta

  # check sim.location
  message("+ Start checking the simulation location")
  check.simulation.location(sim.location)
  message("+ Finished checking the simulation location")
  message("")

  # check filtering parameters
  message("+ Start filtering the data")
  filt.params <- check.filtering.parameters(filt.params)
  # filter features
  feat.filt <- filter.data(feat=feat.original, meta=meta.original,
                           filt.params=filt.params)
  message("+ Finished filtering the data")
  message("")

  # check simulation parameters
  message("+ Start checking the simulation parameters")
  sim.params <- check.simulation.parameters(sim.method, sim.type,
                                            sim.params, meta.original)
  message("+ Finished checking the simulation parameters")
  message("")

  # check putative confounders for associations based on variance analysis
  # confounder.check(feat.filt, meta.original)
  # Keep or remove?
  # TODO

  message("+ Save original data in the h5 file")
  # save original/filtered features in the h5 file
  stopifnot(h5createGroup(sim.location, 'original_data'))
  stopifnot(all(colnames(feat.original) == rownames(meta.original)))
  stopifnot(all(colnames(feat.original) == colnames(feat.filt)))
  suppressMessages(h5write(obj=feat.original, file=sim.location,
                           name='original_data/features'))
  suppressMessages(h5write(obj=rownames(feat.original), file=sim.location,
                           name='original_data/feature_names'))
  h5write(obj=colnames(feat.original), file=sim.location,
          name='original_data/sample_names')
  h5write(obj=feat.filt, file=sim.location,
          name='original_data/filt_features')
  h5write(obj=rownames(feat.filt), file=sim.location,
          name='original_data/filt_feature_names')
  h5write(obj=meta.original, file=sim.location, # TODO
          name='original_data/metadata')

  # save the simulation parameter set in the h5 file
  stopifnot(h5createGroup(sim.location, 'simulation_parameters'))

  h5write(obj=list(sim.method=sim.method, sim.type=sim.type),
          file=sim.location,
          name='simulation_parameters/general_params')
  h5write(obj=filt.params, file=sim.location,
          name='simulation_parameters/filt_params')
  h5write(obj=sim.params, file=sim.location,
          name='simulation_parameters/sim_params')
  message("+ Finished saving original data in the h5 file")
  message("")

  message("+ Starting data generation using the method: ", sim.method)
  # simulate data
  if (sim.method == 'resampling'){
    simulate.resampling(feat.filt, meta.original, sim.location, sim.params)
  } else if (sim.method == 'McMurdie&Holmes'){
    simulate.MMH(feat.filt, meta.original, sim.location, sim.params)
  } else if (sim.method == 'Weiss'){
    simulate.W(feat.filt, meta.original, sim.location, sim.params)
  } else if (sim.method == 'negbin'){
    simulate.negbin(feat.filt, meta.original, sim.location, sim.params)
  } else if (sim.method == 'dirmult'){
    simulate.dirmult(feat.filt, meta.original, sim.location, sim.params)
  } else if (sim.method == 'betabin'){
    simulate.betabin(feat.filt, meta.original, sim.location, sim.params)
  } else if (sim.method=='sparseDOSSA'){
    simulate.sparseDOSSA(feat.filt, meta.original, sim.location, sim.params)
  } else if (sim.method=='pass'){
    label <- as.numeric(as.character(meta.original[,'Label']))
    names(label) <- rownames(meta.original)
    h5write(label, sim.location, paste0('original_data', '/labels'))
    message("++ No artifical/resampled data is created -- labels saved!")
  }
  h5closeAll()
  message("+ Finished data generation")
}
