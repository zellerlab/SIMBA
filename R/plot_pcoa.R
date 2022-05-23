#!/usr/bin/Rscript
### SIMBA
### Simulation of Microbiome data with Biological Accuracy
### Morgan Essex - MDC Berlin
### Jakob Wirbel - EMBL Heidelberg
### 2020 GNU GPL 3.0

#' @title Plot PCoAs for simulated data
#'
#' @description This function will generate a principal coordinate analysis
#' plot comparing the original data and the simulated data for different groups
#'
#' @usage pcoa.plot(sim.location, group)
#'
#' @param sim.location file name for the .h5 file containing the simulations
#'
#' @param group character, name of the group within the .h5 file containing
#' the simulated features
#'
#' @param distance character, name of the distance that should be used to
#' calculate the PCoA, should be either \code{"bray"} or \code{"log-euclidean"}
#'
#' @export
#'
#' @keywords SIMBA pcoa.plot
#'
#' @return This function returns a ggplot2 object with a PCoA plot, coloured by
#' group (either real or simulated data).
#'
pcoa.plot <- function(sim.location, group, distance){
  test.package("ggplot2")

  message("+ Checking parameters")
  log.n0 <- pcoa.check.parameters(sim.location, group, distance)

  # get the data
  # original data
  feat.filt <- h5read(file=sim.location,
                      name='/original_data/filt_features')
  colnames(feat.filt) <- h5read(file=sim.location,
                                name='/original_data/sample_names')
  rownames(feat.filt) <- h5read(file=sim.location,
                                name='/original_data/filt_feature_names')
  meta <- suppressWarnings(h5read(file=sim.location,
                                  name='/original_data/metadata'))
  rownames(meta) <- colnames(feat.filt)
  # simulated data
  feat.sim <- h5read(file=sim.location,
                     name=paste0('/',group, '/features'))
  colnames(feat.sim) <- h5read(file=sim.location,
                               name=paste0('/',group, '/sample_names'))
  rownames(feat.sim) <- h5read(file=sim.location,
                               name=paste0('/',group, '/feature_names'))
  label <- h5read(file=sim.location, name=paste0('/', group, '/labels'))

  if (any(colSums(feat.sim) == 0)){
      idx <- which(colSums(feat.sim)==0)
      feat.sim <- feat.sim[,-idx]
      label <- label[-idx]
  }
  names(label) <- colnames(feat.sim)

  sim.params <- h5read(file=sim.location, name='simulation_parameters')
  conf <- FALSE
  if (!is.null(sim.params$sim_params$conf)){
    if (sim.params$sim_params$conf != 'None'){
      conf <- TRUE
    }
  }
  if (conf){
    # TODO
    if (sim.params$sim_params$conf=='batch'){
      df.meta <- data.frame(group='original',
                            confounder=meta$Study,
                            names=paste0('original_',
                                         seq_len(ncol(feat.filt))))
      df.meta <- rbind(
        df.meta,
        data.frame(group=paste0('simulated-group',
                                ifelse(label==1, 1, 2)),
                   confounder=paste0('Study_',
                                     h5read(file=sim.location,
                                            name=paste0('/',group,
                                                        '/conf_label'))),
                   names=paste0('simulated_', seq_len(ncol(feat.filt)))))
    } else {
      df.meta <- data.frame(group='original',
                            confounder=0,
                            names=paste0('original_',
                                         seq_len(ncol(feat.filt))))
      df.meta <- rbind(df.meta,
                       data.frame(group=paste0('simulated-group',
                                               ifelse(label==1, 1, 2)),
                                  confounder=h5read(file=sim.location,
                                                    name=paste0('/',group,
                                                                '/conf_label')),
                                  names=paste0('simulated_',
                                               seq_len(ncol(feat.filt)))))
    }
  } else {
    df.meta <- data.frame(group='original',
                          names=paste0('original_', seq_len(ncol(feat.filt))))
    df.meta <- rbind(df.meta,
                     data.frame(group=paste0('simulated-group',
                                             ifelse(label==1, 1, 2)),
                                names=paste0('simulated_',
                                             seq_len(ncol(feat.sim)))))
  }

  # check if the feature names match
  x <- rownames(feat.filt)
  y <- rownames(feat.sim)
  stopifnot(all(x == y))
  colnames(feat.filt) <- paste0("original_", seq_len(ncol(feat.filt)))
  colnames(feat.sim) <- paste0("simulated_", seq_len(ncol(feat.sim)))

  # if too large, subsample and combine
  if (ncol(feat.filt) > 200){
    message("++ Large dataset will be subsampled to 200 samples")
    tmp <- cbind(feat.filt[,sample(colnames(feat.filt), 200)],
                 feat.sim[,sample(colnames(feat.sim)[label==1], 100)],
                 feat.sim[,sample(colnames(feat.sim)[label==-1], 100)])
  } else {
    tmp <- cbind(feat.filt, feat.sim)
  }
  # calculate distance matrix
  if (distance=='bray'){
    bc <- as.matrix(vegdist(t(tmp)))
  } else {
    tmp <- log10(tmp + log.n0)
    bc <- as.matrix(vegdist(t(tmp), method='euclidean'))
  }
  # calculate pcoa
  pco.res <- pco(bc)
  df.plot <- as.data.frame(pco.res$points)
  # join with metadata
  df.plot <- cbind(df.plot, df.meta[match(rownames(df.plot), df.meta$names),])

  # plot
  axis.labels <- sprintf(fmt='%.2f',
                         (pco.res$eig[1:2]/sum(
                           pco.res$eig[pco.res$eig > 0])) * 100)

  g <- ggplot(df.plot, aes(x=V1, y=V2, col=group)) +
    theme_classic() +
    scale_colour_manual(values = c('#307FE275', '#FFA30075', '#E4004675'),
                        name='') +
    xlab(paste0('PCo 1 [', axis.labels[1], '%]')) +
    ylab(paste0('PCo 2 [', axis.labels[2], '%]'))
  if ('confounder' %in% colnames(df.meta)){
    if (any(str_detect(df.meta$confounder, 'Study'))){
      g <- g +
        geom_point(aes(shape=as.factor(confounder))) +
        scale_shape_manual(values=c(19, 17), name='Study')
    } else {
      g <- g +
        geom_point(aes(shape=as.factor(confounder))) +
        scale_shape_manual(values=c(19, 3, 17),
                           name=paste0('Confounder: ',
                                       sim.params$sim_params$conf),
                           labels=c('0'='original', '-1'='positive',
                                    '1'='negative')
                           )
    }
  } else {
    g <- g + geom_point()
  }
  return(g)
}

#' @keywords internal
pcoa.check.parameters <- function(sim.location, group, distance){
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
    stop("PCoA plots are not applicable for the 'pass' method!")
  }
  log.n0 <- as.vector(params$filt_params$log.n0)
  params <- params$sim_params
  needed <- setdiff(c('ab.scale', 'repeats'), names(params))
  if (length(needed) > 0){
    stop("Simulation layout is not suitable for checks!\n",
         "Missing the parameters:", paste(needed, collapse = ', '))
  }
  # distance
  if (!is.character(distance)){
    stop("Parameter 'distance' must be a character!")
  }
  if (length(distance) != 1){
    stop("Parmeter 'distance' must be of length 1!")
  }
  if (!distance %in% c("bray", 'log-euclidean')){
    stop("Parmeter 'distance' must be either 'bray' or 'log-euclidean'!")
  }


  all.groups <- h5ls(file=sim.location, recursive = 2)
  if (!group %in% all.groups$name){
    stop("Group '", group, "' not found in the h5 file!")
  }
  red <- all.groups[all.groups$group==paste0('/', group),]
  if (length(setdiff(c('feature_names', 'features', 'labels',
                       'marker_idx', 'sample_names'),red$name)) != 0){
    stop('Could not find all entries for this group!\n',
         'Something seems to have gone wrong before!')
  }
  return(log.n0)
}
