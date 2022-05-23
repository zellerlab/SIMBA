#!/usr/bin/Rscript
### SIMBA
### Simulation of Microbiome data with Biological Accuracy
### Morgan Essex - MDC Berlin
### Jakob Wirbel - EMBL Heidelberg
### 2020 GNU GPL 3.0

# Helper functions to simulate data according to the resampling method

# wrapper for the simulation functions
#' @keywords internal
simulate.resampling <- function(feat, meta, sim.out, sim.params){
  # get parameters
  class.balance <- sim.params$class.balance
  prop.markers <- sim.params$prop.markers
  conf <- sim.params$conf
  conf.params <- sim.params$conf.params
  ab.scale <- sim.params$ab.scale
  prev.scale <- sim.params$prev.scale
  repeats <- sim.params$repeats
  feature.type <- sim.params$feature.type
  if ('balanced' %in% names(sim.params)){
    balanced <- sim.params$balanced
  } else {
    balanced <- TRUE
  }


  if (conf == 'None'){
    conf <- NULL
  }
  no.marker.feat <- round(prop.markers * nrow(feat))

  # compute feature selection probablity
  sel.prob <- select.features(feat, feature.type)
  el.feat.names <- sel.prob$names
  feat.prob <- sel.prob$prob

  # probabilities of class membership (depending on confounder)
  num.sample <- ncol(feat)
  # default is uniform
  prob <- rep(1/num.sample, num.sample)
  names(prob) <- rownames(meta)
  if (!is.null(conf)){
    if (conf == 'artificial'){
      message('++ Implanting artificial confounder with strength ',
              conf.params$bias)
    } else if (conf=='global') {
      message('++ Implanting global confounder with bias ',
              conf.params$bias)
      x <- vegan::vegdist(t(feat))
      pco.res <- labdsv::pco(x)
      pco.res <- as.data.frame(pco.res$points)
      conf.prob <- pco.res$V1
      conf.prob <- (conf.prob - min(conf.prob))/
        (max(conf.prob) - min(conf.prob))
      names(conf.prob) <- rownames(pco.res)
    } else if (conf=='batch') {
      message('++ Implanting confounder ', conf, ' with strength ',
              conf.params$bias)
      temp <- as.numeric(as.factor(meta[,'Study']))
      names(temp) <- rownames(meta)
      if (any(is.na(temp))){
        temp <- temp[!is.na(temp)]
      }
      groups <- unique(temp)
      group.probs <- seq(from=1-conf.params$bias, to=1-(1-conf.params$bias),
                         length.out=length(groups)) * 1/length(groups)
      for (i in seq_along(groups)){
        g <- groups[i]
        group.probs[i] <- group.probs[i] * sum(temp==g)
      }
      group.probs <- group.probs/sum(group.probs)
      prob <- rep(NA, length(temp))
      for (i in seq_along(groups)){
        g <- groups[i]
        prob[temp==g] <- group.probs[i]/sum(temp==g)
      }
      names(prob) <- names(temp)
      num.sample <- length(prob)
      conf.label <- temp
    } else {
      # biased sampling of groups
      message('++ Implanting confounder ', conf, ' with strength ',
              conf.params$bias)
      temp <- as.numeric(as.factor(meta[,conf]))
      names(temp) <- rownames(meta)
      if (any(is.na(temp))){
        temp <- temp[!is.na(temp)]
      }
      groups <- sort(unique(temp))
      group.probs <- seq(from=1-conf.params$bias, to=1-(1-conf.params$bias),
                         length.out=length(groups)) * 1/length(groups)
      for (i in seq_along(groups)){
        g <- groups[i]
        group.probs[i] <- group.probs[i] * sum(temp==g)
      }
      group.probs <- group.probs/sum(group.probs)
      prob <- rep(NA, length(temp))
      for (i in seq_along(groups)){
        g <- groups[i]
        prob[temp==g] <- group.probs[i]/sum(temp==g)
      }
      names(prob) <- names(temp)
      num.sample <- length(prob)
      conf.label <- temp
    }
  }

  # actual implantation
  pb <- progress_bar$new(total = length(ab.scale)*length(prev.scale)*repeats)
  for (a in seq_along(ab.scale)){
    for (b in seq_along(prev.scale)){
      for (r in seq_len(repeats)){

        if (!is.null(conf)){
          if (conf=='artificial'){
            conf.label <- rep(-1, num.sample)
            names(conf.label) <- rownames(meta)
            conf.label[sample(seq_along(conf.label),
                             size=floor(length(conf.label) *
                                          conf.params$prop))] <- 1
            prob <- rep(NA, length(conf.label))
            prob.scale <- c(conf.params$bias*sum(conf.label==1),
                            (1-conf.params$bias)*sum(conf.label==-1))
            prob.scale <- prob.scale/sum(prob.scale)
            prob[conf.label==1] <- prob.scale[1]/sum(conf.label==1)
            prob[conf.label==-1] <- prob.scale[2]/sum(conf.label==-1)
            names(prob) <- names(conf.label)
          } else if (conf=='global'){
            conf.label <- rep(-1, num.sample)
            names(conf.label) <- rownames(meta)
            conf.label[sample(seq_along(conf.label),
                              size=floor(length(conf.label) *
                                           conf.params$prop),
                              prob = conf.prob)] <- 1
            prob <- rep(NA, num.sample)
            prob.scale <- c(conf.params$bias*sum(conf.label==1),
                            (1-conf.params$bias)*sum(conf.label==-1))
            prob.scale <- prob.scale/sum(prob.scale)
            prob[conf.label==1] <- prob.scale[1]/sum(conf.label==1)
            prob[conf.label==-1] <- prob.scale[2]/sum(conf.label==-1)
            names(prob) <- names(conf.label)
          }
        }

        # generate label
        label <- rep(-1, num.sample)
        s <- sample(num.sample,
                    size=round(num.sample*class.balance),
                    prob=prob+1e-20)
        names(label) <- names(prob)
        label[s] <- +1 # switch random half of the samples to positive class

        feat.tmp <- feat[,names(label)]

        # sample markers
        marker.idx <- sample(el.feat.names,
                             size = no.marker.feat,
                             prob = feat.prob)

        # if artificial confounder, select conf.markers
        if (!is.null(conf)){
          if (conf == 'artificial'){
            red.el <- setdiff(rownames(feat), marker.idx)
            conf.marker.idx <- sample(red.el,
                                      size=nrow(feat) * conf.params$feat.prop)
          } else {
            conf.marker.idx <- NULL
          }
        } else {
          conf.marker.idx <- NULL
        }

        # actually simulate
        sim.feat <- create.resampling.simulations(
          feat.tmp=feat.tmp,
          marker.idx=marker.idx,
          label=label,
          conf.label=conf.label,
          prev.shift=prev.scale[b],
          ab.scale=ab.scale[a],
          balanced=balanced,
          conf.marker.idx=conf.marker.idx,
          conf.scale=conf.params$feat.effect)
        
        # gFC = 0 rejection
        sim.feat.rel <- prop.table(sim.feat, 2)
        n.pos <- names(label[label==1])
        n.neg <- names(label[label==-1])
        gFC <- vapply(marker.idx, FUN = function(x){
          g.pos <- quantile(log10(sim.feat.rel[x,n.pos] + 1e-05),
                            prob=seq(from=0.05, to=0.95, by=0.05))
          g.neg <- quantile(log10(sim.feat.rel[x,n.neg] + 1e-05),
                            prob=seq(from=0.05, to=0.95, by=0.05))
          abs(mean(g.pos-g.neg))
        }, FUN.VALUE = double(1))
        if (any(gFC < 0.001)){
          # replace those and remove them from the index
          replace.idx <- names(which(gFC < 0.001))
          sim.feat[replace.idx,] <- feat.tmp[replace.idx,colnames(sim.feat)]
          marker.idx <- setdiff(marker.idx, replace.idx)
        }

        # save data in H5-file
        h5.subdir <- paste0('ab', a, '_prev', b, '_rep', r)
        stopifnot(h5createGroup(sim.out, h5.subdir))
        h5write(label, sim.out, paste0(h5.subdir, '/labels'))
        h5write(sim.feat, sim.out, paste0(h5.subdir, '/features'))
        h5write(marker.idx, sim.out, paste0(h5.subdir, '/marker_idx'))
        if (!is.null(conf.marker.idx)){
          h5write(conf.marker.idx, sim.out, paste0(h5.subdir,
                                                   '/marker_idx_conf'))
        }
        if (!is.null(conf)){
          h5write(conf.label, sim.out, paste0(h5.subdir, '/conf_label'))
        }
        h5write(colnames(sim.feat), sim.out,
                paste0(h5.subdir, '/sample_names'))
        h5write(rownames(sim.feat), sim.out,
                paste0(h5.subdir, '/feature_names'))
        pb$tick()
      }
    }
  }
}

# simulation function for the resampling type of simulations
#' @keywords internal
create.resampling.simulations <- function(feat.tmp,
                                          marker.idx,
                                          label,
                                          conf.label,
                                          prev.shift,
                                          ab.scale,
                                          balanced=TRUE,
                                          conf.marker.idx=NULL,
                                          conf.scale=NULL){
  stopifnot(all(label %in% c(1,-1)))
  stopifnot(all(names(label) %in% colnames(feat.tmp)))

  # set target group for implantation (will switch in each loop below)
  target.group <- 1

  sim.feat <- feat.tmp
  # implant marker features
  for (idx in marker.idx) {

    # switch between target groups to implement features positively AND
    #   negatively in the same (arbitrary group)
    if (balanced) target.group <- -target.group

    # Prevalence shift
    prev <- mean(as.numeric(feat.tmp[idx,]) > 0)
    # swap a few samples with non-zero abundance to create a class-bias
    n.shift <- min(c(round(prev.shift*prev*ncol(feat.tmp)),
                     sum(label==target.group &
                           as.numeric(feat.tmp[idx,]) == 0),
                     sum(label==-target.group &
                           as.numeric(feat.tmp[idx,]) > 0)))
    repeat {
      sp <- sample(which(label==target.group & as.numeric(feat.tmp[idx,]) == 0),
                   size=n.shift)
      sn <- sample(which(label==-target.group & as.numeric(feat.tmp[idx,]) > 0),
                   size=n.shift)
      stopifnot(length(sp) == length(sn))
      if (length(intersect(sp, sn)) == 0)
        break
    }
    tmp <- feat.tmp[idx, sp]
    sim.feat[idx, sp] <- feat.tmp[idx, sn]
    sim.feat[idx, sn] <- tmp

    # scale the abundances in a class-dependent manner
    sim.feat[idx, label==target.group] <-
      round(sim.feat[idx, label==target.group]*ab.scale)
    sim.feat[idx, label==-target.group] <-
      round(sim.feat[idx, label==-target.group]*(1/ab.scale))
  }

  # if applicable, implant confounder markers
  if (!is.null(conf.marker.idx)){
    stopifnot(all(conf.label %in% c(1,-1)))
    target.group <- 1
    if (!is.null(conf.scale)){
      for (idx in conf.marker.idx){
        if (balanced) target.group <- -target.group
        # scale the abundances in a class-dependent manner
        sim.feat[idx, conf.label==target.group] <-
          round(sim.feat[idx, conf.label==target.group]*conf.scale)
        sim.feat[idx, conf.label==-target.group] <-
          round(sim.feat[idx, conf.label==-target.group]*(1/conf.scale))
      }
    } else {
      for (idx in conf.marker.idx){
        if (balanced) target.group <- -target.group
        prev <- mean(as.numeric(feat.tmp[idx,]) > 0)
        n.shift <- min(c(round(prev.shift*prev*ncol(feat.tmp)),
                         sum(conf.label==target.group &
                               as.numeric(feat.tmp[idx,]) == 0),
                         sum(conf.label==-target.group &
                               as.numeric(feat.tmp[idx,]) > 0)))
        repeat {
          sp <- sample(which(conf.label==target.group &
                               as.numeric(feat.tmp[idx,]) == 0),
                       size=n.shift)
          sn <- sample(which(conf.label==-target.group &
                               as.numeric(feat.tmp[idx,]) > 0),
                       size=n.shift)
          stopifnot(length(sp) == length(sn))
          if (length(intersect(sp, sn)) == 0)
            break
        }
        tmp <- feat.tmp[idx, sp]
        sim.feat[idx, sp] <- feat.tmp[idx, sn]
        sim.feat[idx, sn] <- tmp

        # scale the abundances in a class-dependent manner
        sim.feat[idx, conf.label==target.group] <-
          round(sim.feat[idx, conf.label==target.group]*ab.scale)
        sim.feat[idx, conf.label==-target.group] <-
          round(sim.feat[idx, conf.label==-target.group]*(1/ab.scale))
      }
    }


  }
  if (!balanced){
    sim.feat <- t(vegan::rrarefy(t(sim.feat), colSums(feat.tmp)))
  }
  return(sim.feat)
}

# feature selection function
#' @keywords internal
select.features <- function(feat, type){

  if (type == 'all'){
    el.feat.names <- rownames(feat)
    feat.prob <- rep(1, length(el.feat.names))
  } else {
    df.select <- data.frame(
      idx=rownames(feat),
      median=rowMedians(prop.table(feat, 2)),
      log.median=rowMedians(log10(prop.table(feat, 2) + 1e-05)),
      mean=rowMeans(prop.table(feat, 2)),
      log.mean=rowMeans(log10(prop.table(feat, 2) + 1e-05)),
      quant75=rowQuantiles(prop.table(feat, 2), probs=0.75),
      stringsAsFactors = FALSE)
    df.select$type <- ifelse(df.select$median==0 & df.select$quant75==0, 'low',
                             ifelse(df.select$median==0 & df.select$quant75!=0,
                                    'middle','high'))
    if (type %in% c('abundance', 'inverse-abundance')){
      el.feat.names <- rownames(feat)
      if (type=='abundance'){
        feat.prob <- df.select$log.mean
        feat.prob <- feat.prob - min(feat.prob) + 0.05
        feat.prob <- feat.prob/(max(feat.prob) + 0.05*max(feat.prob))
        feat.prob <- feat.prob/sum(feat.prob)
      } else {
        feat.prob <- df.select$log.mean
        feat.prob <- feat.prob - min(feat.prob) + 0.05
        feat.prob <- feat.prob/(max(feat.prob) + 0.05*max(feat.prob))
        feat.prob <- 1 - feat.prob
      }
    } else {
      if (grepl('_', type)){
        type <- strsplit(type, split='_')[[1]]
      }
      df.select <- df.select[(df.select$type %in% type),]
      el.feat.names <- df.select$idx
      feat.prob <- rep(1, length(el.feat.names))
    }
  }

  return(list(names=el.feat.names, prob=feat.prob))
}
