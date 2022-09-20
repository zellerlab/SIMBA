#!/usr/bin/Rscript
### SIMBA
### Simulation of Microbiome data with Biological Accuracy
### Morgan Essex - MDC Berlin
### Jakob Wirbel - EMBL Heidelberg
### 2020 GNU GPL 3.0

# Helper functions to check data and parameters


# check the simulation destination file
#' @keywords internal
check.simulation.location <- function(sim.out){
  c.sim <- class(sim.out)
  if (length(c.sim) !=1){
    stop("Parameter 'sim.location' should be a filename!")
  } else if (c.sim!='character'){
    stop("Parameter 'sim.location' should be a filename!")
  }
  if (length(sim.out) !=1){
    stop("Parameter 'sim.location' should have length 1!")
  }
  if (sim.out == ''){
    stop("Parameter 'sim.location' is empty!")
  }
  if (!endsWith(sim.out, '.h5')){
    stop("Parameter 'sim.location' should end in .h5 for rdhf5 files!")
  }
  # if a file exists in the same location already,
  # check that it is empty
  if (file.exists(sim.out)){
    info <- file.info(sim.out)
    if (info$size > 800){
      stop('Simulation file already exists! Exiting...')
    } else {
      success <- file.remove(sim.out)
      if (success){
        warning("Removed empty simulation file before creating a new one!")
      } else {
        stop("Could not remove old empty simulation file!")
      }
    }
  }
  # check that a rdhf5 file can be created there
  success <- h5createFile(sim.out)
  if (!success){
    stop('Simulation file could not be created!')
  }
}

# check input data
#' @keywords internal
check.original.data <- function(feat, meta, sim.type, sim.method){
  
  # check that feat and meta make sense
  # check the class of the features
  c.feat <- class(feat)
  if (length(c.feat) != 1){
    stop('Your feature table has several classes (should be only one)!')
  }
  if (!c.feat %in% c('list', 'matrix', 'data.frame')){
    stop('Your feature input should be either a matrix, a data.frame',
         ', or a list!')
  }
  if (c.feat == 'list'){
    indiv.classes <- vapply(feat, class, FUN.VALUE = character(1))
    if (any(!indiv.classes %in% c('matrix', 'data.frame'))){
      stop('Your feature input should be a list of matrices or data.frames!')
    }
  }
  # check the class of meta
  c.meta <- class(meta)
  if (length(c.meta) != 1){
    stop('Your metadata information has several classes (should be only one)!')
  }
  if (!c.meta %in% c('list', 'data.frame')){
    stop('Your metadata input should be either a data.frame or a list!')
  }
  if (c.meta == 'list'){
    indiv.classes <- vapply(meta, class, FUN.VALUE = character(1))
    if (any(indiv.classes != 'data.frame')){
      stop('Your metadata input should be a list of data.frames!')
    }
  }
  # check that meta and feat agree and create original.data list
  original.data <- list()
  if (c.feat == 'list'){
    if (c.meta != 'list'){
      stop("Feature and metadata input should agree (either both lists or not)")
    } else {
      if (length(feat) != length(meta)){
        stop("Both feature and metadata list should be of the same length!")
      }
      for (i in seq_along(feat)){
        original.data[[i]] <- list()
        original.data[[i]][['feat']] <- feat[[i]]
        original.data[[i]][['meta']] <- meta[[i]]
      }
    }
  } else {
    if (c.meta=='list'){
      stop("Feature and metadata input should agree (either both lists or not)")
    } else {
      original.data[[1]] <- list()
      original.data[[1]][['feat']] <- feat
      original.data[[1]][['meta']] <- meta
    }
  }
  
  # validate original data
  data.list <- validate.original.data(original.data, sim.method)
  feat <- data.list$feat
  meta <- data.list$meta
  # remove repeated measurements/individuals with single samples
  if (sim.type == 'cross-section'){
    message("++ Removing repeated measurements for individuals")
    # if cross section, remove repeated measurements for indiviuals
    red.rownames <- vapply(unique(meta$Individual_ID), FUN = function(x){
      temp <- meta[meta$Individual_ID==x,,drop=FALSE]
      id <- sample(which(temp$Timepoint == min(temp$Timepoint)), size=1)
      rownames(temp)[id]
    }, FUN.VALUE = character(1))
    if (length(red.rownames) < 50){
      stop("Dataset contains a very limited number of individuals! (n=",
           length(red.rownames), ")")
    } else {
      meta.red <- meta[red.rownames,]
      feat.red <- feat[,red.rownames]
    }
    # sample to same number per study
    if (length(unique(meta.red$Study)) > 1){
      x <- table(meta.red$Study)
      message("++ The distribution of samples across studies is:")
      for (i in seq_along(x)){
        message('\t', names(x)[i], '\t', x[i])
      }
      message("++ Samples will be randomly selected to be balanced ",
              "across studies")
      m <- min(x)
      # subsample
      f.idx <- c()
      for (i in names(x)){
        f.idx <- c(f.idx, sample(which(meta.red$Study==i), size = m))
      }
      meta.red <- meta.red[f.idx,]
      feat.red <- feat.red[,rownames(meta.red)]
    }
    
  } else if (sim.type == 'time-course'){
    message("++ Removing individuals with insufficient number of measurements")
    # if time-series, check that there are enough samples with repeated
    # measuremnts
    # TODO: number of timepoints to be simulated? (currently, we look at 3)
    # more design options?
    indiv.count <- vapply(unique(meta$Individual_ID), FUN=function(x){
      nrow(meta[meta$Individual_ID==x,])
    }, FUN.VALUE = integer(1))
    indiv.count.red <- indiv.count[which(indiv.count >= 3)]
    message("++ There are ", length(indiv.count.red),
            " individuals with at least 3 data points")
    if (length(indiv.count.red) < 35){
      stop("Dataset contains a very limited number of individuals! (n=",
           length(indiv.count.red), ")")
    } else {
      meta.red <- meta[meta$Individual_ID %in% names(indiv.count.red),]
      feat.red <- feat[,rownames(meta.red)]
    }
  }
  
  return(list(feat=feat.red, meta=meta.red))
}

# validate original data
#' @keywords internal
validate.original.data <- function(d.list, sim.method){
  
  feat.final <- list()
  meta.final <- list()
  colnames.meta <- list()
  feat.names <- list()
  for (i in seq_along(d.list)){
    feat.temp <- as.matrix(d.list[[i]]$feat)
    meta.temp <- d.list[[i]]$meta
    # check sample name overlap
    # do we need a cutoff?
    if (is.null(rownames(feat.temp)) | is.null(colnames(feat.temp))){
      stop("All features need to have rownames and column names!")
    }
    if (is.null(rownames(meta.temp)) | is.null(rownames(meta.temp))){
      stop("All metadata need to have rownames and column names!")
    }
    
    # browser()
    # probably need some refactoring for e.g. sim.method params
    # these are good checks, but if data preprocessed they're annoying
    if (sim.method != 'pass') {
      # check that features are count tables
      if (typeof(feat.temp) != 'integer'){
        stop("Features for dataset No. ", i, " contains non-integer values!")
      } else if (any(feat.temp < 0)){
        stop("Features for dataset No. ", i, " contains negative entries!")
      }
      # filter
      removed.samples <- sum(colSums(feat.temp) < 1000)
      message("+ Removing ", removed.samples, 
              ' samples with less than 1000 counts!')
      feat.temp <- feat.temp[,colSums(feat.temp) > 1000, drop=FALSE] 
    }
    
    # check that metadata contains needed entries
    if (any(!c('Individual_ID', "Timepoint") %in% colnames(meta.temp))){
      stop("Metadata for dataset No. ", i,
           " should contain the columns:\n\t",
           paste(setdiff(c('Individual_ID', "Timepoint"),
                         colnames(meta.temp)), collapse=', '))
    } 
    if (!("Study" %in% colnames(meta.temp))) {
      meta.temp$Study <- paste0("Study_", i)
    }
    
    # check overlap
    colnames(feat.temp) <- make.names(colnames(feat.temp))
    rownames(meta.temp) <- make.names(rownames(meta.temp))
    # rownames(meta.temp) <- make.names(meta.temp$Sample_ID)
    x.f <- colnames(feat.temp)
    x.m <- rownames(meta.temp)
    all.samples <- intersect(x.f, x.m)
    if (length(all.samples) < 50){
      if (length(all.samples) == 0){
        stop("There is no overlap between samples in feat and meta",
             " for dataset No. ", i)
      }
      stop("Less than 100 samples for the dataset No. ", i)
    }
    
    feat.temp <- feat.temp[,all.samples]
    meta.temp <- meta.temp[all.samples,]
    colnames.meta[[i]] <- colnames(meta.temp)
    feat.names[[i]] <- rownames(feat.temp)
    
    rownames(meta.temp) <- paste0("Sample_", i, "_", seq_len(nrow(meta.temp)))
    colnames(feat.temp) <- paste0("Sample_", i, "_", seq_len(nrow(meta.temp)))
    
    # adjust the list
    feat.final[[i]] <- feat.temp
    meta.final[[i]] <- meta.temp
  }
  
  # test feature names
  if (length(feat.names) == 1){
    test <- TRUE
  } else if (length(feat.names) == 2){
    test <- all.equal(feat.names[[1]], feat.names[[2]])
  } else {
    test <- TRUE
    for (i in seq_len(length(feat.names))){
      test.sub <- all.equal(feat.names[[1]], feat.names[[i]])
      if (length(test.sub) > 1 | !is.logical(test.sub)){
        test <- FALSE
      }
    }
  }
  if (!test){
    stop("Feature names do not match across datasets!")
  }
  
  # combine everything
  # features
  feat.all <- do.call(cbind, feat.final)
  # metadata
  colnames.final <- reduce(colnames.meta, intersect)
  for (j in seq_along(meta.final)){
    meta.final[[j]] <- meta.final[[j]][,colnames.final]
  }
  meta.all <- do.call(rbind, meta.final)
  meta.all <- clean.meta.data(meta.all)
  feat.all <- feat.all[,rownames(meta.all)]
  rownames(feat.all) <- paste0('bact_otu_', seq_len(nrow(feat.all)))
  return(list(feat=feat.all, meta=meta.all))
}

# clean metadata
#' @keywords internal
clean.meta.data <- function(df.meta){
  must.haves <- c('Individual_ID', 'Timepoint', 'Study')
  
  metadata.categories <- setdiff(colnames(df.meta), must.haves)
  if (length(metadata.categories) == 0){
    return(df.meta)
  } else {
    
    # factorize
    for (x in metadata.categories){
      vec <- df.meta[,x]
      if (toupper(x) == 'BMI') {
        if (is.factor(vec)){
          lvl <- levels(vec)
          if (length(setdiff(lvl, c('Normal', 'Obese',
                                    'Overweight', 'Underweight'))) == 0){
            next()
          }
        }
        message("++ ", x, " is interpreted as body mass index:\n",
                "  will be converted to underweight/normal/overweight/obese")
        vec <- as.numeric(as.character(vec))
        temp <- vapply(vec, FUN=function(x) {
          if (is.na(x)) {return(NA_character_)}
          else if (x < 18.5) {return("Underweight")}
          else if ((x >= 18.5) & (x <= 24.9)) {return("Normal")}
          else if ((x > 24.9) & (x <= 29.9)) {return("Overweight")}
          else if (x > 29.9) {return("Obese")}},
          FUN.VALUE = character(1), USE.NAMES = TRUE)
        vec <- as.factor(temp)
      } else if (mode(vec)=='numeric' & length(unique(vec)) > 5){
        message("++ ", x, " will be converted to a factor via quartiles")
        quart <- quantile(vec, probs = seq(0, 1, 0.25), na.rm = TRUE)
        temp <- cut(vec, unique(quart), include.lowest = TRUE)
        vec <- factor(temp, labels = seq_along(levels(temp)))
      } else {
        vec <- as.factor(vec)
      }
      df.meta[,x] <- vec
    }
    
    # check
    n.levels <- vapply(metadata.categories, FUN=function(x){
      length(levels(df.meta[[x]]))}, FUN.VALUE = integer(1))
    remove.single <- which(n.levels==1)
    if (length(remove.single) > 0){
      message('++ meta-variables:\n\t',
              paste(names(remove.single), collapse = ', '),
              "\n++ have only a single value across all samples",
              " and will be removed")
      df.meta <- df.meta[,-which(colnames(df.meta) %in% names(remove.single))]
    }
    remove.many <- which(n.levels > 0.6*nrow(df.meta))
    if (length(remove.many) > 0){
      message('++ meta-variables:\n\t',
              paste(names(remove.many), collapse = ', '),
              "\n++ have too many values across samples",
              " and will be removed")
      df.meta <- df.meta[,-which(colnames(df.meta) %in% names(remove.many))]
    }
    return(df.meta)
  }
}

# check simulation parameters
#' @keywords internal
check.filtering.parameters <- function(filt.params){
  
  if (is.null(filt.params)){
    message("+ No filtering will be performed!")
  }
  x <- names(filt.params)
  unused.params <- setdiff(x, c('prev.cutoff', 'ab.cutoff', 'log.n0'))
  if (length(unused.params) != 0){
    warning("+ Filtering parameters: ", paste0(unused.params, collapse = ','),
            " will not be used for filtering!")
  }
  prev.cutoff <- filt.params$prev.cutoff
  ab.cutoff <- filt.params$ab.cutoff
  log.n0 <- filt.params$log.n0
  # check filtering parameters
  # check that prevalence,
  if (!is.null(prev.cutoff)){
    if (mode(prev.cutoff)!='numeric'){
      stop("Prevalence cutoff is non-numeric!")
    } else if (prev.cutoff > 0.5){
      stop("Prevalence cutoff is too aggressive!")
    } else if (prev.cutoff < 0){
      warning("Prevalence cutoff is negative and will be ignored")
      prev.cutoff <- 0
    }
  }
  # abundance cutoff,
  if (!is.null(ab.cutoff)){
    if (mode(ab.cutoff)!='numeric'){
      stop("Abundance cutoff is non-numeric!")
    } else if (ab.cutoff > 0.01){
      stop("Abundance cutoff is too aggressive!")
    } else if (ab.cutoff < 0){
      warning("Abundance cutoff is negative and will be ignored")
      ab.cutoff <- 0
    }
  }
  # log.n0
  if (is.null(log.n0)){
    warning("Parameter 'log.n0' is needed and will be set to 1e-05!")
    log.n0 <- 1e-05
  } else {
    if (mode(log.n0) != 'numeric'){
      stop("Paramter 'log.n0' is non-numeric!")
    } else if (log.n0 > 0.5){
      stop("Paramter 'log.n0' is unreasonably large (", log.n0, ")!")
    } else if (log.n0 < 0){
      stop("Paramter 'log.n0' is negative (", log.n0, ")!")
    }
  }
  return(list(ab.cutoff=ab.cutoff, prev.cutoff=prev.cutoff, log.n0=log.n0))
}

# check confounders for microbiome associations
#' @keywords internal
confounder.check <- function(feat, meta){
  # relative abundances
  feat.rel <- prop.table(feat, 2)
  # get meta-variables to test
  voi <- setdiff(colnames(meta), c('Individual_ID', 'Timepoint', 'Study'))
  df.plot.all <- data.frame()
  for (v in voi){
    temp <- meta[[v]]
    names(temp) <- rownames(meta)
    if (any(is.na(temp))){
      temp <- temp[-which(is.na(temp))]
    }
    var.batch <- vapply(rownames(feat), FUN=function(x){
      x <- feat[x,names(temp)]
      x <- rank(x)/length(x)
      ss.tot <- sum((x - mean(x))^2)/length(x)
      ss.o.i <- sum(vapply(levels(temp), function(s){
        sum((x[temp==s] - mean(x[temp==s]))^2)
      }, FUN.VALUE = double(1)))/length(x)
      return(1-ss.o.i/ss.tot)
    }, FUN.VALUE = double(1))
    p.vals <- vapply(rownames(feat), FUN=function(x){
      x <- feat[x,names(temp)]
      t <- kruskal.test(x=x, g=temp, exact=FALSE)
      t$p.value
    }, FUN.VALUE = double(1))
    df.plot <- data.frame(feature=rownames(feat),
                          p.val=p.adjust(p.vals, method='fdr'),
                          effect.size=var.batch,
                          variable=v)
    rownames(df.plot) <- NULL
    df.plot.all <- rbind(df.plot.all, df.plot)
  }
  
  # create plot
  ggplot(df.plot.all, aes(x=effect.size, y=-log10(p.val))) +
    geom_point() +
    facet_grid(~variable)
}

# check simulation parameters
#' @keywords internal
check.simulation.parameters <- function(sim.method,
                                        sim.type,
                                        sim.params,
                                        meta){
  
  if (mode(sim.method) == 'character'){
    allowed.sim.methods <- c('resampling', 'McMurdie&Holmes',
                             "Weiss", "negbin", "betabin",
                             "dirmult", "sparseDOSSA", "pass")
    if (!sim.method %in% allowed.sim.methods){
      stop("Parameter 'sim.method' must be one of those: ",
           paste(allowed.sim.methods, collapse = ', '))
    }
    if (sim.method=='resampling'){
      needed.params <- c('class.balance', 'prop.markers',
                         'ab.scale', 'prev.scale', 'repeats',
                         'feature.type')
    } else if (sim.method=='McMurdie&Holmes'){
      needed.params <- c('class.balance', 'prop.markers',
                         'ab.scale',  'repeats')
    } else if (sim.method=='Weiss'){
      needed.params <- c('class.balance', 'prop.markers',
                         'ab.scale',  'repeats')
    } else if (sim.method=='negbin'){
      needed.params <- c('ab.scale', 'prop.markers', 'class.balance',
                         'repeats', 'correlation')
    } else if (sim.method=='dirmult'){
      needed.params <-  c('class.balance', 'prop.markers',
                          'ab.scale',  'repeats')
    } else if (sim.method=='betabin'){
      needed.params <-  c('class.balance', 'prop.markers',
                          'ab.scale',  'repeats')
    } else if (sim.method=='sparseDOSSA'){
      needed.params <-  c('class.balance', 'prop.markers',
                          'ab.scale',  'repeats')
    } else if (sim.method=='pass'){
      needed.params <- NULL
    }
  } else {
    if (mode(sim.method) == 'function'){
      message("++ Using custom simulation function... all bets are off!")
    }
    needed.parameters <- names(formals(sim.method))
  }
  
  # check general simulation parameters
  allowed.sim.types <- c('cross-section', 'time-course')
  if (!sim.type %in% allowed.sim.types){
    stop("Parameter 'sim.type' must be one of those: ",
         paste(allowed.sim.types, collapse = ', '))
  }
  
  # other parameters
  unused.params <- setdiff(names(sim.params), needed.params)
  opt <- c('conf', 'conf.params', 'balanced')
  if (sim.method=='resampling'){
    unused.params <- setdiff(names(sim.params), c(needed.params, opt))
  }
  missing.params <- setdiff(needed.params, names(sim.params))
  
  if (length(missing.params) > 0){
    stop("Not all needed parameters have been supplied, missing are:\n",
         paste(missing.params, collapse=', '))
  }
  if (length(unused.params) > 0){
    message("++ Not all supplied parameters are needed, superfluous are:\n\t",
            paste(unused.params, collapse = ', '))
  }
  cleaned.params <- list()
  for (x in needed.params){
    cleaned.params[[x]] <- sim.params[[x]]
  }
  if (sim.method=='resampling'){
    if ('conf' %in% names(sim.params)){
      cleaned.params[['conf']] <- sim.params[['conf']]
    }
    if ('conf.params' %in% names(sim.params)){
      cleaned.params[['conf.params']] <- sim.params[['conf.params']]
    }
    if ('balanced' %in% names(sim.params)){
      cleaned.params[['balanced']] <- sim.params[['balanced']]
    }
  }
  
  # check most common parameters
  cleaned.params <- check.common.simulation.parameters(cleaned.params, meta,
                                                       sim.method)
  
  return(cleaned.params)
}

# check simulation parameters
#' @keywords internal
check.common.simulation.parameters <- function(sim.params, meta, sim.method){
  n <- names(sim.params)
  
  if ("prop.markers" %in% n){
    prop.markers <- sim.params$prop.markers
    # more simulation stuff
    # check that prop and class balances are in reasonable ranges
    if (mode(prop.markers) != 'numeric'){
      stop("Parameter 'prop.markers' should be numeric")
    } else if (prop.markers > 0.7){
      stop("Parameter 'prop.markers' seems too high (", prop.markers, ')')
    } else if (prop.markers < 0.01){
      stop("Parameter 'prop.markers' seems too low (", prop.markers, ')')
    }
  }
  if ("class.balance" %in% n){
    class.balance <- sim.params$class.balance
    if (mode(class.balance) != 'numeric'){
      stop("Parameter 'class.balance' should be numeric")
    } else if (class.balance > 0.85){
      stop("Parameter 'class.balance' seems too high (", class.balance, ")")
    } else if (class.balance < 0.45){
      stop("Parameter 'class.balance' seems too low (", class.balance, ")")
    }
  }
  if ("ab.scale" %in% n){
    ab.scale <- sim.params$ab.scale
    # check that ab.scale and prev.scale make sense
    if (mode(ab.scale)!='numeric'){
      stop("Abundance scaling vector should be numeric!")
    } else if (any(ab.scale < 1)){
      stop("Abundance scaling vector should not contain entries smaller than 1")
    }
  }
  if ("prev.scale" %in% n){
    prev.scale <- sim.params$prev.scale
    if (mode(prev.scale)!='numeric'){
      stop("Prevalence scaling vector should be numeric!")
    } else if (any(prev.scale < 0)){
      stop("Prevalence scaling vector should not contain entries smaller than 0")
    } else if (any(prev.scale > 0.99)){
      stop("Prevalence scaling vector should not contain entries bigger than 1")
    }
  }
  if ('repeats' %in% n){
    repeats <- sim.params$repeats
    # check that repeats make sense
    if (mode(repeats)!='numeric'){
      message("++ Repeats are non-numeric, will be set to 100")
      repeats <- 100
    } else if (repeats < 0){
      message("++ Repeats are negative, will be set to 100")
      repeats <- 100
    } else if (repeats > 300){
      stop("Unreasonably high number of repeats")
    } else if (repeats < 5){
      message("++ Unreasonably low number of repeats, will be set to 10")
      repeats <- 10
    }
  }
  if ("feature.type" %in% n){
    feature.type <- sim.params$feature.type
    if (mode(feature.type) != 'character'){
      stop("Parameter 'feature.type' must of type 'character'!")
    }
    if (length(feature.type) != 1){
      stop("Parameter 'feature.type' must of length 1")
    }
    allowed.values <- c('low', 'middle', 'high', 'low_middle',
                        'low_high', 'middle_high', 'all', 'abundance',
                        'inverse-abundance')
    if (!feature.type %in% allowed.values){
      stop("Parameter 'feature.type' must be one of: \n\t",
           paste(allowed.values, collapse = ', '))
    }
  }
  if ("correlation" %in% n){
    correlation <- sim.params$correlation
    if (mode(correlation) != 'logical'){
      stop("Parameter 'correlation' must be of type 'logical'!")
    }
  }
  
  # confounder stuff
  if (sim.method=='resampling'){
    # confounder info
    if (!"conf" %in% n){
      sim.params$conf <- 'None'
    } else if (is.null(sim.params$conf)) {
      sim.params$conf <- 'None'
    } else {
      conf <- sim.params$conf
      if (conf != 'None'){
        if (mode(conf) != 'character'){
          stop("Parameter 'conf' should be of type 'character'")
        } else if (length(conf) > 1){
          stop("Parameter 'conf' should be of length 1!")
        }
        non.allowed.values <- c('Sample_ID', 'Individual_ID',
                                'Timepoint')
        allowed <- c(colnames(meta), 'artificial', 'global', 'batch')
        allowed <- setdiff(allowed, non.allowed.values)
        if (!conf %in% allowed){
          stop("Parameter 'conf' must be present in metadata!")
        }
        # check that the confounder is binary? -- For now
        if (conf=='batch'){
          x <- length(unique(meta$Study))
          if (x != 2){
            stop("Confounding variable needs to be binary (for now)!\n\t",
                 "Number of unique values: ", x)
          }
        }
        if (conf %in% colnames(meta)){
          x <- length(unique(meta[[conf]]))
          if (x != 2){
            stop("Confounding variable needs to be binary (for now)!\n\t",
                 "Number of unique values: ", x)
          }
        }
        
        # confounder parameter list
        if (!"conf.params" %in% n){
          stop("Parameter 'conf.params' needed but not supplied!")
        } else {
          conf.params <- sim.params$conf.params
          if (mode(conf.params) != 'list'){
            stop("Parameter 'conf.params' should be a list of parameters!")
          }
          m <- names(conf.params)
          if (conf %in% c('global', 'batch', 'Study')){
            needed <- c('bias', 'prop')
            opt <- c()
          } else if (conf == 'artificial'){
            needed <- c('bias', 'prop')
            opt <- c('feat.prop', 'feat.effect')
          } else {
            needed <- 'bias'
            opt <- c()
          }
          
          # check that all needed parameters are given
          if (!all(needed %in% m)){
            stop("Not all parameters for confounder implantations",
                 " are supplied, missing are: ", paste(setdiff(needed, m),
                                                       collapse = ', '))
          }
          # check bias
          if (mode(sim.params$conf.params$bias) != 'numeric'){
            stop("Parameter 'conf.params$bias' should be numeric!")
          }
          if (length(sim.params$conf.params$bias) > 1){
            stop("Parameter 'conf.params$bias' should be of length 1!")
          }
          if (sim.params$conf.params$bias > 1 |
              sim.params$conf.params$bias < 0.5){
            stop("Parameter 'conf.params$bias' should be between 1 and 0.5!",
                 "\n\t(Supplied value: ", sim.params$conf.params$bias,")")
          }
          
          # check prop
          if ('prop' %in% needed){
            if (mode(sim.params$conf.params$prop) != 'numeric'){
              stop("Parameter 'conf.params$prop' should be numeric!")
            }
            if (length(sim.params$conf.params$prop) > 1){
              stop("Parameter 'conf.params$prop' should be of length 1!")
            }
            if (sim.params$conf.params$prop > 0.5 |
                sim.params$conf.params$prop < 0){
              stop("Parameter 'conf.params$prop' should be between 0 and 0.5!",
                   "\n\t(Supplied value: ", sim.params$conf.params$prop,")")
            }
          }
          
          if (length(opt) > 0){
            # check feat.prop
            if ('feat.prop' %in% names(sim.params$conf.params)){
              if (mode(sim.params$conf.params$feat.prop) != 'numeric'){
                stop("Parameter 'conf.params$feat.prop' should be numeric!")
              }
              if (length(sim.params$conf.params$feat.prop) > 1){
                stop("Parameter 'conf.params$feat.prop' should be of length 1!")
              }
              if (sim.params$conf.params$feat.prop > 0.5 |
                  sim.params$conf.params$feat.prop < 0.05){
                stop("Parameter 'conf.params$feat.prop' seems extreme.",
                     "\n\t(Supplied value: ",
                     sim.params$conf.params$feat.prop,")")
              }
            } else {
              sim.params$conf.params$feat.prop <- sim.params$prop.markers
            }
            x <- sim.params$prop.markers + sim.params$conf.params$feat.prop
            if (x > 0.8){
              stop("The combination of markers and confounder",
                   " markers would impact ", x, "% of the features!")
            }
            # check feat.effect
            if ('feat.effect' %in% names(sim.params$conf.params)){
              if (mode(sim.params$conf.params$feat.effect) != 'numeric'){
                stop("Parameter 'conf.params$feat.effect' should be numeric!")
              }
              if (length(sim.params$conf.params$feat.effect) > 1){
                stop("Parameter 'conf.params$feat.effect' should be of length 1!")
              }
              if (sim.params$conf.params$feat.effect > 20 |
                  sim.params$conf.params$feat.effect < 2){
                stop("Parameter 'conf.params$feat.effect' should not be",
                     " higher than 20 or lower than 2",
                     "\n\t(Supplied value: ",
                     sim.params$conf.params$feat.effect,")")
              }
            } else {
              message("++++ Confounder markers will be scaled the same",
                      " as label markers!")
            }
          }
        }
      }
    }
  }
  return(sim.params)
}

# filter original data
#' @keywords internal
filter.data <- function(feat, meta, filt.params){
  feat.rel <- prop.table(feat, 2)
  ab <- rowMaxs(feat.rel)
  names(ab) <- rownames(feat.rel)
  prev <- rowMeans(feat.rel!=0)
  
  # get filt.params
  prev.cutoff <- filt.params$prev.cutoff
  ab.cutoff <- filt.params$ab.cutoff
  
  # filter
  if (!is.null(ab.cutoff) & !is.null(prev.cutoff)){
    ids <- intersect(names(which(prev > prev.cutoff)),
                     names(which(ab > ab.cutoff)))
    message("++ Performing prevalence and abundance filtering")
  } else if (!is.null(ab.cutoff)){
    ids <- names(which(ab > ab.cutoff))
    message("++ Performing abundance filtering")
  } else if (!is.null(prev.cutoff)){
    ids <- names(which(prev > prev.cutoff))
    message("++ Performing prevalence filtering")
  } else {
    ids <- rownames(feat.rel)
    message("++ Performing no filtering")
  }
  
  feat.final <- feat[ids,]
  # remove unmapped
  names.unmapped <- c("UNMAPPED", "-1", "X.1", "unmapped",
                      "Unclassified", "Unassigned",
                      "UNCLASSIFIED", "unclassified", "UNASSIGNED",
                      "unassigned")
  if (any(rownames(feat.final) %in% names.unmapped)){
    unm <- intersect(names.unmapped, rownames(feat.final))
    feat.final <- feat.final[-which(rownames(feat.final) == unm),]
    message("++ Removed unmapped reads")
  }
  
  message("++ After filtering, the dataset contains ", nrow(feat.final),
          " features and ", ncol(feat.final), " samples")
  return(feat.final)
}
