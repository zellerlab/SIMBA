#!/usr/bin/Rscript
### SIMBA
### Simulation of Microbiome data with Biological Accuracy
### Morgan Essex - MDC Berlin
### Jakob Wirbel - EMBL Heidelberg
### 2020 GNU GPL 3.0

# Helper functions to simulate data according to McMurdie&Holmes or Weiss et al.

#' # wrapper for the McMurdie&Holmes simulations
#' @keywords internal
simulate.MMH <- function(feat, meta, sim.out, sim.params){

  # get parameters
  ab.scale <- sim.params$ab.scale
  prop.markers <- sim.params$prop.markers
  class.balance <- sim.params$class.balance
  repeats <- sim.params$repeats

  test.package("phyloseq")
  no.marker.feat <- round(prop.markers * nrow(feat))
  el.feat.names <- rownames(feat)

  template <- phyloseq::otu_table(feat, taxa_are_rows = TRUE)

  lib.size <- phyloseq::sample_sums(template)
  n <- sample(lib.size, size = ncol(feat), replace = TRUE)
  message("++ Create simulation template from data")
  sim.feat.mmh <- microbesim(template, J=ncol(feat), n=n)
  message("++ Finished creating simulation template")

  num.sample <- ncol(feat)
  # actual implantation
  pb <- progress_bar$new(total = length(ab.scale)*repeats)
  for (a in seq_along(ab.scale)){


    for (r in seq_len(repeats)){
      # generate label
      label <- rep(-1, ncol(feat))
      s <- sample(num.sample,
                  size=round(num.sample*class.balance))
      names(label) <- colnames(feat)
      label[s] <- +1
      # actually simulate
      sim.feat <- simulate.markers.MMH(sim.feat.mmh,
                                       label,
                                       no.marker.feat,
                                       ab.scale[a])
      sim.feat <- t(sim.feat@.Data)

      # adjust some specialties from the MM&H method
      # column names
      colnames(sim.feat) <- names(label)
      # row names
      marker.idx <- grep('-TP$', rownames(sim.feat))
      rownames(sim.feat) <- gsub(rownames(sim.feat),
                                 pattern = '-TP$',
                                 replacement = '')
      marker.idx <- rownames(sim.feat)[marker.idx]

      # save data in H5-file
      h5.subdir <- paste0('ab', a, '_rep', r)
      stopifnot(h5createGroup(sim.out, h5.subdir))
      h5write(label, sim.out, paste0(h5.subdir, '/labels'))
      h5write(sim.feat, sim.out, paste0(h5.subdir, '/features'))
      h5write(marker.idx, sim.out, paste(h5.subdir, '/marker_idx', sep=''))
      h5write(colnames(sim.feat), sim.out,
              paste(h5.subdir, '/sample_names', sep=''))
      h5write(rownames(sim.feat), sim.out,
              paste(h5.subdir, '/feature_names', sep=''))
      pb$tick()
    }
  }
}

#' # wrapper for the Weiss simulations
#' @keywords internal
simulate.W <- function(feat, meta, sim.out, sim.params){

  # get parameters
  ab.scale <- sim.params$ab.scale
  prop.markers <- sim.params$prop.markers
  class.balance <- sim.params$class.balance
  repeats <- sim.params$repeats

  test.package("phyloseq")
  no.marker.feat <- round(prop.markers * nrow(feat))
  el.feat.names <- rownames(feat)

  template <- phyloseq::otu_table(feat, taxa_are_rows = TRUE)

  lib.size <- phyloseq::sample_sums(template)
  n <- sample(lib.size, size = ncol(feat), replace = TRUE)
  message("++ Create simulation template from data")
  sim.feat.mmh <- microbesim(template, J=ncol(feat), n=n)
  message("++ Finished creating simulation template")
  num.sample <- ncol(feat)
  browser()
  # actual implantation
  pb <- progress_bar$new(total = length(ab.scale)*repeats)
  for (a in seq_along(ab.scale)){


    for (r in seq_len(repeats)){

      # generate label
      label <- rep(-1, ncol(feat))
      s <- sample(num.sample,
                  size=round(num.sample*class.balance))
      names(label) <- colnames(feat)
      label[s] <- +1
      # actually simulate
      sim.feat <- simulate.markers.W(sim.feat.mmh,
                                     label,
                                     no.marker.feat,
                                     ab.scale[a])
      sim.feat <- sim.feat@.Data


      # adjust some specialties from the MM&H method carried over to
      # the method from Weiss et al.
      # column names
      colnames(sim.feat) <- names(sort(label))
      sim.feat <- sim.feat[,names(label)]
      # row names
      marker.idx <- grep('-TP$', rownames(sim.feat))
      rownames(sim.feat) <- gsub(rownames(sim.feat),
                                 pattern = '-TP$',
                                 replacement = '')
      marker.idx <- rownames(sim.feat)[marker.idx]

      # save data in H5-file
      h5.subdir <- paste0('ab', a, '_rep', r)
      stopifnot(h5createGroup(sim.out, h5.subdir))
      h5write(label, sim.out, paste0(h5.subdir, '/labels'))
      h5write(sim.feat, sim.out, paste0(h5.subdir, '/features'))
      h5write(marker.idx, sim.out, paste(h5.subdir, '/marker_idx', sep=''))
      h5write(colnames(sim.feat), sim.out,
              paste(h5.subdir, '/sample_names', sep=''))
      h5write(rownames(sim.feat), sim.out,
              paste(h5.subdir, '/feature_names', sep=''))
      pb$tick()
    }
  }
}

# ##############################################################################
# simulate features according to McMurdie&Holmes or Weiss et al
#   (both use the same simulation approach,
#    but different implantation of features)
microbesim <- function(template, J, n=1E4){
  # Generate `J` simulated microbiomes with `n` total reads each
  # (all the same, or n has length equal to the value of `J`),
  # with subsamples drawn from `template`.
  # simulated samples in downstream code.
  test.package("phyloseq")
  # call the proporitions vector `pi`, similar to nomenclature from DMN
  pi = phyloseq::taxa_sums(template)
  # n must be a scalar (recycled as the number of reads for every simulation)
  # or it can be vector of length equal to J,
  # which is the number of samples being simulated.
  if(length(J) != 1){ stop("Length of J should be 1.") }
  if(length(n) != 1 & length(n) != J){
    stop("n should be length 1, or length J.")
  }
  # Actually create the simulated abundance table
  # simat = mapply(function(i, x, sample.size){
  #   if(FALSE){print(i)} # i is a dummy iterator
  #   phyloseq:::rarefaction_subsample(x, sample.size)
  # }, i=1:J, sample.size=n, MoreArgs=list(x=pi), SIMPLIFY=TRUE)
  
  # ############
  # EDIT by Jakob
  # alternative for KEGG-type profiles?
  # the previous way to simulate fails because of overflow (as is feared in the
  # phyloseq source code... alternatively, we could use the vegan rarefaction 
  # for very similar results?)
  # # UPDATE
  # i checked with PCoA and the results are virtually indistinguishable, 
  # but the vegan version is way faster and more memory-friendly
  # # SECOND UPDATE
  # fails for KEGG profiles since the counts are too large
  while (TRUE){
    pi.2 <- pi
    mode(pi.2) <- "integer"
    if (any(is.na(pi.2))){
      message('+++ Reduce counts in the taxa_sums vector to prevent errors')
      pi <- pi/10
      n <- n/10
    } else {
      pi <- pi.2
      break
    }
  }
  message("+++ Create base simulations")
  simat <- vapply(n, FUN=function(x){vegan::rrarefy(pi, x)}, 
                  FUN.VALUE = double(length(pi)))
  # i checked with PCoA and the results are virtually indistinguishable, 
  # but the vegan version is way faster and more memory-friendly
  # ############
  simat = t(simat)
  # Add the OTU names to the OTU (column) indices
  colnames(simat) <- names(pi)
  # Add new simulated sample_names to the row (sample) indices
  rownames(simat) <- paste(1:nrow(simat), '_sim', sep="")
  # Put simulated abundances together with metadata as a phyloseq object
  OTU = phyloseq::otu_table(simat, taxa_are_rows=FALSE)

  # Return a phyloseq object
  return(phyloseq::phyloseq(OTU))
}

# ##############################################################################
# simulate markers according to McMurdie&Holmes
#' @keywords internal
simulate.markers.MMH <- function(physeq, label, nTP, effectsize){
  # Randomly sample from the available OTU names in physeq and
  # assign them as TP
  TPOTUs = sample(phyloseq::taxa_names(physeq), nTP, replace=FALSE)
  # Define the samples that will have effect applied
  effectsamples = which(label == 1)
  # Apply effect (multiply abundances by effectsize scalar)
  phyloseq::otu_table(physeq)[effectsamples, TPOTUs] <- effectsize *
    phyloseq::otu_table(physeq)[effectsamples, TPOTUs]
  # Rename these new "true positive" OTUs, with 'TP'
  wh.TP = phyloseq::taxa_names(physeq) %in% TPOTUs
  newname = paste0(phyloseq::taxa_names(physeq)[wh.TP], "-TP")
  colnames(physeq@.Data)[wh.TP] <- newname
  rownames(physeq@.Data)[which(label == 1)] <-
    paste0(rownames(physeq@.Data)[which(label == 1)], '-pos')
  return(physeq)
}

# ##############################################################################
# simulate makers according to Weiss et al.
#' @keywords internal
simulate.markers.W <- function(physeq, label, nTP, effectsize){
  # transpose otu table for the Weiss function to work
  # browser()
  phyloseq::otu_table(physeq) <- t(phyloseq::otu_table(physeq))
  # also set the taxa_are_rows variable?
  physeq@taxa_are_rows <- TRUE
  ## Randomly sample from the available OTU names in physeq and
  # assign them as TP
  TPOTUs = sample(phyloseq::taxa_names(physeq), nTP, replace = FALSE)
  physeq.2 = physeq
  # Apply effect (multiply abundances by effectsize scalar)
  #### DIFF: round the otu count to integer after adding effect size
  phyloseq::otu_table(physeq)[TPOTUs, ] <-
    apply(effectsize * phyloseq::otu_table(physeq)[TPOTUs, ], 2, round)

  # Rarely a simulation has a weird value and fails.  Catch these with `try`,
  # and repeat the simulation call if error (it will be a new seed)
  tryAgain = TRUE
  infiniteloopcounter = 1
  while (tryAgain & infiniteloopcounter < 5) {
    lib.size.1 <- phyloseq::sample_sums(physeq)
    n1 <- sample(lib.size.1, size = length(which(label==1)), replace = TRUE)
    lib.size.1 <- phyloseq::sample_sums(physeq.2)
    n2 <- sample(lib.size.1, size = length(which(label==-1)), replace = TRUE)
    sim1 = microbesim(physeq, length(which(label==1)), n1)
    sim2 = microbesim(physeq.2, length(which(label==-1)), n2)
    if (is.null(sim1) | is.null(sim2) | is.null(n1) | is.null(n2) |
        inherits(sim1, "try-error") | inherits(sim2, "try-error")) {
      tryAgain = TRUE
      infiniteloopcounter = infiniteloopcounter + 1
    } else {
      tryAgain = FALSE
    }
  }
  if (infiniteloopcounter >= 5) {
    stop("Consistent error found during simulation. Need to investigate cause.")
  }
  # Merge the two simulated datasets together into one otu_table
  rownames(sim1) <- paste0(rownames(sim1), '-pos')
  sim = phyloseq::otu_table(cbind(t(sim1), t(sim2)), taxa_are_rows = TRUE)

  # Rename the new 'true positive' OTUs, with 'TP'
  wh.TP = phyloseq::taxa_names(sim) %in% TPOTUs
  newname = paste0(phyloseq::taxa_names(sim)[wh.TP], "-TP")
  rownames(sim)[wh.TP] <- newname
  return(sim)
}
