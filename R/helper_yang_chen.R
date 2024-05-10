#!/usr/bin/Rscript
### SIMBA
### Simulation of Microbiome data with Biological Accuracy
### Morgan Essex - MDC Berlin
### Jakob Wirbel - EMBL Heidelberg
### 2020 GNU GPL 3.0

# Helper functions to simulate data according to Yang and Chen
# Based on the SimulateMSeq function in the GUniFrac package

#' # wrapper for the Yang and Chen 
#' @keywords internal
simulate.SimMSeq <- function(feat, meta, sim.out, sim.params){
  
  test.package("GUniFrac")
  
  ab.scale <- sim.params$ab.scale
  prop.markers <- sim.params$prop.markers
  class.balance <- sim.params$class.balance
  repeats <- sim.params$repeats
  
  no.marker.feat <- round(prop.markers * nrow(feat))
  el.feat.names <- rownames(feat)
  num.sample <- ncol(feat)
  
  libSizesOrig <- colSums(feat)
  
  pb <- progress_bar$new(total = length(ab.scale)*repeats)
  for (a in seq_along(ab.scale)){
    for (r in seq_len(repeats)){
      
      # create new dataset with their function
      sim.chen <- SimulateMSeq(feat,
                               nSam=ncol(feat), nOTU=nrow(feat),
                               diff.otu.pct = prop.markers, 
                               diff.otu.mode = 'mix', 
                               diff.otu.direct = 'balanced',
                               covariate.type = 'binary', 
                               grp.ratio = 0.5/class.balance,
                               covariate.eff.mean = a,
                               depth.mu = mean(libSizesOrig))
      # get the info out of sim.chen
      sim.feat <- sim.chen$otu.tab.sim[rownames(feat),]
      marker.idx <- sim.chen$otu.names[which(sim.chen$diff.otu.ind)]
      label <- sim.chen$covariate[,1]
      label <- label[colnames(sim.feat)]
      label[label==0] <- -1
      
      # save data in H5-file
      h5.subdir <- paste0('ab', a, '_rep', r)
      stopifnot(h5createGroup(sim.out, h5.subdir))
      h5write(label, sim.out, paste0(h5.subdir, '/labels'))
      h5write(sim.feat, sim.out, paste0(h5.subdir, '/features'))
      h5write(marker.idx, sim.out, paste0(h5.subdir, '/marker_idx'))
      h5write(colnames(sim.feat), sim.out,
              paste(h5.subdir, '/sample_names', sep=''))
      h5write(rownames(sim.feat), sim.out,
              paste(h5.subdir, '/feature_names', sep=''))
      pb$tick()
    }
  }
}
