#!/usr/bin/Rscript
### SIMBA
### Simulation of Microbiome data with Biological Accuracy
### Morgan Essex - MDC Berlin
### Jakob Wirbel - EMBL Heidelberg
### 2020 GNU GPL 3.0

# Helper functions to simulate data according to Assefa, Vandesompele & Thas
# Based on the SPsimSeq function in the package with the same name

#' # wrapper for the Yang and Chen 
#' @keywords internal
simulate.SPsimSeq <- function(feat, meta, sim.out, sim.params){
  
  test.package("SPsimSeq")
  browser()
  
  ab.scale <- sim.params$ab.scale
  prop.markers <- sim.params$prop.markers
  class.balance <- sim.params$class.balance
  repeats <- sim.params$repeats
  
  no.marker.feat <- round(prop.markers * nrow(feat))
  el.feat.names <- rownames(feat)
  num.sample <- ncol(feat)
  
  sim.assefa <- SPsimSeq(n.sim=repeats, 
                         s.data = feat, 
                         n.genes = nrow(feat),
                         group.config = 1,
                         pDE = prop.markers,
                         cand.DE.genes = list(null.genes=rownames(feat)),
                         result.format = 'list')

  for (x in seq_len(repeats)){
    sim.feat <- sim.assefa[[x]]$counts
    # save data in H5-file
    h5.subdir <- paste0('ab', 0, '_rep', x)
    stopifnot(h5createGroup(sim.out, h5.subdir))
    h5write(sim.feat, sim.out, paste0(h5.subdir, '/features'))
    h5write(colnames(sim.feat), sim.out,
            paste(h5.subdir, '/sample_names', sep=''))
    h5write(rownames(sim.feat), sim.out,
            paste(h5.subdir, '/feature_names', sep=''))
  }
    
}
