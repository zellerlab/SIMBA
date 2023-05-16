#!/usr/bin/Rscript
### SIMBA
### Simulation of Microbiome data with Biological Accuracy
### Morgan Essex - MDC Berlin
### Jakob Wirbel - EMBL Heidelberg
### 2020 GNU GPL 3.0

#' @keywords internal
run.test <- function(data, label, test, conf){

  if (test == 'wilcoxon'){
    p.val <- test.wilcoxon(data=data, label=label, conf=conf)
  } else if (test == 'fisher'){
    if (!is.null(conf)) {
      stop("Test '", test, "' is not a confounder-aware test!") }
    p.val <- test.fisher(data=data, label=label)
  } else if (test == 'DESeq2'){
    if (!is.null(conf)) {
      stop("Test '", test, "' is not a confounder-aware test!") }
    p.val <- test.DESeq(data=data, label=label)
  } else if (test == 'metagenomeSeq'){
    if (!is.null(conf)) {
      stop("Test '", test, "' is not a confounder-aware test!") }
    p.val <- test.MGS(data=data, label=label)
  } else if (test == 'metagenomeSeq2'){
    if (!is.null(conf)) {
      stop("Test '", test, "' is not a confounder-aware test!") }
    p.val <- test.MGS(data=data, label=label, type=2)
  } else if (test == 'mdc-FE'){
    if (is.null(conf)) {
      stop("Test '", test, "' is strictly a confounder-aware test!") }
    p.val <- test.metadeconfoundR(data=data, label=label, type=1, conf=conf)
  } else if (test == 'mdc-RE') {
    if (is.null(conf)) {
      stop("Test '", test, "' is strictly a confounder-aware test!") }
    p.val <- test.metadeconfoundR(data=data, label=label, type=2, conf=conf)
  } else if (test == 'edgeR'){
    p.val <- test.edgeR(data=data, label=label)
  } else if (test == 'ZIBSeq'){
    p.val <- test.ZIBseq(data=data, label=label, conf=conf)
  } else if (test == 'ZIBSeq-sqrt'){
    p.val <- test.ZIBseq_sqrt(data=data, label=label, conf=conf)
  } else if (test == 'ANCOM'){
    p.val <- test.ANCOM.2(data=data, label=label, conf=conf)
  } else if (test == 'ANCOM_old'){
    p.val <- test.ANCOM(data=data, label=label, conf=conf)
  } else if (test == 'ANCOMBC'){
    p.val <- test.ANCOMBC(data=data, label=label, conf=conf)
  } else if (test=='corncob'){
    p.val <- test.corncob(data=data, label=label, conf=conf)
  } else if (test=='limma'){
    p.val <- test.via.limma(data=data, label=label, conf=conf)
  } else if (test=='lm'){
    p.val <- test.lm(data=data, label=label, conf=conf)
  } else if (test=='lme'){
    p.val <- test.lme(data=data, label=label, conf=conf)
  } else if (test=='ALDEx2') {
    p.val <- test.aldex2(data=data, label=label, conf=conf)
  } else if (test=='KS') {
    p.val <- test.ks(data=data, label=label, conf=conf)
  } else if (test=='scde') {
    p.val <- test.scde(data=data, label=label, conf=conf)
  } else if (test=='MAST') {
    p.val <- test.MAST(data=data, label=label, conf=conf)
  } else if (test=='mixMC') {
    p.val <- test.mixMC(data=data, label=label, conf=conf)
  } else if (test=='distinct') {
    p.val <- test.distinct(data=data, label=label, conf=conf)
  } else if (test=='songbird') {
    p.val <- test.songbird(data=data, label=label, conf=conf)
  } else if (test=='ZINQ') {
    p.val <- test.ZINQ(data=data, label=label, conf=conf)
  } else if (test=='gFC') {
    p.val <- test.gFC(data=data, label=label, conf=conf)
  } else {
    stop("Test not found!") }
  return(p.val)
}

#' @keywords internal
test.wilcoxon <- function(data, label, conf){

  stopifnot(ncol(data) == length(label))
  stopifnot(all(colnames(data) %in% names(label)))
  stopifnot(all(sort(unique(label)) == c(-1, 1)))
  p.val = rep(1, nrow(data))
  names(p.val) <- rownames(data)
  if (!is.null(conf)){
    test.package("coin")
    stopifnot(all(names(label) %in% rownames(conf)))
    if (ncol(conf) > 1){
      stop("Blocked wilcoxon test only works with a single confounder (atm)!")
    }
    data <- data[,names(label)]
    conf <- conf[colnames(data),,drop=FALSE]
    # apply blocked Wilcoxon test

    for (f in rownames(data)){
      df.temp <- data.frame(feat=as.numeric(data[f,]),
                            l=as.factor(label),
                            conf=as.factor(conf[,1]))
      t <- wilcox_test(feat~l|conf, data=df.temp)
      p.val[f] <- pvalue(t)
    }
    # browser()
  } else {
    data <- data[,names(label)]
    # apply Wilcoxon test
    for (f in rownames(data)) {
      t = wilcox.test(as.numeric(data[f,label==-1]),
                      as.numeric(data[f,label==+1]), exact=FALSE)
      p.val[f] = t$p.value
    }
  }

  # Nan p.val can occurr for constant samples
  p.val[!is.finite(p.val)] = 1.0
  return(p.val)
}

#' @keywords internal
test.fisher <- function(data, label){

  stopifnot(ncol(data) == length(label))
  stopifnot(all(colnames(data) %in% names(label)))
  data <- data[,names(label)]
  stopifnot(all(sort(unique(label)) == c(-1, 1)))

  p.val = rep(1, nrow(data))
  for (f in 1:nrow(data)) {
    # test for enrichment of non-zero abundances
    p = ifelse(data[f,] > 0, 1, 0)
    counts = matrix(0, 2, 2)
    t = table(p[label==-1])
    counts[1, as.numeric(names(t))+1] = t
    t = table(p[label==+1])
    counts[2, as.numeric(names(t))+1] = t
    t = fisher.test(counts)
    p.val[f] = t$p.value
  }
  return(p.val)
}

#' @keywords internal
test.DESeq <- function(data, label){
  # based on vignette in phyloseq
  # https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html
  test.package('DESeq2')
  test.package('phyloseq')
  stopifnot(ncol(data) == length(label))
  stopifnot(all(sort(unique(label)) == c(-1, 1)))

  # convert to phyloseq object
  temp_otu <- otu_table(data, taxa_are_rows = TRUE)
  temp_sample <- as.matrix(label)
  rownames(temp_sample) <- colnames(data)
  colnames(temp_sample) <- 'label'
  temp_sample <- data.frame(temp_sample)
  temp_sample$label <- factor(temp_sample$label)
  physeq <- phyloseq(otu_table(temp_otu),
                     sample_data(temp_sample))
  ####
  # taken from phyloseq vignette
  diagdds = phyloseq_to_deseq2(physeq, ~ label)
  # calculate geometric means prior to estimate size factors
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(diagdds), 1, gm_mean)
  diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
  diagdds = DESeq(diagdds, fitType="local")
  res = results(diagdds)
  p.val <- res$pvalue
  names(p.val) <- rownames(res)
  return(p.val)
}

#' @keywords internal
test.MGS <- function(data, label, type=1, conf=NULL){
  test.package('metagenomeSeq')

  stopifnot(ncol(data) == length(label))
  stopifnot(all(sort(unique(label)) == c(-1, 1)))

  # first exclude features with low prevalence
  # (as recommended by the authors)
  non.sparse.feat = which(rowMeans(data > 0) >= 0.05)
  trim.feat = data[non.sparse.feat, ]
  non.sparse.samples = which(colSums(trim.feat) > 10)
  trim.feat = trim.feat[,non.sparse.samples]
  label <- label[names(non.sparse.samples)]
  mgs.obj = newMRexperiment(trim.feat)
  p = tryCatch({cumNormStat(mgs.obj, pFlag=FALSE, main="Trimmed data")},
               error=function(err){0.5})
  mgs.obj = cumNorm(mgs.obj, p=p)

  mod = model.matrix(~as.factor(label))
  rownames(mod) <- colnames(trim.feat)
  p.val <- tryCatch({
    if (type==1){
      settings = zigControl(maxit=20, verbose=TRUE)
      fit = fitZig(mgs.obj, mod=mod, control=settings)
    } else if (type==2){
      fit = fitFeatureModel(obj=mgs.obj, mod=mod)
    } else {
      stop("Unsupported type!")
    }

    res <- MRcoefs(fit, number = Inf)
    p.val = rep(1, nrow(data))
    names(p.val) = rownames(data)
    p.val[row.names(res)] = res$pvalues
    p.val},
    error=function(err){
      p.val <- rep(1, nrow(data))
      names(p.val) <- rownames(data)
      p.val
    })

  return(p.val)
}

#' @keywords internal
test.edgeR <- function(data, label){
  
  test.package('edgeR')

  stopifnot(ncol(data) == length(label))
  stopifnot(all(sort(unique(label)) == c(-1, 1)))

  # prepare count data
  nonzero.idx = which(rowSums(data) > 0)
  f.counts = data[nonzero.idx,]
  groups = rep('pos', ncol(f.counts))
  groups[label==-1] = 'neg'
  names(groups) = colnames(data)

  y = DGEList(counts=f.counts, group=groups)
  # test for error message
  t <- exp(rowMeans(log(f.counts)))
  if (any(t > 0)){
    y = edgeR::calcNormFactors(y, method='RLE')
  }
  y = estimateCommonDisp(y)
  y = estimateTagwiseDisp(y)
  et = exactTest(y)
  res = topTags(et, n=length(nonzero.idx), sort.by='none')
  stopifnot(all(names(nonzero.idx) == rownames(res$table)))
  p.val = rep(1.0, nrow(data))
  names(p.val) <- rownames(data)
  p.val[nonzero.idx] = res$table[,'PValue']
  return(p.val)
}

#' @keywords internal
test.ANCOM <- function(data, label, conf=FALSE){
  # This function is based on the scripts that i downloaded from:
  # https://sites.google.com/site/siddharthamandal1985/research
  # on 2019-01-28
  #
  # some minor edits from that original script, since here, we expect
  # data which are not repeated and the formula should not be adjusted
  test.package('exactRankTests')
  test.package('nlme')

  stopifnot(ncol(data) == length(label))
  stopifnot(all(sort(unique(label)) == c(-1, 1)))

  # adjusted ANCOM script:
  n_otu <- nrow(data)
  otu_ids <- rownames(data)
  base.formula <- "lr ~ label"
  fformula  <- formula(base.formula)
  tfun <- exactRankTests::wilcox.exact

  logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
  pb <- progress_bar$new(total=((n_otu^2)/2))
  for(ii in 1:(n_otu-1)){
    for(jj in (ii+1):n_otu){
      data.pair <- data[which(rownames(data)%in%otu_ids[c(ii,jj)]),]
      lr <- log((1+as.numeric(data.pair[1,]))/(1+as.numeric(data.pair[2,])))
      lr_dat <- data.frame(lr=lr, label=label,row.names=NULL )
      logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      pb$tick()
    }
  }
  ind <- lower.tri(logratio.mat)
  logratio.mat[ind] <- t(logratio.mat)[ind]
  logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1

  mc.pval <- t(apply(logratio.mat,1,function(x){
    s <- p.adjust(x, method = "BH")
    return(s)
  }))
  a <- logratio.mat[upper.tri(logratio.mat,diag=FALSE)==TRUE]
  b <- matrix(0,ncol=n_otu,nrow=n_otu)
  b[upper.tri(b)==TRUE] <- p.adjust(a, method = "BH")
  diag(b)  <- NA
  ind.1    <- lower.tri(b)
  b[ind.1] <- t(b)[ind.1]

  W <- apply(mc.pval,1,function(x){
    subp <- length(which(x<0.05))
  })

  # problem: the ANCOM method does not output p values, only effect sizes
  ## there is a recommendation for a cutoff (0.7*(n_otu - 1))
  # generate fake p-values (with a breakpoint at 0.05)
  p.val <- rep(1, length(W))
  names(p.val) <- rownames(data)
  fake.signif <- which(W > 0.7*(n_otu - 1))
  fake.non.signif <- which(W < 0.7*(n_otu - 1))
  if (length(fake.signif) == 0){
    fake.n.signif.pvalue <- seq(from=0.5, to=1,
                                length.out=length(fake.non.signif))
    p.val[fake.non.signif[order(W[fake.non.signif])]] <- fake.n.signif.pvalue
  } else {
    fake.signif.pvalue <- seq(from=0, to=0.04, length.out=length(fake.signif))
    p.val[fake.signif[order(W[fake.signif])]] <- fake.signif.pvalue
    if (length(fake.non.signif) != 0){
      fake.n.signif.pvalue <- seq(from=0.5, to=1,
                                  length.out=length(fake.non.signif))
      p.val[fake.non.signif[order(W[fake.non.signif])]] <- fake.n.signif.pvalue
    }
  }
  return(p.val)
}

#' @keywords internal
test.ANCOM.2 <- function(data, label, conf){
  test.package('dplyr')
  test.package('nlme')

  if (is.null(conf)){
    meta.dat <- data.frame('Sample'=names(label),
                           'label'=label)
    adj_formula <- NULL
  } else {
    meta.dat <- cbind(data.frame('Sample'=names(label),
                                 'label'=label), conf)
    adj_formula <- paste(colnames(conf), collapse='+')
  }
  preprocess <- feature_table_pre_process(data, meta.dat,
                                          sample_var='Sample',
                                          group_var=NULL,
                                          lib_cut = 100,
                                          neg_lb = FALSE)
  p.val.res <- tryCatch({
    p.val <- rep(1, nrow(data))
    names(p.val) <- rownames(data)
    res <- ANCOM(preprocess$feature_table, meta_data = preprocess$meta_data,
                 struc_zero = preprocess$structure_zeros, adj_formula=adj_formula,
                 main_var='label', p_adj_method = 'fdr')

    # problem: the ANCOM method does not output p values, only effect sizes
    # same problem here for ANCOM-2
    # solution: create fake p-value for evaluation
    fake.signif <- rownames(preprocess$feature_table)[
      which(res$out$detected_0.7)]
    W <- res$out$W
    names(W) <- res$out$taxa_id
    if (length(fake.signif) != 0){
      non.signif <- setdiff(res$out$taxa_id, fake.signif)
      # significant ones
      fake.signif.pvalue <- seq(from=0, to=0.04,
                                length.out=length(fake.signif))
      W.signif <- W[fake.signif]
      W.signif <- W.signif[order(W.signif)]
      p.val[names(W.signif)] <- fake.signif.pvalue

      # others
      fake.pvalue <- seq(from=0.5, to=0.99, length.out=length(non.signif))
      W.nonsignif <- W[non.signif]
      W.nonsignif <- W.nonsignif[order(W.nonsignif)]
      p.val[names(W.nonsignif)] <- fake.pvalue
      p.val
    } else {
      fake.pvalue <- seq(from=0.5, to=0.99, length.out=length(W))
      W <- W[order(W)]
      p.val[names(W)] <- fake.pvalue
      p.val
    }}, error=function(err){
      p.val <- rep(1, nrow(data))
      names(p.val) <- rownames(data)
      p.val})
  return(p.val.res)
}

#' @keywords internal
test.ANCOMBC <- function(data, label, conf){
  
  test.package('ANCOMBC')
  test.package("phyloseq")
  # deal with nonunique samples from biased idx generation
  names(label) <- paste0(names(label), '_', seq_along(label))
  colnames(data) <- paste0(colnames(data), '_', seq_along(label))
  label <- label + 1
  label <- label/2
  if (is.null(conf)){
    s.data <- data.frame(label=label, dummy=2)
    f.form <- 'label'
  } else {
    s.data <- cbind(data.frame(label=label, dummy=2), conf)
    f.form <- paste0('label+', paste(colnames(conf), collapse = '+'))
  }
  x.phylo <- phyloseq(otu_table = otu_table(data, taxa_are_rows = TRUE),
                      sample_data = sample_data(s.data))
  temp <- ancombc(phyloseq = x.phylo, formula = f.form,
                  p_adj_method = 'fdr', lib_cut = 100)
  p.val <- rep(1, nrow(data))
  names(p.val) <- rownames(data)
  p.val[rownames(temp$res$q_val)] <- temp$res$p_val$label
  return(p.val)
}

#' @keywords internal
test.ZIBseq <- function(data, label, conf, transform=FALSE){
  if (!is.null(conf)){
    stop("'ZIBSeq' is not a confounder-aware test!")
  }
  test.package('ZIBseq')
  stopifnot(ncol(data) == length(label))
  stopifnot(all(sort(unique(label)) == c(-1, 1)))
  # actual code form ZIBSeq

  genodata = t(data)
  Y = label
  useF = which(colSums(genodata) > 2 * dim(genodata)[1])
  X = genodata[, useF]
  ST = rowSums(X)
  P = dim(X)[2]
  beta = matrix(data = NA, P, 2)
  .f_test <- function(a, b){
    bereg <- gamlss(a ~ b, family = BEZI(sigma.link = "identity"),
                    trace = FALSE, control = gamlss.control(n.cyc = 100))
    out = summary(bereg)
    return(out[2, c(1, 4)])
  }
  for (i in 1:P) {
    x.prop = X[, i]/ST
    if (transform == TRUE) {
      x.prop = sqrt(x.prop)
    }
    # need to add a try-catch statement
    out <- tryCatch({
      .f_test(x.prop, Y)
    }, error=function(err){
      c(NA_real_, NA_real_)
    })
    beta[i, ] = out
  }
  pvalues = beta[, 2]
  # some edits to match p values back to original data
  p.val <- rep(1, nrow(data))
  names(p.val) <- rownames(data)
  p.val[names(useF)] <- pvalues
  return(p.val)
}

#' @keywords internal
test.ZIBseq_sqrt <- function(data, label, conf){
  test.ZIBseq(data, label, conf, transform=TRUE)
}

#' @keywords internal
test.corncob <- function(data, label, conf){
  test.package('corncob')
  test.package('phyloseq')
  if (!is.null(conf)){
    stop("'Corncob' is not a confounder-aware test!")
  }
  label <- label + 1
  label <- label/2
  x.phylo <- phyloseq(otu_table = otu_table(data, taxa_are_rows = TRUE),
                      sample_data = sample_data(as.data.frame(label)))
  x <- differentialTest(formula = ~label,
                        phi.formula = ~label,
                        formula_null = ~1,
                        phi.formula_null = ~label,
                        data=x.phylo,
                        test='Wald')

  p.vals <- x$p
  return(p.vals)
}

#' @keywords internal
test.via.limma <- function(data, label, conf){
  test.package('limma')
  if (is.null(conf)){
    fit <- lmFit(data, design = data.frame(intercept = -1, label))
  } else {
    x <- duplicateCorrelation(data, design = data.frame(label), block = conf[colnames(data), 'conf'])
    fit <- lmFit(data, design = data.frame(intercept = -1, label), block = conf[colnames(data), 'conf'], correlation = x[[1]])
  }
  res <- eBayes(fit)
  p.val <- res$p.value[,'label']
  return(p.val)
}

#' @keywords internal
test.lm <- function(data, label, conf){
  if (is.null(conf)){
    p.vals <- vapply(rownames(data), FUN=function(x){
      fit <- lm(as.numeric(data[x,])~label)
      res <- anova(fit)
      return(res$`Pr(>F)`[1])
    }, FUN.VALUE = double(1))
  } else {
    f <- paste0("feat~label+", paste(colnames(conf), collapse = '+'))
    p.vals <- vapply(rownames(data), FUN=function(x){
      df <- cbind(data.frame(feat=data[x,], label=label), conf)
      fit <- lm(data=df, formula = as.formula(f))
      res <- anova(fit)
      return(res$`Pr(>F)`[1])
    }, FUN.VALUE = double(1))
  }
  return(p.vals)
}

#' @keywords internal
test.lme <- function(data, label, conf){
  if (is.null(conf)){
    stop("Test 'LME' needs confounders!")
  } else {
    test.package("lmerTest")
    formula <- paste0("feat~label+",
                paste(paste0('(', paste0('1|', colnames(conf)), ')'),
                      collapse = '+'))
    conf <- conf[colnames(data),,drop=FALSE]
    p.vals <- vapply(rownames(data), FUN=function(x){
      p <- tryCatch({
        df <- cbind(data.frame(feat=data[x,], label=label), conf)
        fit <- suppressMessages(suppressWarnings(
          lmer(data=df, formula = as.formula(formula))))
        res <- coefficients(summary(fit))
        res['label', 'Pr(>|t|)']
      },
      error=function(err){
        p <- NA_real_
        p})
      return(p)
    }, FUN.VALUE = double(1))
  }
  return(p.vals)
}

#' @keywords internal
#' lrtest of nested models; type1=FE, type2=RE
test.metadeconfoundR <- function(data, label, type=1, conf){
  # https://github.com/TillBirkner/metadeconfoundR/blob/
  # 28e860a60815bf3ce5b21a810e037af801826945/R/CheckReducibility.R#L261
  # browser()
  test.package('lmtest')
  if (is.null(conf)){
    stop("Test 'metadeconfoundR' needs confounders!") }
  else {
    conf <- conf[names(label),,drop=FALSE]
    # wilcoxon test as precondition for modeling, usually
    wilcox <- vapply(rownames(data), FUN=function(x) {
      # check sparsity -- return high p-val
      if (var(data[x,]) == 0) return(1.0)
      df <- cbind(data.frame(feat=as.numeric(data[x,]), label=label), conf)
      t = wilcox.test(df[label==-1, 'feat'],
                      df[label==+1, 'feat'], exact=FALSE)
      return(t$p.value) }, FUN.VALUE = double(1))
    # confounder as a fixed effect
    if (type == 1) {
      f.both <- paste0('rank(feat)~label+', paste(colnames(conf),
                                                  collapse = '+'))
      f.conf <- paste0('rank(feat)~', paste(colnames(conf),
                                            collapse = '+'))
      f.lab <- 'rank(feat)~label' }
    # confounder as a random effect
    else if (type == 2) {
      f.both <- paste0('rank(feat)~label+', paste(paste0(
        '(',paste0('1|', colnames(conf)), ')'), collapse = '+'))
      f.conf <- paste0('rank(feat)~1+', paste(paste0(
        '(', paste0('1|', colnames(conf)), ')'), collapse = '+'))
      f.lab <- 'rank(feat)~label' }
    else { stop('Model configuration not supported!') }

    # significance of label beyond confounder, i.e.
    # signal reducible to label
    sig.lab.conf <- vapply(rownames(data), FUN=function(x){
      p <- tryCatch({
        df <- cbind(data.frame(feat=as.numeric(data[x,]), label=label), conf)
        fit.both <- suppressMessages(suppressWarnings(
          lm(data = df,
             formula = as.formula(f.both),
             REML = FALSE)))
        fit.conf <- suppressMessages(suppressWarnings(
          lm(data = df,
             formula = as.formula(f.conf),
             REML = FALSE)))
        lmtest::lrtest(fit.both, fit.conf)$'Pr(>Chisq)'[2] },
        error = function(err){
          p <- NA_real_
          p })
      return(p) }, FUN.VALUE = double(1))

    # significance of confounder beyond label, i.e.
    # signal reducible to covariate
    sig.conf.lab <- vapply(rownames(data), FUN=function(x){
      p <- tryCatch({
        df <- cbind(data.frame(feat=as.numeric(data[x,]), label=label), conf)
        fit.both <- suppressMessages(suppressWarnings(
          lm(data = df,
             formula = as.formula(f.both),
             REML = FALSE)))
        fit.label <- suppressMessages(suppressWarnings(
          lm(data = df,
             formula = as.formula(f.lab),
             REML = FALSE)))
        lmtest::lrtest(fit.both, fit.label)$'Pr(>Chisq)'[2] },
        error = function(err){
          p <- NA_real_
          p })
      return(p) }, FUN.VALUE = double(1))

    # status determination -- update significances accordingly
    p.vals <- cbind(wilcox, sig.lab.conf, sig.conf.lab)
    corrected <- data.frame(FDR = p.adjust(p.vals[,'wilcox'], method = 'fdr'),
                            A = sig.conf.lab,
                            B = sig.lab.conf)
    status <- vapply(rownames(corrected), FUN=function(x) {
      # ignore lrtest if wilcoxon not significant (match mdc behavior)
      if (corrected[x,'FDR'] >= 0.05) s <- corrected[x,'FDR']
      # corrected pval significant -- check lrtest
      else if (corrected[x,'FDR'] < 0.05) {
        # conf significant, label not
        if (corrected[x,'A'] < 0.05 && corrected[x,'B'] >= 0.05) {
          # s <- 'conf' -- take the insignificant p-val
          s <- corrected[x,'B'] }
        # neither significant
        else if (corrected[x,'A'] >= 0.05 && corrected[x,'B'] >= 0.05) {
          # s <- 'LD' -- take either insignificant p-val
          s <- max(corrected[x,'A'], corrected[x,'B']) }
        else if (corrected[x,'A'] >= 0.05 && corrected[x,'B'] < 0.05) {
          # s <- 'SD' -- take the significant p-val
          s <- corrected[x,'B'] }
        else if (corrected[x,'A'] < 0.05 && corrected[x,'B'] < 0.05) {
          # both 'deconfounded' -- take either significant p-val
          s <- min(corrected[x,'A'], corrected[x,'B']) }
        else s <- NA }
      else s <- NA
      return(s) }, FUN.VALUE = double(1))

    return(status)
  }
}

#' @keywords internal
test.aldex2 <- function(data, label, conf){
  test.package("ALDEx2")
  data <- vapply(colnames(data), FUN=function(x){round(data[,x])},
                 FUN.VALUE = double(nrow(data)))
  if (is.null(conf)){
    temp <- aldex(reads=data, conditions=label)
    p.val <- rep(1, nrow(data))
    names(p.val) <- rownames(data)
    p.val[rownames(temp)] <- temp$wi.ep
  } else if (nrow(conf)==1){
    design <- model.matrix(~label + conf[,1])
    x <- aldex.clr(reads=data, conds = design)
    glm.test <- aldex.glm(x, design)
    glm.effect <- aldex.glm.effect(x)
    p.val <- rep(1, nrow(data))
    names(p.val) <- rownames(data)
    p.val[rownames(glm.test)] <- glm.test[['model.label Pr(>|t|)']]
  } else {
    stop("Not yet implemented for multiple confounders!")
  }
  return(p.val)
}

#' @keywords internal
test.ks <- function(data, label, conf){

  if (is.null(conf)){
    data <- data[,names(label)]
    # apply Kolmogorov-Smirnoff test
    p.val = rep(1, nrow(data))
	names(p.val) <- rownames(data)
    for (f in rownames(data)) {
      t = ks.test(as.numeric(data[f,label==-1]),
                  as.numeric(data[f,label==+1]), exact=FALSE)
      p.val[f] = t$p.value
    }
  } else {
    stop("Kolmogorov-Smirnoff test is not implemented for confounders!")
  }
  return(p.val)
}

#' @keywords internal
test.scde <- function(data, label, conf){
  test.package('scde')
  data <- data[,names(label)]
  cd <- clean.counts(data, min.lib.size=10,
                     min.reads = 1,
                     min.detected = 1)
  cd <- apply(cd, 2, function(x){storage.mode(x) <- 'integer'; x})
  o.ifm <- scde.error.models(counts = cd, groups = label[colnames(cd)],
                             threshold.segmentation = TRUE,
                             save.crossfit.plots = FALSE,
                             min.size.entries = 10,
                             save.model.plots = FALSE)
  valid.cells <- o.ifm$corr.a > 0

  o.ifm <- o.ifm[valid.cells, ]
  o.prior <- scde.expression.prior(models = o.ifm, counts = cd,
                                   length.out = 400, show.plot = FALSE)
  if (is.null(conf)){
    ediff <- scde.expression.difference(o.ifm, cd, o.prior,
                                        groups=as.factor(label[colnames(cd)]),
                                        n.randomizations=100,
                                        n.cores=1, verbose=1)
    ediff$cZ <- abs(ediff$cZ)
    p.val <- pnorm(ediff$cZ, lower.tail = FALSE)
    names(p.val) <- rownames(ediff)
  } else if (nrow(conf)==1) {
    ediff <- scde.expression.difference(o.ifm, cd, o.prior,
                                        groups=as.factor(label[colnames(cd)]),
                                        n.randomizations=100,
                                        batch=as.factor(conf[,1]),
                                        n.cores=1, verbose=1)
    p.val <- abs(ediff$batch.adjusted$cZ)
    p.val <- pnorm(p.val, lower.tail = FALSE)
    names(p.val) <- rownames(ediff$batch.adjusted)
  } else {
    stop("Not yet implemented for multiple confounders!")
  }

  return(p.val)
}

#' @keywords internal
test.MAST <- function(data, label, conf){
  test.package('MAST')

  p.val <- rep(1, nrow(data))
  names(p.val) <- rownames(data)

  # transform into single cell data experiment
  x <- FromMatrix(data, cData = as.data.frame(cbind(label, conf)),
                  check_sanity = FALSE)

  ## differential expression
  if (is.null(conf)){
    formula <- as.formula(~label)
  } else {
    formula <- as.formula(paste0('~label+', paste(colnames(conf),
                                                  collapse='+')))
  }
  fit <- zlm(formula, sca = x, ebayes = TRUE)
  coefs <- MAST::summary(fit, logFC=FALSE, doLRT=TRUE)$datatable
  mast.results <- coefs[coefs$component=='H' & coefs$contrast == 'label',]

  p.val[mast.results$primerid] <- mast.results$`Pr(>Chisq)`
  p.val <- p.adjust(p.val, method = 'fdr')
  return(p.val)
}

#' @keywords internal
test.mixMC <- function(data, label, conf){
  test.package('mixOmics')
  if (!is.null(conf)){
    stop("Method 'mixMC' is not compatible with confounders!")
  }
  data <- 10^data
  splsda.tune <- tune.splsda(X=t(data),
                             Y = label,
                             ncomp = 1,
                             multilevel = NULL,
                             logratio = 'CLR',
                             test.keepX = c(seq(5,150, 5)),
                             validation = c('Mfold'),
                             folds = 5,
                             dist = 'max.dist',
                             nrepeat = 10)
  res.splsda <- splsda(X=t(data),
                       Y = label, ncomp=1,
                       keepX = splsda.tune$choice.keepX[1],
                       logratio = 'CLR')

  p.val <- rep(1, nrow(data))
  names(p.val) <- rownames(data)

  vars <- selectVar(res.splsda)
  temp <- abs(vars$value$value.var)
  fake.p.val <- seq(from=0, to=0.049, length.out=length(temp))
  fake.p.val <- fake.p.val[order(-temp)]
  p.val[vars$name] <- fake.p.val

  return(p.val)
}

#' @keywords internal
test.distinct <- function(data, label, conf){
  test.package('distinct')
  test.package('SingleCellExperiment')
  # transform into single cell data experiment
  x <- SingleCellExperiment(list('norm_data'=data))
  colData(x)$sample_id <- colnames(data)
  colData(x)$cluster_id <- 'all'
  colData(x)$group <- label
  if (is.null(conf)){
    design <- model.matrix(~label)
    rownames(design) <- colnames(data)
  } else if (nrow(conf) == 1){
    design <- model.matrix(~label + conf[,1])
    rownames(design) <- colnames(data)
  } else {
    stop("Not yet implemented for mulitple confounders!")
  }
  res <- distinct_test(x, name_assays_expression = 'norm_data',
                       name_sample = 'sample_id', name_cluster = 'cluster_id',
                       design = design)
  p.val <- rep(1, nrow(data))
  names(p.val) <- rownames(data)
  p.val[res$gene] <- res$p_val
  return(p.val)
}

#' @keywords internal
test.gFC <- function(data, label, conf){
  n.pos <- names(which(label==1))
  n.neg <- names(which(label==-1))
  val.gFC <- vapply(rownames(data), FUN = function(x){
    g.pos <- quantile(data[x,n.pos], prob=seq(from=0.05, to=0.95, by=0.05))
    g.neg <- quantile(data[x,n.neg], prob=seq(from=0.05, to=0.95, by=0.05))
    abs(mean(g.pos-g.neg))
  }, FUN.VALUE = double(1))
  return(val.gFC)
}

#' @keywords internal
test.songbird <- function(data, label, conf){
  browser()

  return(p.val)
}

#' @keywords internal
test.ZINQ <- function(data, label, conf){
  test.package('ZINQ')
  p.val <- rep(1, nrow(data))
  names(p.val) <- rownames(data)

  df.test <- as.data.frame(t(data))
  df.test$group <- label
  if (!is.null(conf)){
    df.test <- cbind(df.test, as.data.frame(conf))
    cov <- c('group', colnames(conf))
  } else {
    cov <- 'group'
  }
  for (sp in rownames(data)){
    temp <- df.test[,c(sp, cov)]
    if (ncol(temp) == 2){
      colnames(temp) <- c('Y', 'X')
      p.sp <- tryCatch({
        res <- ZINQ_tests(Y~X, Y~X, C='X', data=temp)
        ZINQ_combination(res)
      }, error=function(err){1})
      p.val[sp] <- p.sp
    } else {
      colnames(temp)[1:2] <- c('Y', 'X')
      formula <- as.formula(paste0('Y~X+', paste0(colnames(temp)[-c(1,2)],
                                                  collapse = '+')))
      p.sp <- tryCatch({
        res <- ZINQ_tests(formula, formula, C='X', data=temp)
        ZINQ_combination(res)
      }, error=function(err){1})
      p.val[sp] <- p.sp
    }
  }
  return(p.val)
}

#' @keywords internal
test.package <- function(x){
  suppressWarnings(
    suppressPackageStartupMessages(
      success <- library(x, character.only = TRUE, logical.return = TRUE,
                         verbose = FALSE, quietly = FALSE)))
  if (!success){
    stop("Package ", x, " is not yet installed!")
  }
}
