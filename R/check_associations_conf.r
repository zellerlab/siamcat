check.associations.conf <- function(siamcat, fn.plot, color.scheme="RdYlBu",
                                    alpha=0.05, mult.corr="fdr", sort.by="fc", probs.fc=seq(.1, .9, .05),
                                    detect.lim=NULL, pr.cutoff=10^-6, max.show=50,
                                    plot.type="quantile.box", panels=c("fc", "auroc"), verbose=2){
  
  meta.mat <- as.matrix(siamcat@phyloseq@sam_data, nrow=nrow(siamcat@phyloseq@sam_data, ncol=ncol(siamcat@phyloseq@sam_data)))
  feat <- matrix(siamcat@phyloseq@otu_table,nrow=nrow(siamcat@phyloseq@otu_table), ncol=ncol(siamcat@phyloseq@otu_table),
                 dimnames = list(rownames(siamcat@phyloseq@otu_table), colnames(siamcat@phyloseq@otu_table)))
  col <- check.color.scheme(color.scheme, siamcat@label)
  
  
  if(verbose>2) cat("+++ analysing features\n")
  result.list <- marker.analysis.with.metadata(feat=feat, label=siamcat@label, meta=meta.mat, detect.lim=detect.lim, colors=col,
                                       pr.cutoff=pr.cutoff, mult.corr=mult.corr, alpha=alpha,
                                       max.show=max.show, sort.by=sort.by, probs.fc=probs.fc,verbose=verbose)
}


### performs pairwise association tests while blocking for metadata;
### additionally performs other differential abundance methods from standard 
### association testing module
marker.analysis.with.metadata <- function(feat, label, meta, detect.lim, colors,
                                          pr.cutoff, mult.corr, alpha,
                                          max.show, sort.by, probs.fc, verbose) {
  
  if(verbose>1) cat("+ starting marker analysis with metadata\n")
  s.time <- proc.time()[3]
  meta.colnames <- colnames(meta)

  # define set of matrices that have the indeces and "description" of all positively and negatively labeled training examples.
  p.val.mat   <- matrix(NA, nrow=nrow(feat), ncol=length(meta.colnames)+1, dimnames=list(row.names(feat), c(meta.colnames, "no block")))
  p.adj       <- matrix(NA, nrow=nrow(feat), ncol=length(meta.colnames)+1, dimnames=list(row.names(feat), c(meta.colnames, "no block")))
  fc          <- matrix(NA, nrow=nrow(feat), ncol=length(meta.colnames), dimnames=list(row.names(feat), meta.colnames))
  pr.shift    <- matrix(NA, nrow=nrow(feat), ncol=length(meta.colnames), dimnames=list(row.names(feat), meta.colnames))
  aucs        <- vector('list', nrow(feat))
  
  if(verbose>1) cat('+++ calculating effect size for each feature.\n')
  if(verbose) pb = txtProgressBar(max=nrow(feat), style=3)

  # loop over features
  for (f in 1:nrow(feat)) {
    
    # loop over metadata columns
    for (c in meta.colnames) { 
      meta.curr <- meta[,c]
      continuous <- FALSE
      
      # split continuous metadata - matrix
      if (length(unique(meta.curr[!is.na(meta.curr)])) > 5) {
        continuous <- TRUE
        quartiles <- quantile(meta.curr, probs=seq(0,1,0.25), na.rm=TRUE) 
        meta.curr.disc <- cut(meta.curr, breaks=quartiles, labels=1:4, include.lowest=TRUE)}
      
      # set up blocking variable
      if (continuous == TRUE) {
        block <- as.factor(meta.curr.disc)} else {
          block <- as.factor(meta.curr)}
      
      # keep track of NA values to exclude later
      na.samples <- vector('character')
      for (n in 1:length(block)) {
        if (is.na(block[n])) {
          na.samples <- c(na.samples, rownames(meta.curr[n]))}
      }
      block <- block[!names(label@label) %in% na.samples]
      
      # wilcoxon test blocking for this metadata category
      d <- data.frame(y=feat[f, which(!names(label@label) %in% na.samples)], 
                      x=as.factor(label@label[!(names(label@label) %in% na.samples)]), 
                      block, row.names=names(label@label[!(names(label@label) %in% na.samples)]))
      p.val.mat[f,c] <- pvalue(wilcox_test(y ~ x | block, data=d))
    }
    
    # redefine x & y for unstratified differential abundance methods
    x <- feat[f, label@p.idx] #cases
    y <- feat[f, label@n.idx] #controls
    
    # wilcoxon
    p.val.mat[f,"no block"] <- wilcox.test(y, x, exact = FALSE)$p.value
    
    # AUC
    #temp <- roc(predictor=feat[f,], response=siamcat@label@label, ci=TRUE, direction='<')
    temp <- roc(controls=y, cases=x, direction='<', ci=TRUE, auc=TRUE)
    aucs[[f]] <- c(temp$ci)
    
    # Prevalence-shift
    temp.n <- sum(y >= pr.cutoff)/length(y)
    temp.p <- sum(x >= pr.cutoff)/length(x)
    pr.shift[f] <- temp.p-temp.n
    
    # FC
    q.p <- quantile(log10(x+detect.lim), probs.fc)
    q.n <- quantile(log10(y+detect.lim), probs.fc)
    fc[f] <- sum(q.p - q.n)/length(q.p)
    if(verbose) setTxtProgressBar(pb, (pb$getVal()+1))
  } 
  if(verbose>1) cat('+++ finished looping over features\n')
  
  # apply multi-hypothesis testing correction
  if(!tolower(mult.corr) %in% c('none','bonferroni','holm','fdr','bhy')) {
    stop("! Unknown multiple testing correction method:', mult.corr,' Stopping!\n  Must of one of c('none','bonferroni', 'holm','fdr','bhy')")}
  if (mult.corr == 'none') {
    warning('WARNING: No multiple hypothesis testing performed.')
    p.adj <- p.val.mat} else {
      if(verbose>1) cat('+++ correcting for multiple testing\n')
      for (c in colnames(p.val.mat)) {
        p.adj[,c] <- p.adjust(p.val.mat[,c], method=tolower(mult.corr))}}
  
  if(verbose>1) cat('+++ found', sum(p.adj[,"no block"] < alpha, na.rm=TRUE), 'significant associations at a significance level <', alpha, '\n')
  
  # order species according to unstratified pvalues for fc, aucs, pr.shift...
  idx <- which(p.adj[,"no block"] < alpha)
  
  if (length(idx) == 0){stop('No significant associations found. Stopping.\n')}
  
  idx <- idx[order(p.adj[idx,"no block"], decreasing=TRUE)]
  
  e.time <- proc.time()[3]
  if(verbose>1) cat("+ finished analysing markers in",e.time-s.time,"s\n")
  return(list("p.val" = p.val.mat[idx,],
                 "fc"=fc[idx],
                 "aucs"=aucs[idx],
                 "pr.shift"=pr.shift[idx],
                 "p.adj"=p.adj[idx,]))
}


### anova on ranks for each feature to partition variance between label and metadata
feature.wise.anova.with.metadata <- function(feat, label, meta, colors, pr.cutoff,
                                             alpha, sort.by, verbose) {
  
  anova.heatmap <- matrix(NA, nrow=nrow(feat), ncol=length(meta.colnames)+1, 
                          dimnames=list(row.names(feat), c("disease status", meta.colnames)))
  
  # total sum of squares
  ss.total <- data.frame(cbind(apply(feat, 1, FUN=function(x) {
    rank.x <- rank(x)/length(x)
    return(sum((rank.x - mean(rank.x))^2)/length(rank.x))})))
  
  # label sum of squares - add this in with other groups later!
  ss.label <- data.frame(cbind(apply(feat, 1, FUN=function(x, label) {
    rank.p <- rank(x[label@p.idx])/length(label@p.idx)
    rank.n <- rank(x[label@n.idx])/length(label@n.idx)
    return(sum(sum((rank.p - mean(rank.p))^2), sum((rank.n - mean(rank.n))^2))/
      length(c(label@p.idx, label@n.idx)))
    }, label=label)))
  
  # metadata sum of squares - fix this
  if(verbose) pb = txtProgressBar(max=nrow(feat), style=3)
  ss.groups <- data.frame(cbind(apply(feat, 1, FUN=function(x, meta, label) {
    group.means <- list();
    
    for (c in unique(colnames(meta))) {
      if (length(unique(meta[,c][!is.na(meta[,c])])) > 5) {
        quartiles <- quantile(meta[,c], probs=seq(0,1,0.25), na.rm=TRUE) 
        meta.disc <- cut(meta[,c], breaks=quartiles, labels=1:4, include.lowest=TRUE)} else {
          meta.disc <- as.factor(meta[,c])}
      for (m in meta.disc) {
        if (is.na(m)) {next} else{
          idx.m <- which(meta.disc == m)
          rank.m <- rank(x[idx.m])/length(idx.m)
          index <- paste(c,m)
          group.means[[index]] <- mean(rank.m)}}}
    
    # need to do this for each metadata category... not all at once
    # to be continued
    
    # rank.p <- rank(x[label@p.idx])/length(label@p.idx)
    # rank.n <- rank(x[label@n.idx])/length(label@n.idx)
    # group.means[['cases']] <- mean(rank.p)
    # group.means[['controls']] <- mean(rank.n)
    grand.mean <- mean(rank(x)/length(x))
    if(verbose) setTxtProgressBar(pb, (pb$getVal()+1))
    squares <- lapply(group.means, FUN=function(y, grand.mean) {
      (y-grand.mean)^2}, grand.mean=grand.mean)
    return(do.call(sum, squares)/length(squares))
  }, meta=meta, label=label)))


}

check.color.scheme <- function(color.scheme, label, meta.studies=NULL,  verbose=1){
  if(verbose>2) cat("+ starting check.color.scheme\n")
  n.classes = ifelse(label@info$type == 'BINARY', 2, length(unique(label@label)))
  
  if (length(color.scheme) == 1 && class(color.scheme) == 'character'){
    if (n.classes == 2){
      # if color scheme and binary label, make colors as before
      if (!color.scheme %in% row.names(brewer.pal.info)){
        warning("Not a valid RColorBrewer palette name, defaulting to RdBu.\nSee brewer.pal.info for more information about RColorBrewer palettes.")
        color.scheme <- 'RdYlBu'
      }
      colors <- rev(colorRampPalette(brewer.pal(brewer.pal.info[color.scheme,'maxcolors'], color.scheme))(2))
    } else {
      # if color scheme and multiclass label, make colors either directly out of the palette (if n.classes smaller than maxcolors) or like before
      if (!color.scheme %in% row.names(brewer.pal.info)){
        warning("Not a valid RColorBrewer palette name, defaulting to Set3.\n  See brewer.pal.info for more information about RColorBrewer palettes.")
        color.scheme <- 'Set3'
      }
      # if color scheme and multiclass label, check that the palette is not divergent or sequential, but qualitative. Only issue warning.
      if (brewer.pal.info[color.scheme,'category'] != 'qual'){warning("Using a divergent or sequential color palette for multiclass data.")}
      if (n.classes <= brewer.pal.info[color.scheme, 'maxcolors']){
        colors <- brewer.pal(n.classes, color.scheme)
      } else {
        warning("The data contains more classes than the color.palette provides.")
        colors <- rev(colorRampPalette(brewer.pal(brewer.pal.info[color.scheme,'maxcolors'], color.scheme))(n.classes))
      }
    }
  } else if (length(color.scheme == n.classes) && all(is.color(color.scheme))){
    # if colors, check that all strings are real colors and check that the same length as n classes
    # convert color names to hex representation
    colors <- sapply(color.scheme, FUN=function(x){rgb(t(col2rgb(x)), maxColorValue = 255)}, USE.NAMES = FALSE)
  } else {
    stop("Supplied colors do not match the number of classes or are no valid colors")
  }
  # add transparency
  colors <- sapply(colors, FUN=function(x){paste0(x, '85')}, USE.NAMES = FALSE)
  if(verbose>2) cat("+ finished check.color.scheme\n")
  return(colors)
}