###
# SIAMCAT -  Statistical Inference of Associations between Microbial Communities And host phenoTypes
# RScript flavor
#
# written by Georg Zeller
# with additions by Nicolai Karcher and Konrad Zych
# EMBL Heidelberg 2012-2017
#
# version 0.2.0
# file last updated: 26.06.2017
# GNU GPL 3.0
###

#' @title Check for potential confounders in the metadata
#' @description This function checks for associations between class labels and potential confounders (e.g. age, sex, or BMI) that are present in the metadata.
#' Statistical testing is performed with Fisher's exact test or Wilcoxon test, while associations are visualized either as barplot or Q-Q plot, depending on the type of metadata.
#' @param meta metadata object
#' @param label labels object
#' @param fn.plot string, filename for the pdf-plot
#' @keywords SIAMCAT confounder.check
#' @export
#' @return Does not return anything, but produces a single plot for each metadata category
confounder.check <- function(meta, label, fn.plot){

  pdf(fn.plot, onefile=TRUE)
  par(xpd=FALSE)
  for (m in 1:ncol(meta)) {
    mname <- gsub('[_.-]', ' ', colnames(meta)[m])
    mname <- paste(toupper(substring(mname, 1, 1)), substring(mname, 2), sep="")
    cat('checking', mname, 'as a potential confounder...\n')

    mvar  <- as.numeric(meta[,m])
    u.val <- unique(mvar)
    u.val <- u.val[!is.na(u.val)]
    colors <- brewer.pal(5,"Dark2")

    if (length(u.val) == 1) {
      cat('  skipped because all subjects have the same value\n')
    } else if (length(u.val) <= 5) {
      cat('  using a bar plot\n')
      par(mar=c(6.1,4.1,4.1,4.1))
      ct     <- matrix(NA, nrow=2, ncol=length(u.val))

      for (i in 1:length(u.val)) {
        ct[1,i] = sum(mvar[label$n.idx] == u.val[i], na.rm=TRUE)
        ct[2,i] = sum(mvar[label$p.idx] == u.val[i], na.rm=TRUE)
      }

      freq   <- t(ct)

      for (i in 1:dim(freq)[2]) {
        freq[,i] <- freq[,i] / sum(freq[,i])
      }

      barplot(freq, ylim=c(0,1), main=mname, names.arg=c(label$n.lab, label$p.lab), col=colors)
      #legend("right", legend=u.val, col=colors)

      p.val  <- fisher.test(ct)$p.value
      mtext(paste('Fisher test p-value:', format(p.val, digits=4)), side=1, line=3, at=1, adj=0)
    } else {
      cat('  using a Q-Q plot\n')
      par(mar=c(5.1,4.1,4.1,4.1))
      ax.int <- c(min(mvar, na.rm=TRUE), max(mvar, na.rm=TRUE))

      qqplot(mvar[label$n.idx], mvar[label$p.idx], xlim=ax.int, ylim=ax.int, pch=16, cex=0.6,
             xlab=label$n.lab, ylab=label$p.lab, main=paste('Q-Q plot for', mname))
      abline(0, 1, lty=3)

      p.val  <- wilcox.test(mvar[label$n.idx], mvar[label$p.idx], exact=FALSE)$p.value
      text(ax.int[1]+0.9*(ax.int[2]-ax.int[1]), ax.int[1]+0.1*(ax.int[2]-ax.int[1]),
           paste('MWW test p-value:', format(p.val, digits=4)), pos=2)
    }
  }
  tmp <- dev.off()
}
