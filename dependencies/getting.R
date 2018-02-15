package.list <- c("RColorBrewer",
                  "beanplot",
                  "glmnet",
                  "LiblineaR",
                  "pROC",
                  "optparse",
                  "colorRamps",
                  "gelnet",
                  "mlr",
                  "PRROC",
                  "knitr",
                  "testthat",
                  "rmarkdown",
                  "viridis",
                  "gridExtra",
                  "gridBase",
                  "grid",
                  "phyloseq",
                  "ade4",
                  "ape",
                  "Biobase",
                  "BiocGenerics",
                  "biomformat",
                  "Biostrings",
                  "cluster",
                  "data.table",
                  "foreach",
                  "ggplot2",
                  "igraph",
                  "methods",
                  "multtest",
                  "plyr",
                  "reshape2",
                  "scales",
                  "vegan")

download.packages(package.list, destdir="/Users/zych/Documents/git/SIAMCAT/dependencies/", 
                  type="source")

#stopifnot(require("tools")) ## load tools

instPkgPlusDeps <- function(pkg, install = FALSE,
                            which = c("Depends", "Imports", "LinkingTo"),
                            inc.pkg = TRUE) {
  print(pkg)
  ap <- available.packages() ## takes a minute on first use
  ## get dependencies for pkg recursively through all dependencies
  deps <- package_dependencies(pkg, db = ap, which = which, recursive = TRUE)
  ## the next line can generate warnings; I think these are harmless
  ## returns the Priority field. `NA` indicates not Base or Recommended
  pri <- sapply(deps[[1]], packageDescription, fields = "Priority")
  ## filter out Base & Recommended pkgs - we want the `NA` entries
  deps <- deps[[1]][is.na(pri)]
  ## install pkg too?
  if (inc.pkg) {
    deps = c(pkg, deps)
  }
  ## are we installing?
  download.packages(deps, destdir="/Users/zych/Documents/git/SIAMCAT/dependencies/", 
                    type="source")
  deps ## return dependencies
}

#sapply(package.list,instPkgPlusDeps,FALSE,"Depends")
