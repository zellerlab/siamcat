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
                  "rmarkdown","viridis", "gridExtra", "gridBase",
                  "energy", "corrplot", "coin")

notInst      <- which(!package.list%in%installed.packages())
if(length(notInst)>0) install.packages(package.list[notInst], repos="http://cran.uni-muenster.de")
biocLite('phyloseq')