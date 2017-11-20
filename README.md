### Overview

SIAMCAT is a pipeline for Statistical Inference of Associations between Microbial
Communities And host phenoTypes. A primary goal of analyzing microbiome data is to
determine changes in community composition that are associated with environmental factors.
In particular, linking human microbiome composition to host phenotypes such as diseases
has become an area of intense research. For this, robust statistical modeling and
biomarker extraction toolkits are crucially needed. SIAMCAT provides a full pipeline
supporting data preprocessing, statistical association testing, statistical modeling
(LASSO logistic regression) including tools for evaluation and interpretation of these
models (such as cross validation, parameter selection, ROC analysis and diagnostic model
plots). SIAMCAT is available as Galaxy web server or can be downloaded and run as a
command line tool, as detailed below.


### Input data format

All files are in tab-separated column format

**Label data**:
                     First row is expected to be
                     #BINARY:1=[label for cases];-1=[label for controls]
                     Second row should contain the sample identifiers as tab-separated list
                     (consistent with feature and metadata).
                     Third row is expected to contain the actual class labels (tab-separated)
                     1 for each case and -1 for each control.
                     Note: Labels can take other numeric values (but not characters or strings);
                     importantly, the label for cases has to be greater than the one for controls.

**Feature matrix**:      features (in rows) x samples (in columns)
                     First row should contain sample labels (consistent with label data).
                     First column should contain feature labels (e.g. taxonomic identifiers).
                     The remaining entries are expected to be real values >= 0
                     that quantify the abundance of each feature in each sample.

**Metadata (optional)**: samples (in rows) x metadata (in columns)
                     Metadata needs to be converted to numerical values by the user
                     (This is necessary for heatmap displays)!



### Galaxy interface:

http://siamcat.embl.de/

#### Galaxy in brief


Left panel:      TOOLS lists available analysis modules,
                 click to choose which ones you'd like to run

Right panel:     HISTORY keeps track of every analysis step you've been performing
                 click on the "eye" icon to view data, or
                 click on the "floppy disk" icon to download results
                 if something goes wrong the "i" icon can provide useful details
                 the "circular arrows" icon allows to rerun a job
                 you can delete analysis steps from your history using the "x" icon

center panel:    ANALYZE DATA allows to specify input data sets and parameters for each
                 analysis module

Additional info: https://usegalaxy.org/ (in particular the Help menu) and
                 https://wiki.galaxyproject.org/Learn


#### Getting started


Start by uploading your data (see above for input data formats) using the
DATA IMPORT / Import Data module / Upload File

Then procede by executing all SIAMCAT modules in order (from A to I).
See example history / Workflow as well as each module's description for specific information on input and output data


#### Stable version:
https://git.embl.de/grp-zeller/SIAMCAT

#### Development version (the most up to date but not always stable):
https://git.embl.de/grp-zeller/SIAMCAT/tree/development


R packages required to run SIAMCAT:
```r
install.packages(c("RColorBrewer","beanplot","glmnet","LiblineaR","pROC","optparse","colorRamps","gelnet","mlr"))
```


### Commandline version (bash script calling modules implemented in R)
https://git.embl.de/grp-zeller/SIAMCAT/tree/master/Rscript_flavor

### License

Distributed under the GNU GPLv3 License. See accompanying file [LICENSE](https://git.embl.de/grp-zeller/SIAMCAT/blob/development/LICENSE) or copy at http://www.gnu.org/licenses/gpl-3.0.html.

### Support


Google user group for support:

https://groups.google.com/d/forum/siamcat-users
mailto:zeller@embl.de
