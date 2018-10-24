# SIAMCAT - Statistical Inference of Associations between Microbial Communities And host phenoType

## Overview
<<<<<<< HEAD
`SIAMCAT` is part of the suite of computational microbiome analysis tools
hosted at [EMBL](https://www.embl.org) by the groups of
[Peer Bork](https://www.embl.de/research/units/scb/bork/index.html) and
[Georg Zeller](https://www.embl.de/research/units/scb/zeller/index.html). Find
out more at [EMBL-microbiome tools](http://microbiome-tools.embl.de/).

## Starting with SIAMCAT
In order to start with SIAMCAT, you need to install it from Bioconductor:
```
source("https://bioconductor.org/biocLite.R")
biocLite("SIAMCAT")
```

There are few manuals that will kick-start you and help you analyse your
data with SIAMCAT:
```
browseVignettes("SIAMCAT")
```

## Contact 

If you run into any issue:
- mail Georg Zeller (mailto: zeller@embl.de)
or
- create an issue in this repository
or
- ask at the [SIAMCAT support group](https://groups.google.com/forum/#!forum/siamcat-users)

## Citation

If you use `SIAMCAT`, please cite us by using

```
citation("SIAMCAT")
```

or by

> Zych K, Wirbel J, Essex M, Breuer K, Karcher N, Costea PI, Sunagawa S, Bork P,
Zeller G (2018). _SIAMCAT: Statistical Inference of Associations between Microbial
Communities And host phenoTypes_. doi: 10.18129/B9.bioc.SIAMCAT (URL:
http://doi.org/10.18129/B9.bioc.SIAMCAT), R package version 1.0.1, <URL:
https://bioconductor.org/packages/SIAMCAT/>.
=======

SIAMCAT is a pipeline for Statistical Inference of Associations between Microbial Communities And host phenoTypes.  
A primary goal of analyzing microbiome data is to determine changes in community composition that are associated with 
environmental factors. In particular, linking human microbiome composition to host phenotypes such as diseases has 
become an area of intense research. For this, robust statistical modeling and biomarker extraction toolkits are crucially 
needed.  
SIAMCAT provides a full pipeline supporting data preprocessing, statistical association testing, statistical 
modeling (LASSO logistic regression) including tools for evaluation and interpretation of these models (such as cross 
validation, parameter selection, ROC analysis and diagnostic model plots).  
SIAMCAT is available in three different flavors: 
 + Galaxy web server
 + command line tool
 + R package

Please see the Support Section if you run into problems when using SIAMCAT.

## Input data format

The input data should be organized in the same way for every version of SIAMCAT. All files are in tab-separated column format

 + **Label data**: First row is expected to be `#BINARY:1=[label for cases];-1=[label for controls]`  
    Second row should contain the sample identifiers as tab-separated list (consistent with feature and metadata). 
    Third row is expected to contain the actual class labels (tab-separated), e.g. `1` for each case and `-1` for each 
    control.  
    Note: Labels can take other numeric values (but not characters or strings); importantly, the label for cases has to 
    be greater than the one for controls.

 + **Feature matrix**: _features_ (in rows) _x samples_ (in columns)  
    First row should contain sample labels (consistent with label data), while the first column should contain feature 
    labels (e.g. taxonomic identifiers). The remaining entries are expected to be real values >= 0 that quantify the 
    abundance of each feature in each sample.

 + **Metadata** (optional): _samples_ (in rows) _x metadata_ (in columns)  
    Metadata needs to be converted to numerical values by the user (This is necessary for heatmap displays)!


## Galaxy interface:

The Galaxy interface can be found here: http://siamcat.embl.de/

### Galaxy in brief


 + **Left panel**: `TOOLS` lists available analysis modules.  
    Click to choose which ones you'd like to run.
 + **Right panel**: `HISTORY` keeps track of every analysis step you have perfomed.  
  + Click on the "eye" icon to view data, or click on the "floppy disk" icon to download results  
  + If something goes wrong the "i" icon can provide useful details
  + the "circular arrows" icon allows to rerun a job
  + you can delete analysis steps from your history using the "x" icon

 + **Central panel**: `ANALYZE DATA` allows to specify input data sets and parameters for each analysis module


Additional info: https://usegalaxy.org/ (in particular the Help menu) and
                 https://wiki.galaxyproject.org/Learn


### Getting started with Galaxy

Start by uploading your data (see above for input data formats) using the `DATA IMPORT / Import Data module / Upload File`

Then procede by executing all SIAMCAT modules in order (from A to I). 
See example history / Workflow as well as each module's description for specific information on input and output data


## Commandline version 

The commandline version are a collection of modules implemented in R which are called via a bash script.

+ Stable version: https://github.com/gezel/siamcat/

+ Developmental version (only available inside the EMBL intranet): `beta:/g/bork4/zeller/dev/siamcat`

```
# type
git clone beta:/g/bork4/zeller/dev/siamcat
# in the folder in which you'd like to clone the siamcat repository
```

R packages required to run SIAMCAT:
```R
install.packages('optparse')
install.packages('LiblineaR')
install.packages('pROC')
install.packages('colorRamps')
install.packages('RColorBrewer')
install.packages('beanplot')
```

### Using the Commandline version

...COMING SOON...

## R package

The SIAMCAT R package ...COMING SOON...

### Using the R package

...COMING SOON...

## Support

Google user group for support:

https://groups.google.com/d/forum/siamcat-users


### Known issues

Examples are weighted differently between classes (a remnant of our colorectal cancer 
microbiome study). Fixed in Galaxy, will be pushed to GitHub soon.

Class labels are somehow swapped in the LASSO module, so that prediction scores are 1 - p
instead of p (posterior probability), consequently precision-recall curves are incorrect,
but ROC-curves are unaffected. Appears to only occur in a recent version of R and/or the
LiblineaR package; will be fixed with high priority.

### Contact 

Please let me know if you run into any issues (mailto: zeller@embl.de)
>>>>>>> upstream/master
