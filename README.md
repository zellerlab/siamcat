# SIAMCAT <img src="man/figures/logo.png" align="right" width="120" />

## Overview
`SIAMCAT` is an R package for Statistical Inference of Associations
between Microbial Communities And host phenoType. It is part of the suite of 
computational microbiome analysis tools hosted at [EMBL](https://www.embl.org)
by the groups of
[Peer Bork](https://www.embl.de/research/units/scb/bork/index.html) and
[Georg Zeller](https://www.embl.de/research/units/scb/zeller/index.html). Find
out more at [EMBL-microbiome tools](http://microbiome-tools.embl.de/).

## Starting with SIAMCAT

### Installation

In order to start with SIAMCAT, you need to install it from Bioconductor:
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SIAMCAT", version = "3.8")
```

Alternatively, you can install the current development version via `devtools`:
```R
require("devtools")
devtools::install_github(repo = 'zellerlab/siamcat')
```

### Quick start

There are a few manuals that will kick-start you and help you analyse your
data with SIAMCAT. You can find links to those on the
[Bioconductor website of SIAMCAT](https://bioconductor.org/packages/release/bioc/html/SIAMCAT.html)
or you can type into `R`:
```R
browseVignettes("SIAMCAT")
```
> Please note:  
`browseVignettes` only works if `SIAMCAT` has been installed via Bioconductor

## Contact

If you run into any issue:
- create an
[issue in this repository](https://github.com/zellerlab/siamcat/issues/new) or
- mail Georg Zeller (mailto: zeller@embl.de) or
- ask at the
[SIAMCAT support group](https://groups.google.com/forum/#!forum/siamcat-users)

## Citation

If you use `SIAMCAT`, please cite us by using

```R
citation("SIAMCAT")
```

or by

> Zych K, Wirbel J, Essex M, Breuer K, Karcher N, Costea PI, Sunagawa S, Bork P,
Zeller G (2018). _SIAMCAT: Statistical Inference of Associations between Microbial
Communities And host phenoTypes_. doi: 10.18129/B9.bioc.SIAMCAT (URL:
http://doi.org/10.18129/B9.bioc.SIAMCAT), R package version 1.0.1, (URL:
https://bioconductor.org/packages/SIAMCAT/).
