# SIAMCAT <img src="man/figures/logo.png" align="right" width="120" />

[![Build Status](https://travis-ci.com/zellerlab/siamcat.svg?branch=master)](https://travis-ci.com/zellerlab/siamcat)

## Overview
`SIAMCAT` is a pipeline for Statistical Inference of Associations between
Microbial Communities And host phenoTypes. A primary goal of analyzing
microbiome data is to determine changes in community composition that are
associated with environmental factors. In particular, linking human microbiome
composition to host phenotypes such as diseases has become an area of intense
research. For this, robust statistical modeling and biomarker extraction
toolkits are crucially needed. `SIAMCAT` provides a full pipeline supporting
data preprocessing, statistical association testing, statistical modeling
(LASSO logistic regression) including tools for evaluation and interpretation
of these models (such as cross validation, parameter selection, ROC analysis
and diagnostic model plots).

<a href='https://microbiome-tools.embl.de'> <img src="man/figures/embl_microbiome_tools_logo.png" align="right" width="200"> </a>

`SIAMCAT` is developed in the
[Zeller group](https://www.embl.de/research/units/scb/zeller/index.html)
and is part of the suite of computational microbiome analysis tools hosted at
[EMBL](https://www.embl.org/).

## Starting with SIAMCAT

### Installation

In order to start with `SIAMCAT`, you need to install it from Bioconductor:
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
data with `SIAMCAT`. You can find links to those on the
[Bioconductor website of SIAMCAT](https://bioconductor.org/packages/release/bioc/html/SIAMCAT.html)
or you can type into `R`:
```R
browseVignettes("SIAMCAT")
# Please Note:
# `browseVignettes` only works if `SIAMCAT` has been installed via Bioconductor
```

## Contact

If you run into any issue:
- create an
[issue in this repository](https://github.com/zellerlab/siamcat/issues/new) or
- mail [Georg Zeller](mailto:zeller@embl.de) or
- ask at the
[SIAMCAT support group](https://groups.google.com/forum/#!forum/siamcat-users)

If you found `SIAMCAT` useful, please consider giving us
[feedback](https://www.surveymonkey.de/r/denbi-service?sc=hd-hub&tool=siamcat).

## Citation

If you use `SIAMCAT`, please cite us by using

```R
citation("SIAMCAT")
```

or by

> Wirbel J, Zych K, Essex M, Karcher N, Kartal E, Salazar G, Bork P,
Sunagawa S, Zeller G (2020) _SIAMCAT: user-friendly and versatile machine
learning workflows for statistically rigorous microbiome analyses_
bioRxiv 2020.02.06.931808; doi: https://doi.org/10.1101/2020.02.06.931808

## Examples of primary package output

To give you a small preview about the primary package output, here are some
example plots taking from the main `SIAMCAT` vignette.

In this vignette, we use an example dataset which is also included in
the `SIAMCAT` package. The dataset is taken from the publication of
[Zeller et al](http://europepmc.org/abstract/MED/25432777), which demonstrated
the potential of microbial species in fecal samples to distinguish patients
with colorectal cancer (CRC) from healthy controls.

### Association testing

The result of the `check.associations` function is an association plot.
For significantly associated microbial features, the plot shows:
- the abundances of the features across the two different classes (CRC vs.
controls)
- the significance of the enrichment calculated by a Wilcoxon test (after
multiple hypothesis testing correction)
- the generalized fold change of each feature
- the prevalence shift between the two classes, and
- the Area Under the Receiver Operating Characteristics Curve (AU-ROC) as
non-parametric effect size measure.

![Association testing](man/figures/associations_plot.png)


### Model interpretation plot

After statistical models have been trained to distinguish cancer cases
from controls, the models can be investigated by the function
`model.interpretation.plot`. The plots shows:
- the median relative feature weight for selected features (barplot on the left)
- the robustness of features (i.e. in how many of the models the specific
feature has been selected)
- the distribution of selected features across samples (central heatmap)
- which proportion of the weight of all different models are shown in the plot
(boxplot on the right), and
- distribution of metadata across samples (heatmap below).

![Model interpretation plot](man/figures/interpretation_plot.png)

## Where SIAMCAT has been used already

Several publications already used `SIAMCAT` (or previous versions thereof).

- __[Potential of fecal microbiota for early-stage detection of colorectal cancer](http://europepmc.org/abstract/MED/25432777)__  
_Zeller G,  Tap J,  Voigt AY,  Sunagawa S,  Kultima JR,  Costea PI,  Amiot A,
Böhm J,  Brunetti F,  Habermann N,  Hercog R,  Koch M,  Luciani A,  Mende DR,
Schneider MA,  Schrotz-King P,  Tournigand C,  Tran Van Nhieu J,  Yamada T,
Zimmermann J,  Benes V,  Kloor M,  Ulrich CM,  von Knebel Doeberitz M,
Sobhani I,  Bork P_  
Molecular Systems Biology, (__2014__) 10, 766  
>Original Publication that inspired `SIAMCAT`

- __[Gut Microbiota Linked to Sexual Preference and HIV Infection](https://doi.org/10.1016/j.ebiom.2016.01.032)__  
_Noguera-Julian M, Rocafort M, Guillén Y, Rivera J, Casadellà M, Nowak P,
Hildebrand F, Zeller G, Parera M, Bellido R, Rodríguez C,Carrillo J, Mothe B,
Coll J, Bravo I, Estany C, Herrero C, Saz J, Sirera G, Torrela A, Navarro J,
Crespo M, Brander C, Negredo E, Blanco J, Guarner F, Calle ML, Bork P,
Sönnerborgd A, Clotet B, Paredes R_  
EBioMedicine 5 (__2016__) 135-146
>See Figure 5

- __[Extensive transmission of microbes along the gastrointestinal tract](https://elifesciences.org/articles/42693)__  
_Schmidt TSB, Hayward MR, Coelho LP, Li SS, Costea PI, Voigt AY, Wirbel J,
Maistrenko OM, Alves RJC, Bergsten E, de Beaufort C, Sobhani I,
Heintz-Buschart A, Sunagawa S, Zeller G, Wilmes P, Bork P_  
eLife, (__2019__) 8:e42693  
> See Figure 3 - figure supplement 1

- __[Meta-analysis of fecal metagenomes reveals global microbial signatures
that are specific for colorectal cancer](https://www.nature.com/articles/s41591-019-0406-6)__  
_Wirbel J, Pyl PT, Kartal E, Zych K, Kashani A, Milanese A, Fleck JS, Voigt AY,
Palleja A, Ponnudurai R, Sunagawa S, Coelho LP, Schrotz-King P, Vogtmann E,
Habermann N, Niméus E, Thomas AM, Manghi P, Gandini S, Serrano D, Mizutani S,
Shiroma H, Shiba S, Shibata T, Yachida S, Yamada T, Waldron L, Naccarati A,
Segata N, Sinha R, Ulrich CM, Brenner H, Arumugam M, Bork P, Zeller G_  
Nature Medicine, (__2019__) [Epub ahead of print]  
> In this publication, `SIAMCAT` is used extensively for holdout testing

If you used `SIAMCAT` in your publication,
[let us know](mailto:zeller@embl.de)!
