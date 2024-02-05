[![R-CMD-check](https://github.com/bhklab/consensusOV/workflows/R-CMD-check/badge.svg)](https://github.com/bhklab/consensusOV/actions)

**Bioc-Release**: ![Bioconductor RELEASE](http://bioconductor.org/shields/build/release/bioc/consensusOV.svg)

**Bioc-Devel**: ![Bioconductor DEVEL](http://bioconductor.org/shields/build/devel/bioc/consensusOV.svg)


Overview
--------

This package implements four major subtype classifiers for high-grade serous (HGS) ovarian cancer as described by Helland et al. (PLoS One, 2011), Bentink et al. (PLoS One, 2012), Verhaak et al. (J Clin Invest, 2013), and Konecny et al. (J Natl Cancer Inst, 2014). In addition, the package implements a consensus classifier, which consolidates and improves on the robustness of the proposed subtype classifiers, thereby providing reliable stratification of patients with HGS ovarian tumors of clearly defined subtype.

Author: Gregory M Chen, Lavanya Kannan, Ludwig Geistlinger, *Victor Kofia*, Levi Waldron, *Benjamin Haibe-Kains*

Installation
------------

``` r
# Installing the development version from GitHub:
# install.packages("devtools")
devtools::install_github("bhklab/consensusOV")

# Installing the release version from Bioconductor:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("consensusOV")
```

Usage
-----
Extensive usage examples are provided in the consensusOV vignette on Bioconductor: [consensusOV.pdf](https://bioconductor.org/packages/release/bioc/manuals/consensusOV/man/consensusOV.pdf)

Getting help
------------

Contact us by filing an issue in the consensusOV [issues](https://github.com/bhklab/consensusOV/issues) page.
