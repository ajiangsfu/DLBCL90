DLBCL90: R Package for Molecular Classification of DLBCL
================

## Overview

**DLBCL90** is an R package that performs robust molecular
classification of diffuse large B-cell lymphoma (DLBCL) samples using
NanoString gene expression data. The primary function for users is
`DLBCL90calls()`, which supports input from either:

- NanoString RCC files (packaged as a `.zip`/`.ZIP` archive)
- Raw expression matrices in `.csv` format (directly extracted from
  RCCs, no preprocessing)

The `DLBCL90calls()` function produces three **separate** classification
results for each sample, based on:

- **Cell-of-Origin (COO)** classifier [(Wright et al.,
  2003)](#references)
- **Double-hit signature (DHITsig / DZsig)** classifier [(Ennishi et
  al., 2019)](#references)
- **PMBL-like** classifier [(Sha et al., 2021)](#references)

Each classification is executed independently and returned in the final
output.

> **Note:** This package is specifically designed for the *original
> DLBCL90 NanoString codeset*. If you are using a different or
> customized codeset, please contact <aijiang@bccrc.ca> for a calibrated
> version.

## Installation

This package is not yet on CRAN. You can install it from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("ajiangsfu/DLBCL90")
```

Note: Make sure the dependency **nanostringr** is installed beforehand:

``` r
devtools::install_github("NCICB/nanostringr")
```

## Example Usage

You can use either a `.zip` file containing NanoString RCC files, or a
`.csv` file with raw expression values:

``` r
library(DLBCL90)

# Example 1: Using a ZIP or zip file containing RCC data
zip_file <- "path/to/your_data.zip"
results_from_zip <- DLBCL90calls(zip_file)

# Example 2: Using a CSV file with expression data
csv_file <- "path/to/your_expression_data.csv"
results_from_csv <- DLBCL90calls(csv_file)

# View the classification results (for example from the ZIP file)
head(results_from_zip)

# If you like to modify default nHeading and geomeanCut parameter in the function, please feel free to do so
```

## Authors

This package was developed by:

- Aixiang Jiang, Stacy Hung, George Wright, David Scott

## References

- Wright G, Tan B, Rosenwald A, Hurt EH, Wiestner A, Staudt LM. (2003).
  *A gene expression-based method to diagnose clinically distinct
  subgroups of diffuse large B cell lymphoma*. **Proc Natl Acad Sci U S
  A**, 100(17):9991–6. [PMID:
  12900505](https://pubmed.ncbi.nlm.nih.gov/12900505/)

- Ennishi D, Jiang A, Boyle M, Collinge B, Grande BM, Ben-Neriah S,
  Rushton C, Tang J, Thomas N, Slack GW, Steidl C, Gascoyne RD, Mungall
  AJ, Marra MA, Kridel R. (2019). *Double-Hit Gene Expression Signature
  Defines a Distinct Subgroup of Germinal Center B-Cell-Like Diffuse
  Large B-Cell Lymphoma*. **Cell Reports**, 25(5):1446–1456.e6. [PMID:
  30523716](https://pubmed.ncbi.nlm.nih.gov/30523716/)

- Sha C, Barrans S, Cucco F, Bentley MA, Care MA, Cummin T, Kennedy H,
  Thompson JS, Uddin R, Wang J, Ahmad S, Worrillow L, Wright DW,
  Rosenwald A, Ott G, Campo E, Gascoyne RD, Staudt LM, Alizadeh AA, Dave
  SS, Scott DW, Smeland EB, Rosenquist R, Youngson JH, Pileri SA, Kluin
  PM, Holte H, Chan WC, Jaffe ES, et al. (2021). *Molecular High-Grade
  B-cell Lymphoma: Defining a Poor-Risk Group That Requires Different
  Approaches to Therapy*. **J Clin Oncol**, 39(8):802–813. [PMID:
  33684939](https://pubmed.ncbi.nlm.nih.gov/33684939/)
