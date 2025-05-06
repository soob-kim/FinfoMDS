
# *FinfoMDS*: Multidimensional scaling informed by *F*-statistic

*F*-informed MDS is a new multidimensional scaling-based ordination
method that configures data distribution based on the *F*-statistic
(i.e., the ratio of dispersion between groups with shared or differing
labels). An R package, `FinfoMDS`, for computing the *F*-informed MDS is
currently being incorporated into Bioconductor. A preprint describing
the method in full is available at:

- H Kim⋆, S Kim⋆, JA Kimbrel, MM Morris, X Mayali and CR Buie (2025).
  Multidimensional scaling informed by *F*-statistic: Visualizing
  grouped microbiome data with inference, *arXiv*.
  (<https://arxiv.org/abs/2308.00354v2>).

## Installation

### GitHub

A development version can be installed from [GitHub
repository](https://github.com/soob-kim/fmds) by entering:

``` r
devtools::install_github("soob-kim/FinfoMDS")
```

### Bioconductor

In the future, the official released version can be installed from
Bioconductor by entering:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("FinfoMDS")
```

## Implementation

We outline steps for users to implement `FinfoMDS` package to a
microbiome dataset and obtain 2D representation of the microbiome. Let’s
take an algal-associated bacterial community for example (Kim et al.,
2022). First, load the data by typing

``` r
data("microbiome", package = "FinfoMDS")
```

Next, compute the *F*-informed MDS by running:

``` r
result <- fmds(lambda = 0.3, D = microbiome$dist, y = microbiome$host)
```

.. where it will iterate until the 2D distributions converge or for 100
times by default—whichever occurs first. While we have observed that
setting `lambda` between 0.3 and 0.5 typically yields optimal results,
this hyperparameter can be adjusted as long as it does not exceed 1.

Finally, the 2D representation of the community dataset is returned as a
matrix `result` and can be viewed by typing:

``` r
plot(result, pch=microbiome$host)
```
