
# *FinfoMDS*: Multidimensional scaling informed by *F*-statistic

*F*-informed MDS is a new multidimensional scaling method that
configures data distribution based on the *F*-statistic (i.e., the ratio
of dispersion between groups with shared or differing labels). An R
package, `FinfoMDS`, for computing the *F*-informed MDS is currently
under review at Bioconductor
([link](https://github.com/Bioconductor/Contributions/issues/3811)). A
preprint describing the method in full is available at:

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
take an algal-associated bacterial community for example ([Kim et al.,
2022](https://doi.org/10.1038/s41396-021-01147-x)). First, load the data
by typing

``` r
data("microbiome", package = "FinfoMDS")
```

Next, compute the weighted UniFrac distance from this dataset and obtain
its label set:

``` r
D <- distance(microbiome, method = 'wunifrac') # requires phyloseq package
y <- sample_data(microbiome)$Treatment
```

Then, compute the *F*-informed MDS by running:

``` r
result <- fmds(D = D, y = y, lambda = 0.3, threshold_p = 0.05)
```

This procedure will iterate until the 2D distributions converge, as long
as the *p*-value does not deviate by more than `threshold_p`, or until
reaching the default maximum of 100 iterations, whichever occurs first.
While lambda between 0.3 and 0.5 has typically yielded optimal results,
it can be adjusted as long as it does not exceed 1.

The `fmds()` function returns a two-column matrix representing the
community dataset, which can be visualized by typing:

``` r
plot(result, pch=microbiome$host)
```
