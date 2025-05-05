
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

### GitHub development version

A development version can be installed from [GitHub
repository](https://github.com/soob-kim/fmds) by entering:

``` r
devtools::install_github("soob-kim/FinfoMDS")
```

### Bioconductor official release

In the future, the official released version can be installed from
Bioconductor by entering:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("FinfoMDS")
```

## Example

Here we illustrate steps how you can use `FinfoMDS` and obtain 2D
representation of your microbiome dataset. As an example, we take an
algal-associated bacterial community (Kim et al., 2022). First, load the
data by typing

``` r
data("microbiome", package = "FinfoMDS")
```

Next, compute the *F*-informed MDS by running:

``` r
result <- fmds(lambda = 0.5, D=microbiome$dist, y=microbiome$host)
```

.. where it will iterate for 100 times by default, or until the 2D
distributions converge—whichever occurs first. While we have observed
that setting `lambda = 0.5` typically yields optimal results, this
hyperparameter can be adjusted as long as it does not exceed 1.

Next, compute the *F*-informed MDS by running

- `D`: Pairwise distance, matrix

- y: Label or group set, vector
