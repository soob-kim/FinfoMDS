#' Microbiome dataset
#'
#' The dataset consists of 36 balanced microbial community samples cultured with or without the microbial host (alga Phaeodactylum tricornutum). The community comprised 72 bacterial taxa as identified by amplicon sequence variants (ASV) of the 16S rRNA gene, PCR-amplified. ASV counts were obtained through Illumina MiSeq system and subsequently converted into relative abundances using cumulative sum scaling (CSS). A phylogeny-based metric (weighted Unifrac) and compositional data were used to calculate distance matrix.
#'
#' @format A list of length 2. Each element corresponds to:
#' \describe{
#'   \item{dist}{A numeric matrix representing pairwise weighted Unifrac distances between samples at the site.}
#'   \item{host}{A character vector indicating the presence ("+") or absence ("-") of the host for each sample.}
#' }
#' @source Data derived from an experimental study on host-microbe interactions.
#' @references
#' Kim H., Kimbrel J.A., Vaiana C.A., Wollard J.R., Mayali X., Buie C.R. (2022).
#' Bacterial response to spatial gradients of algal-derived nutrients in a porous microplate.
#' \emph{The ISME Journal}, 16(4), 1036–1045.
#' \doi{10.1038/s41396-021-01163-4}
#' @usage data(microbiome)
"microbiome"
