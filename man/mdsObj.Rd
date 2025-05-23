% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fmds.R
\name{mdsObj}
\alias{mdsObj}
\title{Objective term of MDS}
\usage{
mdsObj(D, z)
}
\arguments{
\item{D}{Original distance matrix}

\item{z}{Lower dimension representation}
}
\value{
Scalar of objective function value of MDS
}
\description{
Objective term of MDS
}
\examples{
z <- cmdscale(d = microbiome$dist)
mdsObj(D = microbiome$dist, z = z)
}
