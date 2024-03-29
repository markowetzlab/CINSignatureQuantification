% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/calculateFeatures.R
\docType{methods}
\name{calculateFeatures}
\alias{calculateFeatures}
\alias{calculateFeatures,CNQuant-method}
\title{calculateFeatures}
\usage{
calculateFeatures(
  object,
  method = "drews",
  smooth.diploid = TRUE,
  cores = 1,
  DCIN = 20
)

\S4method{calculateFeatures}{CNQuant}(
  object,
  method = NULL,
  smooth.diploid = TRUE,
  cores = 1,
  DCIN = 20
)
}
\arguments{
\item{object}{CNQuant object}

\item{method}{Method to extract copy number features. Default is "drews".}

\item{smooth.diploid}{Binary variable indicating whether segments close to 2
should be collapsed to 2 and merged together. Default is TRUE.}

\item{cores}{Number of CPU threads/cores to utilise via doParallel. Default
is 1. Maximum number is equal to the number of features to extract (drews &
mac methods = 6 features).}

\item{DCIN}{Threshold for required number of non-diploid segments to compute
copy number features (and subsequently copy number signatures) using method
"drews". Default is 20. This parameter should not need to be changed and
will affect feature values and signature activity.}
}
\value{
A CNQuant class object with extracted features stored in the
"featData" slot
}
\description{
Calculates and returns copy number features from copy number profiles in a
\code{CNQuant} class object. The output from this function is stored within the
\code{featData} slot in the returned \code{CNQuant} class object.
}
\examples{
  data(TCGA_478_Samples_SNP6_GOLD)
  cnobj <- createCNQuant(TCGA_478_Samples_SNP6_GOLD)
  cnobj <- calculateFeatures(cnobj[1:10],method="drews")
}
\seealso{
\code{\link[=getFeatures]{getFeatures()}}
}
