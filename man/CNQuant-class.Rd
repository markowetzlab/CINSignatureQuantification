% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{CNQuant-class}
\alias{CNQuant-class}
\alias{CNQuant}
\title{CNQuant class}
\description{
The \code{CNQuant} class is a structured S4 class object designed to contain copy
number and copy number signature-related data. This is initial class prior to
computation of copy number signatures which extends the \code{CNquant} class to
\code{SigQuant} class.
}
\section{Slots}{

\describe{
\item{\code{segments}}{Object of type \code{"list"} containing data.frame objects for each
sample segment table. Access with \code{\link[=getSegments]{getSegments()}}.}

\item{\code{featData}}{Object of type \code{"list"} containing data.frame objects for each
copy number feature extracted by \code{\link[=calculateFeatures]{calculateFeatures()}}. Access with
\code{\link[=getFeatures]{getFeatures()}}.}

\item{\code{featFitting}}{Object of type \code{"list"} containing objects feature
component fitting, including the sample by component matrix and model
parameters, calculated by \code{\link[=calculateSampleByComponentMatrix]{calculateSampleByComponentMatrix()}}. Access with
\code{\link[=getSampleByComponent]{getSampleByComponent()}}.}

\item{\code{samplefeatData}}{Object of type \code{"data.frame"} containing sample-level
information related to samples provided on \code{CNQuant} class initialisation,
including segment counts and sample ploidy. Access with
\code{\link[=getSamplefeatures]{getSamplefeatures()}} and add additional sample features using \link{addsampleFeatures}}

\item{\code{ExpData}}{An object of class \code{"ExpQuant"}, see \linkS4class{ExpQuant}.}
}}

\seealso{
\linkS4class{SigQuant}

\code{\link[=createCNQuant]{createCNQuant()}}
}
