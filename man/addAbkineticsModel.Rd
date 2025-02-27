% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_model.R
\name{addAbkineticsModel}
\alias{addAbkineticsModel}
\title{addAbkineticsModel}
\usage{
addAbkineticsModel(id, biomarker, exposureType, pars, funcForm)
}
\arguments{
\item{id}{The name of the biomarker.}

\item{biomarker}{The name of the biomarker.}

\item{exposureType}{The name of the exposure type.}

\item{pars}{The parameters of the model.}

\item{funcForm}{The antibody kinetics function.}
}
\value{
A list with the biomarker name, the exposure name, whether the exposure is inferred, the parameters and the antibody kinetics function.
}
\description{
This function adds an antibody kinetics model to the model definition.
}
