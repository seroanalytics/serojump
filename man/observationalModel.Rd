% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_data.R
\name{observationalModel}
\alias{observationalModel}
\title{Define observational model}
\format{
A list with the following columns:
\describe{
\item{names}{Unique identifier for each biomarker.}
\item{model}{Description of the model using the \code{makeModel} and \code{addObservationalModel} functions.}
\item{prior}{Description of the prior distribution using the \code{addPrior} function.}
}
}
\description{
This data format represents information about the observation model
}
\details{
Add more data about this model here.
}
\examples{
# Example usage. This describes the observation model for a SARS-CoV-2 delta wave using IgG data. First define the log likelihood function, which is cauchy, with a LOD at a titre value of log10(40):

obsFunction = function(ll, titre_val, titre_est, pars) {
   if (titre_val <= log10(40)) {
       ll <- ll + pcauchy(log10(40), titre_est, pars[1], log.p = TRUE)
   } else {
       ll <- ll + dcauchy(titre_val, titre_est, pars[1], log = TRUE)
   }
   ll
}
# Now define the observation model in the format required for the rjmc package:

observationalModel <- list(
    names = c("IgG"),
    model = makeModel(addObservationalModel("IgG", c("sigma"), obsFunction)),
    prior = addPrior("sigma", 0.0001, 4, "unif", 0.0001, 4) # observational model,
)


}
\seealso{
\code{\link{makeModel}}, \code{\link{addObservationalModel}}, \code{\link{addPrior}},  for related functions.
}
\author{
David Hodgson
}
\keyword{data}
\keyword{format}
\keyword{model}
\keyword{observation}
