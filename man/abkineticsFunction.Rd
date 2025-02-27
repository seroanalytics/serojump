% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_model.R
\name{abkineticsFunction}
\alias{abkineticsFunction}
\title{Define the function form of th antibody kinetics}
\arguments{
\item{titre_est}{The model-estimated titre value at the observational time}

\item{inf_status}{How long since the last exposur event occured}

\item{pars}{The fitted parameters needed to calculate the estimate titre value. These are defined in the \code{prior} entry of the abkineticsModel list.}
}
\value{
A function that returns estimate titre value
}
\description{
This function is defined by the user and describes antibody kinetics given the model.
}
\details{
Add information here.
}
\examples{
# Example usage: we define two functional form examples, one describing waning until infection, and then the other describing kinetics post infection. 
noInfSerumKinetics <- function(titre_est, timeSince, pars) {
   titre_est <- titre_est - pars[1] * (timeSince)
   titre_est
}

infSerumKinetics <- function(titre_est, timeSince, pars) {
   a <- pars[1] 
   b <- pars[2] 
   c <- pars[3]
   if (timeSince < 14) {
       titre_est <- titre_est + (log(exp(a) + exp(c)) * (timeSince) / 14) * 1;
   } else {
       titre_est <- titre_est + (log(exp(a) * exp(-b/10 * (timeSince - 14))) + exp(c)) * 1;
   }
   titre_est
}

}
\author{
Your Name
}
\keyword{data}
\keyword{export}
\keyword{format}
\keyword{function}
