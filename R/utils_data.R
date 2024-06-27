#' @title Input the serological data
#' @name data_titre
#'
#' @description This data format represents the serological data used by the rjmcmc model
#'
#' @format A list with the following columns:
#' \describe{
#'   \item{id}{Unique identifier for each biomarker.}
#'   \item{time}{Description of the model using the \code{makeModel} and \code{addObservationalModel} functions.}
#'   \item{titre}{Description of the prior distribution using the \code{add_par_df} function.}
#'   \item{biomarker}{Description of the prior distribution using the \code{add_par_df} function.}
#' }
#' 
#' @details Add more data about this model here.
#'
#' @seealso Add related functions
#' 
#' @examples
#' # Example usage. For example we can make some summy data for IgG SARS-CoV-2
#' 
#' data_titre <- data.frame(
#'      id = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
#'     time = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'    titre = c(1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4),
#'      biomarker = rep("IgG", 10)
#' )
#'
#' 
#' @author David Hodgson
#' @keywords data format data_titre
#' @export
NULL


#' @title Input the known exposures
#' @name data_known_exposures
#'
#' @description This data format represents the serological data used by the rjmcmc model
#'
#' @format A list with the following columns:
#' \describe{
#'   \item{id}{Unique identifier for each biomarker.}
#'   \item{time}{Time at which the known exposure occurred}
#'   \item{exposure_type}{Key for the type of exposure experienced}
#' }
#' 
#' @details Add more data about this model here.
#'
#' @seealso Add related functions
#' 
#' @examples
#' # Example usage. For example if you have some information on known vaccinatoins, you can define:
#' 
#' data_titre <- data.frame(
#'      id = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
#'      time = c(10, 4, 100, 40, 20, 23, 26, 70, 40, 10),
#'      exposure_type = rep("vax", 10)
#' )
#'
#' 
#' @author David Hodgson
#' @keywords data format data_known_exposures
#' @export
NULL


#' @title Input the known exposures
#' @name data_known_exposures
#'
#' @description This data format represents the serological data used by the rjmcmc model
#'
#' @format A list with the following columns:
#' \describe{
#'   \item{id}{Unique identifier for each biomarker.}
#'   \item{time}{Time at which the known exposure occurred}
#'   \item{exposure_type}{Key for the type of exposure experienced}
#' }
#' 
#' @details Add more data about this model here.
#'
#' @seealso Add related functions
#' 
#' @examples
#' # Example usage. For example if you have some information on known vaccinatoins, you can define:
#' 
#' data_titre <- data.frame(
#'      id = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
#'      time = c(10, 4, 100, 40, 20, 23, 26, 70, 40, 10),
#'      exposure_type = rep("vax", 10)
#' )
#'
#' 
#' @author David Hodgson
#' @keywords data format data_known_exposures
#' @export
NULL





#' @title Define antibody kinetics model
#' @name abkineticsModel
#'
#' @description This data format represents information about the antibody kinetics
#'
#' @format A list with the following columns:
#' \describe{
#'   \item{names}{Unique identifier for each exposure type.}
#'   \item{model}{Description of the model using the \code{makeModel} and \code{addAbkineticsModel} functions.}
#'   \item{prior}{Description of the prior distribution using the \code{add_par_df} function.}
#' }
#' 
#' @details Add more data about this model here.
#'
#' @seealso \code{\link{makeModel}}, \code{\link{addAbkineticsModel}}, \code{\link{add_par_df}},  for related functions.
#' 
#' @examples
#' # Example usage. This describes the antibody kinetics for a SARS-CoV-2 delta wave using IgG data. First define the antibody kinetics function:
#' 
#' infSerumKinetics <- function(titre_est, timeSince, pars) {
#'    a <- pars[1] 
#'    b <- pars[2] 
#'    c <- pars[3]
#'    if (timeSince < 14) {
#'        titre_est <- titre_est + (log(exp(a) + exp(c)) * (timeSince) / 14);
#'    } else {
#'        titre_est <- titre_est + (log(exp(a) * exp(-b/10 * (timeSince - 14))) + exp(c));
#'    }
#'    titre_est
#'}
#' 
#' # Now define the antibody kinetics function in the format required for the rjmc package:
#' 
#' abkineticsModel <- list(
#'    names = c("delta"),
#'    model = makeModel(
#'            addAbkineticsModel("IgG", "delta", TRUE, c("a_d", "b_d", "c_d"), infSerumKinetics),
#'        ),
#'    prior = bind_rows(
#'        add_par_df("a_d", -2, 2, "norm",  0, 1), # ab kinetics
#'        add_par_df("b_d", 0, 1, "norm",  0.3, 0.05), # ab kinetics
#'        add_par_df("c_d", 0, 4, "unif", 0,  4), # ab kinetics
#'    )
#' 
#' 
#' @author David Hodgson
#' @keywords data format antibody kinetics model
NULL

#' @title Define observational model
#' @name observationalModel
#'
#' @description This data format represents information about the observation model
#'
#' @format A list with the following columns:
#' \describe{
#'   \item{names}{Unique identifier for each biomarker.}
#'   \item{model}{Description of the model using the \code{makeModel} and \code{addObservationalModel} functions.}
#'   \item{prior}{Description of the prior distribution using the \code{add_par_df} function.}
#' }
#' 
#' @details Add more data about this model here.
#'
#' @seealso \code{\link{makeModel}}, \code{\link{addObservationalModel}}, \code{\link{add_par_df}},  for related functions.
#' 
#' @examples
#' # Example usage. This describes the observation model for a SARS-CoV-2 delta wave using IgG data. First define the log likelihood function, which is cauchy, with a LOD at a titre value of log10(40):
#' 
#' obsFunction = function(ll, titre_val, titre_est, pars) {
#'    if (titre_val <= log10(40)) {
#'        ll <- ll + pcauchy(log10(40), titre_est, pars[1], log.p = TRUE)
#'    } else {
#'        ll <- ll + dcauchy(titre_val, titre_est, pars[1], log = TRUE)
#'    }
#'    ll
#'}
#' # Now define the observation model in the format required for the rjmc package:
#' 
#' observationalModel <- list(
#'     names = c("IgG"),
#'     model = makeModel(addObservationalModel("IgG", c("sigma"), obsFunction)),
#'     prior = add_par_df("sigma", 0.0001, 4, "unif", 0.0001, 4) # observational model,
#' )
#' 
#' 
#' @author David Hodgson
#' @keywords data format observation model
#' 
NULL

#' @title Define correlate of protection model
#' @name copModel
#'
#' @description This data format represents information about the model for the correlate of protection (probability of infection given exposure) given the titre of a biomarker for a predefined exposure type
#'
#' @format A list with the following columns:
#' \describe{
#'   \item{names}{Unique identifier for each biomarker.}
#'   \item{model}{Description of the model using the \code{makeModel} and \code{addObservationalModel} functions.}
#'   \item{prior}{Description of the prior distribution using the \code{add_par_df} function.}
#' }
#' 
#' @details Add more data about this model here.
#'
#' @seealso \code{\link{makeModel}}, \code{\link{addObservationalModel}}, \code{\link{add_par_df}},  for related functions.
#' 
#' @examples
#' # Example usage. This describes the observation model for a SARS-CoV-2 delta wave using IgG data. First define the COP log likelihood function, for example if the COP is a logistic function of the titre of the biomarker IgG (see \ref{obsLogLikelihood} for format):
#' 
#' copLogLikelihood <- function(ll, inf_status, esttitreExp, params) {
#'    # COP parameters
#'    beta0 <- params[1]
#'    beta1 <- params[2]
#'    p <- 1.0 / (1.0 + exp(- (beta0 + beta1 * esttitreExp) ) )
#'    ll <- inf_status * log(p) + (1 - inf_status) * log(1 - p)
#'    ll
#'}
#'
#' 
#' @author David Hodgson
#' @keywords data format COP model
#' @export
NULL