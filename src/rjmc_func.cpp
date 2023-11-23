#include <Rcpp.h>
#include <RcppEigen.h>

#include "./headers/mvn.hpp"
#include "./headers/rjmc.hpp"
#include "./headers/rjmc_sero.hpp"

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins("cpp14")]]


// [[Rcpp::export]]
List run_rjmc(Rcpp::List model, Rcpp::RObject dataList, Rcpp::List settings, bool update_ind, Rcpp::List RJMCpar, int i)
{
  rjmc::RJMC_D RJMC; List output_full;
  MatrixXd output;
  MatrixXd jump;
  rjmc::init_samplePriorDistributions(&RJMC, model["samplePriorDistributions"]);
  rjmc::init_evaluateLogPrior(&RJMC, model["evaluateLogPrior"]);
  rjmc::init_evaluateLogLikelihood(&RJMC, model["evaluateLogLikelihood"]);
  rjmc::init_jumpSampling(&RJMC, model["jumpSampling"]);
  rjmc::init_initialiseJump(&RJMC, model["initialiseJump"]);
  rjmc::init_evaluateAcceptanceRatioR(&RJMC, model["evaluateAcceptanceRatioR"]);

  if (update_ind) {
    RJMC.updateClass(settings, dataList, RJMCpar);
  } else {  
    RJMC.initialiseClass(settings, dataList, i);
  }

  output_full = RJMC.runRJMCC();

  RJMCpar = RJMC.saveRJMCpar();
  output = output_full[0];
  jump = output_full[1];
  return Rcpp::List::create(_["output"] = output, _["jump"] = jump, _["RJMCpar"] = RJMCpar);

}

// [[Rcpp::export]]
List run_rjmc_sero(Rcpp::List model, Rcpp::RObject dataList, Rcpp::List settings, bool update_ind, Rcpp::List RJMCpar, int i)
{
  rjmc_sero::RJMC_SERO_D RJMC_SERO; List output_full;
  MatrixXd output;
  MatrixXd jump;
  rjmc_sero::init_samplePriorDistributions(&RJMC_SERO, model["samplePriorDistributions"]);
  rjmc_sero::init_evaluateLogPrior(&RJMC_SERO, model["evaluateLogPrior"]);
  rjmc_sero::init_initialiseJump(&RJMC_SERO, model["initialiseJump"]);
  rjmc_sero::init_evaluateLogLikelihood(&RJMC_SERO, model["evaluateLogLikelihood"]);

  rjmc_sero::init_exposureFunctionSample(&RJMC_SERO, model["exposureFunctionSample"]);
  rjmc_sero::init_exposureFunctionDensity(&RJMC_SERO, model["exposureFunctionDensity"]);


  if (update_ind) {
    RJMC_SERO.updateClass(settings, dataList, RJMCpar);
  } else {  
    RJMC_SERO.initialiseClass(settings, dataList, i);
  }

  output_full = RJMC_SERO.runRJMCC();

  RJMCpar = RJMC_SERO.saveRJMCpar();
  output = output_full[0];
  jump = output_full[1];
  return Rcpp::List::create(_["output"] = output, _["jump"] = jump, _["RJMCpar"] = RJMCpar);

}
