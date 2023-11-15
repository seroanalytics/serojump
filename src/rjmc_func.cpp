#include <Rcpp.h>
#include <RcppEigen.h>

#include "./headers/mvn.hpp"
#include "./headers/rjmc.hpp"

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins("cpp14")]]

////////////////////////////////////////
////////////////////////////////////////
//////////// CONTINUOUS PTMC //////////
///////////////////////////////////
///////////////////////////////////

/*
void init_samplePriorDistributions(RPTMC* model, Rcpp::Function samplePriorDistributions) {
  auto func = [samplePriorDistributions]() {
    PutRNGstate();
    auto rData = samplePriorDistributions();
    GetRNGstate();
    return Rcpp::as<VectorXd>(rData);
  };
  model->samplePriorDistributions = func;
}

void init_evaluateLogPrior(RPTMC* model, Rcpp::Function evaluateLogPrior) {
  auto func = [evaluateLogPrior](VectorXd params) {
    PutRNGstate();
    auto rData = evaluateLogPrior(params);
    GetRNGstate();
    return Rcpp::as<double>(rData);
  };
  model->evaluateLogPrior = func;
}

void init_evaluateLogLikelihood(RPTMC* model, Rcpp::Function evaluateLogLikelihood) {
  auto func = [evaluateLogLikelihood](VectorXd params, MatrixXd covariance, RObject dataList) {
    PutRNGstate();
    auto rData = evaluateLogLikelihood(params, covariance, dataList);
    GetRNGstate();
    return Rcpp::as<double>(rData);
  };
  model->evaluateLogLikelihood = func;
}
*/


////////////////////////////////////////
////////////////////////////////////////
//////////// DISCRETE PTMC //////////
///////////////////////////////////
///////////////////////////////////
/*
void init_samplePriorDistributions_discrete(ptmc_discrete::PTMC_D* model, Rcpp::Function samplePriorDistributions) {
  auto func = [samplePriorDistributions](S4 dataList) {
    PutRNGstate();
    auto rData = samplePriorDistributions(dataList);
    GetRNGstate();
    return Rcpp::as<VectorXd>(rData);
  };
  model->samplePriorDistributions = func;
}

void init_evaluateLogPrior_discrete(ptmc_discrete::PTMC_D* model, Rcpp::Function evaluateLogPrior) {
  auto func = [evaluateLogPrior](VectorXd params, VectorXi discrete, S4 dataList) {
    PutRNGstate();
    auto rData = evaluateLogPrior(params, discrete, dataList);
    GetRNGstate();
    return Rcpp::as<double>(rData);
  };
  model->evaluateLogPrior = func;
}

void init_evaluateLogLikelihood_discrete(ptmc_discrete::PTMC_D* model, Rcpp::Function evaluateLogLikelihood) {
  auto func = [evaluateLogLikelihood](VectorXd params, VectorXi discrete, MatrixXd covariance, S4 dataList) {
    PutRNGstate();
    auto rData = evaluateLogLikelihood(params, discrete, covariance, dataList);
    GetRNGstate();
    return Rcpp::as<double>(rData);
  };
  model->evaluateLogLikelihood = func;
}

void init_initialiseDiscrete_discrete(ptmc_discrete::PTMC_D* model, Rcpp::Function initialiseDiscrete) {
  auto func = [initialiseDiscrete](S4 dataList) {
    PutRNGstate();
    auto rData = initialiseDiscrete(dataList);
    GetRNGstate();
    return Rcpp::as<VectorXi>(rData);
  };
  model->initialiseDiscrete = func;
}

void init_discreteSampling_discrete(ptmc_discrete::PTMC_D* model, Rcpp::Function discreteSampling) {
  auto func = [discreteSampling](VectorXi discrete, S4 dataList) {
    PutRNGstate();
    auto rData = discreteSampling(discrete, dataList);
    GetRNGstate();
    return Rcpp::as<VectorXi>(rData);
  };
  model->discreteSampling = func;
}
*/
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
