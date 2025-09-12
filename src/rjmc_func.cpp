
#include "./headers/rjmc_inf.hpp"

void set_rng_seed(unsigned int seed) {
    rng.seed(seed);
}

// [[Rcpp::export]]
List run_rjmc_sero(Rcpp::List model, Rcpp::RObject dataList, Rcpp::List settings, bool update_ind, Rcpp::List RJMCpar, int i, int seed = -1)
{

  // Define a random number generator
  if (seed > -1) {
    rng.seed(seed);
  }
  double a = uniformContinuousDist(0, 1);

  List exposureInfo = model["infoModel"];
  List observationalModel = model["observationalModel"];
  List abkineticsModel = model["abkineticsModel"];
  List copModel = model["copModel"];


  auto SeroJumpRunInst = SeroJumpRun::create(exposureInfo, observationalModel, abkineticsModel, copModel);

  List output_full;
  MatrixXd output, jump, inf;
  RObject titreexp, obstitre, obsloglik;

  init_samplePriorDistributions(SeroJumpRunInst.get(), model["samplePriorDistributions"]);
  init_evaluateLogPrior(SeroJumpRunInst.get(), model["evaluateLogPrior"]);
  init_evaluateLogPriorInfExp(SeroJumpRunInst.get(), model["evaluateLogPriorInfExp"]);
  init_initialiseJump(SeroJumpRunInst.get(), model["initialiseJump"]);
  init_exposureFunctionSample(SeroJumpRunInst.get(), model["exposureFunctionSample"]);
  init_exposureFunctionDensity(SeroJumpRunInst.get(), model["exposureFunctionDensity"]);

  
  output_full = SeroJumpRunInst->runRJMCC(settings, dataList, i);
  RJMCpar = SeroJumpRunInst->saveRJMCpar();
  output = output_full[0];
  jump = output_full[1];
  titreexp = output_full[2];
  obstitre = output_full[3];
  obsloglik = output_full[4];


  return Rcpp::List::create(_["output"] = output, _["jump"] = jump, _["titreexp"] = titreexp, 
      _["obstitre"] = obstitre, _["obsloglik"] = obsloglik, _["RJMCpar"] = RJMCpar);

}