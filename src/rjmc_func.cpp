
#include "./headers/rjmc_inf.hpp"
//#include "./headers/rjmc_full.hpp"
//#include "./headers/rjmc_pp.hpp"



/*List run_rjmc_pp(Rcpp::List model, Rcpp::RObject dataList, Rcpp::List settings, bool update_ind, Rcpp::List RJMCpar, int i)
{
  Rcpp::Rcout << "Start: run_rjmc_full" << std::endl;

  Rcpp::Rcout << "Initiate class" << std::endl;
  rjmc_pp::RJMC_PP_D RJMC_PP;
  Rcpp::Rcout << "End: Initiate class" << std::endl;

  List output_full;
  MatrixXd output;
  List jump;
  // Priors 
  Rcpp::Rcout << "Start: initiate functions" << std::endl;
  rjmc_pp::init_sampleInitPrior(&RJMC_PP, model["sampleInitPrior"]);
  rjmc_pp::init_sampleInitJump(&RJMC_PP, model["sampleInitJump"]);
  rjmc_pp::init_evaluateLogPrior(&RJMC_PP, model["evaluateLogPrior"]);
  rjmc_pp::init_evaluateLogLikelihood(&RJMC_PP, model["evaluateLogLikelihood"]);
  rjmc_pp::init_sampleBirthProposal(&RJMC_PP, model["sampleBirthProposal"]);
  rjmc_pp::init_sampleDeathProposal(&RJMC_PP, model["sampleDeathProposal"]);
  rjmc_pp::init_evaluateBirthProposal(&RJMC_PP, model["evaluateBirthProposal"]);
  rjmc_pp::init_evaluateDeathProposal(&RJMC_PP, model["evaluateDeathProposal"]);
  rjmc_pp::init_sampleJump(&RJMC_PP, model["sampleJump"]);
  rjmc_pp::init_sampleProposal(&RJMC_PP, model["sampleProposal"]);

  Rcpp::Rcout << "End: initiate functions" << std::endl;

  if (update_ind) {
    RJMC_PP.updateClass(settings, dataList, RJMCpar);
  } else {  
    RJMC_PP.initialiseClass(settings, dataList, i);
  }

  output_full = RJMC_PP.runRJMCC();

  RJMCpar = RJMC_PP.saveRJMCpar();
  output = output_full[0];
  jump = output_full[1];

  Rcpp::Rcout << "End: run_rjmc_pp" << std::endl;

  return Rcpp::List::create(_["output"] = output, _["jump"] = jump, _["RJMCpar"] = RJMCpar);

}*/


// [[Rcpp::export]]
List run_rjmc_sero(Rcpp::List model, Rcpp::RObject dataList, Rcpp::List settings, bool update_ind, Rcpp::List RJMCpar, int i)
{

  List exposureInfo = model["infoModel"];
  List observationalModel = model["observationalModel"];
  List abkineticsModel = model["abkineticsModel"];
  List copModel = model["copModel"];

  Rcpp::Rcout << "Initiate class" << std::endl;
  Rcpp::Rcout << "Initiate class2" << std::endl;
    Rcpp::Rcout << "Initiate class2: " << i << std::endl;

  auto SeroJumpRunInst = SeroJumpRun::create(exposureInfo, observationalModel, abkineticsModel, copModel);
  Rcpp::Rcout << "End: Initiate class" << std::endl;

  List output_full;
  MatrixXd output, jump, inf;
  RObject titreexp, obstitre;

  Rcpp::Rcout << "Start: initiate functions" << std::endl;
  init_samplePriorDistributions(SeroJumpRunInst.get(), model["samplePriorDistributions"]);
  init_evaluateLogPrior(SeroJumpRunInst.get(), model["evaluateLogPrior"]);
  init_evaluateLogPriorInfExp(SeroJumpRunInst.get(), model["evaluateLogPriorInfExp"]);
  init_initialiseJump(SeroJumpRunInst.get(), model["initialiseJump"]);
  init_exposureFunctionSample(SeroJumpRunInst.get(), model["exposureFunctionSample"]);
  init_exposureFunctionDensity(SeroJumpRunInst.get(), model["exposureFunctionDensity"]);
  Rcpp::Rcout << "End: initiate functions" << std::endl;

  
  output_full = SeroJumpRunInst->runRJMCC(settings, dataList, i);
  RJMCpar = SeroJumpRunInst->saveRJMCpar();
  output = output_full[0];
  jump = output_full[1];
  titreexp = output_full[2];
  obstitre = output_full[3];

  Rcpp::Rcout << "End: run_rjmc_sero" << std::endl;

  return Rcpp::List::create(_["output"] = output, _["jump"] = jump, _["titreexp"] = titreexp, _["obstitre"] = obstitre, _["RJMCpar"] = RJMCpar);

}


/*List run_rjmc_full(Rcpp::List model, Rcpp::RObject dataList, Rcpp::List settings, bool update_ind, Rcpp::List RJMCpar, int i)
{
  Rcpp::Rcout << "Start: run_rjmc_full" << std::endl;

  List exposureInfo = model["infoModel"];
  List observationalModel = model["observationalModel"];
  List abkineticsModel = model["abkineticsModel"];
  List copModel = model["copModel"];

  Rcpp::Rcout << "Initiate class" << std::endl;
  rjmc_full::RJMC_FULL_D RJMC_FULL(exposureInfo, observationalModel, abkineticsModel, copModel);
  Rcpp::Rcout << "End: Initiate class" << std::endl;

  List output_full;
  MatrixXd output, jump, inf;
  RObject titreexp, obstitre;
  // Priors 
  Rcpp::Rcout << "Start: initiate functions" << std::endl;
  rjmc_full::init_samplePriorDistributions(&RJMC_FULL, model["samplePriorDistributions"]);
  rjmc_full::init_evaluateLogPrior(&RJMC_FULL, model["evaluateLogPrior"]);
  rjmc_full::init_initialiseJump(&RJMC_FULL, model["initialiseJump"]);
  rjmc_full::init_exposureFunctionSample(&RJMC_FULL, model["exposureFunctionSample"]);
  rjmc_full::init_exposureFunctionDensity(&RJMC_FULL, model["exposureFunctionDensity"]);
  Rcpp::Rcout << "End: initiate functions" << std::endl;

  if (update_ind) {
    RJMC_FULL.updateClass(settings, dataList, RJMCpar);
  } else {  
    RJMC_FULL.initialiseClass(settings, dataList, i);
  }

  output_full = RJMC_FULL.runRJMCC();

  RJMCpar = RJMC_FULL.saveRJMCpar();
  output = output_full[0];
  jump = output_full[1];
  inf = output_full[2];
  titreexp = output_full[3];
  obstitre = output_full[4];

  Rcpp::Rcout << "End: run_rjmc_full" << std::endl;

  return Rcpp::List::create(_["output"] = output, _["jump"] = jump, _["inf"] = inf, _["titreexp"] = titreexp, _["obstitre"] = obstitre, _["RJMCpar"] = RJMCpar);

}*/
