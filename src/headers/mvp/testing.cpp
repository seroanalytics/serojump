
#include "./rjmc_inf.hpp"
//#include "./headers/rjmc_full.hpp"


// [[Rcpp::export]]
void runSeroJump() {
    // Dummy Rcpp::List to simulate model
    Rcpp::List model = Rcpp::List::create(
        Rcpp::Named("infoModel") = Rcpp::List::create(),
        Rcpp::Named("observationalModel") = Rcpp::List::create(),
        Rcpp::Named("abkineticsModel") = Rcpp::List::create()
    );

    Rcpp::List exposureInfo = model["infoModel"];
    Rcpp::List observationalModel = model["observationalModel"];
    Rcpp::List abkineticsModel = model["abkineticsModel"];

    Rcpp::Rcout << "Initiate class" << std::endl;
    auto SeroJumpRunInst = SeroJumpRun::create(exposureInfo, observationalModel, abkineticsModel);
    Rcpp::Rcout << "Class run" << std::endl;
}

// Include Rcpp interface to be recognized
RCPP_MODULE(SeroJumpModule) {
    Rcpp::function("runSeroJump", &runSeroJump);
}