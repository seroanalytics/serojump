#ifndef RJMC_F_HPP
#define RJMC_F_HPP

#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/random.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <vector>
#include <random>
#include <math.h>

//[[Rcpp::depends(RcppClock)]]
#include <RcppClock.h>

#define EIGEN_DONT_VECTORIZE
#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT


using namespace Rcpp;
using namespace std;
using namespace Eigen;
using namespace boost::math;
// [[Rcpp::plugins("cpp14")]]

#define PI 3.14159265358979323846

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b)) // define MAX function for use later
#endif

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b)) // define MAX function for use later
#endif



namespace rjmc_full{
    Rcpp::Clock clock;



    struct RJMC_FULL_D
    {
        struct abKineticsInfo {
            Function func;
            std::string name;
            NumericVector params; 

            abKineticsInfo(const std::string& name_i, Function func_i, NumericVector params_i) : name(name_i), func(func_i), params(params_i) {}
            // Define a comparison operator to sort by string values
        };

        Rcpp::List observationalList, abkineticsModelList, copModelList; // List of models and parameters
        Rcpp::List evalLoglikelhoodObs, evalabkineticsFunc, evalLoglikelhoodCOP; // List of functions to evaluate likelihoods and kinetics
        Rcpp::List observationalModel, abkineticsModel, copModel; // List of models and parameters

        StringVector observationalNames, abkineticNames, copNames; // Names of all the models
        List parsAbKinN, parsCOPN, parsObsN; // Lists of all the parameters names for each model
        std::vector<abKineticsInfo> abinfo;
        std::map<std::string, int> abmap;

        RJMC_FULL_D(Rcpp::List observationalModel_in, Rcpp::List abkineticsModel_in, Rcpp::List copModel_in) : observationalList(observationalModel_in), abkineticsModelList(abkineticsModel_in), copModelList(copModel_in) {

            this->observationalNames = observationalList["names"];
            this->abkineticNames = abkineticsModelList["names"];
            this->copNames = copModelList["names"];

            this->observationalModel = observationalList["model"];
            this->abkineticsModel = abkineticsModelList["model"];
            this->copModel = copModelList["model"];


            for (int i = 0; i < this->observationalNames.size(); i++) {
                List temp1 = this->observationalModel[i];
                Function temp2 = temp1["logLikelihood"];
                StringVector temp3 = temp1["pars"];
                this->evalLoglikelhoodObs[as<string>(this->observationalNames[i])] = temp2;
                this->parsObsN[as<string>(this->observationalNames[i])] = temp3;
            }
            for (int i = 0; i < this->abkineticNames.size(); i++) {
                List temp1 = this->abkineticsModel[i];
                Function temp2 = temp1["funcForm"];
                StringVector temp3 = temp1["pars"];
                this->evalabkineticsFunc[as<string>(this->abkineticNames[i])] = temp2;
                this->parsAbKinN[as<string>(this->abkineticNames[i])] = temp3;
            }
            for (int i = 0; i < this->copNames.size(); i++) {
                List temp1 = this->copModel[i];
                Function temp2 = temp1["logLikelihood"];
                StringVector temp3 = temp1["pars"];
                this->evalLoglikelhoodCOP[as<string>(this->copNames[i])] = temp2;
                this->parsCOPN[as<string>(this->copNames[i])] = temp3;
            }
        }


        struct DoubleWithString {
            double value;
            std::string name;

            DoubleWithString(const std::string& nm, double val) : value(val), name(nm) {}

            // Define a comparison operator to sort by string values
            bool operator<(const DoubleWithString& other) const {
                if (value != other.value) {
                    return value < other.value; // Sort by value
                } else {
                    // If values are equal, put "bleed" first
                    return (name != "bleed" && other.name == "bleed");
                }            
            }
        };

        List currentParsCOP, currentParsObs, currentParsAb;
        StringVector exposureNames;
        StringVector fittedParamNames;

        // call each step
        void updateParams(NumericVector paramsN) {
            // parsCOPN is a list of the parameter names for the COP model
            for (int i = 0; i < parsCOPN.size(); ++i) {
                NumericVector currentParsCOP_i;
                StringVector parnams = parsCOPN[i];
                for (int j = 0; j < parnams.size(); ++j) {
                    currentParsCOP_i.push_back(paramsN[as<string>(parnams[j])]);
                }
                this->currentParsCOP[as<string>(this->copNames[i])] = currentParsCOP_i;
            }
            for (int i = 0; i < parsObsN.size(); ++i) {
                NumericVector currentParsObs_i;
                StringVector parnams = parsObsN[i];
                for (int j = 0; j < parnams.size(); ++j) {
                    currentParsObs_i.push_back(paramsN[as<string>(parnams[j])]);
                }
                this->currentParsObs[as<string>(this->observationalNames[i])] = currentParsObs_i;

            }
            for (int i = 0; i < parsAbKinN.size(); ++i) {
                NumericVector currentParsAb_i;
                StringVector parnams = parsAbKinN[i];
                for (int j = 0; j < parnams.size(); ++j) {
                    currentParsAb_i.push_back(paramsN[as<string>(parnams[j])]);
                }
                this->currentParsAb[as<string>(this->abkineticNames[i])] = currentParsAb_i;

            }
        }

        bool conPropIn = true;
        bool disPropIn = false;

        VectorXd lowerParBounds, upperParBounds;
        double nonadaptiveScalar, adaptiveScalar;
        MatrixXd nonadaptiveCovarianceMat, adaptiveCovarianceMat;
        MatrixXd currentSample, currentSampleMean;

        // Outputs for posterior
        MatrixXd posteriorOut;
        MatrixXd posteriorJump;
        MatrixXd posteriorInf;
        MatrixXd posteriorTitreExp;
        MatrixXd posteriorObsTitre;

        MatrixXd currentCovarianceMatrix;
        double currentLogPosterior;
        VectorXd proposalSample;
        
        VectorXd currentJump;
        VectorXd proposalJump;

        VectorXd currentInf;
        VectorXd proposalInf;

        VectorXd currentTitreExp;
        VectorXd proposalTitreExp;

        VectorXd currentObsTitre;
        VectorXd proposalObsTitre;

        VectorXd historicJump;

        int iPosterior;
        int iterations, posteriorSamplesLength, thin, burninPosterior, burninAdaptiveCov, consoleUpdates, updatesAdaptiveCov, chainNumber;
        int lengthJumpVec;
        int numberFittedPar;
        int workingIteration;
        bool onDebug, onAdaptiveCov, profile;
        bool isSampleAccepted, isProposalAdaptive;

        RObject dataList;
        List dataListCPP;

        int counterFuncEval, counterAccepted, counterPosterior ,counterAdaptive;
        int counterNonAdaptive;
        double proposedLogPosterior, alpha, covarMaxVal, covarInitVal, covarInitValAdapt;

        std::function<VectorXd(RObject)> samplePriorDistributions;
        std::function<VectorXd(RObject)> initialiseJump;
        std::function<double(VectorXd, VectorXd, RObject)> evaluateLogPrior;

        std::function<double()> exposureFunctionSample;
        std::function<double(double)> exposureFunctionDensity;

        // Functions for internal 
        Mvn Mvn_sampler;
        vector<vector<DoubleWithString> > currentEventsFull, currentTitreFull;
        vector<vector<DoubleWithString> > proposalEventsFull, proposalTitreFull;

        void updateLists() {
            if (onDebug) Rcpp::Rcout << "In: updateLists" << std::endl;
            int i_idx;
            //vector<DoubleWithString> df_order_exp_i;
            for (int i = 0; i < this->diff_ind.size(); i++) {
                i_idx = this->diff_ind[i];
                std::vector<DoubleWithString> df_order_exp_i = this->sortevents(i_idx, this->proposalJump, this->proposalInf);
                this->proposalEventsFull[i_idx] = df_order_exp_i;
                this->proposalTitreFull[i_idx] = this->calculateTitre(df_order_exp_i, i_idx);
            }
            if (onDebug) Rcpp::Rcout << "End: updateLists" << std::endl;
        }

        double stepSizeRobbinsMonro;
        double evalLogPosterior(const VectorXd& param, const VectorXd& jump, const VectorXd& jumpInf, const MatrixXd& covariance, const RObject& dataList, bool init = false)
        {

            if (onDebug) Rcpp::Rcout << "In: evalLogPosterior" << std::endl;
            NumericVector paramsN = this->createNamedParam(param);
            updateParams(paramsN);
            define_abs();
            updateLists();

            // call this to update the events and titre space

            double logPrior = this->evaluateLogPrior(param, jump, dataList);

            if (isinf(logPrior)) {
                return log(0);
            }
            // need to be a converstion here as observationalModel is ObsFuncTemplate, but function call is Function

            double logLikelihood_ab = this->evaluateLogLikelihoodCOP_cpp(param, jump, jumpInf, init) + this->evaluateLogLikelihoodObs_cpp(param, jump, jumpInf, init);

            double logLikelihood_time = 0; 
            for (int i = 0; i < this->N; i++) {
                if (jump[i] > -1) {
                    logLikelihood_time += exposureFunctionDensity(jump[i]);
                }
            }
            
            // Prior distribution on number of infected persons given number of exposed people, this is the evlaution of a betabinomial(E, 1, 1) distribution 
            double logPriorInf;
            if (this->propInferredExpN == 0) {
                logPriorInf = -1000;
            }  else {
                logPriorInf = log(1.0 / this->propInferredExpN);
            }
            double logPriorExp = log(1.0 / this->N);
    

            return logPrior + logPriorExp + logPriorInf + logLikelihood_ab + logLikelihood_time;
        }
        
        // "A handy approximation for the error function and its inverse" by Sergei Winitzki.
        double ErfInv(float x){
            double tt1, tt2, lnx, sgn;
            sgn = (x < 0) ? -1.0 : 1.0;
            
            x = (1 - x)*(1 + x);
            lnx = logf(x);
            
            double pi = atan(1)*4;    

            tt1 = 2/(pi*0.147) + 0.5 * lnx;
            tt2 = 1/(0.147) * lnx;
            
            return(sgn*sqrtf(-tt1 + sqrtf(tt1*tt1 - tt2)));
        }

        // Important information from data
        int N, N_data;
        int knownInfsN;
        VectorXd knownExpVec;
        bool knownExpInd = false;
        VectorXd knownInfsVec, knownInfsTimeVec;
        int currInferredInfN;
        int propInferredInfN;
        int currInferredExpN;
        int propInferredExpN;
        int currJumpType = 1;
        int currJumpIdx;
        VectorXd counterAdaptiveGibbs;
        VectorXd adaptiveGibbsSD;
        int gibbsIdx = 0;
        int noGibbsSteps;
        VectorXd initialTitreTime;
        VectorXd initialTitreValue;
        VectorXd endTitreTime;
        VectorXd titre_full, times_full, id_full;
        List titre_list, times_list;

        StringVector exposurePeriods;
        string exposureNameInf;


         std::vector<DoubleWithString> sortevents(int i_idx, const VectorXd& jump_inf, const VectorXd& jump) {
            NumericVector titre_full_i = this->titre_list[i_idx];
            NumericVector times_full_i = this->times_list[i_idx];
            std::vector<DoubleWithString> df_order_exp;

            // Add all bleed times
            for (int i = 0; i < titre_full_i.size(); i++) {
                df_order_exp.emplace_back("bleed", times_full_i[i]); // End event
            }
            // Add all exposure events
            for (int i = 0; i < this->abkineticNames.size(); i++) {
                List here = this->abkineticsModel[as<string>(this->abkineticNames[i])]; 
                bool inferred = here["inferred"];
                if (!inferred){
                    NumericVector knowninf = here["known_inf"];
                    if (knowninf[i_idx] > -1) {
                        df_order_exp.emplace_back(as<string>(this->abkineticNames[i]), knowninf[i_idx]); // "pre-delta" event
                    }
                } else {
                    if (jump_inf[i_idx] == 1) {
                        df_order_exp.emplace_back(as<string>(this->abkineticNames[i]), jump[i_idx]); // "delta" event
                    }
                }
            }

            std::sort(df_order_exp.begin(), df_order_exp.end());
            return df_order_exp;

        }

        std::vector<DoubleWithString> calculateTitre(std::vector<DoubleWithString>& orderedEvents, double i_idx)  {
            if (onDebug) Rcpp::Rcout << "In: calculateTitre" << std::endl;

            double initTitre = this->initialTitreValue[i_idx];
            double initTime = this->initialTitreTime[i_idx];
            std::vector<DoubleWithString> df_order_titre;

            double time_since, titre_obs;
            double anchor_titre = initTitre;
            double curr_time = initTime;
            string anchor_func = orderedEvents[0].name;
            double anchor_time = curr_time;

            df_order_titre.emplace_back(anchor_func, anchor_titre);

            for (int i = 1; i < orderedEvents.size(); i++) {
                    time_since = orderedEvents[i].value - anchor_time;
                    titre_obs = calTitre(anchor_titre, anchor_func, time_since );
                    df_order_titre.emplace_back(orderedEvents[i].name, titre_obs);
                if (orderedEvents[i].name != "bleed") {
                    anchor_titre = titre_obs;
                    anchor_time = orderedEvents[i].value;
                    anchor_func = orderedEvents[i].name;
                }
            }
            return df_order_titre;

        }

        void initialiseClass(List settings, RObject dataList, int i)
        {
            Rcpp::Rcout << "In: initialiseClass" << std::endl;
            Rcpp::Rcout << "In: initialiseClass" << std::endl;
            if (onDebug) Rcpp::Rcout << "In: initialiseClass" << std::endl;

            if (onDebug) Rcpp::Rcout << "In: Extract data from dataList" << std::endl;

            this->dataList = dataList;            
            this->dataListCPP = as<List>(dataList);
            this->fittedParamNames = this->dataListCPP["par_names"];

            // Useful for internal functions
            this->N = this->dataListCPP["N"];
            this->N_data = this->dataListCPP["N_data"];

            this->titre_list = this->dataListCPP["titre_list"];
            this->times_list = this->dataListCPP["times_list"];
            
            this->titre_full = this->dataListCPP["titre_full"];
            this->times_full = this->dataListCPP["times_full"];
            this->id_full = this->dataListCPP["id_full"];

            this->knownExpVec = this->dataListCPP["knownExpVec"];
            if (this->knownExpVec.size() > 1) {
                this->knownExpInd = true;
            }
            this->knownInfsVec = this->dataListCPP["knownInfsVec"];
            this->knownInfsTimeVec = this->dataListCPP["knownInfsTimeVec"];
            this->knownInfsN = this->dataListCPP["knownInfsN"];

            this->currInferredInfN = this->knownInfsN;
            this->propInferredInfN = this->knownInfsN;
            if (this->knownExpInd) {
                this->currInferredExpN = this->knownExpVec.sum();
                this->propInferredExpN = this->knownExpVec.sum();
            } else {
                this->currInferredExpN = this->knownInfsN;
                this->propInferredExpN = this->knownInfsN;    
            }

            this->N_data = this->dataListCPP["N_data"];

            this->initialTitreTime = this->dataListCPP["initialTitreTime"];
            this->endTitreTime = this->dataListCPP["endTitreTime"];
            this->initialTitreValue = this->dataListCPP["initialTitreValue"];
            
            if (onDebug) Rcpp::Rcout << "Out: Extract data from dataList" << std::endl;

            if (onDebug) Rcpp::Rcout << "In: Initialise rjmcmc" << std::endl;


            this->exposurePeriods = this->abkineticsModel["exposurePeriods"];
            this->exposureNames = this->abkineticsModel["exposureNames"];
            this->exposureNameInf = as<string>(this->abkineticsModel["exposureNameInf"]);

            if (onDebug) Rcpp::Rcout << "In: Check 1" << std::endl;

            this->adaptiveGibbsSD = VectorXd::Constant(this->N, 0);
            this->counterAdaptiveGibbs = VectorXd::Zero(this->N);

            this->chainNumber = i;
            if (onDebug) Rcpp::Rcout << "In: Check 2" << std::endl;

            this->noGibbsSteps = settings["noGibbsSteps"];
            this->numberFittedPar = settings["numberFittedPar"];
            this->iterations = settings["iterations"];
            this->thin = settings["thin"];
            this->burninPosterior = settings["burninPosterior"];
            this->burninAdaptiveCov = settings["burninAdaptiveCov"];
            this->consoleUpdates = settings["consoleUpdates"];
            this->onAdaptiveCov = settings["onAdaptiveCov"];
            this->updatesAdaptiveCov = settings["updatesAdaptiveCov"];
            this->onDebug = settings["onDebug"];
            this->profile = settings["profile"];
            this->lengthJumpVec = settings["lengthJumpVec"];
            if (onDebug) Rcpp::Rcout << "In: Check 3" << std::endl;

            this->lowerParBounds = settings["lowerParBounds"];
            this->upperParBounds = settings["upperParBounds"];

            this->covarInitVal = settings["covarInitVal"];
            this->covarInitValAdapt = settings["covarInitValAdapt"];
            this->covarMaxVal = settings["covarMaxVal"];

            this->counterFuncEval = 0;
            this->counterAccepted = 0;
            this->counterPosterior = 0;
            this->counterAdaptive = 0;
            this->counterNonAdaptive = 0;
            
            this->iPosterior = 0;
            if (onDebug) Rcpp::Rcout << "In: Check 4" << std::endl;

            this->posteriorSamplesLength = (this->iterations-this->burninPosterior)/(this->thin);
            this->posteriorOut = MatrixXd::Zero(this->posteriorSamplesLength, this->numberFittedPar + 2);
            this->posteriorJump = MatrixXd::Zero(this->posteriorSamplesLength, this->lengthJumpVec);
            this->posteriorInf = MatrixXd::Zero(this->posteriorSamplesLength, this->lengthJumpVec);
            this->posteriorTitreExp = MatrixXd::Zero(this->posteriorSamplesLength, this->lengthJumpVec);
            this->posteriorObsTitre = MatrixXd::Zero(this->posteriorSamplesLength, this->lengthJumpVec);

            this->currentLogPosterior = 0;
            this->currentSampleMean = VectorXd::Zero(this->numberFittedPar);
            this->currentSample = VectorXd::Zero(this->numberFittedPar);

            this->proposalSample = VectorXd::Zero(this->numberFittedPar);

            this->nonadaptiveScalar = 0;
            this->adaptiveScalar = 0;
            
            this->currentCovarianceMatrix = MatrixXd::Zero(this->numberFittedPar, this->numberFittedPar );
            this->nonadaptiveCovarianceMat = MatrixXd::Zero(this->numberFittedPar, this->numberFittedPar );
            this->adaptiveCovarianceMat = MatrixXd::Zero(this->numberFittedPar ,this->numberFittedPar );
            if (onDebug) Rcpp::Rcout << "In: Check 5" << std::endl;

            VectorXd initialSample;
            VectorXd initialJump;
            VectorXd initialInf;

            double initialLogLikelihood;
            MatrixXd initialCovarianceMatrix;
            
            for(int parNum = 0; parNum < this->numberFittedPar ; parNum++){
                this->nonadaptiveCovarianceMat(parNum,parNum) = this->covarInitVal*(this->upperParBounds(parNum) - this->lowerParBounds(parNum));
                this->adaptiveCovarianceMat(parNum,parNum) = this->covarInitValAdapt*(this->upperParBounds(parNum) - this->lowerParBounds(parNum));
            }
            
            this->nonadaptiveScalar = log(0.1*0.1/(double)this->numberFittedPar);
            this->adaptiveScalar = log(2.382*2.382/(double)this->numberFittedPar);
            if (onDebug) Rcpp::Rcout << "In: Check 6" << std::endl;

            // Get initial conditions for the proposal sample
            initialSample = this->samplePriorDistributions(this->dataList);

            this->historicJump = VectorXd::Constant(this->N, -1);

            // Get initial conditions for the exposure states
            if (knownExpInd) {
                initialJump = this->knownExpVec;
            } else {
                initialJump = this->knownInfsTimeVec;
                /*for (int i = 0; i < this->N; i ++) {
                    if ((initialJump[i] > -1) & (this->knownInfsVec(i) != 1)) {
                        double u = uniformContinuousDist(0, 1);
                        if (u > 0.5) {
                            initialJump[i] = this->exposureFunctionSample();
                        }
                    }
                }*/ 
            }
            this->currentJump = initialJump;
            // Get initial conditions for the infection states
            initialInf = initInfFunc(); 
            this->currentInf = initialInf;
            if (onDebug) Rcpp::Rcout << "In: Check 7" << std::endl;


            NumericVector paramsN = this->createNamedParam(initialSample);
            updateParams(paramsN);
            define_abs();
            // Get initial conditions for the Loglikelihood probabilities
            std::vector<DoubleWithString> df_order_exp_i, df_order_titre_i;
            for (int i_idx = 0; i_idx < this->N; i_idx++) {

                df_order_exp_i = this->sortevents(i_idx, initialJump, initialInf);
                int l = df_order_exp_i.size();
                this->currentEventsFull.push_back(df_order_exp_i);

                if (onDebug) Rcpp::Rcout << "In: Check 7iii: " << i_idx << std::endl;

                df_order_titre_i = this->calculateTitre(df_order_exp_i, i_idx);
                if (onDebug) Rcpp::Rcout << "In: Check 7iv: " << i_idx << std::endl;
                //this->currentTitreFull.insert(this->currentTitreFull.begin() + i_idx, df_order_titre_i);
                this->currentTitreFull.push_back(df_order_titre_i);
                if (onDebug) Rcpp::Rcout << "In: Check 7iv: " << i_idx << std::endl;
                this->proposalEventsFull.push_back(df_order_exp_i);
                this->proposalTitreFull.push_back(df_order_titre_i);
            }
            this->currentCovarianceMatrix = this->nonadaptiveScalar*this->nonadaptiveCovarianceMat;
            initialLogLikelihood = this->evalLogPosterior(initialSample, initialJump, initialInf, this->currentCovarianceMatrix, this->dataList, true);

            if (onDebug) Rcpp::Rcout << "In: Check 8" << std::endl;

            while(isinf(initialLogLikelihood) || isnan(initialLogLikelihood)){
                initialSample = this->samplePriorDistributions(this->dataList);
                initialLogLikelihood = this->evalLogPosterior(initialSample, initialJump, initialInf, this->currentCovarianceMatrix, this->dataList, true);
            }
            if (onDebug) Rcpp::Rcout << "In: Check 9" << std::endl;
            this->currentSample = initialSample;
            this->currentSampleMean = initialSample;
            this->currentLogPosterior = initialLogLikelihood;

            ///this->proposalEventsFull = this->currentEventsFull;
            //this->proposalTitreFull = this->currentTitreFull;

            
            double alphaMVN = -sqrt(2)*ErfInv(0.234-1);
            this->stepSizeRobbinsMonro = (1.0-1.0/(double)this->numberFittedPar)*(pow(2*3.141, 0.5)*exp(alphaMVN*alphaMVN*0.5))/(2*alphaMVN) + 1.0/(this->numberFittedPar*0.234*(1-0.234));

            if (onDebug) Rcpp::Rcout << "End: Initialise rjmcmc" << std::endl;
            if (onDebug) Rcpp::Rcout << "End: initialiseClass" << std::endl;
        }

        VectorXd initInfFunc() {
            VectorXd initialInf = VectorXd::Zero(this->N);
            for (int i = 0; i < this->N; i++) {
                if ((this->currentJump(i) > -1) && (this->knownInfsVec(i) != 1)){
                    boost::random::beta_distribution<> b(1, 1); 
                    double p = b(rng);
                    boost::random::bernoulli_distribution<> l(p); 
                    initialInf(i) = l(rng);
                } else if ((this->currentJump(i) > -1) && (this->knownInfsVec(i) == 1)) {
                    initialInf(i) = 1;
                }
            }
            return initialInf;
        }

        void updateClass(List settings, RObject dataList, List PTMCpar)
        {
            this->dataList = dataList;


            this->numberFittedPar = settings["numberFittedPar"];
            this->iterations = settings["iterations"];
            this->thin = settings["thin"];
            this->burninPosterior = settings["burninPosterior"];
            this->burninAdaptiveCov = settings["burninAdaptiveCov"];
            this->consoleUpdates = settings["consoleUpdates"];
            this->onAdaptiveCov = settings["onAdaptiveCov"];
            this->updatesAdaptiveCov = settings["updatesAdaptiveCov"];
            this->onDebug = settings["onDebug"];
            
            this->lowerParBounds = settings["lowerParBounds"];
            this->upperParBounds = settings["upperParBounds"];

            // Counters 
            this->counterFuncEval = PTMCpar["counterFuncEval"];
            this->counterAccepted = PTMCpar["counterAccepted"];
            this->counterPosterior = PTMCpar["counterPosterior"];
            this->counterAdaptive = PTMCpar["counterAdaptive"];
            this->counterNonAdaptive = PTMCpar["counterNonAdaptive"];
            
            this->posteriorSamplesLength = (int)PTMCpar["posteriorSamplesLength"] + (this->iterations-this->burninPosterior)/(this->thin);
            
            this->posteriorOut = MatrixXd::Zero(this->posteriorSamplesLength, this->numberFittedPar+2);
            this->posteriorJump = MatrixXd::Zero(this->posteriorSamplesLength, this->lengthJumpVec);
            this->posteriorInf = MatrixXd::Zero(this->posteriorSamplesLength, this->lengthJumpVec);
            this->posteriorTitreExp = MatrixXd::Zero(this->posteriorSamplesLength, this->lengthJumpVec);
            this->posteriorObsTitre = MatrixXd::Zero(this->posteriorSamplesLength, this->lengthJumpVec);

            MatrixXd posteriorOutPrev = PTMCpar["posteriorOut"];
            for (int i = 0; i < (int)PTMCpar["posteriorSamplesLength"]; i++)
                this->posteriorOut.row(i) = posteriorOutPrev.row(i);

            this->currentLogPosterior = PTMCpar["currentLogPosterior"];
            this->currentSampleMean = PTMCpar["currentSampleMean"];
            this->currentSample = PTMCpar["currentSample"];
            
            this->proposalSample = PTMCpar["proposalSample"];
            this->nonadaptiveScalar = PTMCpar["nonadaptiveScalar"];
            this->adaptiveScalar = PTMCpar["adaptiveScalar"];
            
            this->currentCovarianceMatrix = PTMCpar["currentCovarianceMatrix"];
            this->nonadaptiveCovarianceMat =PTMCpar["nonadaptiveCovarianceMat"];
            this->adaptiveCovarianceMat = PTMCpar["adaptiveCovarianceMat"];
            
            double alphaMVN = -sqrt(2)*ErfInv(0.234-1);
            this->stepSizeRobbinsMonro = (1.0-1.0/(double)this->numberFittedPar)*(pow(2*3.141, 0.5)*exp(alphaMVN*alphaMVN*0.5))/(2*alphaMVN) + 1.0/(this->numberFittedPar*0.234*(1-0.234));
        }
        
        List saveRJMCpar() 
        {
            List RJMCpar = 
                Rcpp::List::create(
                    _["counterFuncEval"] = this->counterFuncEval,
                    _["counterAccepted"] = this->counterAccepted,
                    _["counterPosterior"] = this->counterPosterior,
                    _["counterAdaptive"] = this->counterAdaptive,
                    _["counterNonAdaptive"] = this->counterNonAdaptive,

                    _["currentLogPosterior"] = this->currentLogPosterior,
                    _["currentSampleMean"] = this->currentSampleMean,
                    _["currentSample"] = this->currentSample,

                    _["proposalSample"] = this->proposalSample,
                    _["nonadaptiveScalar"] = this->nonadaptiveScalar,
                    _["adaptiveScalar"] = this->adaptiveScalar,

                    _["currentCovarianceMatrix"] = this->currentCovarianceMatrix,
                    _["nonadaptiveCovarianceMat"] = this->nonadaptiveCovarianceMat,
                    _["adaptiveCovarianceMat"] = this->adaptiveCovarianceMat,

                    _["posteriorSamplesLength"] = this->posteriorSamplesLength,
                    _["posteriorOut"] = this->posteriorOut
            );
            return RJMCpar;
        }
        
        List runRJMCC()
        {
            if (onDebug) Rcpp::Rcout << "Pre: iterations" << std::endl;
            for (int i = 0; i < this->iterations; i++){

                this->workingIteration = i;
                if (onDebug) Rcpp::Rcout << "Pre: updateAllChains" << std::endl;
                if (this->profile) clock.tick("Astep");
                updateAllChains();
                if (this->profile) clock.tock("Astep");
            }

            List out = Rcpp::List::create(
                _["pars"] = this->posteriorOut,
                _["jump"] = this->posteriorJump,
                _["inf"] = this->posteriorInf,
                _["titreexp"] = this->posteriorTitreExp,
                _["obstitre"] = this->posteriorObsTitre
            );
            clock.stop("profiletimes");

            return out;
        }
        
        void updateAllChains()
        {
            if (onDebug) Rcpp::Rcout << "Pre: updateRJMC" << std::endl;
            updateRJMC();            
            if (onDebug) Rcpp::Rcout << "Pre: updateGibbsTiming" << std::endl;
            if (onDebug) Rcpp::Rcout << "currInferredExpN: " << currInferredExpN << std::endl;
            if (onDebug) Rcpp::Rcout << "knownInfsN: " << knownInfsN << std::endl;

            if ((this->currInferredExpN != this->knownInfsN) & !this->knownExpInd) { 
                updateGibbsTiming();
            }

           this->counterFuncEval++;

            if (onDebug) Rcpp::Rcout << "Pre: updateOutputPosterior" << std::endl;
            updateOutputPosterior();
            if (onDebug) Rcpp::Rcout << "Pre: updateProposal" << std::endl;
        
            consoleUpdatefunction();
        }
        
        void updateRJMC() {
            if (onDebug) Rcpp::Rcout << "Pre: getAcceptanceRate" << std::endl;
            getAcceptanceRate();
        //    Rcpp::Rcout << "alpha (kin): " << this->alpha << std::endl;

            if (onDebug) Rcpp::Rcout << "Pre: updateSampleAndLogPosterior" << std::endl;
            updateSampleAndLogPosterior();
            if(this->conPropIn & (this->currJumpType == 1)) {
                updateProposal();
            }
        }
    

        void getAcceptanceRate()
        {        
            this->isSampleAccepted = false;
           // Evaluate the ab kinetics parameters and do rjmcmc
            if (onDebug) Rcpp::Rcout << "Pre: selectProposalDist" << std::endl;
            selectProposalDist();
            if (onDebug) Rcpp::Rcout << "Pre: JumpProposalDist" << std::endl;
            JumpProposalDist();
            if (onDebug) Rcpp::Rcout << "Pre: evalLogPosterior" << std::endl;
            // CALCULATE TITREATINFECTION
            this->proposedLogPosterior = this->evalLogPosterior(this->proposalSample, this->proposalJump, this->proposalInf, this->currentCovarianceMatrix, this->dataList);
            if (onDebug) Rcpp::Rcout << "Pre: evaluateMetropolisRatio" << std::endl;
            evaluateMetropolisRatio();
        }

        // Sample an individual who is exposed in current chain
        int sampleExposed() {
            if (onDebug) Rcpp::Rcout << "In: sampleExposed" << std::endl;

            int s = uniformDiscreteDist(0, this->N - 1); 
            while ((this->knownInfsVec(s) == 1) || (this->currentJump(s) == -1) ) {
                s = uniformDiscreteDist(0, this->N - 1); 
            }
            return s;
        }

        // Sample an individual who is not exposed in current chain
        int sampleNotExposed() {
            int s = uniformDiscreteDist(0, this->N - 1); 
            while ((this->knownInfsVec(s) == 1) || (this->currentJump(s) != -1 )) {
                s = uniformDiscreteDist(0, this->N - 1); 
            }
            return s;
        }

        // Sample an new infection status for an exposed individual
        void proposeInfection(int t) {
            if (onDebug) Rcpp::Rcout << "In: proposeInfection" << std::endl;

            if (this->propInferredInfN == 0) {
                this->proposalInf(t) = 1;
            } else {
                boost::random::bernoulli_distribution<> l(0.5); 
                this->proposalInf(t) = l(rng);
            }
            bookKeepingSize();
        }

        // Resample the exposure time fo and individual, does not resample the infection rate
        void resampleTime(int t) {
            if (onDebug) Rcpp::Rcout << "In: resampleTime" << std::endl;
            double r = uniformContinuousDist(0, 1);
            double r_inf = uniformContinuousDist(0, 1);

            // Completely new sample time is drawn 5% of the time for an individual.
            if (r < 0.05) {
                this->proposalJump(t) = this->exposureFunctionSample(); // this->exposureFunctionSample(); //runif(1, 1, 150) + 121 //
                while ((this->proposalJump(t) >= this->endTitreTime(t) - 7) || (this->proposalJump(t) < this->initialTitreTime(t))) {
                    // check this
                    this->proposalJump(t) = this->exposureFunctionSample(); // this->exposureFunctionSample(); //runif(1, 1, 150) + 121 //
                }
                proposeInfection(t);
            // An exisiting sample time is updated 95% of the time for an individual.
            } else {
                double temp = this->currentJump(t) + normalDistSample(0, exp(this->adaptiveGibbsSD(t)));
                if (temp >= this->endTitreTime(t) - 7) {
                    temp = this->endTitreTime(t) - 7; // this->exposureFunctionSample(); //runif(1, 1, 150) + 121 //
                } else if (temp < this->initialTitreTime(t)) {
                    temp = this->initialTitreTime(t); // this->exposureFunctionSample(); //runif(1, 1, 150) + 121 //
                }
                this->proposalJump(t) = temp;
            }
            bookKeepingSize();
        }

        // 'birth' process of rjmcmc regime. Samples a new exposure from list unexposed individuals.
        void sampleNewTime(int t) {
            if (onDebug) Rcpp::Rcout << "In: sampleNewTime" << std::endl;

            double r = uniformContinuousDist(0, 1);


            if ((this->historicJump(t) < 0) | (r < 0.05)) { 
                // New exposure time is resampled from scratch, only do if first time individual exposed in mcmc chain or 5% of time thereafter.
                    this->proposalJump(t) = this->exposureFunctionSample(); 
                    while ((this->proposalJump(t) >= this->endTitreTime(t) - 7) || (this->proposalJump(t) < this->initialTitreTime(t))) {
                        this->proposalJump(t) = this->exposureFunctionSample(); 
                    }
                    bookKeepingSize();
                    proposeInfection(t);
            } else { 
                // New exposure time is taken from previous exposure time in mcmc chain (95% of the time)
                    this->proposalJump(t) = this->historicJump(t);
                    bookKeepingSize();
            }
            bookKeepingSize();
        }

        void bookKeepingSize() {
            int count_exp = 0;
            int count_inf = 0;
            for (int i = 0; i < this->N; i++) {
                if (this->proposalJump[i] > -1) {
                    count_exp++;
                    if (this->proposalInf[i] == 1) {
                        count_inf++;
                    }
                }
            }            
            this->propInferredExpN = count_exp;
            this->propInferredInfN = count_inf;
        }
        
        void bookKeepingSizeInf() {
            int count_inf = 0;
            for (int i = 0; i < this->N; i++) {
                if (this->proposalJump[i] >= -1) {
                    if (this->proposalInf[i] == 1) {
                        count_inf++;
                    }
                }
            }            
            this->propInferredInfN = count_inf;
        }

        void JumpProposalDist() {
            if (onDebug) Rcpp::Rcout << "Redfine prop = curr" << std::endl;
            this->proposalJump = this->currentJump;
            this->proposalInf = this->currentInf;

            double q = uniformContinuousDist(0, 1);
          //  double r;
            VectorXd q_prob(3);
            if (onDebug) Rcpp::Rcout << "get q_prob" << std::endl;
            if (onDebug) Rcpp::Rcout << "this->knownInfsN: " << this->knownInfsN << std::endl;


            if (this->knownExpInd) {
                q_prob << 0, 1, 1;
            } else {
                if (this->currInferredExpN == this->knownInfsN) {
                    q_prob << 0, 0.67, 1.0;
                } else if (this->currInferredExpN == this->N) {
                    q_prob << 0.33, 1.0, 0;
                } else {
                    q_prob << 0.33, 0.67, 1.0;
                }
            }

            // Substract a value
            if (onDebug) Rcpp::Rcout << "Find jumping" << std::endl;
            if (q < q_prob(0)) {
             //   Rcpp::Rcout << "In deletion" << std::endl;
                this->currJumpIdx = this->sampleExposed(); // length is n_ - 114
                this->historicJump(this->currJumpIdx) = this->currentJump(this->currJumpIdx);
                this->proposalJump(this->currJumpIdx) = -1; // length now n_ - 144 - 1
                this->currJumpType = 0;
                bookKeepingSize();
                //this->propInferredExpN = this->currInferredExpN - 1;
            // Stay same
            } else if (q < q_prob(1)) {
               // Rcpp::Rcout << "In stay same" << std::endl;
                this->propInferredExpN = this->currInferredExpN;
                if(this->currInferredExpN != this->knownInfsN) {
                    this->currJumpIdx = this->sampleExposed(); // length is n_ - 114
                    // RESAMPLE INFECTION
                    proposeInfection(this->currJumpIdx);
                }
                this->currJumpType = 1;
            // Add an infection same
            } else if (q < q_prob(2)) {
              //  Rcpp::Rcout << "In addition same" << std::endl;
                this->currJumpIdx = this->sampleNotExposed(); // length is 220 - n_
                sampleNewTime(this->currJumpIdx);
                this->currJumpType = 2;
                //this->propInferredExpN = this->currInferredExpN + 1;
            }
            if (onDebug) Rcpp::Rcout << "END jump sample" << std::endl;
        }

        void updateGibbsTiming() {
            if (onDebug) Rcpp::Rcout << "In: updateGibbsTiming" << std::endl;

            this->currJumpType = 1;
            for (int i = 0; i < this->noGibbsSteps; i++) { 
                this->proposalJump = this->currentJump;

                this->gibbsIdx = this->sampleExposed(); // stuck here
                resampleTime(this->gibbsIdx);

                selectProposalDist(false);
                this->proposedLogPosterior = this->evalLogPosterior(this->proposalSample, this->proposalJump, this->proposalInf, this->currentCovarianceMatrix, this->dataList);

                evaluateMetropolisRatio();
                updateJumpHistoric();
                updateProposalGibbs();
            }
        }
        
        void selectProposalDist(bool update = true){
          //  if (onDebug) Rcpp::Rcout << "Post: selectProposalDist" << std::endl;

            if (this->workingIteration < this->burninAdaptiveCov || uniformContinuousDist(0, 1) < 0.05 || !this->onAdaptiveCov ){
                generateSampleFromNonAdaptiveProposalDist(update);
            }
            else{
                generateSampleFromAdaptiveProposalDist(update);
            }
        }
        
        void generateSampleFromNonAdaptiveProposalDist(bool update = true)
        {
            double s;
            s = exp(this->nonadaptiveScalar);
            if(update) {
                this->counterNonAdaptive++; this->isProposalAdaptive = false;
            }
            this->currentCovarianceMatrix = s*this->nonadaptiveCovarianceMat;

            Mvn_sampler.updateCholesky(this->currentSample, this->currentCovarianceMatrix);
            this->proposalSample = Mvn_sampler.sampleTrunc(this->lowerParBounds, this->upperParBounds, 10, this->onDebug);
            
            errorCheckVectorValid(this->proposalSample);
        }
        
        void generateSampleFromAdaptiveProposalDist(bool update = true)
        {
            double s;
            s = exp(this->adaptiveScalar);
            if(update) {
                this->counterAdaptive++; this->isProposalAdaptive = true;
            }
            this->currentCovarianceMatrix = s*this->adaptiveCovarianceMat;
            Mvn_sampler.updateCholesky(this->currentSample, this->currentCovarianceMatrix);
            this->proposalSample = Mvn_sampler.sampleTrunc(this->lowerParBounds, this->upperParBounds, 10, this->onDebug);
            
            errorCheckVectorValid(this->proposalSample);
        }
        
        void evaluateMetropolisRatio()
        {
            double rjadjustmentFactor = 0;
            if(std::isnan(this->proposedLogPosterior) || std::isinf(this->proposedLogPosterior)) {
                this->alpha = 0;
            } else {
                if (this->currJumpType == 0) {
                    Function evalLoglikelhoodCOP_i = this->evalLoglikelhoodCOP[as<string>(this->copNames[0])];
                    NumericVector pars = this->currentParsCOP[as<string>(this->copNames[0])];

                    rjadjustmentFactor = log(this->currInferredExpN) - log((this->N - this->currInferredExpN + 1) ) + 
                        this->exposureFunctionDensity(this->currentJump(this->currJumpIdx)) +
                        as<double>(evalLoglikelhoodCOP_i(this->currentInf(this->currJumpIdx), this->currentTitreExp(this->currJumpIdx), pars) ); 
                } else if (this->currJumpType == 1) {
                    rjadjustmentFactor = 0;
                } else if (this->currJumpType == 2) {
                    Function evalLoglikelhoodCOP_i = this->evalLoglikelhoodCOP[as<string>(this->copNames[0])];
                    NumericVector pars = this->currentParsCOP[as<string>(this->copNames[0])];

                    rjadjustmentFactor = log((this->N - this->currInferredExpN)) - log(this->currInferredExpN + 1) - 
                        this->exposureFunctionDensity(this->proposalJump(this->currJumpIdx)) - 
                        as<double>(evalLoglikelhoodCOP_i(this->proposalInf(this->currJumpIdx), this->proposalTitreExp(this->currJumpIdx), pars) ); // dlnorm(proposalJump[currJumpIdx], 3.395, 0.5961)
                }
                this->alpha = min(1.0, exp((this->proposedLogPosterior - this->currentLogPosterior + rjadjustmentFactor)));
            }
        }

        std::vector<double> diff_ind;

        void finddifferInf() {
            this->diff_ind.clear();
            int count = 0;
            for (int i = 0; i < this->N; i++) {
                if (this->currentInf(i) != this->proposalInf(i) || this->currentJump(i) != this->proposalJump(i)) {
                    this->diff_ind.push_back(i);
                }
            }
        }
        
        void updateSampleAndLogPosterior()
        {
            if (uniformContinuousDist(0, 1) < this->alpha) {
                if (onDebug) Rcpp::Rcout << "In: ACCEPTED" << std::endl;
                this->finddifferInf();
                this->isSampleAccepted = true; this->counterAccepted++;

                this->currentSample = this->proposalSample;
                this->currentJump = this->proposalJump;
                this->currentInf = this->proposalInf;

                this->currentTitreExp = this->proposalTitreExp;
                this->currentObsTitre = this->proposalObsTitre;

                this->currentLogPosterior = this->proposedLogPosterior;
                this->currInferredExpN = this->propInferredExpN;
                this->currInferredInfN = this->propInferredInfN;
            } else {
                this->proposalSample = this->currentSample;
                this->proposalJump = this->currentJump;
                this->proposalInf = this->currentInf;
                this->proposalTitreExp = this->currentTitreExp;
                this->proposalObsTitre = this->currentObsTitre;

                this->proposedLogPosterior = this->currentLogPosterior;
                this->propInferredExpN = this->currInferredExpN;
                this->propInferredInfN = this->currInferredInfN;                
            }
        }
        void updateJumpHistoric()
        {
            if (uniformContinuousDist(0, 1) < this->alpha) {
                if (onDebug) Rcpp::Rcout << "In: ACCEPTED JUMP" << std::endl;
                //Rcpp::Rcout << "ACCEPTED (GIBBS): " << this->alpha << std::endl;
                this->currentSample = this->proposalSample;
                this->currentJump = this->proposalJump;
                this->currentInf = this->proposalInf;
                this->currentTitreExp = this->proposalTitreExp;
                this->currentObsTitre = this->proposalObsTitre;
                this->finddifferInf();

                this->currentLogPosterior = this->proposedLogPosterior;
                this->currInferredExpN = this->propInferredExpN;
                this->currInferredInfN = this->propInferredInfN;
            } else {
                this->proposalSample = this->currentSample;
                this->proposalJump = this->currentJump;
                this->proposalInf = this->currentInf;
                this->proposalTitreExp = this->currentTitreExp;
                this->proposalObsTitre = this->currentObsTitre;

                this->proposedLogPosterior = this->currentLogPosterior;
                this->propInferredExpN = this->currInferredExpN;
                this->propInferredInfN = this->currInferredInfN;     
            }
        }
        void updateOutputPosterior()
        {
            if ((this->workingIteration > (this->burninPosterior-1)) && 
                (this->workingIteration%thin == 0) ) {

                for (int p = 0; p < this->numberFittedPar; p++)
                    this->posteriorOut(this->counterPosterior, p) = this->currentSample(p);
                
                this->posteriorOut(this->counterPosterior, this->numberFittedPar) = this->currentLogPosterior;
                this->posteriorOut(this->counterPosterior, this->numberFittedPar+1) = (double)this->counterAccepted/(double)this->counterFuncEval;


                for (int j = 0; j < this->lengthJumpVec; j++) {
                    this->posteriorJump(this->counterPosterior, j) = this->currentJump(j);
                    this->posteriorInf(this->counterPosterior, j) = this->currentInf(j);
                    this->posteriorTitreExp(this->counterPosterior, j) = this->currentTitreExp(j);
                    this->posteriorObsTitre(this->counterPosterior, j) = this->currentObsTitre(j);
                }
                this->counterPosterior++;            
            }
        }
        
        void updateProposal()
        {       
            if (this->isProposalAdaptive){
                this->adaptiveScalar += this->stepSizeRobbinsMonro * pow(1+this->counterAdaptive,-0.5)*(this->alpha - 0.234);
                errorCheckNumberValid(this->adaptiveScalar);
                trimAdaptiveValues(this->adaptiveScalar);
            }
            else{   
                this->nonadaptiveScalar += this->stepSizeRobbinsMonro * pow(1+this->counterNonAdaptive,-0.5)*(this->alpha - 0.234);
                errorCheckNumberValid(this->nonadaptiveScalar);
                trimNonAdaptiveValues(this->nonadaptiveScalar);
            }
    
            // Update adaptive proposal stuff
            if(((this->workingIteration) % (this->updatesAdaptiveCov) == 0) && (this->workingIteration > this->burninAdaptiveCov)){

                //int iPosterior = (this->workingIteration-this->burninAdaptiveCov);
                this->iPosterior++;
                double gainFactor = pow(1+iPosterior, -0.5);
               // if (iPosterior == this->updatesAdaptiveCov){
                if (iPosterior == 1){
                    this->currentSampleMean = this->currentSample;
                    this->adaptiveScalar = this->nonadaptiveScalar;
                }
                else{
                    
                    this->currentSampleMean = this->currentSampleMean + gainFactor*(this->currentSample-this->currentSampleMean);
                    errorCheckVectorValid(this->currentSampleMean);
                    this->adaptiveCovarianceMat = this->adaptiveCovarianceMat + gainFactor*((this->currentSample-this->currentSampleMean)*((this->currentSample)-(this->currentSampleMean)).transpose()) - gainFactor*this->adaptiveCovarianceMat;
                    errorCheckMatrixValid(this->adaptiveCovarianceMat);
                }
            }
        }

        void updateProposalGibbs()
        {
            // Update adaptive proposal stuff
            this->counterAdaptiveGibbs(this->gibbsIdx)++;
            double gainFactor = pow(this->counterAdaptiveGibbs(this->gibbsIdx), -0.5);
            this->adaptiveGibbsSD(this->gibbsIdx) += gainFactor*(this->alpha - 0.44);

            errorCheckNumberValid(this->adaptiveGibbsSD(this->gibbsIdx));
            trimAdaptiveValues(this->adaptiveGibbsSD(this->gibbsIdx));
        }
        
        void trimNonAdaptiveValues(double value)
        {
            this->nonadaptiveScalar = log(MAX(MIN(exp(value),  this->covarMaxVal), 1e-25));
        }
        
        void trimAdaptiveValues(double value)
        {
            this->adaptiveScalar = log(MAX(MIN(exp(value), this->covarMaxVal), 1e-25));
        }

        void errorCheckNumberValid(double value)
        {
            if (isinf(value)||isnan(value)){
                Rcout << "The number is not finite." << endl;
                stop("Value: ", value);
            }
        }
        
        void errorCheckVectorValid(const VectorXd& proposalSample)
        {
            for (int i = 0; i < this->numberFittedPar; i++){
                if (isinf(proposalSample(i))||isnan(proposalSample(i))){
                    Rcout << "The proposed vector is not finite." << endl;
                    stop("Value: ", proposalSample(i));
                }
            }
        }
        
        void errorCheckMatrixValid(const MatrixXd& covarianceMatrix)
        {
            for (int i = 0; i < this->numberFittedPar; i++){
                for (int j = 0; j < this->numberFittedPar; j++){
                    if (isinf(covarianceMatrix(i,j)) || isnan(covarianceMatrix(i,j))){
                        Rcout << "The proposed matrix is not finite." << endl;
                        stop("Value: ", covarianceMatrix(i,j));
                    }
                }
            }
        }
                
        double uniformContinuousDist(double minValue, double maxValue)
        {
            boost::random::uniform_real_distribution<> u(minValue,maxValue); return u(rng);
        }
        
        double uniformDiscreteDist(int minValue, int maxValue)
        {
            boost::random::uniform_int_distribution<> u(minValue, maxValue); return u(rng);
        }

        double normalDistSample(double mean, double sigma)
        {
            boost::random::normal_distribution<> n(mean, sigma); return n(rng);
        }
            
        void consoleUpdatefunction()
        {
            int i = this->workingIteration;
            if(i%this->consoleUpdates == 0) {
                Rcpp::Rcout << "Running MCMC-PT iteration number: " << this->workingIteration << " of " <<  this->iterations << ". Chain number " << this->chainNumber << ". Current logpost: " << this->currentLogPosterior << ". No exp/inf: " <<  this->currInferredExpN << "/" << this->currInferredInfN << " .          \r";
                if (this->onDebug) {
                    Rcpp::Rcout << "\n Current values: " << this->currentSample << std::endl;
                }
            }
        }


        std::vector<DoubleWithString> orderExposureEvents(
                const VectorXd& jump_inf, 
                const VectorXd& jump, double initialtime, int i_idx, double time) {

                std::vector<DoubleWithString> df_order_exp;
                        
                df_order_exp.emplace_back("none", initialtime); // Initial event
                df_order_exp.emplace_back("bleed", time); // End event

                for (int i = 0; i < this->exposureNames.size(); i++) {
                    List here = this->abkineticsModel[as<string>(this->exposureNames[i])]; 
                    bool inferred = here["inferred"];
                    if (!inferred){
                        NumericVector knowninf = here["known_inf"];
                        if (knowninf[i_idx] > -1 && knowninf[i_idx] < time) {
                            df_order_exp.emplace_back(as<string>(this->exposureNames[i]), knowninf[i_idx]); // "pre-delta" event
                        }
                    } else {
                        if (jump_inf[i_idx] == 1 && jump[i_idx] < time) {
                            df_order_exp.emplace_back(as<string>(this->exposureNames[i]), jump[i_idx]); // "delta" event
                        }
                    }
                }
                std::sort(df_order_exp.begin(), df_order_exp.end());

                return df_order_exp;
        }


        NumericVector createNamedParam(const VectorXd& params) {
            NumericVector paramsName(this->fittedParamNames.size());
            for (int i = 0; i < this->fittedParamNames.size(); ++i) {
                paramsName[i] = params[i];
            }
            paramsName.attr("names") = this->fittedParamNames;
            return paramsName;
        }


        double calTitre(double& titre_est, string& exposureType_i, double& timeSince) {
            if (onDebug) Rcpp::Rcout << "In: calTitre" << std::endl;
            int abkey = this->abmap[exposureType_i];
            titre_est = Rcpp::as<double>(this->abinfo[abkey].func(titre_est, timeSince, this->abinfo[abkey].params));
            return titre_est;
        }

        void define_abs() {
            abinfo.clear();
            abmap.clear();
            for (int i = 0; i < this->abkineticNames.size(); i++) {
                Function thisfunc = this->evalabkineticsFunc[as<string>(this->abkineticNames[i])];
                NumericVector params = this->currentParsAb[as<string>(this->abkineticNames[i])];
                this->abinfo.emplace_back(as<string>(abkineticNames[i]), thisfunc, params); // Initial event
            }

            for (int i = 0; i < this->abkineticNames.size(); i++) {
                this->abmap[as<string>(this->abkineticNames[i])] = i;
            }
        }

        double evaluateLogLikelihoodCOP_cpp(VectorXd params, VectorXd jump, VectorXd jumpinf, bool init) {
            if (this->profile)  clock.tick("cop_ll");
            if (onDebug) Rcpp::Rcout << "In: evaluateLogLikelihoodCOP_cpp" << std::endl;

            VectorXd titreExp(this->N);

            double titre_est, initialtime, time_until;
            Function evalLoglikelhoodCOP_i = this->evalLoglikelhoodCOP[as<string>(this->copNames[0])];
            NumericVector pars = this->currentParsCOP[as<string>(this->copNames[0])];

            double ll = 0;
            for (int i_idx = 0; i_idx < this->N; ++i_idx) {
                if (jump[i_idx] == -1) {
                    titreExp[i_idx] = -1;
                } else {
                    // Add custom func in here which can be used to input pre-determined pre-exposure titre
                    std::vector<DoubleWithString> proposalTitreFull_i = this->proposalTitreFull[i_idx];

                    for (int j = 0; j < proposalTitreFull_i.size(); j++) {
                        if(proposalTitreFull_i[j].name == this->exposureNameInf) {
                            titre_est = proposalTitreFull_i[j].value;
                        }
                    }
                    titreExp[i_idx] = titre_est;

                    ll += as<double>(evalLoglikelhoodCOP_i(jumpinf[i_idx], titreExp[i_idx], pars ) );
                }
            }
            if (init) {
                this->currentTitreExp = titreExp;
            } else {
                this->proposalTitreExp = titreExp;
            }
            if (this->profile) clock.tock("cop_ll");
            return ll;
        }

        double evaluateLogLikelihoodObs_cpp(const VectorXd& params, const VectorXd& jump, const VectorXd& jump_inf, bool init) {
            if (this->profile) clock.tick("obs_ll");
            if (onDebug) Rcpp::Rcout << "In: evaluateLogLikelihoodObs_cpp" << std::endl;

            NumericVector paramsN = this->createNamedParam(params);
            std::vector<DoubleWithString> df_order_exp; 
            VectorXd obsTitre(this->N_data);
        
            double time_until, titre_est, time, initialtime, initialtitre ;
            double titre_val, sigma;
            int i_idx;
            double ll = 0;

            Function evalLoglikelhoodObs_i = this->evalLoglikelhoodObs[as<string>(this->observationalNames[0])];
            NumericVector pars = this->currentParsObs[as<string>(this->observationalNames[0])];
            

            define_abs();

            // For each observation
            int k = 0;
            for (int i_idx = 0; i_idx < this->N; ++i_idx) {

                std::vector<DoubleWithString> proposalTitreFull_i = proposalTitreFull[i_idx];
                NumericVector titre_val_i = this->titre_list[i_idx];

                for (int j = 0; j < titre_val_i.size(); j++) {

                    if(proposalTitreFull_i[j].name == "bleed") {
                        titre_est = proposalTitreFull_i[j].value;
                        obsTitre[k] = titre_est; k++;
                        titre_val = titre_val_i[j];
                        ll += as<double>(evalLoglikelhoodObs_i(titre_val, titre_est, pars) );
                    }
                }
            }
            if (init) {
                this->currentObsTitre = obsTitre;
            } else {
                this->proposalObsTitre = obsTitre;
            }
            if (this->profile) clock.tock("obs_ll");
            return ll;
        }
    };

    void init_samplePriorDistributions(rjmc_full::RJMC_FULL_D* model, Rcpp::Function samplePriorDistributions) {
        auto func = [samplePriorDistributions](RObject dataList) {
            PutRNGstate();
            auto rData = samplePriorDistributions(dataList);
            GetRNGstate();
            return Rcpp::as<VectorXd>(rData);
        };
        model->samplePriorDistributions = func;
    }

    void init_evaluateLogPrior(rjmc_full::RJMC_FULL_D* model, Rcpp::Function evaluateLogPrior) {
        auto func = [evaluateLogPrior](VectorXd params, VectorXd jump, RObject dataList) {
            PutRNGstate();
            auto rData = evaluateLogPrior(params, jump, dataList);
            GetRNGstate();
            return Rcpp::as<double>(rData);
        };
        model->evaluateLogPrior = func;
    }

    void init_initialiseJump(rjmc_full::RJMC_FULL_D* model, Rcpp::Function initialiseJump) {
        auto func = [initialiseJump](RObject dataList) {
            PutRNGstate();
            auto rData = initialiseJump(dataList);
            GetRNGstate();
            return Rcpp::as<VectorXd>(rData);
        };
        model->initialiseJump = func;
    }

    void init_exposureFunctionSample(rjmc_full::RJMC_FULL_D* model, Rcpp::Function exposureFunctionSample) {
        auto func = [exposureFunctionSample]() {
            PutRNGstate();
            auto rData = exposureFunctionSample();
            GetRNGstate();
            return Rcpp::as<double>(rData);
        };
        model->exposureFunctionSample = func;
    }

    void init_exposureFunctionDensity(rjmc_full::RJMC_FULL_D* model, Rcpp::Function exposureFunctionDensity) {
        auto func = [exposureFunctionDensity](double time) {
            PutRNGstate();
            auto rData = exposureFunctionDensity(time);
            GetRNGstate();
            return Rcpp::as<double>(rData);
        };
        model->exposureFunctionDensity = func;
    }
};
// namespace rjmc_full
#endif
