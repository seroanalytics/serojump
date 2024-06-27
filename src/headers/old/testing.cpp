#pragma clang diagnostic ignored "-Wunknown-pragmas"
#pragma clang diagnostic ignored "-Wlanguage-extension-token"
#pragma clang diagnostic ignored "-Wunused-but-set-variable"

#include <RcppCommon.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Core>
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
#include <memory>
#include <iostream>

// Define Rcpp dependicies and plugins
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins("cpp14")]]

// Makes for cleaner code
using namespace Rcpp;
using namespace std;
using namespace Eigen;
using namespace boost::math;

// Define a random number generator
static thread_local std::mt19937 rng(std::random_device{}());

#define EIGEN_DONT_VECTORIZE
#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT

// Define the value of pi
#define PI 3.14159265358979323846

// Define a MIN and MAX function
#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b)) // define MAX function for use later
#endif

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b)) // define MAX function for use later
#endif

#include "../utils.hpp"
#include "../mvn.hpp"

struct SeroJumpBase : public std::enable_shared_from_this<SeroJumpBase>
{
  /*  Rcpp::List evalLoglikelhoodObs, evalabkineticsFunc; // List of functions to evaluate likelihoods and kinetics
    Rcpp::List infoModel, observationalModel, abkineticsModel; // List of models and parameters

    SeroJumpBase(Rcpp::List infoModel_in, Rcpp::List observationalModel_in, Rcpp::List abkineticsModel_in) 
        : infoModel(infoModel_in), observationalModel(observationalModel_in), abkineticsModel(abkineticsModel_in) {
    }
    virtual ~SeroJumpBase() = default;*/

    Rcpp::List evalLoglikelhoodObs, evalabkineticsFunc; // List of functions to evaluate likelihoods and kinetics
    Rcpp::List infoModel, observationalModel, abkineticsModel; // List of models and parameters

    List parsAbKinN, parsObsN; // Lists of all the parameters names for each model
    std::vector<abKineticsInfo> abinfo; // Vector of all the abkinetics information
    std::map<std::string, int> abmap, mapOfObs, mapOfExp; // Maps of the abkinetics and observational models

    std::map<std::pair<std::string, std::string>, int> mapOfAbkinetics; // Map of the abkinetics model
    std::vector<string> biomarkers, exposureType, abID; // Vector of biomarkers and exposure types
    string exposureFitted; // The exposure type which is fitted
    List exposureInfo, knownInf; // List of exposure information and known infections
    int B; // Number of biomarkers

    /** 
     * @brief Constructor for the RJMC base class
     * @param infoModel_in The information model
     * @param observationalModel_in The observational model
     * @param abkineticsModel_in The antibody kinetics model
     * 
     */
    SeroJumpBase(Rcpp::List infoModel_in, Rcpp::List observationalModel_in, Rcpp::List abkineticsModel_in) : infoModel(infoModel_in), observationalModel(observationalModel_in), abkineticsModel(abkineticsModel_in) {
        Rcpp::Rcout << "Testing in" << std::endl;

    }        
    
    // Virtual destructor
    virtual ~SeroJumpBase() = default;

    List currentParsObs, currentParsAb;
    StringVector exposureNames, fittedParamNames;

    bool conPropIn = true;
    bool disPropIn = false;

    // Variables for the RJMC algorithm
    double nonadaptiveScalar, adaptiveScalar;
    double currentLogPosterior;

    VectorXd lowerParBounds, upperParBounds;

    MatrixXd nonadaptiveCovarianceMat, adaptiveCovarianceMat;
    MatrixXd currentSample, currentSampleMean;
    VectorXd proposalSample;
    MatrixXd currentCovarianceMatrix;
    VectorXd currentJump, proposalJump;
    // VectorXd currentInf, proposalInf;
    MatrixXd currentTitreExp, proposalTitreExp;
    MatrixXd currentObsTitre, proposalObsTitre;

    VectorXd historicJump;

    vector<vector<DoubleWithString> >  currentEventsFull, proposalEventsFull;
    vector<vector<vector<DoubleWithString> > > currentTitreFull, proposalTitreFull;

    // Outputs for posterior
    MatrixXd posteriorOut, posteriorJump, posteriorInf;
    std::vector<MatrixXd> posteriorTitreExp, posteriorObsTitre;

    // Information extracted from the settingd
    int iPosterior;
    int iterations, posteriorSamplesLength, thin, burninPosterior, burninAdaptiveCov, consoleUpdates, updatesAdaptiveCov, chainNumber;
    int lengthJumpVec;
    int numberFittedPar;
    int workingIteration;
    bool onDebug, onAdaptiveCov, profile;
    bool isSampleAccepted, isProposalAdaptive;
    double stepSizeRobbinsMonro;

    // Variables to be filled
    int counterFuncEval, counterAccepted, counterPosterior ,counterAdaptive;
    int counterNonAdaptive;
    double proposedLogPosterior, alpha, covarMaxVal, covarInitVal, covarInitValAdapt;

    // Functions for the RJMC algorithm
    std::function<VectorXd(RObject)> samplePriorDistributions;
    std::function<VectorXd(RObject)> initialiseJump;
    std::function<double(VectorXd, VectorXd, RObject)> evaluateLogPrior;

    std::function<double()> exposureFunctionSample;
    std::function<double(double)> exposureFunctionDensity;

    // Functions for internal 
    Mvn Mvn_sampler;

    // Important information from data
     RObject dataList;
    List dataListCPP;
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
    MatrixXd initialTitreValue;
    VectorXd endTitreTime;
    MatrixXd endTitreValue;
    VectorXd times_full, id_full;
    MatrixXd titre_full;
    List titre_list, times_list;
    //StringVector exposurePeriods;
    //string exposureNameInf;


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

    void initialiseBaseClass(List settings, RObject dataList, int i)
    {
        this->onDebug = settings["onDebug"];

        // Extract Data from the dataList
        this->dataList = dataList;            
        this->dataListCPP = as<List>(dataList);
        this->fittedParamNames = this->dataListCPP["par_names"];

        if (this->onDebug) Rcpp::Rcout << "In: initialiseClass" << std::endl;
        if (this->onDebug) Rcpp::Rcout << "In: Extract data from dataList" << std::endl;


        this->N = this->dataListCPP["N"];
        this->N_data = this->dataListCPP["N_data"];
        this->titre_list = this->dataListCPP["titre_list"];
        this->times_list = this->dataListCPP["times_list"];
        this->titre_full = this->dataListCPP["titre_full"];
        this->times_full = this->dataListCPP["times_full"];
        this->id_full = this->dataListCPP["id_full"];
        //if (this->onDebug) Rcpp::Rcout << "In: Extract data from dataList2" << std::endl;

        this->knownExpVec = this->dataListCPP["knownExpVec"];
        if (this->knownExpVec.size() > 1) {
            this->knownExpInd = true;
        }
        this->knownInfsVec = this->dataListCPP["knownInfsVec"];
        this->knownInfsTimeVec = this->dataListCPP["knownInfsTimeVec"];
        this->knownInfsN = this->dataListCPP["knownInfsN"];
       // if (this->onDebug) Rcpp::Rcout << "In: Extract data from dataList3" << std::endl;



      //  if (this->onDebug) Rcpp::Rcout << "In: Extract data from dataList4" << std::endl;

        this->initialTitreTime = this->dataListCPP["initialTitreTime"];
        this->endTitreTime = this->dataListCPP["endTitreTime"];
        this->initialTitreValue = this->dataListCPP["initialTitreValue"];
        this->endTitreValue = this->dataListCPP["endTitreValue"];

       // if (this->onDebug) Rcpp::Rcout << "In: Check 1" << std::endl;
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
    //    this->profile = settings["profile"];
        this->lengthJumpVec = settings["lengthJumpVec"];
        if (this->onDebug) Rcpp::Rcout << "In: Check 3" << std::endl;

        this->lowerParBounds = settings["lowerParBounds"];
        this->upperParBounds = settings["upperParBounds"];

        this->covarInitVal = settings["covarInitVal"];
        this->covarInitValAdapt = settings["covarInitValAdapt"];
        this->covarMaxVal = settings["covarMaxVal"];

        // Initialise fixed variables 
        double alphaMVN = -sqrt(2)*ErfInv(0.234-1);
        this->stepSizeRobbinsMonro = (1.0-1.0/(double)this->numberFittedPar)*(pow(2*3.141, 0.5)*exp(alphaMVN*alphaMVN*0.5))/(2*alphaMVN) + 1.0/(this->numberFittedPar*0.234*(1-0.234));
        this->chainNumber = i;


        // Initialise dynamic variables 
        // // Counters
        this->counterFuncEval = 0;
        this->counterAccepted = 0;
        this->counterPosterior = 0;
        this->counterAdaptive = 0;
        this->counterNonAdaptive = 0;
        this->iPosterior = 0;

        // Adaptice mcmc variables 
        this->nonadaptiveScalar = 0;
        this->adaptiveScalar = 0;
        this->nonadaptiveScalar = log(0.1*0.1/(double)this->numberFittedPar);
        this->adaptiveScalar = log(2.382*2.382/(double)this->numberFittedPar);
        this->adaptiveGibbsSD = VectorXd::Constant(this->N, 0);
        this->counterAdaptiveGibbs = VectorXd::Zero(this->N);

        // Initialise markov chains 
        this->currentLogPosterior = 0;

        this->currInferredInfN = this->knownInfsN;
        this->propInferredInfN = this->knownInfsN;
        if (this->knownExpInd) {
            this->currInferredExpN = this->knownExpVec.sum();
            this->propInferredExpN = this->knownExpVec.sum();
        } else {
            this->currInferredExpN = this->knownInfsN;
            this->propInferredExpN = this->knownInfsN;    
        }

        if (this->onDebug) Rcpp::Rcout << "In: Check 4" << std::endl;

        this->posteriorSamplesLength = (this->iterations-this->burninPosterior)/(this->thin);
        this->posteriorOut = MatrixXd::Zero(this->posteriorSamplesLength, this->numberFittedPar + 2);
        this->posteriorJump = MatrixXd::Zero(this->posteriorSamplesLength, this->lengthJumpVec);

        this->currentSample = VectorXd::Zero(this->numberFittedPar);
        this->currentSampleMean = VectorXd::Zero(this->numberFittedPar);
        this->proposalSample = VectorXd::Zero(this->numberFittedPar);
        
        this->currentCovarianceMatrix = MatrixXd::Zero(this->numberFittedPar, this->numberFittedPar );
        this->nonadaptiveCovarianceMat = MatrixXd::Zero(this->numberFittedPar, this->numberFittedPar );
        this->adaptiveCovarianceMat = MatrixXd::Zero(this->numberFittedPar ,this->numberFittedPar );
      //  if (this->onDebug) Rcpp::Rcout << "In: Check 5" << std::endl;
        
        for(int parNum = 0; parNum < this->numberFittedPar ; parNum++){
            this->nonadaptiveCovarianceMat(parNum,parNum) = this->covarInitVal*(this->upperParBounds(parNum) - this->lowerParBounds(parNum));
            this->adaptiveCovarianceMat(parNum,parNum) = this->covarInitValAdapt*(this->upperParBounds(parNum) - this->lowerParBounds(parNum));
        }
        this->historicJump = VectorXd::Constant(this->N, -1);

    }



    // This is currently BROKEN DO NOT USE
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
    //     this->posteriorTitreExp = MatrixXd::Zero(this->posteriorSamplesLength, this->lengthJumpVec);
    //    this->posteriorObsTitre = MatrixXd::Zero(this->posteriorSamplesLength, this->lengthJumpVec);

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
    
    // SIMILARLY THIS IS NOT WORKING
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

    ///////////////////////////////////////////////////////////////////////////
    //// SUPPORT FUNCTIONS FOR THE RJMC ALGORITHM
    ///////////////////////////////////////////////////////////////////////////
    /**
     * @brief Ensure nonadaptiveScalar are within bounds (1e-25, his->covarMaxVal), if not truncats at min/max value
     */
    void trimNonAdaptiveValues(double value)
    {
        this->nonadaptiveScalar = log(MAX(MIN(exp(value),  this->covarMaxVal), 1e-25));
    }
    /**
     * @brief Ensure adaptiveScalar are within bounds (1e-25, his->covarMaxVal), if not truncats at min/max value
     */
    void trimAdaptiveValues(double value)
    {
        this->adaptiveScalar = log(MAX(MIN(exp(value), this->covarMaxVal), 1e-25));
    }

    /**
     * @brief Checks to make sure number is valid, stops analysis if not
     */
    void errorCheckNumberValid(double value)
    {
        if (isinf(value)||isnan(value)){
            Rcout << "The number is not finite." << endl;
            stop("Value: ", value);
        }
    }

    /**
     * @brief Checks to make sure vector is valid, stops analysis if not
     */
    void errorCheckVectorValid(const VectorXd& proposalSample)
    {
        for (int i = 0; i < this->numberFittedPar; i++){
            if (isinf(proposalSample(i))||isnan(proposalSample(i))){
                Rcout << "The proposed vector is not finite." << endl;
                stop("Value: ", proposalSample(i));
            }
        }
    }

    /**
     * @brief Checks to make sure matrix is valid, stops analysis if not
     */
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

};


// Placeholder classes for RecalculateTitres and EvaluateLogLikelihood
class RecalculateTitres {
public:
    /**
     * @brief Construct a new RecalculateTitres object
     * 
     * @param basePtr A shared pointer to the SeroJumpBase object
     */
    RecalculateTitres( std::shared_ptr<SeroJumpBase> basePtr)
        : parent(basePtr) {
    }

    /** 
     * @brief Recalculate the titre values for all individuals given the new exposure times and fixed parameters
     * @param params The new parameter values
     * 
     * This function updates order of events for each individual, updates the kinetic parameters values and recalculates the titre values for all individuals at the event times
     * 
     */
    void RecalculateTitresAll(VectorXd& params) {
        this->updateEventsFull();
        this->updateAbKineticParams(params);
        this->recalculateTitreAll();
    }

    /**
     * @brief Update the kinetic parameters values
     * @param params The new parameter values
     */
    void updateAbKineticParams(VectorXd& params) {
        NumericVector paramsN = this->createNamedParam(params);
        this->updatePars(paramsN);
        this->initAbInfo();
    }

    /** 
     * @brief Update the order of events for each individual
     * 
     */
    void updateEventsFull() {
        if (parent->onDebug) Rcpp::Rcout << "In: updateLists" << std::endl;
        std::vector<double> diff_ind = finddifferInf();
        std::vector<DoubleWithString> df_order_exp_i;
        int i_idx;
        for (int i = 0; i < diff_ind.size(); i++) {
            i_idx = diff_ind[i];
            df_order_exp_i = this->sortevents(i_idx, parent->proposalJump);
            parent->proposalEventsFull[i_idx] = df_order_exp_i;
        }
        if (parent->onDebug) Rcpp::Rcout << "End: updateLists" << std::endl;
    }


    /**
     * @brief Recalculate the titre values for all individuals at the event times
     * 
     * Now that the events have been updated, and the parameter values names and updated, we can recalculate the titre values for all individuals at the event times
     */
    void recalculateTitreAll() {
        std::vector<DoubleWithString> df_order_exp_i;
        for (int i_idx = 0; i_idx < parent->N; i_idx++) {
            df_order_exp_i = parent->proposalEventsFull[i_idx];
            parent->proposalTitreFull[i_idx] = this->recalculateTitre_i(df_order_exp_i, i_idx);
        }
    }

   
private:


    /**
     * @brief Find the individuals who exposure times have changes bewtween markov steps
     * @return A vector of ids of individuals having different exposure times
     */
    std::vector<double> finddifferInf() {
        std::vector<double> diff_ind; // numeric vector of ids of individuals having different exposure times
        for (int i = 0; i < parent->N; i++) {
            if (parent->currentJump(i) != parent->proposalJump(i)) {
                diff_ind.push_back(i);
            }
        }
        return diff_ind;
    }


    /** 
     * @brief Sort the events for individual i_idx
     * @param i_idx The individual id
     * @param jump The new exposure times
     * @return A vector of events sorted by time
     */
    std::vector<DoubleWithString> sortevents(int i_idx, const VectorXd& jump) {
        NumericVector times_full_i = parent->times_list[i_idx];
        std::vector<DoubleWithString> df_order_exp; // class with time (double) and event (string)
        // Append all bleed times
        for (int i = 0; i < times_full_i.size(); i++) {
            df_order_exp.emplace_back("bleed", times_full_i[i]); // End event
        }
        // Append all exposure times
        for (int i = 0; i < parent->exposureType.size(); i++) {
            string exposureType_i = parent->exposureType[i];

            if (exposureType_i == parent->exposureFitted){
                // Append an exposute time of the fitted exposure
                if (jump[i_idx] > -1) {
                    df_order_exp.emplace_back(exposureType_i, jump[i_idx]);
                }
            } else {
                // Append known exposure times for non-fitted exposures
                NumericVector known_inf_i = parent->knownInf[exposureType_i]; 
                if (known_inf_i[i_idx] > -1) {
                    df_order_exp.emplace_back(exposureType_i, known_inf_i[i_idx]); 
                }
            }
        }
        // Sort the events by time
        std::sort(df_order_exp.begin(), df_order_exp.end());
        return df_order_exp;
    }


    /** 
     * @brief Create a named parameter vector
     * @param params The parameter values
     * 
     * This assigns names to the parameter values
     * 
     * @return A named parameter vector
     */
    NumericVector createNamedParam(const VectorXd& params) {
        NumericVector paramsName(parent->fittedParamNames.size());
        for (int i = 0; i < parent->fittedParamNames.size(); ++i) {
            paramsName[i] = params[i];
        }
        paramsName.attr("names") = parent->fittedParamNames;
        return paramsName;
    }

    /** 
     * @brief Update the kinetic parameters values, 
     * @param paramsN The new names parameter vector
     * 
     * Updates currentParsAb and currentParsObs, the parameters of the Ab and observational model 
     * 
     */
    void updatePars(NumericVector paramsN) {
        for (int i = 0; i < parent->B; ++i) {
            NumericVector currentParsObs_i;
            int i_b = parent->mapOfObs[parent->biomarkers[i]];
            StringVector parnams = parent->parsObsN[i_b];
            for (int j = 0; j < parnams.size(); ++j) {
                currentParsObs_i.push_back(paramsN[as<string>(parnams[j])]);
            }
            parent->currentParsObs[parent->biomarkers[i]] = currentParsObs_i;
        }
        for (int i = 0; i < parent->abID.size(); ++i) {
            NumericVector currentParsAb_i;    
        // int abkey = mapOfAbkinetics[{biomarkerAB[i], exposureTypeAB[i]}];
            StringVector parnams = parent->parsAbKinN[i];

            for (int j = 0; j < parnams.size(); ++j) {
                currentParsAb_i.push_back(paramsN[as<string>(parnams[j])]);
            }
            parent->currentParsAb[parent->abID[i]] = currentParsAb_i;
        }
    }

    /**
     * @brief Define the abinfo vector
     * 
     * This function initialises the abinfo vector with the current values of the antibody kinetics parameters
     * 
     */
    void initAbInfo() {
        parent->abinfo.clear();
        for (int i = 0; i < parent->abID.size(); i++) {
            Function thisfunc = parent->evalabkineticsFunc[parent->abID[i]];
            NumericVector params = parent->currentParsAb[parent->abID[i]];
            parent->abinfo.emplace_back(parent->abID[i], thisfunc, params); // Initial event
        }
    }

    /** 
     * @brief Recalculate the titre values for individual i_idx at the event times
     * @param orderedEvents The ordered events for individual i_idx
     * @param i_idx The individual id
     * @return A vector of titre values for individual i_idx
     * 
     * 
     * 
     */
    std::vector<std::vector<DoubleWithString> > recalculateTitre_i(std::vector<DoubleWithString>& orderedEvents, int i_idx)  {
        if (parent->onDebug) Rcpp::Rcout << "In: calculateTitre" << std::endl;

        std::vector<std::vector<DoubleWithString> > df_order_titre;
        std::vector<DoubleWithString> df_order_titre_b;

        double initTime = parent->initialTitreTime[i_idx];
        string bio;
        double time_since, titre_obs;

        // Across all biomarkers
        for (int b = 0; b < parent->B; b++) {
            double anchor_time = initTime; // intiial time
            string anchor_func = orderedEvents[0].name; // initial event function
            // Extract 
            bio = parent->biomarkers[b]; // biomarker
            double initTitre = parent->initialTitreValue(i_idx, b); // start titre
            double anchor_titre = initTitre; //Anchot this

            // Define state with anchor func and start titre
            df_order_titre_b.clear();
            df_order_titre_b.emplace_back(anchor_func, anchor_titre);

            // For all events
            for (int i = 1; i < orderedEvents.size(); i++) {
                // Update time since last event
                time_since = orderedEvents[i].value - anchor_time;
                titre_obs = abkineticsFunction(anchor_titre, bio, anchor_func, time_since ); // Calcualte titre at new event
                df_order_titre_b.emplace_back(orderedEvents[i].name, titre_obs); // add to vector

                // Only update the anchors if the event is not a bleed, as bleeds don't affect titre trajectory
                if (orderedEvents[i].name != "bleed") {
                    // Update the anchor time, titre and function if the event is the fitted exposure
                    if (orderedEvents[i].name == parent->exposureFitted){
                        anchor_titre = titre_obs;
                        anchor_time = orderedEvents[i].value;
                        anchor_func = orderedEvents[i].name;    
                    } else {
                        // Update the anchor time, titre and function if the event is a known exposure (but not the fitted one)
                        NumericVector known_inf_i = parent->knownInf[orderedEvents[i].name]; 
                        if (known_inf_i[i_idx] > -1) {
                            anchor_titre = titre_obs;
                            anchor_time = orderedEvents[i].value;
                            anchor_func = orderedEvents[i].name;  
                        }
                    }
                }
            }
            // Add to the vector of titre values
            df_order_titre.push_back(df_order_titre_b);
        }
        if (parent->onDebug) Rcpp::Rcout << "Out: calculateTitre" << std::endl;
        return df_order_titre;

    }

    /** 
     * @brief Calculate the titre value at a given time since the last event
     * @param titre_est The titre value at the last event
     * @param biomarker The biomarker
     * @param exposureType_i The type of exposure
     * @param timeSince The time since the last event
     * @return The titre value at the new event
     * 
     * This function takes the string inputs and uses the mapOfAbkinetics to find the correct antibody kinetics function to use
     * 
     */
    double abkineticsFunction(double titre_est, string biomarker, string exposureType_i, double timeSince) {
        if (parent->onDebug) Rcpp::Rcout << "In: calTitre" << std::endl;
        int abkey = parent->mapOfAbkinetics[{biomarker, exposureType_i}];
        string abID_i = parent->abID[abkey] ;
        titre_est = Rcpp::as<double>(parent->abinfo[abkey].func(titre_est, timeSince, parent->abinfo[abkey].params));
        if (parent->onDebug) Rcpp::Rcout << "Out: calTitre" << std::endl;
        return titre_est;
    }

    // A shared pointer to the SeroJumpBase object
    std::shared_ptr<SeroJumpBase> parent;

};

class EvaluateLogLikelihood {
public:
 /**
     * @brief Construct a new RecalculateTitres object
     * 
     * @param basePtr A shared pointer to the SeroJumpBase object
     */
    EvaluateLogLikelihood( std::shared_ptr<SeroJumpBase> basePtr)
        : parent(basePtr) {
    };

    /** 
     * @brief Evaluate the log posterior of the model
     * @param param The parameter vector for the fitted parameters
     * @param jump The jump vector for the exposure times
     * @param covariance The covariance matrix for the fitted parameters
     * @param dataList The data list
     * @param init A boolean to indicate if the function is being called for the initialisation step
     * @return The log posterior of the model
     * 
     * This function evaluates the log posterior of the model given the parameter vector, jump vector, covariance matrix and data list. There are four main components to the log posterior:
     * 1. The log prior of the model
     * 2. The log likelihood of the model given the exposure times
     * 3. The log likelihood of the model given the titre values
     * 4. The log likelihood of the model given the exposure times
     * 
     */
    double evalLogPosterior(const VectorXd& param, const VectorXd& jump, const MatrixXd& covariance, const RObject& dataList, bool init = false) {
        if (parent->onDebug) Rcpp::Rcout << "In: evalLogPosterior" << std::endl;
    //  this->finddifferInf();
    //  this->updateLists();
        long double logPriorPars, logPriorJump, logLikelihood_ab, logPriorExpTime;

        // Evaluate the log prior of the model, CHECK WHERE THIS IS DEFINED
        logPriorPars = parent->evaluateLogPrior(param, jump, dataList);// + log(1.0 / (this->propInferredExpN + 1));

        // If the log prior is infinite, return a log probability of zero
        if (isinf(logPriorPars)) {
            return log(0);
        }

        // Prior distribution on number of infected persons given number of exposed people, this is the evlaution of a betabinomial(E, 1, 1) distribution 
        double long M_adj, k_adj, k2_adj, i_adj, N_adj;
        VectorXd x_vec(3);
        VectorXd alpha_vec(3);
        alpha_vec << 1.1, 1.1, 1.0;
        if (init) {
            if (parent->currInferredExpN == 0) {
                //logPriorInf = -1e10;
                //logPriorExp = -1e10;
                logPriorJump = -1e10;
            }  else {
                // exposure priors
                M_adj = parent->N - parent->currInferredInfN;
                k_adj = parent->currInferredExpN - parent->currInferredInfN;      
                // infection priors
                k2_adj = parent->currInferredExpN - parent->currInferredInfN ;
                i_adj = parent->currInferredInfN - parent->knownInfsN ;
                N_adj = parent->N - parent->currInferredInfN ;//- this->knownInfsN ;
    
                logPriorJump =  betaPDFlog(k2_adj / N_adj, 1.1, 1.1) + logFactorial(k2_adj) + logFactorial(N_adj - k2_adj) - logFactorial(N_adj );// +
                //  log(1 / (k2_adj + 1)) + logFactorial(i_adj) + logFactorial(k2_adj - i_adj) - logFactorial(k2_adj );
                //dirichlet_pdf(x_vec, alpha_vec) + 
            //    Rcpp::Rcout << "logPriorExpInf curr: " << logPriorExpInf << std::endl;
            //  logPriorExpInf = logFactorial(k_adj) + logFactorial(M_adj - k_adj) - logFactorial(M_adj );
            }
        } else { // || this->propInferredExpN == this->propInferredInfN
            if (parent->propInferredExpN == 0 ) {
               // logPriorInf = -1e10;
                //logPriorExp = -1e10;
                logPriorJump = -1e10;
            }  else {
                // exposure priors
                M_adj = parent->N - parent->propInferredInfN;
                k_adj = parent->propInferredExpN - parent->propInferredInfN;
                    // infection priors
                k2_adj = parent->propInferredExpN - parent->propInferredInfN ;
                i_adj = parent->propInferredInfN - parent->knownInfsN ;
                N_adj = parent->N - parent->propInferredInfN ;//- this->knownInfsN ;
                /*Rcpp::Rcout << "M_adj : " << M_adj << std::endl;
                Rcpp::Rcout << "k_adj : " << k_adj << std::endl;
                Rcpp::Rcout << "k2_adj : " << k2_adj << std::endl;
                Rcpp::Rcout << "i_adj : " << i_adj << std::endl;*/
            //  x_vec << (this->N - this->propInferredExpN) / static_cast<double>(this->N),  this->propInferredExpN / static_cast<double>(this->N);
                x_vec << (M_adj - k_adj) / N_adj, k_adj / N_adj, i_adj / N_adj;


            //   Rcpp::Rcout << "alpha_vec1: " << alpha_vec[0] << std::endl;
            //    Rcpp::Rcout << "alpha_vec2: " << alpha_vec[1] << std::endl;

            /*  Rcpp::Rcout << "x_vec1: " << x_vec[0] << std::endl;
                Rcpp::Rcout << "x_vec2: " << x_vec[1] << std::endl;
                Rcpp::Rcout << "x_vec2: " << x_vec[2] << std::endl;

                Rcpp::Rcout << "logFactorial(i_adj): " <<  logFactorial(i_adj) << std::endl;
                Rcpp::Rcout << "logFactorial(k_adj): " << logFactorial(k_adj) << std::endl;
                Rcpp::Rcout << "logFactorial(M_adj - k_adj): " << logFactorial(M_adj - k_adj) << std::endl;
                Rcpp::Rcout << "logFactorial(N_adj ): " << logFactorial(N_adj ) << std::endl;

                Rcpp::Rcout << "dirichlet_pdf(x_vec, alpha_vec): " << dirichlet_pdf(x_vec, alpha_vec) << std::endl;*/
                //logPriorExpInf = 0;
            // logPriorExpInf = logFactorial(k_adj) + logFactorial(M_adj - k_adj) - logFactorial(M_adj);
            // log(1.0 / (N_adj + 1))  +  log(1.0 / (k2_adj + 1))
            //  beta_logpdf(k2_adj / N_adj, 1.1, 1.1)
                logPriorJump = betaPDFlog(k2_adj / N_adj, 1.1, 1.1) + logFactorial(k2_adj) + logFactorial(N_adj - k2_adj) - logFactorial(N_adj ) ;//+
                //  log(1 / (k2_adj + 1)) + logFactorial(i_adj) + logFactorial(k2_adj - i_adj) - logFactorial(k2_adj );
        //     Rcpp::Rcout << "logPriorExpInf prob: " << logPriorExpInf << std::endl;

            //logPriorExpInf = log(1.0 / N_adj) + log(1.0 / k2_adj) + logFactorial(i_adj) + logFactorial(k_adj) + logFactorial(M_adj - k_adj) - logFactorial(N_adj ); // dirichlet_pdf(x_vec, alpha_vec) +
            //  Rcpp::Rcout << "logPriorExpInf: " << logPriorExpInf << std::endl;
            }
        }
               // If the log prior is infinite, return a log probability of zero
        if (isinf(logPriorJump) || isnan(logPriorJump)) {
            return log(0);
        }

        // Evaluate the priors of the model given the exposure times
        logPriorExpTime = 0; 
        for (int i = 0; i < parent->N; i++) {
            if (jump[i] > -1) {
                logPriorExpTime += parent->exposureFunctionDensity(jump[i]);
            }
        }
            
        // Evaluate the log likelihood of the model given the exposure times and titre values
        logLikelihood_ab = this->evaluateLogLikelihoodCOP_cpp(jump, init) + this->evaluateLogLikelihoodObs_cpp(param, jump, init);


        // Evaluate the log likelihood of the model given the titre values
        return logPriorPars + logPriorJump + logLikelihood_ab + logPriorExpTime;
    }

private:

/**
 * @brief Evaluate the log likelihood of the COP part of the model
 * @param jump The jump vector for the exposure times
 * @param init A boolean to indicate if the function is being called for the initialisation step
 * @return The log likelihood of the COP part of the model
 * 
 */
    double evaluateLogLikelihoodCOP_cpp( VectorXd jump, bool init) {
        if (parent->onDebug) Rcpp::Rcout << "In: evaluateLogLikelihoodCOP_cpp" << std::endl;

        // Evaluate the titre values for all individuals given the new exposure times
        MatrixXd titreExp(parent->N, parent->B);
        double titre_est;

        double ll = 0;
        // Evaluate the titre values for all individuals given the new exposure times
        for (int i_idx = 0; i_idx < parent->N; ++i_idx) {
            // If no exposure time is provided, set the titre value to -1
            if (jump[i_idx] == -1) {
                for (int bio = 0; bio < parent->B; bio++) {
                    titreExp(i_idx, bio) = -1;
                }
            } else {
                // Extract the titre values for the individual across all their event times
                //List titre_list_i = parent->titre_list[i_idx];
                std::vector<std::vector<DoubleWithString> > proposalTitreFull_i = parent->proposalTitreFull[i_idx]; 
                // For all biomarkers, extract the titre value for the individual at the new exposure time
                for (int bio = 0; bio < parent->B; bio++) {
                    std::vector<DoubleWithString> proposalTitreFull_i_b = proposalTitreFull_i[bio];
                    for (int j = 0; j < proposalTitreFull_i_b.size(); j++) {
                        // If the titre value is the one at the fitted exposure type, store it
                        if(proposalTitreFull_i_b[j].name == parent->exposureFitted) {
                            titre_est = proposalTitreFull_i_b[j].value;
                        }
                    }
                    // Store the titre value for the individual at the new exposure time
                    titreExp(i_idx, bio) = titre_est;
                }
            }
        }
        if (init) {
            parent->currentTitreExp = titreExp;
        } else {
            parent->proposalTitreExp = titreExp;
        }
        return ll;
    }

/** 
 * @brief Evaluate the log likelihood of the observation model
 * @param params The parameter vector for the fitted parameters
 * @param jump The jump vector for the exposure times
 * @param init A boolean to indicate if the function is being called for the initialisation step
 * @return The log likelihood of the observation model
 */
    double evaluateLogLikelihoodObs_cpp(const VectorXd& params, const VectorXd& jump, bool init) {
        if (parent->onDebug) Rcpp::Rcout << "In: evaluateLogLikelihoodObs_cpp" << std::endl;

        std::vector<DoubleWithString> df_order_exp; 
        MatrixXd obsTitre(parent->N_data, parent->B);
    
        double titre_est ;
        double titre_val;
        double ll = 0;


        int k_idx = 0, k_count = 0;
        for (int i_idx = 0; i_idx < parent->N; ++i_idx) {
            List titre_list_i = parent->titre_list[i_idx];
            std::vector<std::vector<DoubleWithString> > proposalTitreFull_i = parent->proposalTitreFull[i_idx];

            for (int bio = 0; bio < parent->B; bio++) {
                k_count = 0;
                string biomarker_b = parent->biomarkers[bio];
                Function evalLoglikelhoodObs_i = parent->evalLoglikelhoodObs[biomarker_b];
                NumericVector pars = parent->currentParsObs[biomarker_b];

                NumericVector titre_val_i_b = titre_list_i[bio]; 
                int j_data = 0;
                std::vector<DoubleWithString> proposalTitreFull_i_b = proposalTitreFull_i[bio];
                for (int j = 0; j < proposalTitreFull_i_b.size(); j++) {

                    if(proposalTitreFull_i_b[j].name == "bleed") {
                        titre_est = proposalTitreFull_i_b[j].value;
  
                        k_count++;
                        titre_val = titre_val_i_b[j_data];
                        j_data ++;
                        ll += as<double>(evalLoglikelhoodObs_i(titre_val, titre_est, pars) );
                    }
                }
            }
            k_idx = k_idx + k_count;        
        }
        if (init) {
            parent->currentObsTitre = obsTitre;
        } else {
            parent->proposalObsTitre = obsTitre;
        }
        return ll;
    }

    // A shared pointer to the SeroJumpBase object
    std::shared_ptr<SeroJumpBase> parent;
};


class SeroJumpRun : public SeroJumpBase {
/*public:
    std::shared_ptr<SeroJumpBase> base;
    std::unique_ptr<RecalculateTitres> recalInit;
    std::unique_ptr<EvaluateLogLikelihood> loglikInit;

    using SeroJumpBase::SeroJumpBase;

    static std::shared_ptr<SeroJumpRun> create(Rcpp::List infoModel, Rcpp::List observationalModel, Rcpp::List abkineticsModel) {
        // Use shared_ptr to create the instance
        // Rcpp::Rcout << "Hello from create" << std::endl;
        std::shared_ptr<SeroJumpRun> instance(new SeroJumpRun(infoModel, observationalModel, abkineticsModel));
        instance->initialize(); // Pass the shared_ptr instance to initialize
        // Rcpp::Rcout << "Hello from create" << std::endl;

        return instance;
    }

    void initialize() {
        Rcpp::Rcout << "Instance of the previous classes" << std::endl;
        auto base = shared_from_this();  // Safe to use here
        recalInit = std::make_unique<RecalculateTitres>(base);
        loglikInit = std::make_unique<EvaluateLogLikelihood>(base);
    }*/



   private:
   
    using SeroJumpBase::SeroJumpBase;

// Define the shared pointer for the previous classes
    std::shared_ptr<SeroJumpBase> base;
    std::unique_ptr<RecalculateTitres> recalInit;
    std::unique_ptr<EvaluateLogLikelihood> loglikInit;
public:
    /**
     * @brief Construct a new SeroJumpRun object
     * @param infoModel The information model
     * @param observationalModel The observational model
     * @param abkineticsModel The abkinetics model
     * 
    */
  /*  SeroJumpRun(List infoModel, List observationalModel, List abkineticsModel)
        : SeroJumpBase(infoModel, observationalModel, abkineticsModel) {
            Rcpp::Rcout << "Instance of the previous classes" << std::endl;
          //  base = shared_from_this();
           // recalInit = std::make_unique<RecalculateTitres>(base);
            //loglikInit = std::make_unique<EvaluateLogLikelihood>(base);
    }*/

    static std::shared_ptr<SeroJumpRun> create(Rcpp::List infoModel, Rcpp::List observationalModel, Rcpp::List abkineticsModel) {
        // Use shared_ptr to create the instance
        Rcpp::Rcout << "Instance of the previous classes1" << std::endl;
        std::shared_ptr<SeroJumpRun> instance(new SeroJumpRun(infoModel, observationalModel, abkineticsModel));
        instance->initialize(instance);
        Rcpp::Rcout << "Instance of the previous classes2" << std::endl;
        return instance;
    }

    void initialize(std::shared_ptr<SeroJumpRun> self) {
        Rcpp::Rcout << "Instance of the previous classes" << std::endl;
        auto base = this->shared_from_this();  // Safe to use here
        recalInit = std::make_unique<RecalculateTitres>(base);
        loglikInit = std::make_unique<EvaluateLogLikelihood>(base);
        // recalInit = std::make_unique<RecalculateTitres>(base);
        // loglikInit = std::make_unique<EvaluateLogLikelihood>(base);
    }
    /**
     * @brief Initialise the run class
     * @param settings The settings list
     * @param dataList The data list
     * @param i The iteration number
     * 
     */
    void initialiseRunClass(List settings, RObject dataList, int i) {
        // Call the base class initialiser
        this->initialiseBaseClass(settings, dataList, i);

        base = shared_from_this();
        recalInit = std::make_unique<RecalculateTitres>(base);
        loglikInit = std::make_unique<EvaluateLogLikelihood>(base);

        // Define parameters to be defined in function
        VectorXd initialSample, initialJump, initialInf;
        double initialLogLikelihood;
        MatrixXd initialCovarianceMatrix;
        
        if (this->onDebug) Rcpp::Rcout << "In: Check 6" << std::endl;

        // Get parameters values 
        initialSample = this->samplePriorDistributions(this->dataList);

        // Get initial Jump values, ONCE FUNCTION
        if (knownExpInd) {
            initialJump = this->knownExpVec;
        } else {
            if (this->onDebug) Rcpp::Rcout << "In: Check 7i" << std::endl;
            // this gets all the known infections
            initialJump = this->knownInfsTimeVec;
        }

        // Get initial Inf values
        initialInf = this->knownInfsVec;

        if (this->onDebug) Rcpp::Rcout << "In: Check 7ii" << std::endl;
        for (int i = 0; i < this->N; i++) {
            if (this->knownInfsVec(i) == 0) {
                if (this->endTitreValue(i, 0) - this->initialTitreValue(i, 0) > 1) {
                    if (!knownExpInd) {
                        initialJump(i) = this->exposureFunctionSample();
                    }
                    initialInf(i) = 1;
                } else {
                    /*double u = uniformContinuousDist(0, 1);
                    if (u < 0.5) {
                        Rcpp::Rcout << "initialJump(i): randomly added exposure (not infection)" << initialJump << std::endl;
                        if (!knownExpInd) {
                            initialJump(i) = this->exposureFunctionSample();
                        }
                        initialInf(i) = 0;
                    } else {
                        if (!knownExpInd) {
                            initialJump(i) = -1;
                        }
                        initialInf(i) = 0;
                    }*/
                }
            }
        }

        this->currentJump = initialJump;
    //      this->currentInf = initialInf;
        // Get initial conditions for the infection states
        //initialInf = initInfFunc(); 
        if (this->onDebug) Rcpp::Rcout << "In: Check 7" << std::endl;

        // MAYBE MOVE THIS TO INIT CLASS
        //  updateAbKineticParams(initialSample);
       // std::shared_ptr<SeroJumpBase> base = shared_from_this();
       // RecalculateTitres recalInit(base);

        recalInit->updateEventsFull();
        recalInit->updateAbKineticParams(initialSample);
        recalInit->recalculateTitreAll();
        this->currentTitreFull = this->proposalTitreFull;
        this->currentEventsFull = this->proposalEventsFull;

        // Get initial conditions for the Loglikelihood probabilities
        /*   std::vector<DoubleWithString> df_order_exp_i;
        std::vector<std::vector<DoubleWithString> > df_order_titre_i;
        for (int i_idx = 0; i_idx < this->N; i_idx++) {

            df_order_exp_i = this->sortevents(i_idx, initialJump);
            int l = df_order_exp_i.size();
            this->currentEventsFull.push_back(df_order_exp_i);
            this->proposalEventsFull.push_back(df_order_exp_i);

            df_order_titre_i = this->recalculateTitre_i(df_order_exp_i, i_idx);
            //this->currentTitreFull.insert(this->currentTitreFull.begin() + i_idx, df_order_titre_i);
            this->currentTitreFull.push_back(df_order_titre_i);
            this->proposalTitreFull.push_back(df_order_titre_i);
        }*/

        this->proposalJump = this->currentJump;
        this->bookKeepingSize();
        this->currInferredExpN = this->propInferredExpN;
        
        this->currentCovarianceMatrix = this->nonadaptiveScalar*this->nonadaptiveCovarianceMat;
      //  updateAbKineticParams(initialSample);

        initialLogLikelihood = loglikInit->evalLogPosterior(initialSample, initialJump, this->currentCovarianceMatrix, this->dataList, true);

        if (this->onDebug) Rcpp::Rcout << "In: Check 8" << std::endl;

        while(isinf(initialLogLikelihood) || isnan(initialLogLikelihood)){
            if (this->onDebug) Rcpp::Rcout << "In: STUCK 8" << std::endl;

            initialSample = this->samplePriorDistributions(this->dataList);
            recalInit->updateAbKineticParams(initialSample);
            recalInit->recalculateTitreAll();

            initialLogLikelihood = loglikInit->evalLogPosterior(initialSample, initialJump, this->currentCovarianceMatrix, this->dataList, true);
        }
        if (this->onDebug) Rcpp::Rcout << "In: Check 9" << std::endl;
        this->currentSample = initialSample;
        this->currentSampleMean = initialSample;
        this->currentLogPosterior = initialLogLikelihood;

        this->bookKeepingSize();

    }

    /** 
     * @brief Run the RJMC chain
     * @return A list of the posterior values
     */
    List runRJMCC(List settings, RObject dataList, int i)
    {
        initialiseRunClass(settings, dataList, i);
        // initalise FUNCTION
        if (this->onDebug) Rcpp::Rcout << "Pre: iterations" << std::endl;
        for (int i = 0; i < this->iterations; i++){
            
            this->workingIteration = i;
            //if (this->onDebug) Rcpp::Rcout << "Pre: updateAllChains" << std::endl;
            runIteration();
        }

        // List out of values )thinned)
        List out = Rcpp::List::create(
            _["pars"] = this->posteriorOut,
            _["jump"] = this->posteriorJump,
            _["titreexp"] = this->posteriorTitreExp,
            _["obstitre"] = this->posteriorObsTitre
        );

        return out;
    }

    /**
     * @brief Update the markov chain
     */
    void runIteration()
    {
        if (this->onDebug) Rcpp::Rcout << "Pre: updateRJMC" << std::endl;
        updateMarkovChain();            
        if ((this->currInferredExpN != this->knownInfsN) & !this->knownExpInd) { 
            this->updateTiming();
        }

        this->counterFuncEval++;

        if (this->onDebug) Rcpp::Rcout << "Pre: updateOutputPosterior" << std::endl;
        this->updateOutputPosterior();
        if (this->onDebug) Rcpp::Rcout << "Pre: updateProposal" << std::endl;
    
        // Outputs
        this->consoleUpdatefunction();
    }
    
    /** 
     * @brief Update the RJMC chain
     */
    void updateMarkovChain() {
        if (this->onDebug) Rcpp::Rcout << "Pre: getAcceptanceRate" << std::endl;
        getAcceptanceRate();

        if (this->onDebug) Rcpp::Rcout << "Pre: updateSampleAndLogPosterior" << std::endl;
        this->updateSampleAndLogPosterior();
        if(this->conPropIn & (this->currJumpType == 1)) {
            this->updateProposal();
        }
    }    

    /**
     * @brief Get the acceptance rate
     * 
     */
    void getAcceptanceRate()
    {        
        this->isSampleAccepted = false;
        // Evaluate the ab kinetics parameters and do rjmcmc
        if (this->onDebug) Rcpp::Rcout << "Pre: selectProposalDist" << std::endl;
        this->fixedProposalDist();
        if (this->onDebug) Rcpp::Rcout << "Pre: JumpProposalDist" << std::endl;
        this->JumpProposalDist();
        if (this->onDebug) Rcpp::Rcout << "Pre: evalLogPosterior" << std::endl;

        // Calculate the observation titres + titre at infection
        // CALL CLASS HERE
        recalInit->updateEventsFull();
        recalInit->updateAbKineticParams(this->proposalSample);
        recalInit->recalculateTitreAll();      

        // CALL CLASS HERE
        this->proposedLogPosterior = loglikInit->evalLogPosterior(this->proposalSample, this->proposalJump, this->currentCovarianceMatrix, this->dataList);
        if (this->onDebug) Rcpp::Rcout << "Pre: evaluateMetropolisRatio" << std::endl;

        evaluateMetropolisRatio();
    }

    /**
     * @brief Generate a sample from the proposal distribution
     * 
     */
    void fixedProposalDist(bool update = true){
        //  if (this->onDebug) Rcpp::Rcout << "Post: selectProposalDist" << std::endl;

        if (this->workingIteration < this->burninAdaptiveCov || uniformContinuousDist(0, 1) < 0.05 || !this->onAdaptiveCov ){
            generateSampleFromNonAdaptiveProposalDist(update);
        }
        else{
            generateSampleFromAdaptiveProposalDist(update);
        }
    }
    
    /**
     * @brief Generate a sample from the non-adaptive proposal distribution
     * 
     * update: counterNonAdaptive, currentCovarianceMatrix, proposalSample
     */
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
    
    /**
     * @brief Generate a sample from the adaptive proposal distribution
     * 
     * update: counterAdaptive, currentCovarianceMatrix, proposalSample
     */
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
    
    /**
     * @brief Evaluate the Metropolis ratio
     * 
     * Update: alpha
     */
    void evaluateMetropolisRatio()
    {
        double rjadjustmentFactor = 0;
      //  double rjadjustmentFactor2 = 0;
        if(std::isnan(this->proposedLogPosterior) || std::isinf(this->proposedLogPosterior)) {
            this->alpha = 0;
        } else {
            if (this->currJumpType == 0) {
                // in DEATH want to evaluate current state to the proposed parameters
                double N_adj = this->N - this->currInferredInfN;
                double curr_adj = this->currInferredExpN - this->currInferredInfN;

                rjadjustmentFactor = log(curr_adj) - log((N_adj - curr_adj + 1) ) +
                    this->exposureFunctionDensity(this->currentJump(this->currJumpIdx)) + log(0.5); 
                    //  rjadjustmentFactor2;
            } else if (this->currJumpType == 1) {
                rjadjustmentFactor = 0;
            } else if (this->currJumpType == 2) {
                // in birth want to evaluate the Proposed jump state given the current parameters
                double N_adj = this->N - this->currInferredInfN;
                double curr_adj = this->currInferredExpN - this->currInferredInfN;// - this->knownInfsN;

                rjadjustmentFactor = log(N_adj - curr_adj) - log(curr_adj + 1) - 
                    this->exposureFunctionDensity(this->proposalJump(this->currJumpIdx)) - log(0.5);// - 
                    // rjadjustmentFactor2; 
            }
            this->alpha = min(1.0, exp((this->proposedLogPosterior - this->currentLogPosterior + rjadjustmentFactor)));
        }
    }

    /**
     * @brief Update the sample and log posterior
     */
    void updateSampleAndLogPosterior()
    {
        if (uniformContinuousDist(0, 1) < this->alpha) {
            if (this->onDebug) Rcpp::Rcout << "In: ACCEPTED" << std::endl;
            //  this->finddifferInf();
            this->isSampleAccepted = true; this->counterAccepted++;

            this->currentSample = this->proposalSample;
            this->currentJump = this->proposalJump;

            this->currentTitreExp = this->proposalTitreExp;
            this->currentObsTitre = this->proposalObsTitre;

            this->currentEventsFull = this->proposalEventsFull;
            this->currentTitreFull = this->proposalTitreFull;

            this->currentLogPosterior = this->proposedLogPosterior;
            this->currInferredExpN = this->propInferredExpN;
            this->currInferredInfN = this->propInferredInfN;


        } else {
            this->proposalSample = this->currentSample;
            this->proposalJump = this->currentJump;

            this->proposalTitreExp = this->currentTitreExp;
            this->proposalObsTitre = this->currentObsTitre;

            this->proposalEventsFull = this->currentEventsFull;
            this->proposalTitreFull = this->currentTitreFull;

            this->proposedLogPosterior = this->currentLogPosterior;
            this->propInferredExpN = this->currInferredExpN;
            this->propInferredInfN = this->currInferredInfN;                
        }
    }

    /**
     * @brief Update the timing of the RJMC chain
     */
    void updateTiming() {
        if (this->onDebug) Rcpp::Rcout << "In: updateTiming" << std::endl;

            this->currJumpType = 1;
        
            this->proposalJump = this->currentJump;

            for (int i = 0; i < this->noGibbsSteps; i++) { 
                this->gibbsIdx = this->sampleExposed(); // stuck here
                resampleTime(this->gibbsIdx);
                this->gibbsIdx = this->sampleExposed(); // stuck here
                //   proposeInfection(this->gibbsIdx);
            }


            fixedProposalDist(false);

            recalInit->updateEventsFull();
            recalInit->updateAbKineticParams(this->proposalSample);
            recalInit->recalculateTitreAll();    

            this->proposedLogPosterior = loglikInit->evalLogPosterior(this->proposalSample, this->proposalJump, this->currentCovarianceMatrix, this->dataList);

            evaluateMetropolisRatio();
            updateJumpHistoric();
            updateProposalGibbs();
    }

    /**
     * @brief Update the proposal for the Gibbs sampler
     */
    void consoleUpdatefunction()
    {
        int i = this->workingIteration;
        if(i%this->consoleUpdates == 0) {
            Rcpp::Rcout << "Running MCMC-PT iteration number: " << this->workingIteration << " of " <<  this->iterations << ". Chain number " << this->chainNumber << ". Current logpost: " << this->currentLogPosterior << ". No exp/inf: " <<  this->currInferredExpN << "/" << this->currInferredInfN << ". \n";
            if (this->onDebug) {
                Rcpp::Rcout << "\n Current values: " << this->currentSample << std::endl;
            }
        }
    }

    /**
     * @brief Update the jump historic
     */
    void updateJumpHistoric()
    {
        if (uniformContinuousDist(0, 1) < this->alpha) {
            if (this->onDebug) Rcpp::Rcout << "In: ACCEPTED JUMP" << std::endl;
            //Rcpp::Rcout << "ACCEPTED (GIBBS): " << this->alpha << std::endl;
        //     this->finddifferInf();
            this->currentSample = this->proposalSample;
            this->currentJump = this->proposalJump;

            this->currentTitreExp = this->proposalTitreExp;
            this->currentObsTitre = this->proposalObsTitre;

            this->currentEventsFull = this->proposalEventsFull;
            this->currentTitreFull = this->proposalTitreFull;

            this->currentLogPosterior = this->proposedLogPosterior;
            this->currInferredExpN = this->propInferredExpN;
            this->currInferredInfN = this->propInferredInfN;

        } else {
            this->proposalSample = this->currentSample;
            this->proposalJump = this->currentJump;

            this->proposalTitreExp = this->currentTitreExp;
            this->proposalObsTitre = this->currentObsTitre;

            this->proposalEventsFull = this->currentEventsFull;
            this->proposalTitreFull = this->currentTitreFull;

            this->proposedLogPosterior = this->currentLogPosterior;
            this->propInferredExpN = this->currInferredExpN;
            this->propInferredInfN = this->currInferredInfN;     
        }
    }

    /** 
     * @brief Update the output posterior
     * 
     */
    void updateOutputPosterior()
    {
        // Only save once after the burnin and if interation is divisible by the thinning
        if ((this->workingIteration > (this->burninPosterior-1)) && 
            (this->workingIteration%thin == 0) ) {

            // Get parameters
            for (int p = 0; p < this->numberFittedPar; p++)
                this->posteriorOut(this->counterPosterior, p) = this->currentSample(p);
            
            // Get the log posterior
            this->posteriorOut(this->counterPosterior, this->numberFittedPar) = this->currentLogPosterior;
            this->posteriorOut(this->counterPosterior, this->numberFittedPar+1) = (double)this->counterAccepted/(double)this->counterFuncEval;

            // Get the jump values
            for (int j = 0; j < this->lengthJumpVec; j++) {
                this->posteriorJump(this->counterPosterior, j) = this->currentJump(j);
            }
            // Get the titre at exposure values and observational values
            this->posteriorTitreExp.push_back(this->currentTitreExp);
            this->posteriorObsTitre.push_back(this->currentObsTitre);

            this->counterPosterior++;            
        }
    }
    
    /**
     * @brief Update the proposal distributions of the adaptive covariance matrix
     */
    void updateProposal()
    {       
        // Update adaptive proposal stuff   
        if (this->isProposalAdaptive){
            this->adaptiveScalar += this->stepSizeRobbinsMonro * pow(1+this->counterAdaptive,-0.5)*(this->alpha - 0.234);
            // Checks on value
            errorCheckNumberValid(this->adaptiveScalar);
            trimAdaptiveValues(this->adaptiveScalar);
        }
        // Update non-adaptive proposal stuff
        else{   
            this->nonadaptiveScalar += this->stepSizeRobbinsMonro * pow(1+this->counterNonAdaptive,-0.5)*(this->alpha - 0.234);
            // Checks on value 
            errorCheckNumberValid(this->nonadaptiveScalar);
            trimNonAdaptiveValues(this->nonadaptiveScalar);
        }

        // Update mean and covariance matric
        if(((this->workingIteration) % (this->updatesAdaptiveCov) == 0) && (this->workingIteration > this->burninAdaptiveCov)){

            this->iPosterior++;
            double gainFactor = pow(1+iPosterior, -0.5);
            if (this->iPosterior == 1){
                this->currentSampleMean = this->currentSample;
                this->adaptiveScalar = this->nonadaptiveScalar;
            }
            else{
                this->currentSampleMean = this->currentSampleMean + gainFactor*(this->currentSample-this->currentSampleMean);
                errorCheckVectorValid(this->currentSampleMean); // Check if the vector is valid
                this->adaptiveCovarianceMat = this->adaptiveCovarianceMat + gainFactor*((this->currentSample-this->currentSampleMean)*((this->currentSample)-(this->currentSampleMean)).transpose()) - gainFactor*this->adaptiveCovarianceMat;
                errorCheckMatrixValid(this->adaptiveCovarianceMat); // Check if the matrix is valid
            }
        }
    }

    /**
     * @brief Update the proposal distributions associated with the exposure times
     */
    void updateProposalGibbs()
    {
        // Update adaptive proposal stuff
        this->counterAdaptiveGibbs(this->gibbsIdx)++;
        double gainFactor = pow(this->counterAdaptiveGibbs(this->gibbsIdx), -0.5);
        this->adaptiveGibbsSD(this->gibbsIdx) += gainFactor*(this->alpha - 0.44);

        // Checks on value
        errorCheckNumberValid(this->adaptiveGibbsSD(this->gibbsIdx));
        trimAdaptiveValues(this->adaptiveGibbsSD(this->gibbsIdx));
    }
    
    /** 
     * @brief Sample the exposure function
     * 
     * update : currJumpType
     * 
     */
    void JumpProposalDist() {
        if (this->onDebug) Rcpp::Rcout << "Redfine prop = curr" << std::endl;

        this->proposalJump = this->currentJump;
        //  double r;
        VectorXd q_prob(3);
        if (this->onDebug) Rcpp::Rcout << "get q_prob" << std::endl;
        if (this->onDebug) Rcpp::Rcout << "this->knownInfsN: " << this->knownInfsN << std::endl;

        // Definie probability of jumping to death, birth or stay same
        if (this->knownExpInd) {
            q_prob << 0, 1, 1;
        } else {
            if (this->currInferredExpN == this->currInferredInfN ) {
                q_prob << 0, 0.67, 1.0;
                //q_prob << 0.0, 0.5, 0.75, 1.00;
            } else if (this->currInferredExpN == (this->N ) ) {
                //q_prob << 0.25, 1.0, 1.0, 1.0;
                q_prob << 0.33, 1.0, 0;
            } else {
                q_prob << 0.33, 0.67, 1.0;
                // q_prob << 0.25, 0.5, 0.75, 1.0;
            }
        }

        // Substract a exposure 
        if (this->onDebug) Rcpp::Rcout << "Find jumping" << std::endl;
        double q = uniformContinuousDist(0, 1);
        if (q < q_prob(0)) {
            //   Rcpp::Rcout << "In deletion" << std::endl;
            this->currJumpIdx = this->sampleExposedNotInf(); // length is n_ - 114
            this->historicJump(this->currJumpIdx) = this->currentJump(this->currJumpIdx);
            this->proposalJump(this->currJumpIdx) = -1; // length now n_ - 144 - 1
            this->currJumpType = 0;
            this->bookKeepingSize();
            //this->propInferredExpN = this->currInferredExpN - 1;
        // Stay same
        } else if (q < q_prob(1)) {
            // Rcpp::Rcout << "In stay same" << std::endl;
            // this->propInferredExpN = this->currInferredExpN;
            if(this->currInferredExpN != this->knownInfsN) {
                // this->currJumpIdx = uniformDiscreteDist(0, this->N - 1);
                // this->currJumpIdx = this->sampleExposed(); // length is n_ - 114
                // RESAMPLE INFECTION
                // this->proposeInfection(this->currJumpIdx); // always rejecting
            }
            this->currJumpType = 1;
        // Add an infection same
        } else if (q < q_prob(2)) {
            this->currJumpIdx = this->sampleNotExposed(); // length is 220 - n_
            this->sampleNewTime(this->currJumpIdx, 0);
            this->currJumpType = 2;
        }
        // Add an exposure
        else if (q < q_prob(3)) {
            //  Rcpp::Rcout << "In addition same" << std::endl;
            this->currJumpIdx = this->sampleNotExposed(); // length is 220 - n_
            this->sampleNewTime(this->currJumpIdx, 1);
            // proposalInf(this->currJumpIdx) = 0;
            this->currJumpType = 2;
            //this->propInferredExpN = this->currInferredExpN + 1;
        }
        if (this->onDebug) Rcpp::Rcout << "END jump sample" << std::endl;
    }

    /**
     * @brief Sample an individual who is exposed in current chain
     * @return The individual index
     */
    int sampleExposed() {
        if (this->onDebug) Rcpp::Rcout << "In: sampleExposed" << std::endl;

        int s = uniformDiscreteDist(0, this->N - 1); 
        while ((this->knownInfsVec(s) == 1) || (this->currentJump(s) == -1) ) {
            s = uniformDiscreteDist(0, this->N - 1); 
        }
        return s;
    }

    /**
     * @brief Sample an individual who is exposed in current chain
     * @return The individual index
     */
    int sampleExposedNotInf() {
        if (this->onDebug) Rcpp::Rcout << "In: sampleExposed" << std::endl;

        int s = uniformDiscreteDist(0, this->N - 1); 
        while (this->currentJump(s) == -1) {
            s = uniformDiscreteDist(0, this->N - 1); 
        }
        return s;
    }

    /**
     * @brief Sample an individual who is not exposed in current chain
     * @return The individual index
     */
    int sampleNotExposed() {
        int s = uniformDiscreteDist(0, this->N - 1); 
        // Must be not exposured or infectio not known
        while ((this->knownInfsVec(s) == 1) || (this->currentJump(s) != -1 )) {
            /* if (this->knownInfsVec(s) == 1) {
                Rcpp::Rcout << "this->knownInfsVec(s): " << s << std::endl;
                Rcpp::Rcout << "this->currentJump(s): " << this->currentJump(s) << std::endl;
            }*/
            s = uniformDiscreteDist(0, this->N - 1); 
        }
        return s;
    }

    /**
     * @brief Propose an infection
     * @param t The individual index
     * 
     */
    void proposeInfection(int t) {
        if (this->onDebug) Rcpp::Rcout << "In: proposeInfection" << std::endl;

        //  if (this->currInferredInfN <= 1) {
        //      this->proposalInf(t) = 1;
        //   } else if (this->currInferredInfN >= (this->currInferredExpN - 1)) {
            //  boost::random::bernoulli_distribution<> l(0.5); 
            //  this->proposalInf(t) = l(rng);
        //   } else {
            boost::random::bernoulli_distribution<> l(0.5); 
            //  this->proposalInf(t) = l(rng);
        //   }
        bookKeepingSize();
    }

    /**
     * @brief Resample the exposure time
     * @param t The individual index
     * 
     * Updates: proposalJump(t)
     */
    void resampleTime(int t) {
        if (this->onDebug) Rcpp::Rcout << "In: resampleTime" << std::endl;
        double r = uniformContinuousDist(0, 1);

        // Completely new sample time is drawn 5% of the time for an individual.
        if (r < 0.05) {
            this->proposalJump(t) = this->exposureFunctionSample();
            while ((this->proposalJump(t) >= this->endTitreTime(t) - 7) || (this->proposalJump(t) < this->initialTitreTime(t))) {
                // check this
                this->proposalJump(t) = this->exposureFunctionSample(); 
            }
            // proposeInfection(t);

        // An exisiting sample time is updated 95% of the time for an individual.
        } else {
            double temp = this->currentJump(t) + normalDistSample(0, exp(this->adaptiveGibbsSD(t)));
            if (temp >= this->endTitreTime(t) - 7) {
                temp = this->endTitreTime(t) - 7; 
            } else if (temp < this->initialTitreTime(t)) {
                temp = this->initialTitreTime(t); 
            }
            this->proposalJump(t) = temp;
        }
        bookKeepingSize();
    }

    /** 
     * @brief Sample a new exposure time 
     * @param t The individual index
     * @param prop The proposal index
     * 
     * Updates: proposalJump(t)
     * 
     */
    void sampleNewTime(int t, int prop) {
        if (this->onDebug) Rcpp::Rcout << "In: sampleNewTime" << std::endl;

        double r = uniformContinuousDist(0, 1);

        if ((this->historicJump(t) < 0) | (r < 0.05)) { 
            // New exposure time is resampled from scratch, only do if first time individual exposed in mcmc chain or 5% of time thereafter.
                this->proposalJump(t) = this->exposureFunctionSample(); 
                while ((this->proposalJump(t) >= this->endTitreTime(t) - 7) || (this->proposalJump(t) < this->initialTitreTime(t))) {
                    this->proposalJump(t) = this->exposureFunctionSample(); 
                }
                bookKeepingSize();
                //   proposeInfection(t);
        } else { 
            // New exposure time is taken from previous exposure time in mcmc chain (95% of the time)
                this->proposalJump(t) = this->historicJump(t);
                bookKeepingSize();

                //  if (r < 0.5) {
                    //  boost::random::beta_distribution<> b(1, 1); 
                    //  double p = b(rng);
                    //   double p0 = binomialCoefficient(this->propInferredExpN, this->currInferredInfN);
                    //  double p1 = binomialCoefficient(this->propInferredExpN, this->currInferredInfN + 1);
                    //  boost::random::bernoulli_distribution<> l(0.5); 
                    // this->proposalInf(t) = l(rng);
                // }
                bookKeepingSize();
        }
        bookKeepingSize();
    }

   /**
    * @brief Book keeping size of the exposure and infection states
    * 
    * Updates: propInferredExpN
    */
    void bookKeepingSize() {
        int count_exp = 0;
       // int count_inf = 0;
        for (int i = 0; i < this->N; i++) {
            if (this->proposalJump[i] > -1) {
                count_exp++;
            }
        }            
        this->propInferredExpN = count_exp;
    }




};

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