#ifndef RJMC_INIT_HPP
#define RJMC_INIT_HPP

/** 
 * @brief The base structure for the RJMC implementation
 * 
 * This case class defines the base class for the RJMC implementation. It contains all the necessary functions and variables to run the RJMC algorithm. Mainly extracting data from the model and data list, and initialising the RJMC algorithm.
 * 
 */
struct SeroJumpBase : public std::enable_shared_from_this<SeroJumpBase>
{
    // Vartiables which are initialise through the constructor

    Rcpp::List evalLoglikelhoodObs, evalLoglikelhoodCOP, evalabkineticsFunc; // List of functions to evaluate likelihoods and kinetics
    Rcpp::List infoModel, copModel, observationalModel, abkineticsModel; // List of models and parameters

    List parsAbKinN, parsCOPN, parsObsN; // Lists of all the parameters names for each model
    std::vector<abKineticsInfo> abinfo; // Vector of all the abkinetics information
    std::map<std::string, int> abmap, mapOfObs, mapOfCOP, mapOfExp; // Maps of the abkinetics and observational models

    std::map<std::pair<std::string, std::string>, int> mapOfAbkinetics; // Map of the abkinetics model
    std::vector<string> biomarkers, exposureType, abID; // Vector of biomarkers and exposure types
    string exposureFitted; // The exposure type which is fitted
    List exposureInfo, knownInf; // List of exposure information and known infections
    int B; // Number of biomarkers
    bool priorPredFlag = false;
    bool copFlag = false;
    /** 
     * @brief Constructor for the RJMC base class
     * @param infoModel_in The information model
     * @param observationalModel_in The observational model
     * @param abkineticsModel_in The antibody kinetics model
     * 
     */
    SeroJumpBase(Rcpp::List infoModel_in, Rcpp::List observationalModel_in, Rcpp::List abkineticsModel_in, Rcpp::List copModel_in) : infoModel(infoModel_in), observationalModel(observationalModel_in), abkineticsModel(abkineticsModel_in), copModel(copModel_in) {
        Rcpp::Rcout << "Testing in" << std::endl;
    }    
    
    // Virtual destructor
    virtual ~SeroJumpBase() = default;

    List currentParsObs, currentParsCOP, proposalParsCOP, currentParsAb;
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
    MatrixXd currentCORPars, proposalCORPars;

    VectorXd historicJump;
    NumericVector max_titre;

    vector<vector<DoubleWithString> >  currentEventsFull, proposalEventsFull;
    vector<vector<vector<DoubleWithString> > > currentTitreFull, proposalTitreFull;

    // Outputs for posterior
    MatrixXd posteriorOut, posteriorJump, posteriorInf;
    std::vector<MatrixXd> posteriorTitreExp, posteriorObsTitre, posteriorCORPars;

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
    std::function<double(double, double, double, double)> evaluateLogPriorInfExp;

    std::function<double(int)> exposureFunctionSample;
    std::function<double(double, int)> exposureFunctionDensity;

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
    double counterResampleNo, adaptiveResampleNo;
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

    // **
    //  * @brief Function to initialise the model which was in constructor but was taken out due to issues with pointers
    //  *
    //  * This function extracts the information from the infoModel, mainly biomarkers, fitted exposure and exposure information
    //  */
    void initialiseModel() {
        // Extract the information from the infoModel, mainly biomarkers, fitted exposure and exposure information
        this->biomarkers = as<vector<string> >(this->infoModel["biomarkers"]);
        this->B = this->biomarkers.size();
        this->exposureFitted = as<string>(this->infoModel["exposureFitted"]);
        this->exposureInfo = this->infoModel["exposureInfo"];

        // Extract the information from the exposureInfo, mainly the exposure type, known exposures and the map of the exposure
        for (int i = 0; i < this->exposureInfo.size(); i++) {
            List exposureInfo_i = exposureInfo[i];
            string exposureType_i = exposureInfo_i["exposureType"];
            this->exposureType.push_back(exposureType_i);
            NumericVector known_inf_i = exposureInfo_i["known_inf"];
            this->knownInf[exposureType_i] = known_inf_i;
            this->mapOfExp[exposureType_i] = i;
        }

        // Extract the information from the observational model, mainly the biomarker, logliklihood and the parameters for the loglikelihood and the map
        for (int i = 0; i < this->observationalModel.size(); i++) {
            List temp1 = this->observationalModel[i];
            string temp2 = temp1["biomarker"];
            Function temp3 = temp1["logLikelihood"];
            StringVector temp4 = temp1["pars"];
            this->mapOfObs[temp2] = i;
            this->evalLoglikelhoodObs[temp2] = temp3;
            this->parsObsN[temp2] = temp4;
        }
        Rcpp::Rcout << "Testing in: " << this->copModel.size() << std::endl;
       /* if (this->copModel.size() == 0) {
            this->copFlag = false;
        } else {
            this->copFlag = true;
        }*/

        if (this->copFlag) {
            for (int i = 0; i < this->copModel.size(); i++) {
                List temp1 = this->copModel[i];
                string temp2 = temp1["biomarker"];
                Function temp3 = temp1["logLikelihood"];
                StringVector temp4 = temp1["pars"];
                this->mapOfCOP[temp2] = i;
                this->evalLoglikelhoodCOP[temp2] = temp3;
                this->parsCOPN[temp2] = temp4;
            }
        }

        // Extract the information from the abKinetics model, mainly the biomarker, logliklihood and the parameters for the loglikelihood and the map
        for (int i = 0; i < this->abkineticsModel.size(); i++) {
            List temp1 = this->abkineticsModel[i];
            string temp0 = temp1["id"];
            abID.push_back(temp0);
            Function temp2 = temp1["funcForm"];
            StringVector temp3 = temp1["pars"];
            string temp4 = temp1["biomarker"];
            string temp5 = temp1["exposureType"];
            this->mapOfAbkinetics[std::make_pair(temp4, temp5)] = i;

            this->evalabkineticsFunc[temp0] = temp2;
            this->parsAbKinN[temp0] = temp3;
        }
    }    

    void initialiseBaseClass(List settings, RObject dataList, int i)
    {

        initialiseModel();
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
        this->max_titre = this->dataListCPP["max_titre"];

       // if (this->onDebug) Rcpp::Rcout << "In: Extract data from dataList3" << std::endl;


        this->priorPredFlag = this->dataListCPP["priorPredFlag"];
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
        this->counterResampleNo = 0;
        this->adaptiveResampleNo = 0;
        // Initialise markov chains 
        this->currentLogPosterior = 0;

        this->currInferredInfN = this->knownInfsN;
        this->propInferredInfN = this->knownInfsN;
        Rcpp::Rcout << "this->knownExpInd:" << this->knownExpInd << std::endl;
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


// Wrapping functions to put into the RJMCMC algorithm 
namespace {
    void init_samplePriorDistributions(SeroJumpBase* model, Rcpp::Function samplePriorDistributions) {
        auto func = [samplePriorDistributions](RObject dataList) {
            PutRNGstate();
            auto rData = samplePriorDistributions(dataList);
            GetRNGstate();
            return Rcpp::as<VectorXd>(rData);
        };
        model->samplePriorDistributions = func;
    }

    void init_evaluateLogPrior(SeroJumpBase* model, Rcpp::Function evaluateLogPrior) {
        auto func = [evaluateLogPrior](VectorXd params, VectorXd jump, RObject dataList) {
            PutRNGstate();
            auto rData = evaluateLogPrior(params, jump, dataList);
            GetRNGstate();
            return Rcpp::as<double>(rData);
        };
        model->evaluateLogPrior = func;
    }

    void init_evaluateLogPriorInfExp(SeroJumpBase* model, Rcpp::Function evaluateLogPriorInfExp) {
        auto func = [evaluateLogPriorInfExp](double N, double E, double I, double K) {
            PutRNGstate();
            auto rData = evaluateLogPriorInfExp(N, E, I, K);
            GetRNGstate();
            return Rcpp::as<double>(rData);
        };
        model->evaluateLogPriorInfExp = func;
    }

    void init_initialiseJump(SeroJumpBase* model, Rcpp::Function initialiseJump) {
        auto func = [initialiseJump](RObject dataList) {
            PutRNGstate();
            auto rData = initialiseJump(dataList);
            GetRNGstate();
            return Rcpp::as<VectorXd>(rData);
        };
        model->initialiseJump = func;
    }

    void init_exposureFunctionSample(SeroJumpBase* model, Rcpp::Function exposureFunctionSample) {
        auto func = [exposureFunctionSample](int t) {
            PutRNGstate();
            auto rData = exposureFunctionSample(t);
            GetRNGstate();
            return Rcpp::as<double>(rData);
        };
        model->exposureFunctionSample = func;
    }

    void init_exposureFunctionDensity(SeroJumpBase* model, Rcpp::Function exposureFunctionDensity) {
        auto func = [exposureFunctionDensity](double time, int t) {
            PutRNGstate();
            auto rData = exposureFunctionDensity(time, t);
            GetRNGstate();
            return Rcpp::as<double>(rData);
        };
        model->exposureFunctionDensity = func;
    }
}

#endif