#ifndef RJMC_RUN_HPP
#define RJMC_RUN_HPP

// Load previous headers
#include "./rjmc_inf_update.hpp"

/** 
 * @brief A class to represent the RJMC run, inherits SeroJumpBase
 * 
 */
class SeroJumpRun: public SeroJumpBase {
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


#endif