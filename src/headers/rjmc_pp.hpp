#ifndef RJMC_PP_HPP
#define RJMC_PP_HPP

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



namespace rjmc_pp{


    struct RJMC_PP_D
    {
        
       
        RJMC_PP_D() { }

        bool conPropIn = true;
        VectorXd lowerParBounds, upperParBounds;
        double nonadaptiveScalar, adaptiveScalar;
        MatrixXd nonadaptiveCovarianceMat, adaptiveCovarianceMat;
        MatrixXd currentSample, currentSampleMean;

        // Outputs for posterior
        MatrixXd posteriorOut;
        std::vector<MatrixXd > posteriorJump;

        MatrixXd currentCovarianceMatrix;
        double currentLogPosterior;
        VectorXd proposalSample;
        
        MatrixXd currentJump;
        MatrixXd proposalJump;

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


        std::function<VectorXd(RObject)> sampleInitPrior;
        std::function<MatrixXd(VectorXd, RObject)> sampleInitJump;
        std::function<double(VectorXd, MatrixXd, RObject)> evaluateLogPrior;
        std::function<double(VectorXd, MatrixXd, RObject)> evaluateLogLikelihood;
        std::function<MatrixXd(VectorXd, MatrixXd, int, RObject)> sampleBirthProposal;
        std::function<MatrixXd(VectorXd, MatrixXd, int, RObject)> sampleDeathProposal;
        std::function<double(VectorXd, MatrixXd, int, RObject)> evaluateBirthProposal;
        std::function<double(VectorXd, MatrixXd, int, RObject)> evaluateDeathProposal;
        std::function<MatrixXd(VectorXd, MatrixXd, int, RObject)> sampleJump;
        std::function<VectorXd(VectorXd, MatrixXd, RObject)> sampleProposal;

        // Functions for internal 
        Mvn Mvn_sampler;

 
        double stepSizeRobbinsMonro;

        double evalLogPosterior(const VectorXd& param, const MatrixXd& jump, const MatrixXd& covariance, const RObject& dataList, bool init = false)
        {
            if (this->onDebug) Rcpp::Rcout << "In: evalLogPosterior" << std::endl;
          //  this->finddifferInf();
          //  this->updateLists();
        //    Rcpp::Rcout << "LOL1" << std::endl;

            double logPrior = this->evaluateLogPrior(param, jump, dataList);
       //     Rcpp::Rcout << "LOL2" << std::endl;

            if (isinf(logPrior)) {
                return log(0);
            }
            // need to be a converstion here as observationalModel is ObsFuncTemplate, but function call is Function
        //    Rcpp::Rcout << "LOL!" << std::endl;
            double logLikelihood_ab = this->evaluateLogLikelihood(param, jump, dataList);
        //    Rcpp::Rcout << "LOL!1" << std::endl;

            return logPrior + logLikelihood_ab ;
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
        MatrixXd initialTitreValue;
        VectorXd endTitreTime;
        VectorXd times_full, id_full;
        MatrixXd titre_full;
        List titre_list, times_list;

        StringVector exposurePeriods;
        string exposureNameInf;

        void initialiseClass(List settings, RObject dataList, int i)
        {
         //   Rcpp::Rcout << "In: initialiseClass" << std::endl;
            this->onDebug = settings["onDebug"];
         //   Rcpp::Rcout << "In: initialiseClass" << std::endl;
            if (this->onDebug) Rcpp::Rcout << "In: initialiseClass" << std::endl;

            if (this->onDebug) Rcpp::Rcout << "In: Extract data from dataList" << std::endl;

            this->dataList = dataList;            
            this->dataListCPP = as<List>(dataList);
          //  this->fittedParamNames = this->dataListCPP["par_names"];

            if (this->onDebug) Rcpp::Rcout << "In: Extract data from dataList1" << std::endl;

            // Useful for internal functions
            this->N_data = this->dataListCPP["N_data"];
            

            this->adaptiveGibbsSD = VectorXd::Constant(this->N, 0);
            this->counterAdaptiveGibbs = VectorXd::Zero(this->N);

            this->chainNumber = i;

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
            if (this->onDebug) Rcpp::Rcout << "In: Check 3" << std::endl;

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
            if (this->onDebug) Rcpp::Rcout << "In: Check 4" << std::endl;

            this->posteriorSamplesLength = (this->iterations-this->burninPosterior)/(this->thin);
            this->posteriorOut = MatrixXd::Zero(this->posteriorSamplesLength, this->numberFittedPar + 2);
         //   this->posteriorJump;
          
     //       this->posteriorTitreExp = Array::Zero(this->posteriorSamplesLength, this->lengthJumpVec, this->B);
     //  /       this->posteriorObsTitre = Array::Zero((this->posteriorSamplesLength, this->lengthJumpVec, this->B);

            this->currentLogPosterior = 0;
            this->currentSampleMean = VectorXd::Zero(this->numberFittedPar);
            this->currentSample = VectorXd::Zero(this->numberFittedPar);

            this->proposalSample = VectorXd::Zero(this->numberFittedPar);

            this->nonadaptiveScalar = 0;
            this->adaptiveScalar = 0;
            
            this->currentCovarianceMatrix = MatrixXd::Zero(this->numberFittedPar, this->numberFittedPar );
            this->nonadaptiveCovarianceMat = MatrixXd::Zero(this->numberFittedPar, this->numberFittedPar );
            this->adaptiveCovarianceMat = MatrixXd::Zero(this->numberFittedPar ,this->numberFittedPar );
            if (this->onDebug) Rcpp::Rcout << "In: Check 5" << std::endl;

            VectorXd initialSample;
            MatrixXd initialJump;

            double initialLogLikelihood;
            MatrixXd initialCovarianceMatrix;
            
            for(int parNum = 0; parNum < this->numberFittedPar ; parNum++){
                this->nonadaptiveCovarianceMat(parNum,parNum) = this->covarInitVal*(this->upperParBounds(parNum) - this->lowerParBounds(parNum));
                this->adaptiveCovarianceMat(parNum,parNum) = this->covarInitValAdapt*(this->upperParBounds(parNum) - this->lowerParBounds(parNum));
            }
            
            this->nonadaptiveScalar = log(0.1*0.1/(double)this->numberFittedPar);
            this->adaptiveScalar = log(2.382*2.382/(double)this->numberFittedPar);
            if (this->onDebug) Rcpp::Rcout << "In: Check 6" << std::endl;

            // Get initial conditions for the proposal sample
            initialSample = this->sampleInitPrior(this->dataList);

            // Get initial conditions for the exposure states
    
            this->currentJump = this->sampleInitJump(initialSample, this->dataList);; // functionhere
            // Get initial conditions for the infection states
            if (this->onDebug) Rcpp::Rcout << "In: Check 7" << std::endl;

            this->currentCovarianceMatrix = this->nonadaptiveScalar*this->nonadaptiveCovarianceMat;
            initialLogLikelihood = this->evalLogPosterior(initialSample, this->currentJump, this->currentCovarianceMatrix, this->dataList, true);

            if (this->onDebug) Rcpp::Rcout << "In: Check 8" << std::endl;

            while(isinf(initialLogLikelihood) || isnan(initialLogLikelihood)){
                initialSample = this->sampleInitPrior(this->dataList);
                this->currentJump = this->sampleInitJump(initialSample, this->dataList);; // functionhere
                initialLogLikelihood = this->evalLogPosterior(initialSample, this->currentJump, this->currentCovarianceMatrix, this->dataList, true);
            }
            if (this->onDebug) Rcpp::Rcout << "In: Check 9" << std::endl;
            this->currentSample = initialSample;
            this->proposalJump = this->currentJump;
            this->currentSampleMean = initialSample;
            this->currentLogPosterior = initialLogLikelihood;

            ///this->proposalEventsFull = this->currentEventsFull;
            //this->proposalTitreFull = this->currentTitreFull;

            
            double alphaMVN = -sqrt(2)*ErfInv(0.234-1);
            this->stepSizeRobbinsMonro = (1.0-1.0/(double)this->numberFittedPar)*(pow(2*3.141, 0.5)*exp(alphaMVN*alphaMVN*0.5))/(2*alphaMVN) + 1.0/(this->numberFittedPar*0.234*(1-0.234));

            if (this->onDebug) Rcpp::Rcout << "End: Initialise rjmcmc" << std::endl;
            if (this->onDebug) Rcpp::Rcout << "End: initialiseClass" << std::endl;
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
            if (this->onDebug) Rcpp::Rcout << "Pre: iterations" << std::endl;
            for (int i = 0; i < this->iterations; i++){

                this->workingIteration = i;
                if (this->onDebug) Rcpp::Rcout << "Pre: updateAllChains" << std::endl;
               // if (this->profile) 
                updateAllChains();
               // if (this->profile) 
            }

            List out = Rcpp::List::create(
                _["pars"] = this->posteriorOut,
                _["jump"] = this->posteriorJump
            );

            return out;
        }
        
        void updateAllChains()
        {
            if (this->onDebug) Rcpp::Rcout << "Pre: updateRJMC" << std::endl;
            updateRJMC();            
            if (this->onDebug) Rcpp::Rcout << "Pre: updateGibbsTiming" << std::endl;

            updateGibbsTiming();
            

           this->counterFuncEval++;

            if (this->onDebug) Rcpp::Rcout << "Pre: updateOutputPosterior" << std::endl;
            updateOutputPosterior();
            if (this->onDebug) Rcpp::Rcout << "Pre: updateProposal" << std::endl;
        
            consoleUpdatefunction();
        }
        
        void updateRJMC() {
            if (this->onDebug) Rcpp::Rcout << "Pre: getAcceptanceRate" << std::endl;
            getAcceptanceRate();
        //    Rcpp::Rcout << "alpha (kin): " << this->alpha << std::endl;

            if (this->onDebug) Rcpp::Rcout << "Pre: updateSampleAndLogPosterior" << std::endl;
            updateSampleAndLogPosterior();
            if(this->conPropIn & (this->currJumpType == 1)) {
                updateProposal();
            }
        }    

        void getAcceptanceRate()
        {        
            this->isSampleAccepted = false;
           // Evaluate the ab kinetics parameters and do rjmcmc
            if (this->onDebug) Rcpp::Rcout << "Pre: selectProposalDist" << std::endl;
            selectProposalDist();
            if (this->onDebug) Rcpp::Rcout << "Pre: JumpProposalDist" << std::endl;
            JumpProposalDist();
            if (this->onDebug) Rcpp::Rcout << "Pre: evalLogPosterior" << std::endl;
            
            this->proposedLogPosterior = this->evalLogPosterior(this->proposalSample, this->proposalJump, this->currentCovarianceMatrix, this->dataList);
            if (this->onDebug) Rcpp::Rcout << "Pre: evaluateMetropolisRatio" << std::endl;
            evaluateMetropolisRatio();
        }

        // Sample an individual who is exposed in current chain
        int sampleExposed() {
            if (this->onDebug) Rcpp::Rcout << "In: sampleExposed" << std::endl;

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
     
        void JumpProposalDist() {
            if (this->onDebug) Rcpp::Rcout << "Redfine prop = curr" << std::endl;
            this->proposalJump = this->currentJump;

            double q = uniformContinuousDist(0, 1);
          //  double r;
            VectorXd q_prob(3);
            if (this->onDebug) Rcpp::Rcout << "get q_prob" << std::endl;
           // if (this->onDebug) Rcpp::Rcout << "this->knownInfsN: " << this->knownInfsN << std::endl;


            if (this->knownExpInd) {
                q_prob << 0, 1, 1;
            } else {
                q_prob = this->sampleProposal(this->currentSample, this->currentJump, this->dataList);
             //   if (this->currentJump.cols() == 2) {
             //       q_prob << 0, 0.67, 1.0;
            //    } else if (this->currentJump.cols() == 100) {
            ///        q_prob << 0.33, 1.0, 0;
           //     } else {
           //         q_prob << 0.33, 0.67, 1.0;
          //      }
            }
            //Rcpp::Rcout << "q: " << q << std::endl;

            // Substract a value
            if (this->onDebug) Rcpp::Rcout << "Find jumping" << std::endl;
            if (q < q_prob(0)) {
              //  Rcpp::Rcout << "In deletion" << std::endl;
                int N = this->currentJump.cols();
             //   Rcpp::Rcout << "In N: " << N << std::endl;

                this->currJumpIdx = uniformDiscreteDist(1, N); // length is n_ - 114
            //    Rcpp::Rcout << "this->currJumpIdx: " << this->currJumpIdx << std::endl;

                this->proposalJump = sampleDeathProposal(this->currentSample, this->currentJump, this->currJumpIdx, this->dataList); // length now n_ - 144 - 1
                this->currJumpType = 0;
                //this->propInferredExpN = this->currInferredExpN - 1;
            // Stay same
            } else if (q < q_prob(1)) {
              //  Rcpp::Rcout << "In stay same" << std::endl;
                this->propInferredExpN = this->currInferredExpN;
                if(this->currInferredExpN != this->knownInfsN) {
                    this->currJumpIdx = uniformDiscreteDist(1, N); // length is n_ - 114
                }
                this->currJumpType = 1;
            // Add an infection same
            } else if (q < q_prob(2)) {
              //  Rcpp::Rcout << "In addition same" << std::endl;
                int N = this->currentJump.cols();
                this->currJumpIdx = uniformDiscreteDist(1, N); // length is 220 - n_
                this->proposalJump = sampleBirthProposal(this->currentSample, this->currentJump, this->currJumpIdx, this->dataList); // length now n_ - 144 - 1
                this->currJumpType = 2;
                //this->propInferredExpN = this->currInferredExpN + 1;
            }
            if (this->onDebug) Rcpp::Rcout << "END jump sample" << std::endl;
        }

        void updateGibbsTiming() {
            if (this->onDebug) Rcpp::Rcout << "In: updateGibbsTiming" << std::endl;

            this->currJumpType = 1;
         
                this->proposalJump = this->currentJump;

                selectProposalDist(false);

                for (int i = 0; i < this->noGibbsSteps; i++) { 
                    int i_idx = uniformDiscreteDist(1, this->currentJump.cols()); // stuck here
                    this->proposalJump = sampleJump(this->proposalSample, this->currentJump, i_idx, this->dataList);
                }
                this->proposedLogPosterior = this->evalLogPosterior(this->proposalSample, this->proposalJump, this->currentCovarianceMatrix, this->dataList);

                evaluateMetropolisRatio();
                updateJumpHistoric();
               // updateProposalGibbs();
           
        }
        
        void selectProposalDist(bool update = true){
          //  if (this->onDebug) Rcpp::Rcout << "Post: selectProposalDist" << std::endl;

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
            double rjadjustmentFactor2 = 0;
            if(std::isnan(this->proposedLogPosterior) || std::isinf(this->proposedLogPosterior)) {
                this->alpha = 0;
            } else {
                if (this->currJumpType == 0) {
                    rjadjustmentFactor = this->evaluateDeathProposal(this->currentSample, this->currentJump, this->currJumpIdx, this->dataList); 
                } else if (this->currJumpType == 1) {
                    rjadjustmentFactor = 0;
                } else if (this->currJumpType == 2) {
                    rjadjustmentFactor = this->evaluateBirthProposal(this->proposalSample, this->proposalJump, this->currJumpIdx, this->dataList); // dlnorm(proposalJump[currJumpIdx], 3.395, 0.5961)
                }
                this->alpha = min(1.0, exp((this->proposedLogPosterior - this->currentLogPosterior + rjadjustmentFactor)));
            }
        }

        void updateSampleAndLogPosterior()
        {
            if (uniformContinuousDist(0, 1) < this->alpha) {
                if (this->onDebug) Rcpp::Rcout << "In: ACCEPTED" << std::endl;
              //  this->finddifferInf();
                this->isSampleAccepted = true; this->counterAccepted++;

                this->currentSample = this->proposalSample;
                this->currentJump = this->proposalJump;

                this->currentLogPosterior = this->proposedLogPosterior;
            } else {
                this->proposalSample = this->currentSample;
                this->proposalJump = this->currentJump;

                this->proposedLogPosterior = this->currentLogPosterior;
            }
        }
        void updateJumpHistoric()
        {
            if (uniformContinuousDist(0, 1) < this->alpha) {
                if (this->onDebug) Rcpp::Rcout << "In: ACCEPTED JUMP" << std::endl;
                //Rcpp::Rcout << "ACCEPTED (GIBBS): " << this->alpha << std::endl;
           //     this->finddifferInf();
                this->currentSample = this->proposalSample;
                this->currentJump = this->proposalJump;

                this->currentLogPosterior = this->proposedLogPosterior;
            } else {
                this->proposalSample = this->currentSample;
                this->proposalJump = this->currentJump;

                this->proposedLogPosterior = this->currentLogPosterior;
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

                this->posteriorJump.push_back(this->currentJump);

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
        
        void errorCheckVectorValid(const MatrixXd& proposalSample)
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
                        Rcpp::Rcout << "The proposed matrix is not finite." << std::endl;
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
                Rcpp::Rcout << "Running MCMC-PT iteration number: " << this->workingIteration << " of " <<  this->iterations << ". Chain number " << this->chainNumber << ". Current logpost: " << this->currentLogPosterior << ". Length of jump: " <<  this->currentJump.cols() << ".          \r";
                if (this->onDebug) {
                    Rcpp::Rcout << "\n Current values: " << this->currentSample << std::endl;
                }
            }
        }
    };

        //// WRAPPERS FOR R FUNCTIONS

    void init_sampleInitPrior(rjmc_pp::RJMC_PP_D* model, Rcpp::Function sampleInitPrior) {
        auto func = [sampleInitPrior](RObject dataList) {
            PutRNGstate();
            auto rData = sampleInitPrior(dataList);
            GetRNGstate();
            return Rcpp::as<VectorXd>(rData);
        };
        model->sampleInitPrior = func;
    }

    void init_sampleInitJump(rjmc_pp::RJMC_PP_D* model, Rcpp::Function sampleInitJump) {
        auto func = [sampleInitJump](VectorXd params, RObject dataList) {
            PutRNGstate();
            auto rData = sampleInitJump(params, dataList);
            GetRNGstate();
            return Rcpp::as<MatrixXd>(rData);
        };
        model->sampleInitJump = func;
    }

    void init_evaluateLogPrior(rjmc_pp::RJMC_PP_D* model, Rcpp::Function evaluateLogPrior) {
        auto func = [evaluateLogPrior](VectorXd params, MatrixXd jump, RObject dataList) {
            PutRNGstate();
            auto rData = evaluateLogPrior(params, jump, dataList);
            GetRNGstate();
            return Rcpp::as<double>(rData);
        };
        model->evaluateLogPrior = func;
    }

    void init_evaluateLogLikelihood(rjmc_pp::RJMC_PP_D* model, Rcpp::Function evaluateLogLikelihood) {
        auto func = [evaluateLogLikelihood](VectorXd params, MatrixXd jump, RObject dataList) {
            PutRNGstate();
            auto rData = evaluateLogLikelihood(params, jump, dataList);
            GetRNGstate();
            return Rcpp::as<double>(rData);
        };
        model->evaluateLogLikelihood = func;
    }

    void init_sampleBirthProposal(rjmc_pp::RJMC_PP_D* model, Rcpp::Function sampleBirthProposal) {
        auto func = [sampleBirthProposal](VectorXd params, MatrixXd jump, int i_idx, RObject dataList) {
            PutRNGstate();
            auto rData = sampleBirthProposal(params, jump, i_idx, dataList);
            GetRNGstate();
            return Rcpp::as<MatrixXd>(rData);
        };
        model->sampleBirthProposal = func;
    }

    void init_sampleDeathProposal(rjmc_pp::RJMC_PP_D* model, Rcpp::Function sampleDeathProposal) {
        auto func = [sampleDeathProposal](VectorXd params, MatrixXd jump, int i_idx, RObject dataList) {
            PutRNGstate();
            auto rData = sampleDeathProposal(params, jump, i_idx, dataList);
            GetRNGstate();
            return Rcpp::as<MatrixXd>(rData);
        };
        model->sampleDeathProposal = func;
    }

    void init_evaluateBirthProposal(rjmc_pp::RJMC_PP_D* model, Rcpp::Function evaluateBirthProposal) {
        auto func = [evaluateBirthProposal](VectorXd params, MatrixXd jump, int i_idx, RObject dataList) {
            PutRNGstate();
            auto rData = evaluateBirthProposal(params, jump, i_idx, dataList);
            GetRNGstate();
            return Rcpp::as<double>(rData);
        };
        model->evaluateBirthProposal = func;
    }

    void init_evaluateDeathProposal(rjmc_pp::RJMC_PP_D* model, Rcpp::Function evaluateDeathProposal) {
        auto func = [evaluateDeathProposal](VectorXd params, MatrixXd jump, int i_idx, RObject dataList) {
            PutRNGstate();
            auto rData = evaluateDeathProposal(params, jump, i_idx, dataList);
            GetRNGstate();
            return Rcpp::as<double>(rData);
        };
        model->evaluateDeathProposal = func;
    }

    void init_sampleJump(rjmc_pp::RJMC_PP_D* model, Rcpp::Function sampleJump) {
        auto func = [sampleJump](VectorXd params, MatrixXd jump, int i_idx, RObject dataList) {
            PutRNGstate();
            auto rData = sampleJump(params, jump, i_idx, dataList);
            GetRNGstate();
            return Rcpp::as<MatrixXd>(rData);
        };
        model->sampleJump = func;
    }

    void init_sampleProposal(rjmc_pp::RJMC_PP_D* model, Rcpp::Function sampleProposal) {
        auto func = [sampleProposal](VectorXd params, MatrixXd jump, RObject dataList) {
            PutRNGstate();
            auto rData = sampleProposal(params, jump, dataList);
            GetRNGstate();
            return Rcpp::as<VectorXd>(rData);
        };
        model->sampleProposal = func;
    }

};
// namespace RJMC_PP_D
#endif
