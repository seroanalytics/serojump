#ifndef PTMC_D_HPP
#define PTMC_D_HPP

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

namespace rjmc{
    struct RJMC_D
    {
        RJMC_D() {}
        
        bool conPropIn = true;
        bool disPropIn = false;

        VectorXd lowerParBounds, upperParBounds;
        double nonadaptiveScalar, adaptiveScalar;
        MatrixXd nonadaptiveCovarianceMat, adaptiveCovarianceMat;
        MatrixXd currentSample, currentSampleMean;
        MatrixXd posteriorOut;
        MatrixXd posteriorJump;
        MatrixXd currentCovarianceMatrix;
        double currentLogPosterior;
        VectorXd proposalSample;
        
        VectorXd currentJump;
        VectorXd proposalJump;

        int iPosterior;
        int iterations, posteriorSamplesLength, thin, burninPosterior, burninAdaptiveCov, consoleUpdates, updatesAdaptiveCov, chainNumber;
        int lengthJumpVec;
        int numberFittedPar;
        int workingIteration;
        bool onDebug, onAdaptiveCov;
        bool isSampleAccepted, isProposalAdaptive;

        RObject dataList;
        int counterFuncEval, counterAccepted, counterPosterior ,counterAdaptive;
        int counterNonAdaptive;
        
        double proposedLogPosterior, alpha, covarMaxVal, covarInitVal, covarInitValAdapt;

        std::function<VectorXd(RObject)> samplePriorDistributions;
        std::function<VectorXd(RObject)> initialiseJump;
        std::function<double(VectorXd, VectorXd, RObject)> evaluateLogPrior;
        std::function<VectorXd(VectorXd, RObject)> jumpSampling;
        std::function<double(VectorXd, VectorXd, MatrixXd, RObject)> evaluateLogLikelihood;
        std::function<double(VectorXd, double, double, RObject)> evaluateAcceptanceRatioR;

        Mvn Mvn_sampler;

        double stepSizeRobbinsMonro;
        double evalLogPosterior(const VectorXd& param, const VectorXd& jump, const MatrixXd& covariance, const RObject& dataList)
        {
            double logPrior = this->evaluateLogPrior(param, jump, dataList);
            if (isinf(logPrior))
                return log(0);
          
            double logLikelihood = this->evaluateLogLikelihood(param, jump, covariance, dataList);
            return logPrior + logLikelihood;
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
        
        void initialiseClass(List settings, RObject dataList, int i)
        {
            this->dataList = dataList;
            this->chainNumber = i;

            this->numberFittedPar = settings["numberFittedPar"];
            this->iterations = settings["iterations"];
            this->thin = settings["thin"];
            this->burninPosterior = settings["burninPosterior"];
            this->burninAdaptiveCov = settings["burninAdaptiveCov"];
            this->consoleUpdates = settings["consoleUpdates"];
            this->onAdaptiveCov = settings["onAdaptiveCov"];
            this->updatesAdaptiveCov = settings["updatesAdaptiveCov"];
            this->onDebug = settings["onDebug"];
            this->lengthJumpVec = settings["lengthJumpVec"];
            
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

            this->posteriorSamplesLength = (this->iterations-this->burninPosterior)/(this->thin);
            this->posteriorOut = MatrixXd::Zero(this->posteriorSamplesLength, this->numberFittedPar + 2);
            this->posteriorJump = MatrixXd::Zero(this->posteriorSamplesLength, this->lengthJumpVec);

            this->currentLogPosterior = 0;
            this->currentSampleMean = VectorXd::Zero(this->numberFittedPar);
            this->currentSample = VectorXd::Zero(this->numberFittedPar);

            this->proposalSample = VectorXd::Zero(this->numberFittedPar);

            this->nonadaptiveScalar = 0;
            this->adaptiveScalar =0;
            
            this->currentCovarianceMatrix = MatrixXd::Zero(this->numberFittedPar, this->numberFittedPar );
            this->nonadaptiveCovarianceMat = MatrixXd::Zero(this->numberFittedPar, this->numberFittedPar );
            this->adaptiveCovarianceMat = MatrixXd::Zero(this->numberFittedPar ,this->numberFittedPar );

            VectorXd initialSample;
            VectorXd initialJump;
            double initialLogLikelihood;
            MatrixXd initialCovarianceMatrix;
            
            for(int parNum = 0; parNum < this->numberFittedPar ; parNum++){
                this->nonadaptiveCovarianceMat(parNum,parNum) = this->covarInitVal*(this->upperParBounds(parNum) - this->lowerParBounds(parNum));
                this->adaptiveCovarianceMat(parNum,parNum) = this->covarInitValAdapt*(this->upperParBounds(parNum) - this->lowerParBounds(parNum));
            }
            
            this->nonadaptiveScalar = log(0.1*0.1/(double)this->numberFittedPar);
            this->adaptiveScalar = log(2.382*2.382/(double)this->numberFittedPar);
            
            initialSample = this->samplePriorDistributions(this->dataList);
            initialJump = this->initialiseJump(this->dataList);
            this->currentCovarianceMatrix = this->nonadaptiveScalar*this->nonadaptiveCovarianceMat;

            initialLogLikelihood = this->evalLogPosterior(initialSample, initialJump, this->currentCovarianceMatrix, this->dataList);
            while(isinf(initialLogLikelihood) || isnan(initialLogLikelihood)){
                initialSample = this->samplePriorDistributions(this->dataList);
                initialLogLikelihood = this->evalLogPosterior(initialSample, initialJump, this->currentCovarianceMatrix, this->dataList);
            }

            this->currentSample = initialSample;
            this->currentSampleMean = initialSample;
            this->currentJump = initialJump;
            this->currentLogPosterior = initialLogLikelihood;
            
            double alphaMVN = -sqrt(2)*ErfInv(0.234-1);
            this->stepSizeRobbinsMonro = (1.0-1.0/(double)this->numberFittedPar)*(pow(2*3.141, 0.5)*exp(alphaMVN*alphaMVN*0.5))/(2*alphaMVN) + 1.0/(this->numberFittedPar*0.234*(1-0.234));
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
            for (int i = 0; i < this->iterations; i++){
                this->workingIteration = i;
                updateAllChains();
            }
            Rcpp::Rcout << "HERE!!!!!!!" << std::endl;
            List out = Rcpp::List::create(
                _["pars"] = this->posteriorOut,
                _["jump"] = this->posteriorJump
            );
            return out;
        }
        
        void updateAllChains()
        {
            if (onDebug) Rcpp::Rcout << "Pre: getAcceptanceRate" << std::endl;
            getAcceptanceRate();
            if (onDebug) Rcpp::Rcout << "Pre: updateSampleAndLogPosterior" << std::endl;
            updateSampleAndLogPosterior();
            if (onDebug) Rcpp::Rcout << "Pre: updateOutputPosterior" << std::endl;
            updateOutputPosterior();
            if (onDebug) Rcpp::Rcout << "Pre: updateProposal" << std::endl;
            if(this->conPropIn) {
                updateProposal();
            }
        
            consoleUpdatefunction();
        }

        void getAcceptanceRate()
        {
            double p1 = uniformContinuousDist(0, 1);
        
            this->isSampleAccepted = false;
           // if (onDebug) Rcpp::Rcout << "Pre: selectProposalDist" << std::endl;
            selectProposalDist();
            JumpProposalDist();
            this->proposedLogPosterior = this->evalLogPosterior(this->proposalSample, this->proposalJump, this->currentCovarianceMatrix, this->dataList);
            // Jump stuff
            /*double p1 = uniformContinuousDist(0, 1);
            if (p1 < (1.0 / (1.0 + (double)this->updateDiscreteFreq))) {
                this->conPropIn = true; this->disPropIn = false;
                selectProposalDist();
                this->proposedLogPosterior = this->evalLogPosterior(this->proposalSample, this->proposalJump, this->currentCovarianceMatrix, this->dataList);
            }
            else {
                this->conPropIn = false; this->disPropIn = true;
                JumpProposalDist();
                this->proposedLogPosterior = this->evalLogPosterior(this->proposalSample, this->proposalJump, this->currentCovarianceMatrix, this->dataList);
            }*/
            evaluateMetropolisRatio();
            this->counterFuncEval++;
        }

        void JumpProposalDist() {
            this->proposalJump = this->jumpSampling(this->currentJump, this->dataList);
        }
        
        void selectProposalDist(){
          //  if (onDebug) Rcpp::Rcout << "Post: selectProposalDist" << std::endl;

            if (this->workingIteration < this->burninAdaptiveCov || uniformContinuousDist(0, 1) < 0.05 || !this->onAdaptiveCov ){
                generateSampleFromNonAdaptiveProposalDist();
            }
            else{
                generateSampleFromAdaptiveProposalDist();
            }
        }
        
        void generateSampleFromNonAdaptiveProposalDist()
        {
            double s;
            s = exp(this->nonadaptiveScalar);
            this->counterNonAdaptive++; this->isProposalAdaptive = false;
            this->currentCovarianceMatrix = s*this->nonadaptiveCovarianceMat;
         ///   if (onDebug) Rcpp::Rcout << "PRES!!!!!!!!: updateCholesky" << std::endl;

            Mvn_sampler.updateCholesky(this->currentSample, this->currentCovarianceMatrix);
            this->proposalSample = Mvn_sampler.sampleTrunc(this->lowerParBounds, this->upperParBounds, 10, this->onDebug);
            
            errorCheckVectorValid(this->proposalSample);
        }
        
        void generateSampleFromAdaptiveProposalDist()
        {
            double s;
            s = exp(this->adaptiveScalar);
            this->counterAdaptive++; this->isProposalAdaptive = true;
            this->currentCovarianceMatrix = s*this->adaptiveCovarianceMat;
            Mvn_sampler.updateCholesky(this->currentSample, this->currentCovarianceMatrix);
            this->proposalSample = Mvn_sampler.sampleTrunc(this->lowerParBounds, this->upperParBounds, 10, this->onDebug);
            
            errorCheckVectorValid(this->proposalSample);
        }
        
        void evaluateMetropolisRatio()
        {
            if(std::isnan(this->proposedLogPosterior) || std::isinf(this->proposedLogPosterior)) {
                this->alpha = 0;
            }else {
                this->alpha = min(1.0, exp(this->evaluateAcceptanceRatioR(this->proposalJump, this->proposedLogPosterior, this->currentLogPosterior, this->dataList) ) );
                //this->alpha = min(1.0, exp((this->proposedLogPosterior - this->currentLogPosterior)));
            }
        }
        
        void updateSampleAndLogPosterior()
        {
            if (uniformContinuousDist(0, 1) < this->alpha) {
                this->isSampleAccepted = true; this->counterAccepted++;
                this->currentSample = proposalSample;
                this->currentJump = proposalJump;
                this->currentLogPosterior = proposedLogPosterior;
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
                }
                this->counterPosterior++;            
            }
        }
        
        void updateProposal()
        {
            int P = this->numberFittedPar;
            // Update adaptive proposal stuff
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
            
        void consoleUpdatefunction()
        {
            int i = this->workingIteration;
            if(i%this->consoleUpdates == 0) {
                Rcpp::Rcout << "Running MCMC-PT iteration number: " << this->workingIteration << " of " <<  this->iterations << ". Chain number " << this->chainNumber << ". Current logpost: " << this->currentLogPosterior << "\r";
                if (this->onDebug) {
                    Rcpp::Rcout << "\n Current values: " << this->currentSample << std::endl;
                }
            }
        }
    };

    void init_samplePriorDistributions(rjmc::RJMC_D* model, Rcpp::Function samplePriorDistributions) {
        auto func = [samplePriorDistributions](RObject dataList) {
            PutRNGstate();
            auto rData = samplePriorDistributions(dataList);
            GetRNGstate();
            return Rcpp::as<VectorXd>(rData);
        };
        model->samplePriorDistributions = func;
    }

    void init_evaluateLogPrior(rjmc::RJMC_D* model, Rcpp::Function evaluateLogPrior) {
        auto func = [evaluateLogPrior](VectorXd params, VectorXd jump, RObject dataList) {
            PutRNGstate();
            auto rData = evaluateLogPrior(params, jump, dataList);
            GetRNGstate();
            return Rcpp::as<double>(rData);
        };
        model->evaluateLogPrior = func;
    }

    void init_evaluateLogLikelihood(rjmc::RJMC_D* model, Rcpp::Function evaluateLogLikelihood) {
        auto func = [evaluateLogLikelihood](VectorXd params, VectorXd jump, MatrixXd covariance, RObject dataList) {
            PutRNGstate();
            auto rData = evaluateLogLikelihood(params, jump, covariance, dataList);
            GetRNGstate();
            return Rcpp::as<double>(rData);
        };
        model->evaluateLogLikelihood = func;
    }

    void init_initialiseJump(rjmc::RJMC_D* model, Rcpp::Function initialiseJump) {
        auto func = [initialiseJump](RObject dataList) {
            PutRNGstate();
            auto rData = initialiseJump(dataList);
            GetRNGstate();
            return Rcpp::as<VectorXd>(rData);
        };
        model->initialiseJump = func;
    }

    void init_jumpSampling(rjmc::RJMC_D* model, Rcpp::Function jumpSampling) {
        auto func = [jumpSampling](VectorXd discrete, RObject dataList) {
            PutRNGstate();
            auto rData = jumpSampling(discrete, dataList);
            GetRNGstate();
            return Rcpp::as<VectorXd>(rData);
        };
        model->jumpSampling = func;
    }

    void init_evaluateAcceptanceRatioR(rjmc::RJMC_D* model, Rcpp::Function evaluateAcceptanceRatioR) {
        auto func = [evaluateAcceptanceRatioR](VectorXd jumpingSample, double proposed, double current, RObject dataList) {
            PutRNGstate();
            auto rData = evaluateAcceptanceRatioR(jumpingSample, proposed, current, dataList);
            GetRNGstate();
            return Rcpp::as<double>(rData);
        };
        model->evaluateAcceptanceRatioR = func;
    }
};
// namespace rjmc
#endif
