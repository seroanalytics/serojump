#ifndef RJMC_LL_HPP
#define RJMC_LL_HPP

#include "./rjmc_inf_init.hpp"

/** 
 * @brief Class to update the titre values for all individuals given the new exposure times and fixed parameters
 * 
 */
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

# endif