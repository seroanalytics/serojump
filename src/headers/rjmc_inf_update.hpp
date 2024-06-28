#ifndef RJMC_UPDATE_HPP
#define RJMC_UPDATE_HPP

#include "./rjmc_inf_ll.hpp"

/** 
 * @brief Class to update the titre values for all individuals given the new exposure times and fixed parameters
 * 
 */
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
    void updateEventsFull(bool init = FALSE) {
        if (parent->onDebug) Rcpp::Rcout << "In: updateEventsFull" << std::endl;

        if(init) {
            std::vector<DoubleWithString> df_order_exp_i;
            std::vector<std::vector<DoubleWithString> > df_order_titre_i;
            for (int i_idx = 0; i_idx < parent->N; i_idx++) {

                df_order_exp_i = this->sortevents(i_idx, parent->currentJump);
              //  int l = df_order_exp_i.size();
                parent->proposalEventsFull.push_back(df_order_exp_i);

                //df_order_titre_i = this->recalculateTitre_i(df_order_exp_i, i_idx);
                //this->currentTitreFull.insert(this->currentTitreFull.begin() + i_idx, df_order_titre_i);
                parent->proposalTitreFull.push_back(df_order_titre_i);
            }
        } else {
            if (parent->onDebug) Rcpp::Rcout << "In: updateLists1" << std::endl;
            std::vector<double> diff_ind = finddifferInf();
            if (parent->onDebug) Rcpp::Rcout << "In: updateLists2" << std::endl;
            std::vector<DoubleWithString> df_order_exp_i;
            if (parent->onDebug) Rcpp::Rcout << "In: updateLists3" << std::endl;

            int i_idx;
            for (int i = 0; i < diff_ind.size(); i++) {
                i_idx = diff_ind[i];
                df_order_exp_i = this->sortevents(i_idx, parent->proposalJump);
                parent->proposalEventsFull[i_idx] = df_order_exp_i;
            }
                if (parent->onDebug) Rcpp::Rcout << "In: updateLists4" << std::endl;
        }
        if (parent->onDebug) Rcpp::Rcout << "Out: updateEventsFull" << std::endl;
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
        if (parent->onDebug) Rcpp::Rcout << "In: finddifferInf" << std::endl;
        std::vector<double> diff_ind; // numeric vector of ids of individuals having different exposure times
        for (int i = 0; i < parent->N; i++) {
            if (parent->currentJump(i) != parent->proposalJump(i)) {
                diff_ind.push_back(i);
            }
        }
        if (parent->onDebug) Rcpp::Rcout << "Out: finddifferInf" << std::endl;
        return diff_ind;
    }


    /** 
     * @brief Sort the events for individual i_idx
     * @param i_idx The individual id
     * @param jump The new exposure times
     * @return A vector of events sorted by time
     */
    std::vector<DoubleWithString> sortevents(int i_idx, const VectorXd& jump) {
        if (parent->onDebug) Rcpp::Rcout << "In: sortevents" << std::endl;
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
        if (parent->onDebug) Rcpp::Rcout << "Out: sortevents" << std::endl;
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
        if (parent->onDebug) Rcpp::Rcout << "In: calTitre1" << std::endl;
        int abkey = parent->mapOfAbkinetics[{biomarker, exposureType_i}];
        if (parent->onDebug) Rcpp::Rcout << "In: calTitre2" << std::endl;
        string abID_i = parent->abID[abkey] ;
        if (parent->onDebug) Rcpp::Rcout << "In: calTitre3" << std::endl;
      /*  Rcpp::Rcout << "abkey: " << abkey << std::endl;
        Rcpp::Rcout << "titre_est: " << titre_est << std::endl;
        Rcpp::Rcout << "timeSince: " << timeSince << std::endl;
        Rcpp::Rcout << "parent->abinfo[abkey].params: " << parent->abinfo[abkey].params << std::endl;*/

        titre_est = Rcpp::as<double>(parent->abinfo[abkey].func(titre_est, timeSince, parent->abinfo[abkey].params));
        if (parent->onDebug) Rcpp::Rcout << "Out: calTitr4" << std::endl;
        return titre_est;
    }

    // A shared pointer to the SeroJumpBase object
    std::shared_ptr<SeroJumpBase> parent;
};

#endif