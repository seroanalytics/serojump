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
    void updateAbKineticParams(const VectorXd& params) {
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
                
                //Rcpp::Rcout << "i_idx: " << i_idx << std::endl;
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
        if (parent->knownInfsVec(i_idx) == 0) {
            double randomTime = parent->exposureFunctionSample(i_idx + 1);
            df_order_exp.emplace_back("random", randomTime);
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
        if (parent->copFlag) {
            for (int i = 0; i < parent->B; ++i) {
                NumericVector currentParsCOP_i;
                int i_b = parent->mapOfCOP[parent->biomarkers[i]];
                StringVector parnams = parent->parsCOPN[i_b];
                for (int j = 0; j < parnams.size(); ++j) {
                    currentParsCOP_i.push_back(paramsN[as<string>(parnams[j])]);
                } 
                parent->proposalParsCOP[parent->biomarkers[i]] = currentParsCOP_i;
            }
        }

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

           // Rcpp::Rcout << "i_idx: " << i_idx << std::endl;
          // Rcpp::Rcout << "bio: " << bio << std::endl;
           // Rcpp::Rcout << "original orderedEvents[0].name: " << orderedEvents[0].name << std::endl;
           // Rcpp::Rcout << "original orderedEvents[0].value: " << orderedEvents[0].value << std::endl;  
           // Rcpp::Rcout << "original initTitre: " << initTitre << std::endl;  

            // Define state with anchor func and start titre
            df_order_titre_b.clear();
            df_order_titre_b.emplace_back(anchor_func, anchor_titre);

            // For all events
            for (int i = 1; i < orderedEvents.size(); i++) {

                // Update time since last event
                time_since = orderedEvents[i].value - anchor_time;
                titre_obs = abkineticsFunction(anchor_titre, bio, anchor_func, time_since, i_idx ); // Calcualte titre at new event
                df_order_titre_b.emplace_back(orderedEvents[i].name, titre_obs); // add to vector
            /*    if (i == 1) {
                    Rcpp::Rcout << "1" << std::endl;
                    Rcpp::Rcout << "1 orderedEvents[i].name: " << orderedEvents[i].name << std::endl;
                    Rcpp::Rcout << "1 orderedEvents[i].value: " << orderedEvents[i].value << std::endl;
                    Rcpp::Rcout << "1 time_since: " << time_since << std::endl;
                    Rcpp::Rcout << "1 titre_obs: " << titre_obs << std::endl;
                }
                if (i == 2) {
                    Rcpp::Rcout << "2" << std::endl;
                    Rcpp::Rcout << "2 orderedEvents[i].name: " << orderedEvents[i].name << std::endl;
                    Rcpp::Rcout << "2 orderedEvents[i].value: " << orderedEvents[i].value << std::endl;
                    Rcpp::Rcout << "2 time_since: " << time_since << std::endl;
                    Rcpp::Rcout << "2 titre_obs: " << titre_obs << std::endl;
                }
                if (i == 3) {
                    Rcpp::Rcout << "3" << std::endl;
                    Rcpp::Rcout << "3 orderedEvents[i].name: " << orderedEvents[i].name << std::endl;
                    Rcpp::Rcout << "3 orderedEvents[i].value: " << orderedEvents[i].value << std::endl;
                    Rcpp::Rcout << "3 time_since: " << time_since << std::endl;
                    Rcpp::Rcout << "3 titre_obs: " << titre_obs << std::endl;
                }
                if (i == 4) {
                    Rcpp::Rcout << "4" << std::endl;
                    Rcpp::Rcout << "4 orderedEvents[i].name: " << orderedEvents[i].name << std::endl;
                    Rcpp::Rcout << "4 orderedEvents[i].value: " << orderedEvents[i].value << std::endl;
                    Rcpp::Rcout << "4 time_since: " << time_since << std::endl;
                    Rcpp::Rcout << "4 titre_obs: " << titre_obs << std::endl;
                }*/
                // Only update the anchors if the event is not a bleed, as bleeds don't affect titre trajectory
                if ((orderedEvents[i].name != "bleed") && (orderedEvents[i].name != "random")) {
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

    double get_value_par(NumericVector pars_base, StringVector params, string par_name) {
        double value;
        int idx;

        for (int j = 0; j < params.size(); j++) {
            string params_s =  Rcpp::as<std::string>(params[j]);
            if (par_name == params_s) {
                value = pars_base[j];
            }
        }
        //Rcpp::Rcout << "value: " << value << std::endl;
        return value;
    }

    double getArgumentValuesHeir(NumericVector pars_base, StringVector params, string param_name, NumericVector dataHeir, int M, int i_idx) {
        int cov_i = dataHeir[i_idx];

       // Rcpp::Rcout << "cov_i: " << cov_i << std::endl;
        //Rcpp::Rcout << "i_idx: " << i_idx << std::endl;

        double c[3];
        string stringname;
        stringname = param_name;
        c[0] = get_value_par(pars_base, params, stringname);
        //Rcpp::Rcout << "stringname A : " << stringname << std::endl;

        stringname = "z_" + param_name + "_" + to_string(cov_i);
        //Rcpp::Rcout << "stringname B : " << stringname << std::endl;
        //Rcpp::Rcout << "pars_base: " << pars_base << std::endl;
        //Rcpp::Rcout << "params: " << params << std::endl;

        c[1] = get_value_par(pars_base, params, stringname);
        //Rcpp::Rcout << "stringname B : " << stringname << std::endl;

       // Rcpp::Rcout << "cov_i: " << cov_i << std::endl;
      //  Rcpp::Rcout << "c[1]: " << c[1] << std::endl;
        stringname = "sigma_" + param_name;

        c[2] = get_value_par(pars_base, params, stringname);
        //Rcpp::Rcout << "stringname C : " << stringname << std::endl;

        double argument = c[0] + c[1] * c[2];
        //Rcpp::Rcout << "argument: " << i_idx << std::endl;

        return argument;
    }

    NumericVector extractHeirPars(NumericVector pars_base, string abID_i, int i_idx) {
        NumericVector dataHeir = parent->dataAbHeir[abID_i];  // heirarchical data
        int M = parent->dataAbHeirN[abID_i];
        // name of all the parameters
        StringVector params = parent->parsAbKinN[abID_i]; // all the parameter names
        StringVector params_base = parent->parsAbKinNBase[abID_i]; // all the argument parameter names
        StringVector params_heir = parent->parsAbKinNHeir[abID_i]; // all the heirarchical parameter names
        NumericVector arguments;

        for (int i = 0; i < params_base.size(); i++) {
            std::string param_name_A = Rcpp::as<std::string>(params_base[i]);
            if (std::find(params_heir.begin(), params_heir.end(), param_name_A) != params_heir.end()) {
                double argumentvalue = getArgumentValuesHeir(pars_base, params, param_name_A, dataHeir, M, i_idx);
                arguments.push_back(argumentvalue);
            } else {
                //Rcpp::Rcout << "In not here: " << std::endl;
                double argumentvalue = get_value_par(pars_base, params, param_name_A);
                arguments.push_back(argumentvalue);
            }
        }
      
        return arguments;
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
    double abkineticsFunction(double titre_est, string biomarker, string exposureType_i, double timeSince, int i_idx) {
        if (parent->onDebug) Rcpp::Rcout << "In: calTitre1" << std::endl;
        int abkey = parent->mapOfAbkinetics[{biomarker, exposureType_i}];
        if (parent->onDebug) Rcpp::Rcout << "In: calTitre2" << std::endl;
        string abID_i = parent->abID[abkey] ;
        if (parent->onDebug) Rcpp::Rcout << "In: calTitre3" << std::endl;
     
        // VectorXd pars_base = parent->abinfo[abkey].params;
        // VectorXd pars_effects = parent->abinfoheir[abkey].params;
        // findHeirParameters(pars_base, pars_effects);
        // 
        if (!parent->isAbHeir[abkey]) { 
            NumericVector pars_base = parent->abinfo[abkey].params; // the base parameter VALUES
            titre_est = Rcpp::as<double>(parent->abinfo[abkey].func(titre_est, timeSince, pars_base));
        } else {

            NumericVector pars_base = parent->abinfo[abkey].params; // base parameters ALL VALUES
            NumericVector pars_extracted_heir = extractHeirPars(pars_base, abID_i, i_idx);
            //Rcpp::Rcout <<  "pars_extracted_heir: " << pars_extracted_heir << std::endl;
            titre_est = Rcpp::as<double>(parent->abinfo[abkey].func(titre_est, timeSince, pars_extracted_heir));
            //Rcpp::Rcout << "titre_est_end: " << titre_est << std::endl;

        }



        if (parent->onDebug) Rcpp::Rcout << "Out: calTitre4" << std::endl;
        return titre_est;
    }
  //        NumericVector findHeirParameters(NumericVecto pars_base, NumericVector, dataHeir) {

   //     }

    // A shared pointer to the SeroJumpBase object
    std::shared_ptr<SeroJumpBase> parent;
};

#endif