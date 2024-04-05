#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Core>

using namespace Rcpp;
#include "../headers/mvn.hpp"
#include "../headers/rjmc_full.hpp"

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins("cpp14")]]

struct DoubleWithString {
    double value;
    std::string name;

    DoubleWithString(const std::string& nm, double val) : value(val), name(nm) {}

    // Define a comparison operator to sort by string values
    bool operator<(const DoubleWithString& other) const {
        return value < other.value;
    }
};

// IMPORTANT INTERNAL FUNCTIONS

NumericVector createNamedParam(const NumericVector& params) {
    NumericVector paramsName = NumericVector::create(
        _("a_vax") = params[0],
        _("b_vax") = params[1],
        _("c_vax") = params[2],        
        _("a_pd") = params[3],
        _("b_pd") = params[4],
        _("c_pd") = params[5],
        _("a_d") = params[6],
        _("b_d") = params[7],
        _("c_d") = params[8],
        _("tau") = params[9],
        _("wane") = params[10],
        _("sigma_obs_not") = params[11]
    );
    return paramsName;
}

double calTitre(double titre_est, string exposureType, double timeSince, const NumericVector& paramsN) {
    //double titre_est_boost = fmax(0.0, titre_est);
    //double titre_dep = fmax(0.0, 1 - paramsN["tau"] * titre_est_boost );
    double titre_dep = 1;
    if (exposureType == "wane") {
        titre_est -= paramsN["wane"] * (timeSince);
    } else if (exposureType == "pre-delta") {
        double a_pd = paramsN["a_pd"];
        double b_pd = paramsN["b_pd"];
        double c_pd = paramsN["c_pd"];
        if (timeSince < 14) {
            titre_est += (log(exp(a_pd) + exp(c_pd)) * (timeSince) / 14) * titre_dep;
        } else {
            titre_est += (log(exp(a_pd) * exp(-b_pd/10 * (timeSince - 14))) + exp(c_pd)) * titre_dep;
        }
    } else if (exposureType == "vax") {
        double a_vax = paramsN["a_vax"];
        double b_vax = paramsN["b_vax"];
        double c_vax = paramsN["c_vax"];
        if (timeSince < 14) {
            titre_est += (log(exp(a_vax) + exp(c_vax)) * (timeSince) / 14) * titre_dep;
        } else {
            titre_est += (log(exp(a_vax) * exp(-b_vax/10 * (timeSince - 14))) + exp(c_vax)) * titre_dep;
        }
    } else if (exposureType == "delta") {
        double a_d = paramsN["a_d"];
        double b_d = paramsN["b_d"];
        double c_d = paramsN["c_d"];
        if (timeSince < 14) {
            titre_est += (log(exp(a_d) + exp(c_d)) * (timeSince) / 14) * titre_dep;
        } else {
            titre_est += (log(exp(a_d) * exp(-b_d/10 * (timeSince - 14))) + exp(c_d)) * titre_dep;
        }
    }
    titre_est = fmax(0.0, titre_est);
    return titre_est;
}

std::vector<DoubleWithString> orderExposureEvents(
        const NumericVector& inf_pd_time, 
        const NumericVector& vax_time, 
        const NumericVector& jump_inf, 
        const NumericVector& jump, double initialtime, int i_idx, double time) {

        std::vector<DoubleWithString> df_order_exp;
                
        df_order_exp.emplace_back("wane", initialtime); // Initial event
        df_order_exp.emplace_back("bleed", time); // End event

        if (inf_pd_time[i_idx] > -1 || vax_time[i_idx] > -1 || jump_inf[i_idx] == 1) {
            if (inf_pd_time[i_idx] > -1 && inf_pd_time[i_idx] < time) {
                df_order_exp.emplace_back("pre-delta", inf_pd_time[i_idx]); // "pre-delta" event
            }
            if (vax_time[i_idx] > -1 && vax_time[i_idx] < time) {
                df_order_exp.emplace_back("vax", vax_time[i_idx]); // "vax" event
            }
            if (jump_inf[i_idx] == 1 && jump[i_idx] < time) {
                df_order_exp.emplace_back("delta", jump[i_idx]); // "delta" event
            }
            // Sort events by time
            std::sort(df_order_exp.begin(), df_order_exp.end());
        }


        return df_order_exp;
}


double cauchy_pdf(double x, double x0, double gamma) {
    return (1.0 / M_PI) * (gamma / ((x - x0) * (x - x0) + gamma * gamma));
}

double cauchy_cdf(double x, double x0, double gamma) {
    return 0.5 + atan((x - x0) / gamma) / M_PI;
}

double dnorm_cppl(double x, double mean, double sd) {
  return -0.5 * pow((x - mean) / sd, 2) - log(sd * sqrt(2 * M_PI));
}

// Cumulative Distribution Function (CDF) of Normal Distribution
// [[Rcpp::export]]
double pnorm_cppl(double x, double mean, double sd) {
  return log(0.5) + log(1 + erf((x - mean) / (sd * sqrt(2))));
}


// [[Rcpp::export]]
VectorXd calculateTitreExp_cpp(const NumericVector& params, const NumericVector& jump, const NumericVector& jump_inf,
                              NumericMatrix covariance, const List& datalist) {

    int N = as<int>(datalist["N"]);
    IntegerVector knownInfsVec = datalist["knownInfsVec"];
    int knownInfsN = as<int>(datalist["knownInfsN"]);

    NumericVector knownInfsDate = datalist["knownInfsDate"];
    NumericVector endTitreTime = datalist["endTitreTime"];

    NumericVector inf_pd_time = datalist["inf_pd_time"];
    NumericVector vax_time = datalist["vax_time"];
    int N_data = as<int>(datalist["N_data"]);

    NumericVector initialTitreValue = datalist["initialTitreValue"];
    NumericVector initialTitreTime = datalist["initialTitreTime"];
    NumericVector titre_full = datalist["titre_full"];
    NumericVector times_full = datalist["times_full"];
    IntegerVector id_full = datalist["id_full"];

    NumericVector paramsN = createNamedParam(params);

    double ll = 0;
    VectorXd titreExp(N);
    
    std::vector<DoubleWithString> df_order_exp; 

    for (int i_idx = 0; i_idx < N; ++i_idx) {
        if (jump[i_idx] == -1) {
            titreExp[i_idx] = -1;
        } else {

            double titre_est = initialTitreValue[i_idx];
            double initialtime = initialTitreTime[i_idx];

            df_order_exp = orderExposureEvents(inf_pd_time, vax_time, jump_inf, jump, initialtime, i_idx, jump[i_idx]);

            // MODEL-PREDICTED TITRE AT DELTA
            for (int j = 0; j < df_order_exp.size() - 1; ++j) {
                double time_until = df_order_exp[j + 1].value - df_order_exp[j].value;
                if (df_order_exp[j].name == "delta") {
                    break;
                }
                titre_est = calTitre(titre_est, df_order_exp[j].name, time_until, paramsN) ;
            }
            titreExp[i_idx] = titre_est;
        }
    }
    return titreExp;
}

// [[Rcpp::export]]
double evaluateLogLikelihood_cpp(const NumericVector& params, const NumericVector& jump, const NumericVector& jump_inf, const NumericVector& titre_exp,
                              NumericMatrix covariance, const List& datalist) {


    int N = as<int>(datalist["N"]);
    IntegerVector knownInfsVec = datalist["knownInfsVec"];
    int knownInfsN = as<int>(datalist["knownInfsN"]);

    NumericVector knownInfsDate = datalist["knownInfsDate"];
    NumericVector endTitreTime = datalist["endTitreTime"];

    NumericVector inf_pd_time = datalist["inf_pd_time"];
    NumericVector vax_time = datalist["vax_time"];
    int N_data = as<int>(datalist["N_data"]);

    NumericVector initialTitreValue = datalist["initialTitreValue"];
    NumericVector initialTitreTime = datalist["initialTitreTime"];
    NumericVector titre_full = datalist["titre_full"];
    NumericVector times_full = datalist["times_full"];
    IntegerVector id_full = datalist["id_full"];

    NumericVector paramsN = createNamedParam(params);
    std::vector<DoubleWithString> df_order_exp; 

    double ll = 0;

    // For each observation
    for (int i = 0; i < N_data; ++i) {
        double time = times_full[i];
        int i_idx = id_full[i] - 1;

        double initialtime = initialTitreTime[i_idx];    
        double initialtitre = initialTitreValue[i_idx];
        double titre_est = initialtitre;

        if ((time - initialtime) != 0) {
            // Determine order of events 
            df_order_exp = orderExposureEvents(inf_pd_time, vax_time, jump_inf, jump, initialtime, i_idx, time);

            // MODEL-PREDICTED TITRE
            for (int j = 0; j < df_order_exp.size() - 1; ++j) {
                double time_until = df_order_exp[j + 1].value - df_order_exp[j].value;

                titre_est = calTitre(titre_est, df_order_exp[j].name, time_until, paramsN) ;
            }
        }
    
        // OBSERVATIONAL MODEL
        double titre_val = titre_full[i];
        if (titre_val <= log10(40)) {
                // Compute the cumulative density function (CDF) at log10(20)
            ll += log(cauchy_cdf(log10(40), titre_est, paramsN["sigma_obs_not"])); // true for log probability
        } else {
            // Compute the probability density function (PDF) at titre_val
            ll += log(cauchy_pdf(titre_val, titre_est, paramsN["sigma_obs_not"]));  // true for log probability
        }
    }

    return ll;
}


// [[Rcpp::export]]
MatrixXd evaluateLogLikelihood_cpp_trajectories(NumericVector params, NumericVector jump, NumericVector jump_inf, NumericVector titre_exp,
                              NumericMatrix covariance, List datalist) {
    
    int N = as<int>(datalist["N"]);
    IntegerVector knownInfsVec = datalist["knownInfsVec"];
    int knownInfsN = as<int>(datalist["knownInfsN"]);

    NumericVector knownInfsDate = datalist["knownInfsDate"];
    NumericVector endTitreTime = datalist["endTitreTime"];

    NumericVector inf_pd_time = datalist["inf_pd_time"];
    NumericVector vax_time = datalist["vax_time"];
    int N_data = as<int>(datalist["N_data"]);

    NumericVector initialTitreValue = datalist["initialTitreValue"];
    NumericVector initialTitreTime = datalist["initialTitreTime"];
    NumericVector titre_full = datalist["titre_full"];
    NumericVector times_full = datalist["times_full"];
    IntegerVector id_full = datalist["id_full"];

    NumericVector paramsN = createNamedParam(params);

    int T = datalist["T"];
    double ll = 0;

    MatrixXd titreExp(N_data, T);
    std::vector<DoubleWithString> df_order_exp; 

    for (int i_idx = 0; i_idx < N; ++i_idx) {

        double initialtime = initialTitreTime[i_idx];
        double init_titre = initialTitreValue[i_idx];
        df_order_exp = orderExposureEvents(inf_pd_time, vax_time, jump_inf, jump, initialtime, i_idx, T);

        // Get trajectories of titre
        int time_initial = initialTitreTime[i_idx];
        double titre_est = initialTitreValue[i_idx];
        for (int j = 0; j < df_order_exp.size() - 1; ++j) {
            int time_until = df_order_exp[j + 1].value - df_order_exp[j].value;
            for (int t = 0; t < time_until; ++t) {
                titre_est = calTitre(titre_est, df_order_exp[j].name, t, paramsN) ;
                titreExp(i_idx, time_initial) = titre_est;
                time_initial++;
            }
        }    
    }
    return titreExp;
}


      //   Rcpp::Rcout << "titre_est: " <<  titre_est << "\t";
      //  Rcpp::Rcout << "ll: " <<  dnorm_cppl(titre_val, titre_est, sigma_obs) << "\n" << std::endl;

        // Check condition and compute the log-likelihood
        //if ((jump_inf[i_idx] == 1 && jump[i_idx] < time) || (inf_pd_time[i_idx] > -1 && inf_pd_time[i_idx] < time) || 
        //    (vax_time[i_idx] > -1 && vax_time[i_idx] < time)){
        //    if (titre_val < std::log10(20)) {
                // Compute the cumulative density function (CDF) at log10(20)
       //         ll += log(cauchy_cdf(std::log10(20), titre_est, sigma_obs_inf)); // true for log probability
       //     } else {
                // Compute the probability density function (PDF) at titre_val
       /*         ll += log(cauchy_pdf(titre_val, titre_est, sigma_obs_inf));  // true for log probability
            }
        } else {
            if (titre_val < std::log10(20)) {
                // Compute the cumulative density function (CDF) at log10(20)
                ll += log(cauchy_cdf(std::log10(20), titre_est, sigma_obs_not)); // true for log probability
            } else {
                // Compute the probability density function (PDF) at titre_val
                ll += log(cauchy_pdf(titre_val, titre_est, sigma_obs_not));  // true for log probability
            }
        }*/