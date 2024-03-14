#include <Rcpp.h>
using namespace Rcpp;

struct DoubleWithString {
    double value;
    std::string name;

    // Define a comparison operator to sort by string values
    bool operator<(const DoubleWithString& other) const {
        return value < other.value;
    }
};

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
double evaluateLogLikelihood_cpp(NumericVector params, NumericVector jump, NumericVector jump_inf,
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

    double a_vax = params[0];
    double b_vax = params[1];
    double c_vax = params[2];

    double a_pd = params[3];
    double b_pd = params[4];
    double c_pd = params[5];

    double a_d = params[6];
    double b_d = params[7];
    double c_d = params[8];
    
    double wane = params[9];

    double sigma_obs = params[10];

    double ll = 0;

//data_t$inf_pd_time %>% length
//data_t$vax_time %>% length
//jump_inf$vax_time %>% length
//data_t$titre_full
//data_t$times_full
//data_t$id_full

    for (int i = 0; i < N_data; ++i) {
        double titre_val = titre_full[i];
        double time = times_full[i];
        int i_idx = id_full[i] - 1;

        StringVector exposure;
        NumericVector time_exp;

        if (inf_pd_time[i_idx] > -1 && inf_pd_time[i_idx] < time) {
            exposure.push_back("pre-delta");
            time_exp.push_back(inf_pd_time[i_idx]);
        }
        if (vax_time[i_idx] > -1 && vax_time[i_idx] < time) {
            exposure.push_back("vax");
            time_exp.push_back(vax_time[i_idx]);
        }
        if (jump_inf[i_idx] == 1 && jump[i_idx] < time) {
            exposure.push_back("delta");
            time_exp.push_back(jump[i_idx]);
        }

        exposure.push_back("bleed");
        time_exp.push_back(time);

        // Define a map to associate strings with numbers
        std::vector<DoubleWithString> df_order_exp;

        // Populate the map with elements from the vectors
        for (size_t i = 0; i < exposure.size(); ++i) {
            DoubleWithString entry;
            entry.value = time_exp[i];
            entry.name = exposure[i];
            df_order_exp.push_back(entry);
        }

        // Sort the vector of pairs by the number
        std::sort(df_order_exp.begin(), df_order_exp.end());
    //    Rcpp::Rcout << "Individual: " << i_idx ;
        for (int j = 0; j < df_order_exp.size(); ++j) {
     //         Rcpp::Rcout << ". Exposure: " << df_order_exp[j].name << ". Time: " << df_order_exp[j].value << "\t";
        }
    //    Rcpp::Rcout << "titre_val: " <<  titre_val << "\t";

        double titre_est = initialTitreValue[i_idx];
        if (df_order_exp.size() == 1) {

        }
        else {

            for (int j = 0; j < df_order_exp.size() - 1; ++j) {

                double time_until = df_order_exp[j + 1].value - df_order_exp[j].value;
                if (df_order_exp[j].name == "pre-delta") {
                    if (time_until < 14) {
                        titre_est += log(exp(a_pd) + exp(c_pd)) * (time_until) / 14;
                    } else {
                        titre_est += log(exp(a_pd) * exp(-b_pd/10 * (time_until - 14)) + exp(c_pd));
                    }
                } else if (df_order_exp[j].name == "vax") {
                    if (time_until < 14) {
                        titre_est += log(exp(a_vax) + exp(c_vax)) * (time_until) / 14;
                    } else {
                        titre_est += log(exp(a_vax) * exp(-b_vax/10 * (time_until - 14)) + exp(c_vax));
                    }
                } else if (df_order_exp[j].name == "delta") {
                    if (time_until < 14) {
                        titre_est += log(exp(a_d) + exp(c_d)) * (time_until) / 14;
                    } else {
                        titre_est += log(exp(a_d) * exp(-b_d/10 * (time_until - 14)) + exp(c_d));
                    }
                }
            }
        }
      //   Rcpp::Rcout << "titre_est: " <<  titre_est << "\t";
      //  Rcpp::Rcout << "ll: " <<  dnorm_cppl(titre_val, titre_est, sigma_obs) << "\n" << std::endl;

        // Check condition and compute the log-likelihood
        if (titre_val < std::log10(20)) {
            // Compute the cumulative density function (CDF) at log10(20)
            ll += log(cauchy_cdf(std::log10(20), titre_est, sigma_obs)); // true for log probability
        } else {
            // Compute the probability density function (PDF) at titre_val
            ll += log(cauchy_pdf(titre_val, titre_est, sigma_obs));  // true for log probability
        }
    }

    return ll;
}