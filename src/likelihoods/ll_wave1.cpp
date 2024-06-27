#include <Rcpp.h>
using namespace Rcpp;

struct DoubleWithString {
    double value;
    std::string name;

    // Define a comparison operator to sort by string values
    bool operator<(const DoubleWithString& other) const {
        return name < other.name;
    }
};

double cauchy_pdf(double x, double x0, double gamma) {
    return (1.0 / M_PI) * (gamma / ((x - x0) * (x - x0) + gamma * gamma));
}

double cauchy_cdf(double x, double x0, double gamma) {
    return 0.5 + atan((x - x0) / gamma) / M_PI;
}

// [[Rcpp::export]]
double evaluateLogLikelihood_cpp(NumericVector params, NumericVector jump, NumericVector jump_inf,
                              NumericMatrix covariance, List datalist) {


    int N = as<int>(datalist["N"]);
    IntegerVector knownInfsVec = datalist["knownInfsVec"];
    int knownInfsN = as<int>(datalist["knownInfsN"]);

    NumericVector knownInfsDate = datalist["knownInfsDate"];
    NumericVector endTitreTime = datalist["endTitreTime"];

    int N_data = as<int>(datalist["N_data"]);

    NumericVector initialTitreValue = datalist["initialTitreValue"];
    NumericVector initialTitreTime = datalist["initialTitreTime"];
    NumericVector titre_full = datalist["titre_full"];
    NumericVector times_full = datalist["times_full"];
    IntegerVector id_full = datalist["id_full"];

    double a = params[0];
    double b = params[1];
    double c = params[2];

    double sigma_obs = params[3];

    double ll = 0;

    for (int i = 0; i < N_data; ++i) {
        double titre_val = titre_full[i];
        double time = times_full[i];
        int i_idx = id_full[i] - 1;

        double titre_est = initialTitreValue[i_idx];
        if (jump_inf[i_idx] == 0) {
        } else if ((jump_inf[i_idx] == 1) && (time < jump[i_idx])) {
        } else if ((jump_inf[i_idx] == 1) && (time >= jump[i_idx])) {
            if ((time - jump[i_idx]) < 14) {
                titre_est += (log(exp(a) + exp(c)) * (time - jump[i_idx]) / 14);
            } else {
                titre_est += (log(exp(a) * exp(-b / 10 * ((time - jump[i_idx]) - 14)) + exp(c)));
            }
        }

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