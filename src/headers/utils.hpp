#ifndef UTILS_HPP
#define UTILS_HPP

// Define utiliy functions to help with the RJMC implementation

/**
 * @brief Calculate the log of a factorial
 * @param n The number for which to calculate the factorial
 * @return The log of the factorial
 */
double logFactorial(int n) {
    double result = 0.0;
    for (int i = 1; i <= n; ++i) {
        result += log(i);
    }
    return result;
}


/**
 * @brief Calculate the binomialCoefficient
 * @param n The number of trials
 * @param k The number of successes
 * @return The binomial coefficient
 */
unsigned long long binomialCoefficient(int n, int k) {
    if (k == 0 || k == n) return 1;
    unsigned long long result = 1;
    for (int i = 1; i <= k; ++i) {
        result *= n - (k - i);
        result /= i;
    }
    return result;
}

/**
 * @brief Calculate the round down of a double
 * @param value The value to round down
 * @return The rounded down value
 
 */
int roundDown(double value) {
    return static_cast<int>(std::floor(value));
}

/**
 * @brief Calculate the log of the Binomial PMF
 * @param k The number of successes
 * @param n The number of trials
 * @param p The probability of success
 * @return The log of the probability mass function
 */
double binomialPMFlog(int k, int n, double p) {
    // Calculate the binomial coefficient
    double log_coef = std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1);
    
    // Calculate the probability mass function
    double log_pmf = log_coef + k * std::log(p) + (n - k) * std::log(1.0 - p);

    return log_pmf;
}

/**
 * @brief Calculate the log of the beta PDF
 * @param x The value at which to evaluate the PDF
 * @param alpha The alpha parameter of the beta distribution
 * @param beta The beta parameter of the beta distribution
 * @return The log of the PDF
 */
double betaPDFlog(double x, double alpha, double beta) {
    if (x < 0 || x > 1) {
        // x must be in [0, 1]
        return -INFINITY; // Return negative infinity for out-of-range values
    }
    
    double log_numerator = (alpha - 1) * log(x) + (beta - 1) * log(1 - x);
    double log_denominator = lgamma(alpha) + lgamma(beta) - lgamma(alpha + beta);
    
    return log_numerator - log_denominator; // Return the log of the PDF
}



/**
 * @brief Calculate the PDF of the Dirichlet distribution
 * @param x The value at which to evaluate the PDF (vector)
 * @param alpha The alpha parameter of the Dirichlet distribution (vector)
 * @return The PDF
 */
double dirichletPDF(const VectorXd& x, const VectorXd& alpha) {
    // Check if the sizes of x and alpha match
    if (x.size() != alpha.size()) {
        std::cerr << "Error: Size mismatch between x and alpha.\n";
        return 0.0; // Return 0 as an error code
    }

    double sum_alpha = 0.0;
    double sum_log_gamma_alpha = 0.0;
    double sum_log_x_alpha = 0.0;

    for (size_t i = 0; i < x.size(); ++i) {
        sum_alpha += alpha[i];
        sum_log_gamma_alpha += boost::math::lgamma(alpha[i]);
        sum_log_x_alpha += (alpha[i] - 1) * std::log(x[i]);
    }

    double log_pdf = sum_log_gamma_alpha - boost::math::lgamma(sum_alpha);
    log_pdf += sum_log_x_alpha;

    return log_pdf;
}



/**
 * @brief Calculate the beta-Binomial PMF
 * @param k The number of successes
 * @param n The number of trials
 * @param alpha The alpha parameter of the beta distribution
 * @param beta The beta parameter of the beta distribution
 * @return The probability mass function
 */
double betaBinomialPMF(int k, int n, double alpha, double beta) {
    return binomialCoefficient(n, k) * tgamma(alpha + beta) / (tgamma(alpha) * tgamma(beta)) * 
           tgamma(k + alpha) * tgamma(n - k + beta) / tgamma(n + alpha + beta);
}

/**
 * @brief Calculate the log beta-Binomial PMF
 * @param k The number of successes
 * @param n The number of trials
 * @param alpha The alpha parameter of the beta distribution
 * @param beta The beta parameter of the beta distribution
 * @return The probability mass function
 */
double betaBinomialPMFlog(int k, int n, double alpha, double beta) {
    double log_coef = boost::math::lgamma(k + alpha) + boost::math::lgamma(n - k + beta) - 
                      boost::math::lgamma(k + 1) - boost::math::lgamma(n - k + 1) - 
                      boost::math::lgamma(alpha) - boost::math::lgamma(beta) - 
                      boost::math::lgamma(n + 1);
    
    double log_term1 = k * std::log(alpha) + (n - k) * std::log(beta);
    double log_term2 = std::log(n + 1) - std::log(n + alpha + beta);

    return log_coef + log_term1 + log_term2;
}


/** 
* @brief Sample a value from a uniform continuous distribution
* @param minValue The minimum value of the distribution
* @param maxValue The maximum value of the distribution
* @return A sample from the distribution
*/
double uniformContinuousDist(double minValue, double maxValue)
{
    boost::random::uniform_real_distribution<> u(minValue,maxValue); return u(rng);
}

/** 
 * @brief Sample a value from a uniform discrete distribution
 * @param minValue The minimum value of the distribution
 * @param maxValue The maximum value of the distribution
 * @return A sample from the distribution
 */
double uniformDiscreteDist(int minValue, int maxValue)
{
    boost::random::uniform_int_distribution<> u(minValue, maxValue); return u(rng);
}

// [[Rcpp::export]]
NumericVector get_values_by_key(DataFrame df, std::string key) {
    CharacterVector first_col = df[0];  // Extract first column (strings)
    NumericVector col2 = df[1];  // Second column (numeric)
    NumericVector col3 = df[2];  // Third column (numeric)

    for (int i = 0; i < first_col.size(); i++) {
        if (first_col[i] == key) {
            return NumericVector::create(col2[i], col3[i]);  // Return values as vector
        }
    }

    return NumericVector::create(NA_REAL, NA_REAL);  // Return NA if not found
}

double logit_inverse(double x) {
    return exp(x) / (1.0 + exp(x));
}

/** 
 * @brief Sample a value from a normal distribution
 * @param mean The mean of the distribution
 * @param sigma The standard deviation of the distribution
 * @return A sample from the distribution
 */
double normalDistSample(double mean, double sigma)
{
    boost::random::normal_distribution<> n(mean, sigma); return n(rng);
}

 /**
  * @class abKineticsInfo
  * @brief Define a structure to hold the information about the antibody kinetics
  * 
  * This class is used to store the information about the antibody kinetics which has an usual data structure, require a name function and parameters needs to evaluate the function
*/       
struct abKineticsInfo {
    Function func;
    std::string name;
    NumericVector params; 

    abKineticsInfo(const std::string name_i, Function func_i, NumericVector params_i) : func(func_i), name(name_i), params(params_i) {}
};

/**
 * @class DoubleWithString
 * @brief Define a structure to hold a double value and a string
 * 
 * This class is used to store information about the interesting time points (double) for each individual with a string stating the type of event (bleed or exposure)
 */
struct DoubleWithString {
    double value;
    std::string name;

    DoubleWithString(const std::string& nm, double val) : value(val), name(nm) {}

    // Define a comparison operator to sort by string values
    bool operator<(const DoubleWithString& other) const {
        if (value != other.value) {
            return value < other.value; // Sort by value
        } else {
            // If values are equal, put "bleed" first
            return (name != "bleed" && other.name == "bleed");
        }            
    }
};

/** 
 * @brief Calculate the inverse error function
 * @param x The value at which to evaluate the inverse error function
 * 
 * A handy approximation for the error function and its inverse" by Sergei Winitzki.
 * 
 * @return The value of the inverse error function
 */
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

/** 
 * custom implementation of std::make_unique, in case not available in compiler, comment out if it's throwing errors
 
namespace std {
    template<typename T, typename... Args>
    std::unique_ptr<T> make_unique(Args&&... args) {
        return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
    }
}
*/
#endif // UTILS_HPP