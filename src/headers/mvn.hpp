#ifndef MVN_HPP
#define MVN_HPP

/**
 * @brief A class to represent a multivariate normal distribution and sample from it
 */
class Mvn
{
public:
    
    Eigen::VectorXd mean;
    Eigen::MatrixXd sigma;
    Eigen::MatrixXd L;
    
    Mvn() {

    };

    /** 
     * @brief The updateCholesky function updates the mean vector and computes the Cholesky decomposition of a given covariance matrix.
     * @param mu The mean vector
     * @param sigma The covariance matrix
     * @return void
     */
    void updateCholesky(const Eigen::VectorXd& mu, const Eigen::MatrixXd& sigma) {
        mean = mu;
        Eigen::LLT<Eigen::MatrixXd> cholSolver(sigma);
        this->L = cholSolver.matrixL();
        this->sigma = sigma;
    }
    
    /**
     * @brief The sample function generates a sample from the multivariate normal distribution
     * @param nr_iterations The number of iterations to sample
     * @return A sample from the distribution
     */
    Eigen::VectorXd sample(int nr_iterations)
    {
        int n = mean.rows();
        Eigen::VectorXd sum = Eigen::VectorXd::Zero(n);
        // Generate x from the N(0, I) distribution
        std::uniform_real_distribution<double> unif(0.0, 1.0);
        
        for (unsigned int i = 0; i < nr_iterations; i++) {
            Eigen::VectorXd x(n);
            for (int j = 0; j < n; j++) {
                x(j) = unif(rng);
            }
            sum += x;
        }
        // Normalise: Shift and scale the sum to have mean 0 and variance 1
        sum = sum - (static_cast<double>(nr_iterations) / 2) * Eigen::VectorXd::Ones(n);
        Eigen::VectorXd normalized = sum / std::sqrt(static_cast<double>(nr_iterations) / 12.0);
    
        
        // Cholesky Decomposition and Eigen Decomposition
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd > eigen_solver(sigma);  // Computes the eigenvalues and eigenvectors of the covariance matrix sigma
        Eigen::MatrixXd eigenvectors = eigen_solver.eigenvectors().real();   // Extracts the real part of the eigenvectors.
        Eigen::MatrixXd eigenvalues = eigen_solver.eigenvalues().cwiseMax(0).cwiseSqrt().asDiagonal();  // Extracts the eigenvalues, takes their maximum with 0 (to ensure non-negativity), then takes the square root, and forms a diagonal matrix.
        
        //Constructs the transformation matrix Q that will be used to scale the random samples.        
        Eigen::MatrixXd Q = eigenvectors * eigenvalues;
        
        return Q * normalized + mean;
    }

    /** 
     * @brief The sample_chol function generates a sample from the multivariate normal distribution using the Cholesky decomposition
     * @param nr_iterations The number of iterations to sample
     * @return A sample from the distribution
     */
    Eigen::VectorXd sample_chol(int nr_iterations)
    {
        int n = mean.rows();
        Eigen::VectorXd x(n);
        std::uniform_real_distribution<double> unif(0.0, 1.0);

        // Generate a sample with values in [0, 1], then shift to [-1, 1]
        for (int i = 0; i < n; ++i) {
            x(i) = unif(rng);
        }
        x = (x * 2.0) - Eigen::VectorXd::Ones(n);
                
        // Find the transformation matrix        
        return this->L * x + this->mean;
    }
    
    /** 
     * @brief The checkTrunc function checks if the sample is within the truncation bounds
     * @param sampleCheck The sample to check
     * @param lowerBound The lower truncation bound
     *  @param upperBound The upper truncation bound
     * @param parOOS A boolean to indicate if the sample is out of bounds
     * @return void
     */
    void checkTrunc(Eigen::VectorXd sampleCheck, Eigen::VectorXd lowerBound, Eigen::VectorXd upperBound, bool& parOOS)
    {
        for (int i = 0; i < sampleCheck.size(); i++){
            if ((sampleCheck(i) < lowerBound(i)) || (sampleCheck(i) > upperBound(i))){
                parOOS = true;
                return;
            }
        }
        parOOS = false;
    }
    
    /** 
     * @brief The sampleTrunc function generates a sample from the multivariate normal distribution and truncates it if it is out of bounds
     * @param lowerBound The lower truncation bound
     * @param upperBound The upper truncation bound
     * @param nr_iterations The number of iterations to sample
     * @param onDebug A boolean to indicate if debugging is on
     * @return A sample from the distribution
     */
    Eigen::VectorXd sampleTrunc(Eigen::VectorXd lowerBound, Eigen::VectorXd upperBound, int nr_iterations, bool onDebug)
    {
        bool parOOS = FALSE;
        int j = 0;
        Eigen::VectorXd sampleCheck;
        sampleCheck = sample_chol(nr_iterations);
        checkTrunc(sampleCheck, lowerBound, upperBound, parOOS);
        
        while (parOOS){
            j++;
            sampleCheck = sample_chol(nr_iterations);
            if (onDebug) {
                Rcpp::Rcout << sampleCheck << std::endl;
            }
            checkTrunc(sampleCheck, lowerBound, upperBound, parOOS);
            if (onDebug && (j == 1e6)) 
                Rcpp::Rcout << "Struggling to sample points." << std::endl;
        }
        return sampleCheck;
    }
};


#endif
