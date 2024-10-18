#ifndef RJMC_INF_HPP
#define RJMC_INF_HPP

// Load all libraries and headers
// Requires Rcpp, Eigen and boost
#include <RcppCommon.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Core>
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
#include <memory>
#include <iostream>

// Define Rcpp dependicies and plugins
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins("cpp14")]]

// Makes for cleaner code
using namespace Rcpp;
using namespace std;
using namespace Eigen;
using namespace boost::math;

// Define a random number generator
static thread_local std::mt19937 rng(std::random_device{}());

#define EIGEN_DONT_VECTORIZE
#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT

// Define the value of pi
#define PI 3.14159265358979323846

// Define a MIN and MAX function
#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b)) // define MAX function for use later
#endif

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b)) // define MAX function for use later
#endif


// load all headers
#include "./utils.hpp"
#include "./mvn.hpp"
#include "./rjmc_inf_run.hpp" 

#endif