# run_SeroJump.R
library(Rcpp)

# Source the C++ code
sourceCpp(here::here("src", "headers", "old", "testing.cpp"))

# Run the exported function
runSeroJump()