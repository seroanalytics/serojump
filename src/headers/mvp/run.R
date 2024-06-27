# run_SeroJump.R
library(Rcpp)

# Source the C++ code
sourceCpp(here::here("src", "headers", "mvp", "testing.cpp"))

# Run the exported function
runSeroJump()