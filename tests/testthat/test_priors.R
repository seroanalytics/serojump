test_that("Check that addPrior makes correct dataframe", {
    # these take a while, so don't run on CI

        expect_equivalent(addPrior("alpha", 0, 1, "unif", 0, 1), 
            data.frame(par_name = "alpha", lb = 0, ub = 1, dist = "unif", dist_par1 = "0", dist_par2 = "1", part_type = "prior"))
    }
)

test_that("Check that addPrior fails distribtion is unkown", {
        # these take a while, so don't run on CI
        expect_error(addPrior("alpha", 0, 1, "XX", 0, 1), 
            "Error: `dist` = XX does not correspond to a probability density function in the stats package.")
    }
)

test_that("Check that addPrior lower bound is below upper bound", {
        # these take a while, so don't run on CI
        expect_error(addPrior("alpha", 2, 1, "unif", 0, 1), 
           "Error: `lower bound`, 2,  is greater than `upper bound`, 1")
    }
)

test_that("Check that addPrior can't sample from a poorly defined arguments", {
        # these take a while, so don't run on CI
        expect_error(addPrior("alpha", 0, 1, "unif", 5, 1), 
           "Invalid arguments: for a uniform distribution, lower must be less than upper.")
    }
)
