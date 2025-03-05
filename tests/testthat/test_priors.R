test_that("addPrior Check that addPrior makes correct dataframe", {

        expect_equivalent(addPrior("alpha", 0, 1, "unif", 0, 1), 
            data.frame(par_name = "alpha", lb = 0, ub = 1, dist = "unif", dist_par1 = "0", dist_par2 = "1", part_type = "prior"))
    }
)

test_that("addPrior Check that addPrior fails distribtion is unkown", {
        expect_error(addPrior("alpha", 0, 1, "XX", 0, 1), 
            "Error: `dist` = XX does not correspond to a probability density function in the stats package.")
    }
)

test_that("addPrior Check that addPrior lower bound is below upper bound", {
        expect_error(addPrior("alpha", 2, 1, "unif", 0, 1), 
           "Error: `lower bound`, 2,  is greater than `upper bound`, 1")
    }
)

test_that("addPrior Check that addPrior can't sample from a poorly defined arguments", {
        expect_error(addPrior("alpha", 0, 1, "unif", 5, 1), 
           "Invalid arguments: for a uniform distribution, lower must be less than upper.")
    }
)

test_that("addPriorHier Check that boundaries have been given ", {

        expect_error(addPriorHier("alpha", 0, 1, "unif", 1, 5, "exp", 1, NA, 2), 
            "ERROR: logit_scale must be defined")
    }
)

test_that("addPriorHier Check that addPrior makes correct dataframe", {

        expect_equivalent(addPriorHier("alpha", 0, 1, "unif", 1, 5, "exp", 1, NA, 2, c(0, 1)),
    
            bind_rows(
                data.frame(par_name = "alpha", lb = 0, ub = 1, dist = NA, dist_par1 = NA, dist_par2 = NA, part_type = "logit_boundary"),
                data.frame(par_name = "alpha", lb = 0, ub = 1, dist = "unif", dist_par1 = "1", dist_par2 = "5", part_type = "prior"),
                data.frame(par_name = "z_alpha_1", lb = -10, ub = 10, dist = "norm", dist_par1 = "0", dist_par2 = "1", part_type = "prior"),
                data.frame(par_name = "z_alpha_2", lb = -10, ub = 10, dist = "norm", dist_par1 = "0", dist_par2 = "1", part_type = "prior"),
                data.frame(par_name = "sigma_alpha", lb = 0, ub = 5, dist = "exp", dist_par1 = "1", dist_par2 = NA, part_type = "prior")
            )
        )
    }
)

test_that("addPriorHier Check that addPrior fails distribtion is unkown", {

        expect_error(addPriorHier("alpha", 0, 1, "unif", 1, 5, "XX", 1, NA, 2, c(0, 1)), 
            "Error: `dist_sd` = XX does not correspond to a probability density function in the stats package.")
    }
)
