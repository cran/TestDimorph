test_that("raw_gen", {
  # univariate truncated distribution
  set.seed(123)
  testthat::expect_true(round(suppressWarnings(raw_gen(baboon.parms_df))[1, 6][[1]], 2) == 45.39)
  # multivariate distribution
  set.seed(123)
  testthat::expect_true(round(suppressWarnings(raw_gen(baboon.parms_list))[1, 6][[1]], 2) == 58.55)
  set.seed(123)
  testthat::expect_true(round(suppressWarnings(raw_gen(baboon.parms_df, R.res = baboon.parms_R))[1, 6][[1]], 2) == 58.55)
  testthat::expect_error(suppressWarnings(raw_gen(baboon.parms_R)))
  testthat::expect_error(suppressWarnings(raw_gen(baboon.parms_df, Pop = 500)))
  testthat::expect_error(suppressWarnings(raw_gen(baboon.parms_df, Trait = 500)))
  testthat::expect_error(suppressWarnings(raw_gen(baboon.parms_df, R.res = baboon.parms_list)))
  testthat::expect_error(suppressWarnings(raw_gen(baboon.parms_list[1:5])))
})
