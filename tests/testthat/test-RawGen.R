test_that("RawGen", {
  testthat::expect_true(ncol(Cremains_measurements %>% mutate(Pop = rep("A", nrow(.))) %>% RawGen(
    Pop =
      ncol(.)
  )) == 23)
  # univariate log distribution
  set.seed(123)
  testthat::expect_true(round(RawGen(baboon.parms_df, dist = "log")[1, 6][[1]], 2) == 43.96)
  # univariate truncated distribution
  set.seed(123)
  testthat::expect_true(round(RawGen(baboon.parms_df, dist = "trunc")[1, 6][[1]], 2) == 45.39)
  # multivariate distribution
  set.seed(123)
  testthat::expect_true(round(RawGen(baboon.parms_list)[1, 6][[1]], 2) == 58.55)
  set.seed(123)
  testthat::expect_true(round(RawGen(baboon.parms_df, R.res = baboon.parms_R)[1, 6][[1]], 2) == 58.55)
  testthat::expect_error(RawGen(baboon.parms_R))
  testthat::expect_error(RawGen(baboon.parms_df, complete_cases = 95))
  testthat::expect_error(RawGen(baboon.parms_df, Pop = 500))
  testthat::expect_error(RawGen(baboon.parms_df, Trait = 500))
  testthat::expect_true(ncol(RawGen(baboon.parms_df, format = "long")) == 4)
  testthat::expect_error(RawGen(baboon.parms_df, R.res = baboon.parms_list))
  testthat::expect_error(RawGen(baboon.parms_list[1:5]))
})
