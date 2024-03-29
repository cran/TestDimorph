test_that("aov_ss", {
  testthat::expect_true(round(
    aov_ss(
      baboon.parms_df[1:3, ],
      Pop = 2,
      digits = 3,
      letters = TRUE,
      pairwise = TRUE,
    )[[1]]$p.value[1],
    3
  ) == 0.325)
  testthat::expect_true(
    aov_ss(
      baboon.parms_df[1:3, ],
      Pop = 2,
      digits = 3,
      letters = TRUE,
      pairwise = TRUE,
      es_anova = "f"
    )$`Female model`[[8]][1] == 0.028
  )
  testthat::expect_true(
    aov_ss(
      baboon.parms_df[1:3, ],
      Pop = 2,
      digits = 3,
      letters = TRUE,
      pairwise = TRUE,
      es_anova = "eta"
    )$`Female model`[[8]][1] == 0.027
  )
  testthat::expect_error(
    aov_ss(
      baboon.parms_df[1:3, ],
      Pop = 2,
      digits = 3,
      letters = TRUE,
      pairwise = TRUE,
      es_anova = "qq"
    )$`Female model`[[8]][1] == 0.028
  )
  testthat::expect_error(
    aov_ss(
      x = matrix(NA),
      Pop = 500,
      digits = 3,
      es_anova = 900,
      letters = 900,
      pairwise = 900,
      CI = 900
    )
  )
  testthat::expect_true(aov_ss(
    Howells_summary[1:4, -1],
    pairwise = TRUE,
    letters = TRUE,
    es_anova = "eta"
  )[[6]][[2]][[3]] == "c")
})
