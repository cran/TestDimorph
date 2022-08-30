test_that("extract_sum", {
  testthat::expect_true(TestDimorph::extract_sum(Howells, test = "multi")$p.value[3] == 0.0695)
  testthat::expect_true(TestDimorph::extract_sum(Howells,
    test = "tg",
    padjust = "none"
  )$p.value[1] == 0.2708)
  testthat::expect_true(TestDimorph::extract_sum(Howells,
    test = "uni",
    padjust = "none"
  )$p.value[3] ==
    0.7243)
  testthat::expect_true(TestDimorph::extract_sum(Howells, test = "aov")[[1]][[6]][1] == 0)
  testthat::expect_true(TestDimorph::extract_sum(Howells,
    test = "van",
    plot = F
  )[1, 4] == 0.5271)
})
testthat::expect_true(extract_sum(
  NHANES_1999,
  test = "tg",
  letters = TRUE,
  padjust = "bon"
)[[2]][[2]][[3]] == "b")
