#' @title Univariate Analysis Of Sexual Dimorphism
#' @description Calculation and visualization of the differences in degree
#' sexual dimorphism between multiple populations using a modified one way
#' ANOVA and summary statistics as input
#' @param type_anova type of ANOVA test "I","II" or "III", Default:"II".
#' @param interact_anova Logical; if TRUE calculates interaction effect,
#' Default: TRUE.
#' @param es_anova Type of effect size either "f2" for f squared,"eta2" for eta
#' squared, "omega2" for omega squared or "none", Default:"none".
#' @inheritParams t_greene
#' @inheritParams t_test
#' @param lower.tail Logical; if TRUE probabilities are `P[X <= x]`,
#' otherwise, `P[X > x]`., Default: FALSE
#' @param pairwise Logical; if TRUE runs multiple pairwise comparisons on
#' different populations using \link{t_greene} Default: FALSE
#' @param ... Additional arguments that could be passed to the \link{t_greene}
#' function
#' @return  ANOVA table.
#' @details Data is entered as a data frame of summary statistics where
#' the column containing population names is chosen by position (first by
#' default), other columns of summary data should have specific names (case
#' sensitive) similar to \link{baboon.parms_df}
#' @examples
#' #'
#' # See Tables 6 and 8 and from Fidler and Thompson (2001).
#' # The “eta2” and “omega2” CIs match those in Table 8.
#' # See “FT” dataset for Fidler and Thompson (2001) reference
#'
#' # acquiring summary data
#' FT_sum <- extract_sum(FT, test = "uni", run = FALSE)
#' # univariate analysis on summary data
#' univariate(FT_sum, CI = 0.90, es_anova = "eta2", digits = 5)
#' univariate(FT_sum, CI = 0.90, es_anova = "omega2", digits = 5)
#'
#'
#' # Reproduces Table 2 from Shaw and Mitchell-Olds (1993) using their Table 1.
#' # See “SMO” dataset for Shaw and Mitchell-Olds (1993) reference
#' # Note that Table 2 residual df is incorrectly given as 6,
#' # but is correctly given as 7 in Hector et al. (2010)
#'
#' # acquiring summary data
#' univ_SMO <- extract_sum(SMO, test = "uni", run = FALSE)
#' # univariate analysis on summary data
#' print(univariate(univ_SMO, type_anova = "I")[[1]])
#' print(univariate(univ_SMO, type_anova = "II"))
#' univariate(univ_SMO, type_anova = "III")
#'
#' @rdname univariate
#' @export
#' @importFrom stats pf
#' @references
#'
#'   Hector, Andy, Stefanie Von Felten, and Bernhard Schmid. "Analysis of variance
#'   with unbalanced data: an update for ecology & evolution." Journal of animal
#'   ecology 79.2 (2010): 308-316.
#'
univariate <- function(x,
                       Pop = 1,
                       type_anova = "II",
                       interact_anova = TRUE,
                       es_anova = "none",
                       pairwise = FALSE,
                       padjust = "none",
                       ...,
                       lower.tail = FALSE,
                       CI = 0.95,
                       digits = 4) {
  padjust <- match.arg(padjust, choices = p.adjust.methods)
  es_anova <-
    match.arg(es_anova, choices = c("none", "eta2", "f2", "omega2"))

  if (!(is.data.frame(x))) {
    stop("x should be a dataframe")
  }
  if (!all(c("M.mu", "F.mu", "M.sdev", "F.sdev", "m", "f") %in% names(x))) {
    stop(
      "colnames should contain:
            M.mu= Male mean
            F.mu=Female mean
            M.sdev=Male sd
            F.sdev=Female sd
            m= Male sample size
            f=Female sample size
            N.B: colnames are case sensitive"
    )
  }
  if (!(Pop %in% seq_along(x))) {
    stop("Pop should be number from 1 to ncol(x)")
  }
  if (nrow(x) < 2) {
    stop("x should at least have 2 rows")
  }
  if (!is.logical(pairwise)) {
    stop("pairwise should be either TRUE or FALSE")
  }
  if (CI < 0 ||
    CI > 1 || !is.numeric(CI)) {
    stop("CI should be a number between 0 and 1")
  }
  if (!is.logical(interact_anova)) {
    stop("interact_anova should be either TRUE or FALSE")
  }
  if (isFALSE(interact_anova) && type_anova == "III") {
    stop("main effects ANOVA is only available for types (I) and (II)")
  }

  x <- x %>%
    drop_na() %>%
    as.data.frame()
  x$Pop <- x[, Pop]
  x$Pop <- factor(x$Pop)

  if (isFALSE(interact_anova)) {
    out <- switch(type_anova,
      I = anova_main_I(
        x,
        es_anova,
        digits, CI, lower.tail
      ),
      II = anova_main_II(
        x,
        es_anova,
        digits, CI, lower.tail
      )
    )
  } else {
    out <- switch(type_anova,
      I = anova_I(
        x,
        es_anova,
        digits, CI, lower.tail
      ),
      II = anova_II(
        x,
        es_anova,
        digits, CI, lower.tail
      ),
      III = anova_III(
        x,
        es_anova,
        digits, CI, lower.tail
      )
    )
  }

  if (isTRUE(pairwise)) {
    out <-
      list(
        univariate = out,
        pairwise = TestDimorph::t_greene(
          x,
          Pop = Pop,
          padjust = padjust,
          CI = CI,
          ...
        )
      )
    out
  } else {
    out
  }
}
