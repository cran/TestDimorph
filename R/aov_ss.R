#' @title Sex Specific One way ANOVA From Summary statistics
#' @description Calculates sex specific one way ANOVA from summary statistics.
#' @inheritParams t_greene
#' @param pairwise Logical; if TRUE runs multiple pairwise comparisons on
#' different populations using Tukey's post hoc test, Default: TRUE
#' @return Sex specific ANOVA tables and pairwise comparisons in tidy format.
#' @details Data is entered as a tibble/data frame of summary statistics where
#' the column containing population names is chosen by position (first by
#' default), other columns of summary data should have specific names (case
#' sensitive) similar to [baboon.parms_df]
#' @examples
#' \donttest{
#' # Comparisons of femur head diameter in four populations
#' library(TestDimorph)
#' df <- data.frame(
#'   Pop = c("Turkish", "Bulgarian", "Greek", "Portuguese "),
#'   m = c(150.00, 82.00, 36.00, 34.00),
#'   f = c(150.00, 58.00, 34.00, 24.00),
#'   M.mu = c(49.39, 48.33, 46.99, 45.20),
#'   F.mu = c(42.91, 42.89, 42.44, 40.90),
#'   M.sdev = c(3.01, 2.53, 2.47, 2.00),
#'   F.sdev = c(2.90, 2.84, 2.26, 2.90)
#' )
#' aov_ss(x = df)
#' }
#' @rdname aov_ss
#' @export
#' @importFrom stats rnorm aov TukeyHSD
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom multcompView multcompLetters
aov_ss <-
  function(x,
           Pop = 1,
           pairwise = TRUE,
           letters = FALSE,
           es = FALSE,
           digits = 4,
           sig.level = 0.05) {
    # Data preparation --------------------------------------------------------

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
    if (!is.logical(pairwise)) {
      stop("pairwise should be either TRUE or FALSE")
    }
    if (!is.logical(letters)) {
      stop("letters should be either TRUE or FALSE")
    }
    if (!is.logical(es)) {
      stop("es should be either TRUE or FALSE")
    }
    if (sig.level < 0 ||
      sig.level > 1 || !is.numeric(sig.level)) {
      stop("sig.level should be a number between 0 and 1")
    }
    x <- data.frame(x)
    x$Pop <- x[, Pop]
    x$Pop <- factor(x$Pop)
    levels(x$Pop) <- droplevels(x$Pop)
    if (length(unique(x$Pop)) != length(x$Pop[which(!is.na(x$Pop))])) {
      stop("Populations names'should be unique")
    }

    # ANOVA from summary data -------------------------------------------------

    aov_sum <- function(.mu, .sdev, n) {
      N <- length(.mu)
      Pop <- factor(rep(levels(x$Pop), n))
      df <- lapply(seq_len(N), function(i) {
        scale(stats::rnorm(n[i])) * .sdev[i] + .mu[i]
      })
      x <- do.call(rbind, df)
      out <- data.frame(Pop, x)
      return(out)
    }
    av_M <- aov_sum(x$M.mu, x$M.sdev, x$m)
    av_F <- aov_sum(x$F.mu, x$F.sdev, x$f)
    av_M <- stats::aov(x ~ Pop, data = av_M)
    av_F <- stats::aov(x ~ Pop, data = av_F)
    M_model <-
      anova_es(x = av_M, es = es, digits = digits)
    F_model <-
      anova_es(x = av_F, es = es, digits = digits)

    # Pairwise comparisons ----------------------------------------------------

    M_post <-
      tibble::as_tibble(tibble::rownames_to_column(data.frame(
        TukeyHSD(av_M, conf.level = 1 - sig.level)[[1]]
      )))
    colnames(M_post) <-
      c(
        "populations",
        "mean.diff",
        "conf.low",
        "conf.high",
        "p.value"
      )
    p <- M_post$p.value
    M_post$signif <-
      case_when(
        p > 0.05 ~ "ns",
        p < 0.05 & p > 0.01 ~ "*",
        p < 0.01 & p > 0.001 ~ "**",
        p < 0.001 ~ "***"
      )

    M_letters <- M_post$p.value
    names(M_letters) <- M_post$populations
    M_letters <-
      tibble::rownames_to_column(data.frame(
        "letters" = multcompView::multcompLetters(M_letters, threshold = sig.level)[[1]]
      ),
      var = "populations"
      )
    F_post <-
      tibble::as_tibble(tibble::rownames_to_column(data.frame(
        TukeyHSD(av_F, conf.level = 1 - sig.level)[[1]]
      )))
    colnames(F_post) <-
      c(
        "populations",
        "mean.diff",
        "conf.low",
        "conf.high",
        "p.value"
      )
    p <- F_post$p.value
    F_post$signif <-
      case_when(
        p > 0.05 ~ "ns",
        p < 0.05 & p > 0.01 ~ "*",
        p < 0.01 & p > 0.001 ~ "**",
        p < 0.001 ~ "***"
      )

    F_letters <- F_post$p.value
    names(F_letters) <- F_post$populations
    F_letters <-
      tibble::rownames_to_column(data.frame(
        "letters" = multcompView::multcompLetters(F_letters, threshold = sig.level)[[1]]
      ),
      var = "populations"
      )


    if (isTRUE(pairwise)) {
      if (isTRUE(letters)) {
        list(
          "Male model" = M_model,
          "Male posthoc" = M_post,
          "Male letters" = M_letters,
          "Female model" = F_model,
          "Female posthoc" = F_post,
          "Female letters" = F_letters
        )
      } else {
        list(
          "Male model" = M_model,
          "Male posthoc" = M_post,
          "Female model" = F_model,
          "Female posthoc" = F_post
        )
      }
    } else {
      list("Male model" = M_model, "Female model" = F_model)
    }
  }