#' @title Univariate Analysis Of Sexual Dimorphism
#' @description Calculation of the significance of the differences in
#'   interpopulation degree of sexual dimorphism using a modified one-way ANOVA
#'   which utilizes summary statistics as input.
#' @param x Data frame containing summary statistics of both sexes for two or
#'   more populations.
#' @param Pop Number of the column containing populations' names, Default:
#'   \code{1}
#' @param lower_tail Logical; if \code{TRUE} probabilities are \code{P[X <= x]},
#'   otherwise, \code{P[X > x]}., Default: \code{FALSE}
#' @param padjust Method of p value adjustment for multiple comparisons
#'   following \code{p.adjust.methods}, Default: \code{'none'}
#' @param pairwise Logical; if \code{TRUE} runs multiple pairwise comparisons on
#'   different populations using \code{Tg} test, Default: \code{FALSE}
#' @param pairwise_tail Number of t test tails, Default: \code{'two'}
#' @return ANOVA table
#' @details Data is entered in a wide format with each row representing a given
#'   population.\code{Pop}  (first column by default) contains population names,
#'   \code{.mu} and \code{.sdev} contain means and standard deviations with
#'   \code{M} and \code{F} donating males and females respectively. While
#'   \code{m}&\code{f} are the male and female sample sizes. When more than two
#'   populations are tested, \code{p.adjust.methods:
#'   c('holm','hochberg','hommel', 'bonferroni', 'BH','BY','fdr','none')} can be
#'   used for p value adjustment.
#' @examples
#'  # Comparisons of femur head diameter in four populations
#' library(TestDimorph)
#' m <- c(150.00, 82.00, 36.00, 34.00)
#' f <- c(150.00, 58.00, 34.00, 24.00)
#' M.mu <- c(49.39, 48.33, 46.99, 45.20)
#' F.mu <- c(42.91, 42.89, 42.44, 40.90)
#' M.sdev <- c(3.01, 2.53, 2.47, 2.00)
#' F.sdev <- c(2.90, 2.84, 2.26, 2.90)
#' df <- cbind.data.frame(
#'   Pop = c('Turkish', 'Bulgarian', 'Greek', 'Portuguese '),
#'   m,
#'   f,
#'   M.mu,
#'   F.mu,
#'   M.sdev,
#'   F.sdev,
#'   stringsAsFactors = TRUE
#' )
#' univariate(x = df, pairwise = TRUE, padjust = 'bonferroni')
#'
#' @rdname univariate
#' @importFrom  stats pf
#' @importFrom rowr cbind.fill
#' @export
#' @references \insertRef{konigsberg1991historical}{TestDimorph}
#'
#'   \insertRef{timonov2014study}{TestDimorph}
#'
#'   \insertRef{curate2017sex}{TestDimorph}
#'
#'   \insertRef{kranioti2009sex}{TestDimorph}
#'
#'   \insertRef{gulhan2015new}{TestDimorph}
#'
univariate <- function(x, Pop = 1, lower_tail = FALSE, padjust = "none", pairwise = FALSE,
    pairwise_tail = "two") {
    x <- data.frame(x)
    x$Pop <- x[, Pop]
    x$Pop <- factor(x$Pop)
    r <- nrow(x)
    N <- sum(x[, "m"], x[, "f"])
    Df1 <- (r - 1)
    Df2 <- (N - (2 * r))
    x[, "w"] <- (x[, "m"] * x[, "f"])/(x[, "m"] + x[, "f"])
    x[, "d"] <- x[, "M.mu"] - x[, "F.mu"]
    SSE <- sum((x[, "m"] - 1) * (x[, "M.sdev"]^2) + ((x[, "f"] - 1) * (x[, "F.sdev"]^2)))
    SSI <- sum(x[, "w"] * x[, "d"]^2) - (sum((x[, "w"] * x[, "d"]))^2/sum(x[, "w"]))
    within <- SSE/Df2

    between <- SSI/Df1
    f <- between/within
    p <- stats::pf(f, Df1, Df2, lower.tail = lower_tail)
    p <- format.pval(pv = round(p, 4), eps = 0.001)
    SS <- c(SSI, SSE)
    DF <- c(Df1, Df2)
    MS <- c(between, within)

    out <- rowr::cbind.fill(round(DF, 1), round(SS, 4), round(MS, 4), round(f, 4), p,
        fill = NA)
    rownames(out) <- c("Sex", "Residuals")
    colnames(out) <- c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
    if (pairwise == TRUE) {
        out <- list(out, Tg(x, Pop = Pop, tail = pairwise_tail, lower_tail = lower_tail,
            padjust = padjust))
        names(out) <- c("univariate", "pairwise")
        return(out)
    } else {
        return(out)
    }
}
