#' @title Sex-Specific One-way ANOVA From Summary statistics
#' @description Calculates sex specific one-way ANOVA followed by \code{Tukey
#'   HSD} from summary statistics.
#' @param x Data frame containing summary statistics of both sexes for two or
#'   more populations.
#' @param Pop Number of the column containing populations' names, Default:
#'   \code{1}
#' @param pairwise Logical; if \code{TRUE} runs multiple pairwise comparisons on
#'   different populations using \code{Tukey HSD}, Default: \code{TRUE}
#' @return Sex specific ANOVA tables and pairwise comparisons.
#' @details Data is entered in a wide format with each row representing a given
#'   population.\code{Pop}  (first column by default) contains population names,
#'   \code{.mu} and \code{.sdev} contain means and standard deviations with
#'   \code{M} and \code{F} donating males and females respectively. While
#'   \code{m}&\code{f} are the male and female sample sizes.By setting the
#'   option \code{pairwise} to \code{TRUE}, different pairwise combinations of
#'   populations can be compared with \code{Tukey HSD} post hoc test.
#' @examples
#'   # Comparisons of femur head diameter in four populations
#'   library(TestDimorph)
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
#' aovSS(x = df)
#'
#' @rdname aovSS
#' @export
#' @importFrom stats TukeyHSD
#' @importFrom stats rnorm
#' @importFrom stats aov
#'
#' @references \insertRef{timonov2014study}{TestDimorph}
#'
#'   \insertRef{curate2017sex}{TestDimorph}
#'
#'   \insertRef{kranioti2009sex}{TestDimorph}
#'
#'   \insertRef{gulhan2015new}{TestDimorph}
#'
aovSS <- function(x, Pop = 1, pairwise = TRUE) {
    x <- data.frame(x)
    x$Pop <- x[, Pop]
    x$Pop <- factor(x$Pop)
    y <- function(.mu, .sdev, n) {
        N <- length(.mu)
        Pop <- factor(rep(levels(x$Pop), n))
        df <- lapply(1:N, function(i) {
            scale(stats::rnorm(n[i])) * .sdev[i] + .mu[i]
        })
        x <- do.call(rbind, df)
        out <- data.frame(Pop, x)
        out
    }
    Male <- y(x$M.mu, x$M.sdev, x$m)
    Female <- y(x$F.mu, x$F.sdev, x$f)
    av_M <- stats::aov(x ~ Pop, data = Male)
    av_F <- stats::aov(x ~ Pop, data = Female)
    M <- summary(av_M)
    F <- summary(av_F)
    M1 <- stats::TukeyHSD(av_M)
    F1 <- stats::TukeyHSD(av_F)
    if (pairwise == TRUE) {
        out <- list(M, M1, F, F1)
        names(out) <- NULL
        names(out) <- c("Male model", "Male posthoc", "Female model", "Female posthoc")
        return(out)
    } else {
        out <- list(M, F)
        names(out) <- NULL
        names(out) <- c("Male model", "Female model")
        return(out)
    }
}
