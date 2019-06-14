#' @title Tg
#' @description Calculates the significance of the differences in degree sexual
#'   dimorphism between two populations using a modified Greene t-test which
#'   uses summary statistics as input.
#' @param x Data frame containing summary statistics of both sexes for two or
#'   more populations, Default: \code{NULL}
#' @param Pop Number of the column containing populations' names, Default:
#'   \code{1}
#' @param m Number of male sample size in the first population, Default:
#'   \code{NULL}
#' @param m2 Number of male sample size in the second population, Default:
#'   \code{NULL}
#' @param f Number of female sample size in the first population, Default:
#'   \code{NULL}
#' @param f2 Number of female sample size in the second population, Default:
#'   \code{NULL}
#' @param M.mu Means for males in the first population, Default: \code{NULL}
#' @param M.mu2 Means for males in the second population, Default: \code{NULL}
#' @param F.mu Means for females in the first population, Default: \code{NULL}
#' @param F.mu2 Means for females in the second population, Default: \code{NULL}
#' @param M.sdev Standard deviation for males in the first population, Default:
#'   \code{NULL}
#' @param M.sdev2 Standard deviation for males in the second population,
#'   Default: \code{NULL}
#' @param F.sdev Standard deviation for females in the first population,
#'   Default: \code{NULL}
#' @param F.sdev2 Standard deviation for females in the second population,
#'   Default: \code{NULL}
#' @param lower_tail Logical; if \code{TRUE} probabilities are \code{P[X <= x]},
#'   otherwise, \code{P[X > x]}., Default: \code{FALSE}
#' @param tail Number of t test tails, Default: \code{'two'}
#' @param padjust Method of p value adjustment for multiple comparisons
#'   following \code{p.adjust.methods}, Default: \code{'none'}
#' @return Degree of freedom, t and p values for Tg test.
#' @details Summary statistics can be entered directly as arguments in case of
#'   comparing two populations. Alternatively, data can be entered in a wide
#'   format with each row representing a given population.\code{Pop}  (first
#'   column by default) contains population names, \code{.mu} and \code{.sdev}
#'   contain means and standard deviations with \code{M} and \code{F} donating
#'   males and females. While \code{m}&\code{f} are the male and female sample
#'   sizes respectively. When more than two populations are entered as
#'   input,\code{p.adjust.methods:
#'   c('holm','hochberg','hommel','bonferroni','BH','BY','fdr','none')} can be
#'   used for p value adjustment.
#' @examples
#'  # Summary data in a data frame
#' Pop <- c('Turkish', 'Bulgarian')
#' m <- c(150.00, 82.00)
#' f <- c(150.00, 58.00)
#' M.mu <- c(49.39, 48.33)
#' F.mu <- c(42.91, 42.89)
#' M.sdev <- c(3.01, 2.53)
#' F.sdev <- c(2.90, 2.84)
#' df <- cbind.data.frame(
#'   Pop,
#'   m,
#'   f,
#'   M.mu,
#'   F.mu,
#'   M.sdev,
#'   F.sdev,
#'   stringsAsFactors = TRUE
#' )
#' TestDimorph::Tg(x = df)
#'
#' @importFrom stats pt
#' @importFrom stats p.adjust
#' @importFrom utils combn
#' @export
#' @references \insertRef{greene1989comparison}{TestDimorph}
#'
#'   \insertRef{timonov2014study}{TestDimorph}
#'
#'   \insertRef{gulhan2015new}{TestDimorph}
#'
Tg <- function(x = NULL, Pop = 1, m = NULL, m2 = NULL, f = NULL, f2 = NULL, M.mu = NULL, 
    M.mu2 = NULL, F.mu = NULL, F.mu2 = NULL, M.sdev = NULL, M.sdev2 = NULL, F.sdev = NULL, 
    F.sdev2 = NULL, lower_tail = FALSE, tail = "two", padjust = "none") {
    T <- function(m, f, m2, f2, M.mu, F.mu, M.mu2, F.mu2, M.sdev, F.sdev, M.sdev2, F.sdev2) {
        Tg <- ((M.mu - F.mu) - (M.mu2 - F.mu2))/(sqrt(((((m - 1) * M.sdev^2) + ((f - 
            1) * F.sdev^2) + ((m2 - 1) * M.sdev2^2) + ((f2 - 1) * F.sdev2^2)))/(m + f + 
            m2 + f2 - 4)) * sqrt((1/m) + (1/f) + (1/m2) + (1/f2)))
        
        df <- (m + f + m2 + f2 - 4)
        
        
        if (tail == "one") {
            p <- (stats::pt(abs(Tg), df, lower.tail = lower_tail))
        } else {
            p <- (2 * stats::pt(abs(Tg), df, lower.tail = lower_tail))
        }
        if (!is.null(x)) {
            p <- stats::p.adjust(p = p, method = padjust, n = nlevels(x$Pop)^2 - nlevels(x$Pop))
        }
        p <- format.pval(pv = round(p, 4), eps = 0.001)
        return(paste0("t", "(df=", round(df, 1), ")=", round(Tg, 4), "  , P=", p))
    }
    if (is.null(x)) {
        return(T(m, f, m2, f2, M.mu, F.mu, M.mu2, F.mu2, M.sdev, F.sdev, M.sdev2, F.sdev2))
    } else {
        x <- data.frame(x)
        x$Pop <- x[, Pop]
        x$Pop <- factor(x$Pop)
        pairs <- utils::combn(x$Pop, 2, simplify = FALSE)
        
        names(pairs) <- sapply(pairs, paste, collapse = "-")
        
        Tg <- sapply(pairs, function(y) {
            T(m = x[y[1], "m"], f = x[y[1], "f"], m2 = x[y[2], "m"], f2 = x[y[2], "f"], 
                M.mu = x[y[1], "M.mu"], F.mu = x[y[1], "F.mu"], M.mu2 = x[y[2], "M.mu"], 
                F.mu2 = x[y[2], "F.mu"], M.sdev = x[y[1], "M.sdev"], F.sdev = x[y[1], 
                  "F.sdev"], M.sdev2 = x[y[2], "M.sdev"], F.sdev2 = x[y[2], "F.sdev"])
        })
        
        output <- as.data.frame(Tg)
        colnames(output) <- c()
        return(output)
    }
}
