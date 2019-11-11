#' @title Visualization Of t-Greene Pairwise Comparisons
#' @description Returns a graphical or numerical correlational matrix of
#'   p-values for the interpopulation degree of sexual dimorphism as measured by
#'   Greene t-test
#' @param x Data frame containing summary statistics of both sexes for two or
#'   more populations, Default: NULL
#' @param Pop Number of the column containing populations' names, Default: 1
#' @param lower_tail Logical; if \code{TRUE} probabilities are \code{P[X <= x]},
#'   otherwise, \code{P[X > x]}, Default: FALSE
#' @param tail Number of t test tails, Default: 'two'
#' @param padjust padjust Method of p value adjustment for multiple comparisons
#'   following \code{p.adjust.methods}, Default: 'none'
#' @param plot Logical;if \code{TRUE} graphical matrix of p-values, Default:
#'   \code{TRUE}
#' @param ... additional arguments that can be passed to
#'   \link[corrplot]{corrplot} function.
#' @return Graphical or numerical matrix of p-values from Greene t-test pairwise
#'   comparisons.
#' @details Data is entered in a wide format with each row representing a given
#'   population.\code{Pop}  (first column by default) contains population names,
#'   \code{.mu} and \code{.sdev} contain means and standard deviations with
#'   \code{M} and \code{F} donating males and females respectively. While
#'   \code{m}&\code{f} are the male and female sample sizes.When more than two
#'   populations are tested,\code{p.adjust.methods:
#'   c('holm','hochberg','hommel','bonferroni','BH','BY','fdr','none')} can be
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
#' pMatrix(x = df,plot=TRUE,method = 'ellipse', type = 'lower', col = c('#AEB6E5',
#' '#B1A0DB', '#B788CD', '#BC6EB9', '#BC569E', '#B6407D', '#A93154'), is.corr =
#' FALSE, tl.cex = 0.8, tl.col = 'black', sig.level = 0.05,insig =
#' 'label_sig', pch.cex = 2.5, tl.pos = 'ld', win.asp = 1, tl.srt =
#' 0.1,number.cex = 0.5, na.label = 'NA')
#'
#' @rdname pMatrix
#' @export
#' @importFrom stats pt
#' @importFrom stats p.adjust
#' @importFrom plyr adply
#' @importFrom corrplot corrplot
#'
#' @references \insertRef{timonov2014study}{TestDimorph}
#'
#'   \insertRef{curate2017sex}{TestDimorph}
#'
#'   \insertRef{kranioti2009sex}{TestDimorph}
#'
#'   \insertRef{gulhan2015new}{TestDimorph}
#'
pMatrix <- function(x = NULL, Pop = 1, lower_tail = FALSE, tail = "two", padjust = "none",plot=FALSE,...) {
    x$Pop <- x[, Pop]
    x$Pop <- factor(x$Pop)
    T <- function(m, f, m2, f2, M.mu, F.mu, M.mu2, F.mu2, M.sdev, F.sdev, M.sdev2, F.sdev2,
        Pop = 1, lower_tail = lower_tail, tail = tail, padjust = padjust) {
        G <- ((M.mu - F.mu) - (M.mu2 - F.mu2))/(sqrt(((((m - 1) * M.sdev^2) + ((f - 1) *
            F.sdev^2) + ((m2 - 1) * M.sdev2^2) + ((f2 - 1) * F.sdev2^2)))/(m + f + m2 +
            f2 - 4)) * sqrt((1/m) + (1/f) + (1/m2) + (1/f2)))

        df <- (m + f + m2 + f2 - 4)



        if (tail == "one") {
            p <- (stats::pt(abs(G), df, lower.tail = lower_tail))
        } else {
            p <- (2 * stats::pt(abs(G), df, lower.tail = lower_tail))
        }
        p <- stats::p.adjust(p = p, method = padjust, n = ((nlevels(x$Pop)^2 - nlevels(x$Pop))/2))
        return(p)
    }

    ss <- function(x = x, y = x) {
        s <- function(x) {
            T(m = x[1, "m"], f = x[1, "f"], M.mu = x[1, "M.mu"], F.mu = x[1, "F.mu"],
                M.sdev = x[1, "M.sdev"], F.sdev = x[1, "F.sdev"], m2 = y[, "m"], f2 = y[,
                  "f"], M.mu2 = y[, "M.mu"], F.mu2 = y[, "F.mu"], M.sdev2 = y[, "M.sdev"],
                F.sdev2 = y[, "F.sdev"], padjust = padjust, tail = tail, lower_tail = lower_tail)
        }

        q <- plyr::adply(.data = x, .margins = 1, .fun = s)

        q <- as.matrix(q[, -(1:ncol(x))])
        rownames(q) <- levels(x$Pop)
        colnames(q) <- levels(x$Pop)
        return(q)
    }
    if (plot==TRUE) {
       v <- ss(x, x)
       corrplot::corrplot(corr = v,p.mat =v,...)

    }else{
    ss(x, x)
}
}
