#' @title multivariate
#' @description Multivariate extension of Greene t-test \code{Tg}.
#' @param x Data frame or list containing summary statistics of for multiple
#'   parameters measured in both sexes in two or more populations.
#' @param R.res Pooled within correlational matrix, Default: \code{NULL}.
#' @param Trait Number of the column containing names of measured traits,
#'   Default: \code{1}.
#' @param Pop Number of the column containing populations' names, Default:
#'   \code{2}.
#' @param univariate Logical; if \code{TRUE} conducts multiple univariate
#'   analyses on different parameters separately, Default: \code{FALSE}.
#' @param padjust Method of p value adjustment following
#'   \code{p.adjust.methods}, Default: \code{'none'}
#' @param lower_tail Logical; if \code{TRUE}, probabilities are \code{P[X <=
#'   x]}, otherwise, \code{P[X > x]} for multivariate analysis, Default:
#'   \code{FALSE}
#' @return \code{MANOVA} or \code{ANOVA} table(s).
#' @details Data can be entered either as a data frame where the first 2 columns
#'   contain traits to be measured \code{Traits} and the populations' names
#'   \code{Pop} and other columns containing different summary statistics as in
#'   \code{baboon.parms_df}. In that case the pooled within correlational matrix
#'   \code{R.res} should be entered as an argument. Another acceptable format is
#'   a list of matrices containing different summary statistics as well as the
#'   correlational matrix \code{R.res} as in \code{baboon.parms_list}. By
#'   setting the option \code{univariate} to \code{TRUE}, multiple \code{ANOVA}s
#'   can be run on each parameter independently with the required p value
#'   correction \code{padjust}.
#' @examples
#'
#' # x is a data frame with seprate correlational matrix
#' TestDimorph::multivariate(x = baboon.parms_df, R.res = R)
#' # x is a list with the correlational matrix included
#' TestDimorph::multivariate(baboon.parms_list, univariate = TRUE, padjust = 'bonferroni')
#'
#' @importFrom stats pf
#' @importFrom stats p.adjust
#' @importFrom rowr cbind.fill
#' @importFrom reshape2 melt
#' @importFrom purrr map
#' @export
#' @references \insertRef{konigsberg1991historical}{TestDimorph}
#'
multivariate <- function(x, R.res = NULL, Trait = 1, Pop = 2, univariate = FALSE, padjust = "none",
    lower_tail = FALSE) {
    if (is.data.frame(x)) {
        x <- data.frame(x)
        R <- R.res
        i <- list(levels(x[, Pop]), levels(x[, Trait]))
        M <- matrix(data = x$M.mu, nrow = nlevels(x[, Pop]), ncol = nlevels(x[, Trait]),
            dimnames = i)
        F <- matrix(data = x$F.mu, nrow = nlevels(x[, Pop]), ncol = nlevels(x[, Trait]),
            dimnames = i)
        nM <- as.vector(x$m)
        nF <- as.vector(x$f)
        nM <- nM[!is.na(nM)]
        nF <- nF[!is.na(nF)]
        M.sd <- matrix(data = x$M.sdev, nrow = nlevels(x[, Pop]), ncol = nlevels(x[,
            Trait]), dimnames = i)
        F.sd <- matrix(data = x$F.sdev, nrow = nlevels(x[, Pop]), ncol = nlevels(x[,
            Trait]), dimnames = i)
        x <- list(R.res = R, M.mu = M, F.mu = F, m = nM, f = nF, M.sdev = M.sd, F.sdev = F.sd)
    } else {
        R <- x$R.res
        M <- x$M.mu
        F <- x$F.mu
        nM <- x$m
        nF <- x$f
        M.sd <- x$M.sdev
        F.sd <- x$F.sdev
    }
    p <- NROW(R)
    r <- NROW(M)
    o.p <- rep(1, p)
    o.r <- rep(1, r)
    J <- matrix(1, nr <- r, nc <- r)
    D <- M - F
    N <- sum(nM) + sum(nF)
    w <- (nM * nF)/(nM + nF)

    weighted.D <- as.numeric(t(D) %*% w)
    SSCPsex <- weighted.D %o% weighted.D/sum(w)
    SSCPi <- t(D) %*% (w %*% t(o.p) * D) - SSCPsex

    T <- diag(sqrt(apply((nM - 1) %o% o.p * M.sd^2 + (nF - 1) %o% o.p * F.sd^2, 2, sum)))
    SSCPe <- T %*% R %*% T

    Xm <- nM %*% t(o.p) * M
    Xf <- nF %*% t(o.p) * F
    SSCPsamp <- t(Xm) %*% M + t(Xf) %*% F - t(Xm) %*% J %*% Xm/sum(nM) - t(Xf) %*% J %*%
        Xf/sum(nF) - SSCPi

    Lambda <- det(SSCPe)/det(SSCPi + SSCPe)
    Lambda[2] <- det(SSCPi + SSCPe)/det(SSCPsex + SSCPi + SSCPe)
    Lambda[3] <- det(SSCPe)/det(SSCPsamp + SSCPe)

    vh <- r - 1
    vh[2] <- 1
    vh[3] <- vh[1]

    ve <- N - 2 * r
    ve[2] <- N - r - 1
    ve[3] <- ve[1]


    m <- ve + vh - (p + vh + 1)/2
    s <- sqrt(((p * vh)^2 - 4)/(p^2 + vh^2 - 5))

    DF1 <- p * vh
    DF2 <- m * s - p * vh/2 + 1

    Rao.R <- (1 - Lambda^(1/s))/(Lambda^(1/s)) * (m * s - p * vh/2 + 1)/(p * vh)


    p <- stats::pf(Rao.R, DF1, DF2, lower.tail = lower_tail)

    p <- format.pval(pv = round(p, 4), eps = 0.001)


    out <- cbind.data.frame(vh, round(Lambda, 4), round(Rao.R, 4), round(DF1, 1), round(DF2,
        1), p)
    rownames(out) <- c("Sex:Population", "Sex", "Population")
    colnames(out) <- c("Df", "Wilks", "approx F", "num Df", "den DF", "Pr(>F)")
    if (univariate == TRUE) {
        M <- reshape2::melt(data = x$M.mu)
        F <- reshape2::melt(data = x$F.mu)
        m <- as.data.frame(as.numeric(x$m))
        f <- as.data.frame(as.numeric(x$f))
        Msd <- reshape2::melt(data = x$M.sdev)
        Fsd <- reshape2::melt(data = x$F.sdev)
        e <- rowr::cbind.fill(colnames(x$M.mu), rownames(x$M.mu), M$value, F$value, Msd$value,
            Fsd$value)
        colnames(e) <- c("G", "Pop", "M.mu", "F.mu", "M.sdev", "F.sdev")
        e <- rowr::cbind.fill(e, m, f, fill = NA)
        colnames(e) <- c("G", "Pop", "M.mu", "F.mu", "M.sdev", "F.sdev", "f", "m")
        e <- as.data.frame(e)
        q <- split.data.frame(x = e, f = e$G)
        d <- function(x, F = f, M = m) {
            r <- nrow(x)
            N <- sum(M, F)
            DF1 <- (r - 1)
            DF2 <- (N - (2 * r))
            x["w"] <- (M * F)/(M + F)
            x["d"] <- x["M.mu"] - x["F.mu"]
            SSE <- sum((M - 1) * (x["M.sdev"]^2) + ((F - 1) * (x["F.sdev"]^2)))
            SSI <- sum(x["w"] * x["d"]^2) - (sum((x["w"] * x["d"]))^2/sum(x["w"]))
            within <- SSE/DF2

            between <- SSI/DF1
            f <- between/within
            if (is.data.frame(x)) {
                p <- stats::p.adjust(p = stats::pf(f, DF1, DF2, lower.tail = lower_tail),
                  method = padjust, n = nlevels(x[, Trait]))
            } else {
                p <- stats::p.adjust(p = stats::pf(f, DF1, DF2, lower.tail = lower_tail),
                  method = padjust, n = ncol(x[["M.mu"]]))
            }
            p <- format.pval(pv = round(p, 4), eps = 0.001)

            SS <- c(SSI, SSE)
            DF <- c(DF1, DF2)
            MS <- c(between, within)

            out <- rowr::cbind.fill(round(DF, 1), round(SS, 4), round(MS, 4), round(f,
                4), p, fill = NA)
            rownames(out) <- c("Sex", "Residuals")
            colnames(out) <- c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
            return(out)
        }

        out <- list(out, purrr::map(q, .f = d))
        names(out) <- c("multivariate", "univariate")
        return(out)
    } else {
        return(out)
    }
}
