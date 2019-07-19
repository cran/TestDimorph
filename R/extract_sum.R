#' @title Summary Statistics Extraction
#' @description Extract summary data needed for other functions from raw data.
#' @param x Data frame containing raw data.
#' @param Sex Number of the column containing gender (\code{M} for male and
#'   \code{F} for female), Default: \code{1}.
#' @param Pop Number of the column containing populations' names, Default:
#'   \code{2}.
#' @param firstX Number of the first column containing measured parameters,
#'   Default: \code{3}.
#' @param test \code{1} for Greene t-test \code{\link{Tg}}, \code{2} for
#'   \code{\link{univariate}}, \code{3} for sex specific ANOVA
#'   \code{\link{aovSS}}, \code{4} for \code{\link{multivariate}} and \code{5}
#'   for \code{\link{pMatrix}}, Default: \code{1}.
#' @param run Logical;if \code{TRUE} runs the corresponding test after data
#'   extraction, Default: \code{TRUE}.
#' @param pairwise Logical;if TRUE runs multiple comparisons after multi or
#'   univariate analysis, Default: FALSE.
#' @param padjust Method of p value adjustment for multiple comparisons
#'   following \code{p.adjust.methods}, Default: \code{'none'}.
#' @param lower_tail Logical; if TRUE probabilities are \code{P[X <= x]},
#'   otherwise, \code{P[X > x]}, Default: \code{FALSE}.
#' @param pairwise_tail Number of t-test tails, Default: \code{'two'}.
#' @return Input for other functions using raw data.
#' @details Raw data is entered in a data frame similar in format to
#'   \code{\link{Howells}} data set. The first two columns contain gender
#'   \code{Sex} (\code{M} for male and \code{F} for female) (Default: \code{1})
#'   and populations' names \code{Pop} (Default: \code{2}). Starting from
#'   \code{firstX} column (Default: \code{3}), measured parameters are entered
#'   each in a separate column.
#' @examples
#' # for multivariate test
#' library(TestDimorph)
#' extract_sum(Howells)
#' # for univariate test on a specific parameter
#' library(TestDimorph)
#' extract_sum(Howells, test = 2)
#'
#' @rdname extract_sum
#' @importFrom dplyr group_by summarise_all select ungroup %>%
#' @importFrom  Rfast pooled.cov
#' @importFrom stats sd
#' @importFrom plyr ddply
#' @export
#'
#'
extract_sum <- function(x, Sex = 1, Pop = 2, firstX = 3, test = 1, run = TRUE, pairwise = FALSE,
                      padjust = "none", lower_tail = FALSE, pairwise_tail = "two") {
    x <- data.frame(x)
    x$Pop <- x[, Pop]
    x$Sex <- x[, Sex]
    x$Pop <- factor(x$Pop)
    x$Sex <- factor(x$Sex)
    if (test == 4) {
        x <- as.data.frame.list(x)
        sex <- as.numeric(x[, Sex]) - 1
        pop <- as.numeric(x[, Pop])
        pop.names <- names(table(x[, Pop]))
        N.pops <- length(pop.names)
        ina <- pop + N.pops * sex
        X <- x[, -(1:(firstX - 1))]
        Trait.names <- colnames(X)
        V <- Rfast::pooled.cov(as.matrix(X), ina)
        D <- diag(1/sqrt(diag(V)))
        R.res <- D %*% V %*% D
        M.mu <- x %>%  group_by(Pop) %>% filter(Sex=="M") %>% select(firstX:ncol(x))  %>% summarise_all(.funs=mean) %>% ungroup()   %>% data.frame() %>% select(-1) %>% as.matrix()
        row.names(M.mu) <- pop.names
        F.mu <-  x %>%  group_by(Pop) %>% filter(Sex=="F") %>% select(firstX:ncol(x))  %>% summarise_all(.funs=mean) %>% ungroup() %>% data.frame() %>%  select(-1) %>% as.matrix()
        row.names(F.mu) <- pop.names

        m <- table(x[, Pop][x[, Sex] == "M"])
        f <- table(x[, Pop][x[, Sex] == "F"])

        F.sdev <- matrix(NA, nrow = N.pops, ncol = NCOL(x) - firstX + 1)
        for (i in 1:N.pops) {
            F.sdev[i, ] <- apply(X[ina == i, ], 2, stats::sd)
        }

        row.names(F.sdev) <- pop.names
        colnames(F.sdev) <- Trait.names

        M.sdev <- matrix(NA, nrow = N.pops, ncol = NCOL(x) - firstX + 1)
        for (i in 1:N.pops) {
            M.sdev[i, ] <- apply(X[ina == N.pops + i, ], 2, stats::sd)
        }

        row.names(M.sdev) <- pop.names
        colnames(M.sdev) <- Trait.names

        v <- list(R.res = R.res, M.mu = M.mu, F.mu = F.mu, m = m, f = f, M.sdev = M.sdev,
                  F.sdev = F.sdev)
        if (run == TRUE) {
            return(multivariate(v, univariate = pairwise, lower_tail = lower_tail, padjust = padjust))
        } else {
            return(v)
        }
    } else {
        x <- as.data.frame.list(x)

        m <- table(x[, Pop][x[, Sex] == "M"])
        f <- table(x[, Pop][x[, Sex] == "F"])
        v <- function(x) {
            mean(x[, firstX], na.rm = TRUE)
        }
        h <- function(x) {
            stats::sd(x[, firstX], na.rm = TRUE)
        }
        t <- plyr::ddply(x, .variables = c("Pop", "Sex"), .fun = c(m = v, s = h))

        M.mu <- t$m[t$Sex == "M"]
        F.mu <- t$m[t$Sex == "F"]
        M.sdev <- t$s[t$Sex == "M"]
        F.sdev <- t$s[t$Sex == "F"]



        v <- cbind(M.mu, F.mu, M.sdev, F.sdev, m, f)
        v <- data.frame(v)
        v$Pop <- rownames(v)
        rownames(v) <- NULL
        if (test == 2) {
            if (run == TRUE) {
                return(univariate(x = v, Pop = ncol(v), lower_tail = lower_tail, padjust = padjust,
                                  pairwise = pairwise, pairwise_tail = pairwise_tail))
            } else {
                return(v)
            }
        }
        if (test == 3) {
            if (run == TRUE) {
                return(aovSS(x = v, Pop = ncol(v), pairwise = pairwise))
            } else {
                return(v)
            }
        }
        if (test == 1) {
            if (run == TRUE) {
                return(Tg(x = v, Pop = ncol(v), tail = pairwise_tail, padjust = padjust,
                          lower_tail = lower_tail))
            } else {
                return(v)
            }
        }
        if (test == 5) {
            if (run == TRUE) {
                return(pMatrix(x = v, Pop = ncol(v), lower_tail = lower_tail, tail = pairwise_tail,
                               padjust = padjust))
            } else {
                return(v)
            }
        }
    }
}
