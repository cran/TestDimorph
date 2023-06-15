#' @title Raw Data Generation By Normal Or Truncated Normal Distribution
#' @description Generates raw data from summary statistics using
#' uni/multivariate truncated normal distribution
#' @inheritParams multivariate
#' @param lower scalar of lower bounds, Default: -Inf
#' @param upper scalar of upper bounds, Default: Inf
#' @param verbose Logical; if TRUE displays a message with the method used for
#'   generation , Default: FALSE
#' @return a data frame of raw data
#' @details If data generation is desired using multivariate distribution data
#' is entered in the form of a list of summary statistics and pooled within
#' correlation matrix as in \link{baboon.parms_list}, or the summary
#'  statistics are entered separately in the form of a data frame as in
#' \link{baboon.parms_df} with a separate correlation matrix as in
#' \link{baboon.parms_R}. If data frame is entered without a correlation
#' matrix, data generation is carried out using univariate distribution.
#' @examples
#' # Data generation using univariate distributions
#' raw_gen(baboon.parms_df, lower = 0)
#'
#' # another univariate example
#' library(dplyr)
#' data <- Cremains_measurements[1, ] %>% mutate(Pop=c("A")) %>%
#' relocate(Pop,.after=1)
#' raw_gen(data)[, -2]
#'
#' # Data generation using multivariate distribution
#' raw_gen(baboon.parms_list, lower = 0)
#' @rdname raw_gen
#' @export
#' @importFrom truncnorm rtruncnorm
#' @importFrom tidyr drop_na
#' @importFrom stats rlnorm

raw_gen <- function(x,
                    Trait = 1,
                    Pop = 2,
                    R.res = NULL,
                    lower = -Inf,
                    upper = Inf,
                    verbose = FALSE) {
  warning("The user should be aware that this function assumes that for
            univariate generated data the within sex/population variables are
            normally distributed and truncated if `lower` and/or `upper` are
            changed from their defaults. For multivariate generated data the
            within sex/population correlation matrices are all assumed to be
            equal and again truncated if `lower` and/or `upper` are changed
            from their defaults.")
  if (!(is.list(x) || is.data.frame(x))) {
    stop("x should be a list or a dataframe")
  }

  # univariate data generation ----------------------------------------------

  if (is.data.frame(x)) {
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
    if (!(Trait %in% seq_along(x))) {
      stop("Trait should be number from 1 to ncol(x)")
    }
    if (!(Pop %in% seq_along(x))) {
      stop("Pop should be number from 1 to ncol(x)")
    }
    if (is.null(R.res)) {
      x <- x %>%
        drop_na() %>%
        as.data.frame() %>% rename("Trait"=all_of(Trait),"Pop"=all_of(Pop))
      x$Pop <- factor(x$Pop, levels = unique(x$Pop))
      x$Trait <- factor(x$Trait, levels = unique(x$Trait))

      # Data generation --------------------------------------------------


      if (isTRUE(verbose)) {
        message("Data generation was done using univariate truncated distribution")
      }
      gen_m <- function(x) {
        truncnorm::rtruncnorm(
          n = x$m[1],
          a = lower,
          b = upper,
          mean = x$M.mu[1],
          sd = x$M.sdev[1]
        )
      }
      gen_f <- function(x) {
        truncnorm::rtruncnorm(
          n = x$f[1],
          a = lower,
          b = upper,
          mean = x$F.mu[1],
          sd = x$F.sdev[1]
        )
      }

      m_function <- function(x) {
        df <- by(x, list(x$Trait), list)
        df <- lapply(df, gen_m)
        df <- lapply(df, as.data.frame)
        df <- do.call(cbind_fill2, df)
        colnames(df) <- levels(x$Trait)
        df
      }
      f_function <- function(x) {
        df <- by(x, list(x$Trait), list)
        df <- lapply(df, gen_f)
        df <- lapply(df, as.data.frame)
        df <- do.call(cbind_fill2, df)
        colnames(df) <- levels(x$Trait)
        df
      }
      pops <- split.data.frame(x, x$Pop)
      male <- lapply(pops, m_function)
      male <- do.call(rbind.data.frame, male)
      female <- lapply(pops, f_function)
      female <- do.call(rbind.data.frame, female)
      males <- strsplit(rownames(male), split = "\\.")
      females <- strsplit(rownames(female), split = "\\.")
      male$Pop <- as.factor(do.call(rbind.data.frame, males)[, 1])
      male$Sex <- as.factor(rep("M", nrow(male)))
      female$Pop <-
        as.factor(do.call(rbind.data.frame, females)[, 1])
      female$Sex <- as.factor(rep("F", nrow(female)))
      male <-
        male[, c(ncol(male), ncol(male) - 1, seq(nlevels(x$Trait)))]


      female <-
        female[, c(ncol(female), ncol(female) - 1, seq(nlevels(x$Trait)))]

      # Joining both datasets ---------------------------------------------------

      wide <- rbind.data.frame(male, female)
      rownames(wide) <- NULL
      return(wide)
    }

    # multivariate generation with data.frame and correlation matrix -----------

    if (!is.null(R.res)) {
      if (!is.matrix(R.res)) {
        stop("R.res should be a matrix")
      }
      x <- dataframe2list(
        x = x,
        R.res = R.res,
        Trait = Trait,
        Pop = Pop
      )
    }
  }

  # multivariate data generation with list input -------------------------------------

  if (!(is.data.frame(x))) {
    if (!all(c("M.mu", "F.mu", "M.sdev", "F.sdev", "m", "f", "R.res") %in% names(x))) {
      stop(
        "List should have the following named matricies:
            M.mu= Male mean
            F.mu=Female mean
            M.sdev=Male sd
            F.sdev=Female sd
            m= Male sample size
            f=Female sample size
            R.res=Pooled within correlation matrix
            N.B: names are case sensitive"
      )
    }
    if (isTRUE(verbose)) {
      message("Data generation was done using multivariate truncated distribution")
    }
    multi_raw(
      x = x,
      upper = upper,
      lower = lower
    )
  }
}
