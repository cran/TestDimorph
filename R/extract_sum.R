#' @title Summary Statistics Extraction
#' @description Extract summary data needed for other functions from raw data.
#' @param x Data frame of raw data.
#' @param Sex Number of the column containing sex 'M' for male and 'F' for
#' female, Default: 1
#' @param Pop Number of the column containing populations' names, Default: 2
#' @param firstX Number of column containing measured parameters (First of
#' multiple in case of multivariate analysis), Default: 3
#' @param test `tg` for Greene t test \link{t_greene}, `uni` for
#' \link{univariate}, `aov` for sex specific ANOVA \link{aov_ss},
#' `multi` for \link{multivariate}, and `van`  for \link{van_vark},
#' Default: 1
#' @param run Logical; if TRUE runs the corresponding test after data
#' extraction, Default:TRUE
#' @param ... Additional arguments that could be passed to the test of choice
#' @return Input for other functions.
#' @details Raw data is entered in a wide format data frame similar to
#' \link{Howells} data set. The first two columns contain sex `Sex`
#' (`M` for male and `F` for female) (Default: `1`) and populations' names `Pop`
#' (Default:`2`). Starting from `firstX` column (Default: `3`), measured
#' parameters are entered each in a separate column.
#' @examples
#' # for multivariate test
#' \dontrun{
#' extract_sum(Howells, test = "multi")
#' # for univariate test on a specific parameter
#' extract_sum(Howells, test = "uni", firstX = 4)
#' }
#' @rdname extract_sum
#' @export
#' @import dplyr
#' @importFrom stats sd
extract_sum <-
  function(x,
           Sex = 1,
           Pop = 2,
           firstX = 3,
           test = "tg",
           run = TRUE,
           ...) {
    if (!(is.data.frame(x))) {
      stop("x should be a dataframe")
    }
    x <- x %>%
      drop_na() %>%
      as.data.frame()
    if (!(Sex %in% seq_along(x))) {
      stop("Sex should be number from 1 to ncol(x)")
    }
    if (!(Pop %in% seq_along(x))) {
      stop("Pop should be number from 1 to ncol(x)")
    }
    if (!(firstX %in% seq_along(x))) {
      stop("firstX should be number from 1 to ncol(x)")
    }
    if (!is.logical(run)) {
      stop("run should be either TRUE or FALSE")
    }
    if (test %in% c("multi", "van")) {
      if (!all(sapply(x[firstX:ncol(x)], is.numeric))) {
        stop("All selected parameters should be numeric")
      }
    } else {
      if (!is.numeric(x[[firstX]])) {
        stop("The selected parameter should be numeric")
      }
    }
    x <- x %>% rename("Pop"=all_of(Pop),"Sex"=all_of(Sex))
    x$Pop <- factor(x$Pop)
    x$Sex <- factor(x$Sex)
    if (length(unique(x$Sex)) != 2) {
      stop("Sex column should be a factor with only 2 levels `M` and `F`")
    }
    test <- match.arg(test, choices = c("tg", "uni", "aov", "multi", "van"))
    if (test == "multi") {
      nam <- paste(colnames(x)[firstX:ncol(x)], collapse = ",")
      message("The parameters used are ", nam)
      x <- as.data.frame.list(x)
      sex <- as.numeric(x[, Sex]) - 1
      pop <- as.numeric(x[, Pop])
      pop.names <- names(table(x[, Pop]))
      N.pops <- length(pop.names)
      ina <- pop + N.pops * sex
      X <- x[, -(1:(firstX - 1))]
      Trait.names <- colnames(X)
      V <- pooled_cov(as.matrix(X), ina)
      D <- diag(1 / sqrt(diag(V)))
      R.res <- D %*% V %*% D
      colnames(R.res) <- Trait.names
      rownames(R.res) <- Trait.names
      M.mu <-
        x %>%
        filter(Sex == "M") %>%
        select(all_of(Pop), all_of(firstX):ncol(x)) %>%
        group_by(Pop) %>%
        summarise_all(
          .funs =
            mean
        ) %>%
        select(-1) %>%
        as.matrix()
      row.names(M.mu) <- pop.names
      F.mu <-
        x %>%
        filter(Sex == "F") %>%
        select(all_of(Pop), all_of(firstX):ncol(x)) %>%
        group_by(Pop) %>%
        summarise_all(
          .funs =
            mean
        ) %>%
        select(-1) %>%
        as.matrix()
      row.names(F.mu) <- pop.names

      m <- table(x[, Pop][x[, Sex] == "M"])
      f <- table(x[, Pop][x[, Sex] == "F"])

      F.sdev <-
        matrix(NA, nrow = N.pops, ncol = NCOL(x) - firstX + 1)
      for (i in seq_len(N.pops)) {
        F.sdev[i, ] <- apply(X[ina == i, ], 2, stats::sd)
      }

      row.names(F.sdev) <- pop.names
      colnames(F.sdev) <- Trait.names

      M.sdev <-
        matrix(NA, nrow = N.pops, ncol = NCOL(x) - firstX + 1)
      for (i in seq_len(N.pops)) {
        M.sdev[i, ] <- apply(X[ina == N.pops + i, ], 2, stats::sd)
      }

      row.names(M.sdev) <- pop.names
      colnames(M.sdev) <- Trait.names

      multi_list <-
        list(
          R.res = R.res,
          M.mu = M.mu,
          F.mu = F.mu,
          m = m,
          f = f,
          M.sdev = M.sdev,
          F.sdev = F.sdev
        )
      if (isTRUE(run)) {
        return(multivariate(multi_list, ...))
      } else {
        return(multi_list)
      }
    } else {
      M <-
        x %>%
        filter(Sex == "M") %>%
        select(all_of(Pop), all_of(firstX)) %>%
        group_by(Pop) %>%
        summarise_all(list(
          M.mu = mean,
          M.sdev = sd,
          m = length
        ))
      F <-
        x %>%
        filter(Sex == "F") %>%
        select(all_of(Pop), all_of(firstX)) %>%
        group_by(Pop) %>%
        summarise_all(list(
          F.mu = mean,
          F.sdev = sd,
          f = length
        ))
      df <- dplyr::full_join(M, F, by = "Pop")
      if (test == "uni") {
        message("The parameter used is ", colnames(x)[firstX])
        if (isTRUE(run)) {
          return(univariate(x = df, Pop = 1, ...))
        } else {
          return(df)
        }
      }
      if (test == "aov") {
        message("The parameter used is ", colnames(x)[firstX])
        if (isTRUE(run)) {
          return(aov_ss(x = df, Pop = 1, ...))
        } else {
          return(df)
        }
      }
      if (test == "tg") {
        message("The parameter used is ", colnames(x)[firstX])
        if (isTRUE(run)) {
          return(TestDimorph::t_greene(x = df, Pop = 1, ...))
        } else {
          return(df)
        }
      }
    }
    if (test == "van") {
      nam <- paste(colnames(x)[firstX:ncol(x)], collapse = ",")
      message("The parameters used are ", nam)
      if (isTRUE(run)) {
        dat <-
          Van_vark_raw(
            x = x,
            Pop = Pop,
            Sex = Sex,
            firstX = firstX,
            ...
          )
        x <- dat$x
        W <- dat$W
        van_vark(x = x, W = W, ...)
      } else {
        return(Van_vark_raw(
          x = x,
          Pop = Pop,
          Sex = Sex,
          firstX = firstX,
          ...
        ))
      }
    }
  }
