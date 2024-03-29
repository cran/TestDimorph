#' @title Graphical and statistical representation of dimorphism differences.
#' @description Provides testing for differences in patterning of sexual
#'   dimorphism between populations, as well as for evolutionary trends that may
#'   characterize other species. The test is based on the computation of the
#'   first q canonical variates (q=2 by default) or multiple discriminant
#'   functions to develop various tests of sexual dimorphism in any two
#'   populations A and B.
#' @inheritParams multivariate
#' @param x A Data frame of means and sample sizes for different populations or
#' a list of the summary data frame with Pooled within-group variance-covariance
#' matrix.
#' @param W Pooled within-group variance-covariance matrix supplied if x is a
#' dataframe , Default:NULL
#' @param q Number of canonical variates to retain for chi square test, Default:
#'   2
#' @param Trait number of column containing names of traits Default: 1.
#' @param plot Logical; if TRUE returns a graphical representation of dimorphism
#'   differences, Default: TRUE
#' @return The output includes a two-dimensional plot that illustrate the
#'   existing differences between tested populations and a statistical test of
#'   significance for the difference in dimorphism using chi square
#'   distribution.
#' @details Input is a data frame of means and sample sizes similar to
#'   \link{Howells_summary} with the same naming conventions used throughout the
#'   functions but with the standard deviation columns removed.
#' @note For plot labels to be fully visualized, maximizing image size is advised.
#' @examples
#' # selecting means and sample sizes
#' van_vark_data <- Howells_summary[!endsWith(
#'   x = names(Howells_summary),
#'   suffix = "dev"
#' )]
#' # running the function
#' van_vark(van_vark_data, Howells_V)
#' @rdname van_vark
#' @references Van Vark, G. N., et al. "Some multivariate tests for differences
#' in sexual dimorphism between human populations." Annals of human biology
#' 16.4 (1989): 301-310.
#' @import ggplot2
#' @import dplyr
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom stats pchisq
#' @importFrom utils combn
#' @export
van_vark <- function(x,
                     W = NULL,
                     q = 2,
                     Trait = 1,
                     Pop = 2,
                     plot = TRUE,
                     lower.tail = FALSE,
                     digits = 4) {
  if (!(is.list(x) || is.data.frame(x))) {
    stop("x should be a list or a dataframe")
  }
  if (is.data.frame(x) && is.null(W)) {
    stop("If x is a data frame W should be supplied")
  }
  if (is.data.frame(x) && !is.matrix(W)) {
    stop("W should be a matrix")
  }
  if (is.list(x) && !is.data.frame(x)) {
    if (!all(c("W", "x") %in% names(x))) {
      stop(
        "If x is a list it should contain:
            W= Pooled within-group variance-covariance matrix.
            x=Data frame of means and sample sizes for different populations.
            N.B: names are case sensitive"
      )
    }
    if (!is.data.frame(x$x)) {
      stop("x should be a dataframe")
    }
    if (!is.matrix(x$W)) {
      stop("W should be a matrix")
    }
    W <- x$W
    x <- x$x
  }
  if (!(Trait %in% seq_along(x))) {
    stop("Trait should be number from 1 to ncol(x)")
  }
  if (!(Pop %in% seq_along(x))) {
    stop("Pop should be number from 1 to ncol(x)")
  }
  if (!is.logical(plot)) {
    stop("plot should be either TRUE or FALSE")
  }
  if (q == 1 && isTRUE(plot)) {
    warning("plot can't be produced if q = 1")
  }
  Sex <- NULL
  x <- x %>%
    drop_na() %>%
    as.data.frame() %>% rename("Trait" = all_of(Trait), "Pop" = all_of(Pop))
  x$Trait <- gsub(x = x$Trait,
                  pattern = "[^[:alnum:]]",
                  replacement = "_")
  x$Pop <- gsub(x = x$Pop,
                pattern = "[^[:alnum:]]",
                replacement = "_")

  x$Pop <- factor(x$Pop, levels = unique(x$Pop))
  x$Trait <- factor(x$Trait, levels = unique(x$Trait))
  x <- x %>%
    rename("Trait" = Trait) %>%
    rename("Pop" = Pop) %>%
    rename(F =
             "f") %>%
    rename(M = "m") %>%
    as.data.frame()
  x <- x %>% select(!contains("dev"))
  means <-
    x %>%
    pivot_longer(
      cols = c(M.mu, F.mu),
      names_to = "Sex",
      values_to = "no"
    ) %>%
    select(-c(M, F))
  means$Sex <-
    do.call(rbind.data.frame, strsplit(means$Sex, split = "\\."))[, 1]
  size <-
    x %>%
    pivot_longer(cols = c(F, M),
                 names_to = "Sex",
                 values_to = "N") %>%
    select(!contains("mu"))
  x <-
    full_join(means, size, by = c("Trait", "Pop", "Sex")) %>%
    pivot_wider(names_from = Trait, values_from = no) %>%
    mutate(Sex = factor(Sex, levels = c("F", "M"))) %>%
    as.data.frame()
  p <- NCOL(x) - 3
  g <- NROW(x)
  Rank <- min(c(g - 1, p))
  message()
  if (q > Rank) {
    q <- Rank
    warning("q has been decreased to (", Rank, ") match the matrix rank")
  } else {
    message("The maximum possible value of q is (", Rank, ").")
  }

  calc.B0 <- function() {
    y <- as.data.frame(x[, -c(1:3)])
    G.mean <- apply(y, 2, mean)
    B0 <- matrix(0, ncol = p, nrow = p)
    for (i in seq_len(g))
    {
      d <- as.numeric(y[i,]) - G.mean
      B0 <- B0 + (d %*% t(d))
    }
    B0 <- B0 / (g - 1)
    return(B0)
  }

  B0 <- calc.B0()
  eig.struc <- eigen(solve(W) %*% B0)
  e.vecs <- Re(eig.struc$vec[, seq_len(Rank)])
  e.vecs <-
    e.vecs / (rep(1, p) %o% sqrt(diag(t(e.vecs) %*% W %*% e.vecs)))

  Center <- diag(g) - 1 / g
  y <- x[, -c(1:3)]
  X <- Center %*% as.matrix(y)
  CVs <- X %*% e.vecs[, seq_len(min(c(q, Rank)))]
  CVs <- data.frame(CVs)
  names(CVs) <- paste0("x", seq_along(CVs))
  CVs <-
    cbind.data.frame(
      CVs,
      Pop = x$Pop,
      Sex = x$Sex,
      stringsAsFactors = FALSE
    )
  CVs$Sex <- factor(CVs$Sex, levels = c("F", "M"))
  CVs <- arrange(CVs, Pop, Sex)
  means <- CVs %>%
    select(-Pop) %>%
    group_by(Sex) %>%
    summarise_all(mean)
  if (means[means$Sex == "F", "x1"] > means[means$Sex == "M", "x1"]) {
    CVs$x1 <- CVs$x1 * -1
  }
  if (q > 1 &&
      means[means$Sex == "F", "x2"] > means[means$Sex == "M", "x2"]) {
    CVs$x2 <- CVs$x2 * -1
  }
  pairs <- utils::combn(levels(x$Pop), 2, simplify = FALSE)

  names(pairs) <- lapply(pairs, paste, collapse = "-")
  pairs <-
    do.call(cbind.data.frame, c(pairs, stringsAsFactors = FALSE))

  chi <- function(i) {
    q1_M <- as.numeric(subset(CVs, Pop == pairs[1, i] &
                                Sex == "M", seq_len(q))[seq_len(q)])
    q1_F <-
      as.numeric(subset(CVs, Pop == pairs[1, i] &
                          Sex == "F", seq_len(q))[seq_len(q)])
    q2_M <-
      as.numeric(subset(CVs, Pop == pairs[2, i] &
                          Sex == "M", seq_len(q))[seq_len(q)])
    q2_F <-
      as.numeric(subset(CVs, Pop == pairs[2, i] &
                          Sex == "F", seq_len(q))[seq_len(q)])
    N1_M <- as.numeric(subset(x, Pop == pairs[1, i] &
                                Sex == "M", 3)[1])
    N1_F <- as.numeric(subset(x, Pop == pairs[1, i] &
                                Sex == "F", 3)[1])
    N2_M <- as.numeric(subset(x, Pop == pairs[2, i] &
                                Sex == "M", 3)[1])
    N2_F <- as.numeric(subset(x, Pop == pairs[2, i] &
                                Sex == "F", 3)[1])

    dif <-
      q1_M + q2_F - q1_F - q2_M
    X2 <-
      as.numeric(t(dif) %*% dif) / (1 / N1_F + 1 / N1_M + 1 / N2_F + 1 / N2_M)

    X2 <- round(X2, digits)
    DF <- q
    p <-
      round(pchisq(X2, DF, lower.tail = lower.tail), digits)
    out <- data.frame(statistic = X2,
                      df = DF,
                      p.value = p)
    out <- out %>% mutate(signif = with(
      .data,
      case_when(
        p >= 0.05 ~ "ns",
        p < 0.05 & p > 0.01 ~ "*",
        p <= 0.01 & p > 0.001 ~ "**",
        p <= 0.001 ~ "***"
      )
    ))

    out
  }
  chi <-
    Vectorize(FUN = chi,
              vectorize.args = "i",
              SIMPLIFY = FALSE)

  pairs_list <- chi(seq_along(pairs))
  pairs_df <-
    data.frame(populations = names(pairs), do.call(rbind.data.frame, pairs_list))

  if (isTRUE(plot) && q > 1) {
    line1 <-
      vector(mode = "double", length = nlevels(x$Pop))
    line2 <-
      vector(mode = "double", length = nlevels(x$Pop))
    point1 <-
      vector(mode = "double", length = nlevels(x$Pop))
    point2 <-
      vector(mode = "double", length = nlevels(x$Pop))

    for (i in seq(1, (g - 1), 2))
    {
      line1[i] <- CVs[i:(i + 1), 1][1]
      line2[i] <- CVs[i:(i + 1), 2][1]
    }
    for (i in seq(2, g, 2))
    {
      point1[i] <- CVs[i, 1]
      point2[i] <- CVs[i, 2]
    }
    line1 <-
      line1 %>%
      as.data.frame() %>%
      mutate(n = row_number()) %>%
      filter(n %% 2 != 0) %>%
      select(-n) %>%
      drop_na()
    line2 <-
      line2 %>%
      as.data.frame() %>%
      mutate(n = row_number()) %>%
      filter(n %% 2 != 0) %>%
      select(-n) %>%
      drop_na()
    point1 <-
      point1 %>%
      as.data.frame() %>%
      mutate(n = row_number()) %>%
      filter(n %% 2 == 0) %>%
      select(-n) %>%
      drop_na()
    point2 <-
      point2 %>%
      as.data.frame() %>%
      mutate(n = row_number()) %>%
      filter(n %% 2 == 0) %>%
      select(-n) %>%
      drop_na()
    male <- x$Pop[seq(2, g, 2)]
    female <- x$Pop[seq(1, (g - 1), 2)]

    CVs <-
      cbind_fill2(CVs, line1, line2, point1, point2, female, male)

    names(CVs)[(q + 3):ncol(CVs)] <-
      c("line1", "line2", "point1", "point2",
        "female", "male")
    x_all <- as.numeric(c(CVs$x1, CVs$x2))
    min_a <- min(x_all) - 1
    max_a <- 1 + max(x_all)
    q_plot <- ggplot(CVs,
                     aes(
                       x = x1,
                       y = x2,
                       color = Sex,
                       fill = Sex
                     )) +
      scale_x_continuous(limits = c(min_a, max_a)) +
      scale_y_continuous(limits = c(min_a, max_a)) +
      geom_segment(
        aes(
          xend = point1,
          yend = point2,
          x = line1,
          y = line2
        ),
        color = "black",
        na.rm = TRUE
      ) +
      geom_point(aes(x = point1,
                     y = point2),
                 color = "blue",
                 na.rm = TRUE) +
      geom_point(aes(x = line1,
                     y = line2),
                 color = "red",
                 na.rm = TRUE) +
      geom_text(
        na.rm = TRUE,
        aes(x = point1,
            y = point2,
            label = male),
        nudge_y = 0.05,
        vjust = "inward",
        hjust = "inward",
        color = "blue",
        check_overlap = TRUE
      ) +
      geom_text(
        aes(x = line1,
            y = line2,
            label = female),
        check_overlap = TRUE,
        color = "red",
        na.rm = TRUE,
        nudge_y = 0.05,
        vjust = "inward",
        hjust = "inward"
      ) +
      scale_color_manual(values = c("red", "blue")) +
      scale_fill_manual(values = c("red",
                                   "blue")) +
      guides(fill = guide_legend(override.aes = list(color = c("red",
                                                               "blue")))) +
      xlab("First Canonical Axis") +
      ylab("Second Canonical Axis") +
      theme_minimal() +
      theme(legend.title = element_blank())
    plot(q_plot)
    pairs_df
  } else {
    pairs_df
  }
}
