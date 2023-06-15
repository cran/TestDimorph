#' @title Greene t test of Sexual Dimorphism
#' @description Calculation and visualization of the differences in degree
#' sexual dimorphism between two populations using summary statistics as
#' input.
#' @param x A data frame containing summary statistics.
#' @param Pop Number of the column containing populations' names, Default: 1
#' @param plot Logical; if TRUE graphical matrix of p values, Default: FALSE
#' @param colors color palette used in the corrplot
#' @param alternative a character string specifying the alternative
#' hypothesis, must be one of "two.sided", "greater" or "less", Default:"two.sided"
#' @param padjust Method of p.value adjustment for multiple comparisons
#' following \link[stats]{p.adjust}  Default: "none".
#' @param letters Logical; if TRUE returns letters for pairwise comparisons
#' where significantly different populations are given different letters,
#' Default: FALSE'
#' @param digits Number of significant digits, Default: 4
#' @param CI confidence interval coverage takes value from 0 to 1, Default: 0.95.
#' @return data frame of t.test results
#' @details The input is a data frame of summary statistics where the column
#' containing population names is chosen by position (first by default), other
#' columns of summary data should have specific names (case sensitive) similar
#' to \link{baboon.parms_df}.For the visualization of pairwise comparisons using
#' the corrplot, the rounder the image in the plot grid the lower the p-value
#' (see the color scale for similar information). The default colors used in the
#' corrplot are from the "MetBrewer" "Egypt" palette which is listed under the
#' "colorblind_palettes". Different colors palettes can be selected from
#' "RColorBrewer" package.
#' @examples
#' # Comparisons of femur head diameter in four populations
#' df <- data.frame(
#'   Pop = c("Turkish", "Bulgarian", "Greek", "Portuguese"),
#'   m = c(150.00, 82.00, 36.00, 34.00),
#'   f = c(150.00, 58.00, 34.00, 24.00),
#'   M.mu = c(49.39, 48.33, 46.99, 45.20),
#'   F.mu = c(42.91, 42.89, 42.44, 40.90),
#'   M.sdev = c(3.01, 2.53, 2.47, 2.00),
#'   F.sdev = c(2.90, 2.84, 2.26, 2.90)
#' )
#' t_greene(
#'   df,
#'   plot = TRUE,
#'   padjust = "none"
#' )
#' @rdname t_greene
#' @export
#' @importFrom stats qt pt
#' @importFrom utils combn
#' @importFrom multcompView multcompLetters vec2mat
#' @importFrom corrplot corrplot
#' @importFrom dplyr contains
#' @references
#'   # for the t-test
#'
#'   Greene, David Lee. "Comparison of t-tests for differences in sexual
#'   dimorphism between populations." American Journal of Physical Anthropology
#'   79.1 (1989): 121-125.
#'
#'   Relethford, John H., and Denise C. Hodges. "A statistical test for
#'   differences in sexual dimorphism between populations." American Journal of
#'   Physical Anthropology 66.1 (1985): 55-61.
#'
#'   #For the femur head diameter data
#'
#'   F. Curate, C. Umbelino, A. Perinha, C. Nogueira, A.M. Silva, E.
#'   Cunha, Sex determination from the femur in Portuguese populations with
#'   classical and machinelearning classifiers, J. Forensic Leg. Med. (2017) ,
#'   doi:http://dx.doi.org/10.1016/j. jflm.2017.08.011.
#'
#'   O. Gulhan, Skeletal Sexing Standards of Human Remains in Turkey (PhD thesis), Cranfield
#'   University, 2017 [Dataset].
#'
#'   P. Timonov, A. Fasova, D. Radoinova, A.Alexandrov, D. Delev, A study of sexual dimorphism
#'   in the femur among contemporary Bulgarian population, Euras. J. Anthropol. 5 (2014) 46–53.
#'
#'   E.F. Kranioti, N. Vorniotakis, C. Galiatsou, M.Y. Iscan , M.
#'   Michalodimitrakis, Sex identification and software development using
#'   digital femoral head radiographs, Forensic Sci. Int. 189 (2009) 113.e1–7.
#'
t_greene <- function(x,
                     Pop = 1,
                     plot = FALSE,
                     colors = c(
                       "#DD5129", "#985F51", "#536D79", "#0F7BA2", "#208D98", "#319F8E",
                       "#43B284", "#7FB274", "#BCB264", "#FAB255"
                     ),
                     alternative = c("two.sided", "less", "greater"),
                     padjust = "none",
                     letters = FALSE,
                     digits = 4,
                     CI = 0.95) {
  # t-test for a data.frame -------------------------------------------------

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
  if (nrow(x) < 2) {
    stop("x should at least have 2 rows")
  }
  padjust <- match.arg(padjust, choices = p.adjust.methods)
  x <- x %>%
    drop_na() %>%
    as.data.frame() %>% rename(Pop=all_of(Pop))
  if (length(dplyr::contains("-", vars = x$Pop)) != 0) {
    x$Pop <-
      gsub(
        x = x$Pop,
        pattern = "[^[:alnum:]]",
        replacement = "_"
      )
  }
  x$Pop <- factor(x$Pop, levels = x$Pop)
  x$Pop <- droplevels(x$Pop)
  if (length(unique(x$Pop)) != length(which(!is.na(x$Pop)))) {
    warning("Population names are not unique")
  }
  pairs <- utils::combn(x$Pop, 2, simplify = FALSE)

  names(pairs) <- lapply(pairs, paste, collapse = "-")

  tg <-
    lapply(pairs, function(y) {
      t_test(
        m = x[y[1], "m"],
        f = x[y[1], "f"],
        m2 = x[y[2], "m"],
        f2 = x[y[2], "f"],
        M.mu = x[y[1], "M.mu"],
        F.mu = x[y[1], "F.mu"],
        M.mu2 = x[y[2], "M.mu"],
        F.mu2 = x[y[2], "F.mu"],
        M.sdev = x[y[1], "M.sdev"],
        F.sdev = x[
          y[1],
          "F.sdev"
        ],
        M.sdev2 = x[y[2], "M.sdev"],
        F.sdev2 = x[y[2], "F.sdev"],
        alternative = alternative,
        CI = CI,
        digits = digits
      )
    })
  tg <- do.call(rbind.data.frame, tg)
  tg <- rown_col(tg, var = "populations")
  tg$p.value <-
    padjust_n(
      p = tg$p.value,
      method = padjust,
      n = ((nlevels(x$Pop)^2 - nlevels(x$Pop)) / 2)
    )
  tg <- add_sig(tg)

  # Pairwise comparisons and corrplot ---------------------------------------

  pval <- tg$p.value
  names(pval) <- tg$populations
  pmatrix <- multcompView::vec2mat(pval)
  if (!is.logical(letters)) {
    stop("letters should be either TRUE or FALSE")
  }
  if (isTRUE(letters)) {
    tg <-
      list(
        "t.test" = tg,
        "pairwise letters" = rown_col(
          data.frame("letters" = multcompView::multcompLetters(pval,
            threshold = 1 - CI
          )[[1]]),
          var = "populations"
        )
      )
  }
  if (!is.logical(plot)) {
    stop("plot should be either TRUE or FALSE")
  }
  if (isTRUE(plot)) {
    plot_list <- structure(list(
      "t.greene" = tg,
      corrplot(
        corr = pmatrix,
        main = "p-values", method = "ellipse",
        type = "lower",
        mar = c(0, 0, 1, 0),
        col = colors,
        tl.cex = 0.8,
        tl.col = "black",
        insig =
          "blank",
        tl.srt = 0.1,
        pch.cex = 2.5,
        tl.pos = "ld",
        win.asp = 1,
        number.cex = 0.5,
        na.label = "NA",
        diag = FALSE,
        is.corr = FALSE
      ),
      class = "tg"
    ))
    print_tg <- function(x) {
      x[[1]]
    }

    print_tg(plot_list)
  } else {
    tg
  }
}
