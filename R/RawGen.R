#' @title Raw Data Generation By Truncated Distribution
#' @description Generates raw data from summary statistics using left truncated
#'   normal distribution
#' @param v a data frame or a named list containing summary statistics for
#'   different parameters
#' @param format The form of the resultant data frame either : \code{'long'} or
#'   : \code{'wide'}, Default: \code{'wide'}
#' @param complete_cases Logical; if \code{TRUE} rows with missing values will
#'   be removed, Default: \code{FALSE}
#' @details Data can be entered as a data frame,with \code{Pop} (first column by
#'   default) containing population names,\code{Parms} containing names of
#'   tested parameter(s),  .mu and .sdev containing means and standard
#'   deviations
#'   with M and F donating males and females respectively. While m&f are the
#'   male and female sample sizes. Data also can be entered as a list of multiple
#'   data frames with similar structure but different parameters are entered as
#'   names of the list.
#' @examples
#'  # Comparison of two femur parameters in two populations
#'  library(TestDimorph)
#'  Pop <- as.factor(rep(c("Bulgarian","Greek"),2))
#'  Parms <- as.factor(c(rep("MXFL", 2),rep("MLD", 2)))
#'  m <-c(82.0,36.0,82.00,36.00)
#'  M.mu <- c(461.80,440.40,27.67,27.74)
#'  M.sdev <- c(19.9,19.6,2.21,1.79)
#'  f <- c(58.00,34.0,58.00,34.00)
#'  F.mu <- c(411.70,409.80,24.89,26.69)
#'  F.sdev <- c(23.2,21.4,1.78,2.42)
#'  df <- cbind.data.frame(Pop,Parms,m,M.mu,M.sdev,f,F.mu,F.sdev)
#'  RawGen(df)
#' @export
#' @rdname RawGen
#'
#' @references \insertRef{HUSSEIN2019}{TestDimorph}
#'
#' @importFrom truncnorm rtruncnorm
#' @importFrom plyr alply
#' @importFrom stringr str_split
#' @importFrom rowr cbind.fill
#' @importFrom stats complete.cases
#' @importFrom utils head
#' @importFrom purrr map_df
#'
RawGen <-
  function(v,
           format = "wide",
           complete_cases = FALSE) {
    u <- function(x) {

      x <- x[order(x$Pop), ]
      x$Pop <- ordered(x$Pop)
      t <- function(x) {
        truncnorm::rtruncnorm(
          n = x$m[1],
          a = 0,
          b = Inf,
          mean = x$M.mu[1],
          sd = x$M.sdev[1]
        )

      }
      t2 <- function(x) {
        truncnorm::rtruncnorm(
          n = x$f[1],
          a = 0,
          b = Inf,
          mean = x$F.mu[1],
          sd = x$F.sdev[1]
        )

      }
      q <- function(x) {
        data.frame(x)
      }

      v <- plyr::alply (.data = x,
                        .margins = 1,
                        .fun = t)
      v2 <- plyr::alply (.data = x,
                         .margins = 1,
                         .fun = t2)

      o <- sapply(v, data.frame)

      o2 <- sapply(v2, data.frame)

      names(o) <- levels(x$Pop)[1:length(o)]

      names(o2) <- levels(x$Pop) [1:length(o)]

      M <- plyr::adply(.data = o,
                       .margins = 1,
                       .fun = q)
      F <- plyr::adply(.data = o2,
                       .margins = 1,
                       .fun = q)
      r <- nrow(M)
      r2 <- nrow(F)
      Sex1 <- cbind(rep("M", r))
      Sex2 <- cbind(rep("F", r2))
      Sex <- rbind(Sex1, Sex2)
      n <- rbind.data.frame(M, F)
      h <- cbind.data.frame(Sex, n)
      colnames(h) <- c("Sex", "Pop", "No")
      return(h)
    }
    if (is.data.frame (v)) {
      v$Pop <- droplevels(v$Pop)
      if(is.null(v$Parms)){
        return(u(v))
      }else{
        v$Parms <- droplevels(v$Parms)
        v$Parms <- ordered(v$Parms)
        v <- by(v, v$Parms, list)
      }
    }
    if (is.null(names(v))) {
      names(v) <- c(1:length(v))

    }

    n <- lapply(v, u)
    if (length(v)==1) {
      n <-as.data.frame(n)
      colnames(n) <-c("Sex", "Pop", "No")
      return(n)
    }

    q <- function(x) {
      h <- rep(names(v)[1], nrow(x))
      cbind(x, "Parms" = h)

    }

    G <- lapply(n, q)
    B <- lapply(G, data.frame)
    A <- do.call(rbind.data.frame, B)
    e <-
      stringr::str_split(
        string = as.character(rownames(A)),
        pattern = "\\.",
        n = 2,
        simplify = TRUE
      )

    A$Parms <- as.factor(e[, 1])
    rownames(A) <- NULL

    if (format == "long") {
      if (nlevels(A$Parms)==1) {
        A$Parms <- NULL
        return(A)

      }else{
        return(A)
      }
    } else{
      zz <- function(A) {
        Female <- function(A) {
          M <-
            head(sort(c(as.numeric(
              table(A$Sex == "F")
            )[2]), decreasing = TRUE), 1)

          SexM <- rep("F", M)

          M2 <-
            rowr::cbind.fill(A$No[which(A$Sex == "F")], fill = NA)

          cbind("Sex" = SexM, M2)
        }
        Male <- function(A) {
          M <-
            head(sort(c(as.numeric(
              table(A$Sex == "M")
            )[2]), decreasing = TRUE), 1)

          SexM <- rep("M", M)

          M2 <-
            rowr::cbind.fill(A$No[which(A$Sex == "M")], fill = NA)

          cbind("Sex" = SexM, M2)
        }


        vvv <- by(data = A,
                  INDICES = A$Parms,
                  FUN = Male)
        M <- do.call(rowr::cbind.fill, c(vvv, list(fill = NA)))
        vvv2 <- by(data = A,
                   INDICES = A$Parms,
                   FUN = Female)
        F <- do.call(rowr::cbind.fill, c(vvv2, list(fill = NA)))
        M <- M[, seq(2, ncol(M), by = 2)]

        M <- as.data.frame(apply(M, 2, as.numeric))
        M$Sex <- as.factor(rep("M", nrow(M)))
        F <- F[, seq(2, ncol(F), by = 2)]
        F <- as.data.frame(apply(F, 2, as.numeric))
        F$Sex <- as.factor(rep("F", nrow(F)))

        j <- rbind(M, F)

        return(j)

      }

      q <- by(data = A,
              INDICES = A$Pop,
              FUN = zz)

      q <- purrr::map_df(.x = q,
                         .f = ~ as.data.frame(.x),
                         .id = "Pop")
      q$Pop <- as.factor(q$Pop)
      levels(q$Pop) <- levels(A$Pop)
      col_idx <- grep("Sex", names(q))
      q <- q[, c(col_idx, (1:ncol(q))[-col_idx])]
      colnames(q)[3:ncol(q)] <- levels(A$Parms)


      if (complete_cases == TRUE) {
        return(q[complete.cases(q),])

      } else{
        return(q)
      }
    }
  }
