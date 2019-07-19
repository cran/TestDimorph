#' @title Evaluation Of Sex-prediction Accuracy
#' @description Testing the accuracy of different sex prediction models
#'   using the \link[caret]{confusionMatrix} function
#' @param f Formula in the form \code{groups ~ x1 + x2 + ...}. The grouping
#'   factor is placed to the left hand side while the numerical measurements are
#'   placed to the right hand side
#' @param x Data frame to be fitted to the model
#' @param y New data frame to be tested
#' @param byPop Logical;if \code{TRUE} returns the accuracy in different
#'   populations, Default: \code{TRUE}.
#' @param method Different methods of modeling see \code{details} ,
#'   Default:'lda'
#' @param cutoff cutoff value when using logistic regression, Default: 0.5
#' @param ref. reference category in the grouping factor, Default: 'F'
#' @param post. positive category in the grouping factor, Default: 'M'
#' @param ... additional arguments that can passed to modeling and
#'   \link[caret]{confusionMatrix} function.
#' @return Accuracy parameters for the tested model
#' @details Data frames to be entered as input need to be arranged in a similar
#'   manner to \code{\link{Howells}} dataset. Methods used for modeling are:
#'   \describe{ \item{\code{lda}}{linear discriminant analysis}
#'   \item{\code{qda}}{quadratic discriminant analysis}
#'   \item{\code{mda}}{mixture discriminant analysis} \item{\code{fda}}{flexible
#'   discriminant analysis} \item{\code{rda}}{regularized discriminant analysis}
#'   \item{\code{glm}}{binomial logistic regression} }
#'
#' @examples
#' library(TestDimorph)
#' AccuModel(Sex~GOL+NOL+BNL,x = Howells,y = Howells,byPop = FALSE,method = "lda")
#' @export
#' @importFrom dplyr filter
#' @importFrom stats predict binomial relevel
#' @importFrom caret confusionMatrix
#' @importFrom rowr cbind.fill
#' @importFrom rlang .data
#' @importFrom MASS lda qda
#' @importFrom mda mda fda
#' @importFrom klaR rda
AccuModel <-
  function(f,
           x,
           y,
           byPop = TRUE,
           method = "lda",
           cutoff = 0.5,
           ref. = "F",
           post.="M",
           ...) {
    Acc <- function(y) {
      x$Sex <- stats::relevel(x$Sex, ref = ref.)
      y$Sex <- stats::relevel(y$Sex, ref = ref.)
      if (method == "lda") {
        model <- MASS::lda(f, data = x, ...)
        preds <- stats::predict(model, newdata = y)
      }
      if (method == "qda") {
        model <- MASS::qda(f, data = x, ...)
        preds <- stats::predict(model, newdata = y)
      }
      if (method == "mda") {
        model <- mda::mda(f, data = x, ...)
        preds <- stats::predict(model, newdata = y)
        preds <- as.data.frame(preds)
        colnames(preds)[1] <- "class"
      }
      if (method == "fda") {
        model <- mda::fda(f, data = x, ...)
        preds <- stats::predict(model, newdata = y)
        preds <- as.data.frame(preds)
        colnames(preds)[1] <- "class"
      }
      if (method == "rda") {
        model <- klaR::rda(f, data = x, ...)
        preds <- stats::predict(model, newdata = y)
      }
      if (method == "glm") {
        model <-
          stats::glm(
            f,
            family = stats::binomial(link = 'logit'),
            data = x,
            maxit = 100,
            ...
          )
        preds <-
          stats::predict(model, newdata = y, type = 'response')
        class <- ifelse(test = preds > cutoff,
                        yes = "M",
                        no = "F")
        preds <- cbind.data.frame(preds, class)
      }


      xtab <- table(preds$class, y$Sex)

      caret::confusionMatrix(xtab,positive=post.,reference=ref.,...)
    }
    if (byPop == FALSE) {
      Acc(y)


    } else{
      by(y, y$Pop, Acc)
    }
  }
