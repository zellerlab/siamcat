#' @title LiblineaR Classification Learner
#' @name LearnerClassifLiblineaR
#'
#' @details Type of SVC depends on \code{type} argument:
#' \itemize{
#' \item \code{0} L2-regularized logistic regression (primal)
#' \item \code{1} L2-regularized L2-loss support vector classification (dual)
#' \item \code{3} L2-regularized L1-loss support vector classification (dual)
#' \item \code{2} L2-regularized L2-loss support vector classification (primal)
#' \item \code{4} Support vector classification by Crammer and Singer
#' \item \code{5} L1-regularized L2-loss support vector classification
#' \item \code{6} L1-regularized logistic regression
#' \item \code{7} L2-regularized logistic regression (dual)
#' }
#' If number of records > number of features, \code{type = 2} is faster 
#' than  \code{type = 1}
#' (Hsu et al. 2003).
#'
#' Note that probabilistic predictions are only available for 
#' types \code{0}, \code{6}, and \code{7}.
#' The default \code{epsilon} value depends on the \code{type} parameter, 
#' see [LiblineaR::LiblineaR].
#' 
#' @encoding UTF-8
# taken from https://github.com/mlr-org/mlr3extralearners
LearnerClassifLiblineaR <- R6::R6Class("LearnerClassifLiblineaR",  
    inherit = LearnerClassif, public = list(
        #' @description
        #' #' Creates a new instance of this [R6][R6::R6Class] class.
        initialize = function() {
            ps = ps(
                type = p_int(default = 0, lower = 0, upper = 7, tags = "train"),
                cost = p_dbl(default = 1, lower = 0, tags = "train"),
                epsilon = p_dbl(lower = 0, tags = "train"),
                bias = p_dbl(default = 1, tags = "train"),
                cross = p_int(default = 0L, lower = 0L, tags = "train"),
                verbose = p_lgl(default = FALSE, tags = "train"),
                wi = p_uty(default = NULL, tags = "train"),
                findC = p_lgl(default = FALSE, tags = "train"),
                useInitC = p_lgl(default = TRUE, tags = "train")
            )
            # 50 is an arbitrary choice here
            ps$add_dep("findC", "cross", CondAnyOf$new(seq(2:50)))
            ps$add_dep("useInitC", "findC", CondEqual$new(TRUE))
            
            super$initialize(
                id = "classif.liblinear",
                packages = "LiblineaR",
                feature_types = "numeric",
                predict_types = c("response", "prob"),
                param_set = ps,
                properties = c("twoclass", "multiclass"),
            )
        }
    ),
    private = list(
        .train = function(task) {
            pars = self$param_set$get_values(tags = "train")
            data = task$data()
            train = task$data(cols = task$feature_names)
            target = task$truth()
            
            type = ifelse(is.null(pars$type), 0, pars$type)
            pars = pars[names(pars) != "type"]
            
            mlr3misc::invoke(LiblineaR::LiblineaR, data = train, 
                target = target, type = type, .args = pars)
        },
        .predict = function(task) {
            newdata = task$data(cols = task$feature_names)
            
            type = ifelse(is.null(self$param_set$values$type), 0, 
                self$param_set$values$type)
            
            if (!type %in% c(0, 6, 7) && self$predict_type == "prob") {
                stop("'prob' predict_type only possible if",
                    " `type` is `0`, `6`, or `7`.")
            }
            
            if (self$predict_type == "prob") {
                return(list(
                    prob = mlr3misc::invoke(predict, self$model, 
                        newx = newdata, proba = TRUE)$probabilities))
            } else {
                return(list(
                    response = mlr3misc::invoke(predict, self$model, 
                        newx = newdata)$predictions))
            }
        }
    )
)
