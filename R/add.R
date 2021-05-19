#' basic arithmetic
#'
#' @description Return sum of two numbers
#'
#' @param x,y numbers
#'
#'
#' @return \code{add} sum of and y \code{multiply} multiple of x and y
#' @export
#'
#'
#' @examples add(1,2)
add = function(x,y) {
  return(x+y)
}

#' @rdname add
#' @export
#'
#' @examples multiply(1,2)
multiply = function(x,y){
  return(x*y)
}

#' a variant of glmnet
#'
#' @param covariate the model matrix
#' @param response the response vector
#' @param ... other parameter passed to \code{glmnet}
#'
#' @return a glmnet model
#' @export
#' @seealso \link[glmnet]{glmnet}
#'
#' @examples #no examples
my_glmnet = function(covariate,response,...){
  model = glmnet::glmnet(x = covariate,y = response,...)
  return(model)
}
