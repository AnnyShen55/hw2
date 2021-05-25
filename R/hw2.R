

#' Solving Linear System by Gauss-Seidel and Jacobi
#'
#' This function solves a linear system using Gauss-Seidel or Jacobi method,
#' if choose to do parallel, it allows user to specify how many
#' cores to use for parallel implementation.
#'
#' @param X The design matrix
#' @param Y The response variable vector
#' @param algorithm Includes "Gauss-Seidel", "Jocobi",
#' and "parallel Jacobi". Can use integers 1, 2, and 3 respectively, for these
#' three choices.
#' If not specified, it would do "Gauss-Seidel".
#' @importFrom foreach %dopar%
#' @param ncores The number of the cores that the user would like to use for
#' the parallel. If not specified, it would use the number of the system processors.
#' @param max_it Max number of iteration. Default setting is 10^7.
#' @param tolerance  Return when difference between two consecutive steps of beta are
#' smaller then this value. Default setting is 10^(-4).
#' @return It returns the estimator vector of the parameters.
#' It would pop up error messages when not convergence,
#' shape of the design matrix and response are different,
#' and unknown algorithm name.
#' @export
#'
#' @examples n = 100
#' D = diag(rep(3, n))
#' U = rbind(cbind(rep(0, n - 1), diag(rep(-1, n - 1))), rep(0, n))
#' L = t(U)
#' v = rep(c(1, 0), as.integer(n / 2))
#' print(L+D+U)
#' b = (L + D + U) %*% v
#' X = L + D + U
#' solve_ols(X,b,3)
solve_ols = function(X, Y,
                     algorithm = "Gauss-Seidel",
                     ncores = as.numeric(Sys.getenv("NUMBER_OF_PROCESSORS", "2")), max_it = 10^7, tolerance = 10^(-4)) {
  X = as.matrix(X)
  Y = as.matrix(Y)
  if(nrow(X) != ncol(X)){

    return(print("Design matrix should be a square matrix."))
  }

  if (nrow(X) != nrow(Y) && nrow(X) != ncol(Y)){
    return(print("The dimensions of the two inputs don't match."))
  }


  n = nrow(X)
  #get D, U and L
  D = diag(diag(X))
  #print(D)
  U = ramify::triu(X, diag = FALSE)
  #print(U)
  L = ramify::tril(X, diag = FALSE)
  #print(L)
  R_Gauss = -solve(D) %*% (L + U)
  norm_Gauss = norm(R_Gauss, type = "2")
  if(norm_Gauss >= 1){
    return(print("Not converge."))
  }
  beta_update = rep(0, n)

  if (algorithm == "Gauss-Seidel" | algorithm == 1) {

    for (j in 1:max_it) {
      beta_update_old = beta_update
      beta_update = solve(L + D) %*% (Y - U %*% beta_update)

      if (norm(beta_update_old - beta_update, type = "2") < tolerance) {
        break
      }

    }
  }else if(algorithm == "Jocobi" | algorithm == 2) {
    for (j in 1:max_it) {
      beta_update_old = beta_update
      beta_update = solve(D) %*% (Y - (L + U) %*% beta_update)

      if (norm(beta_update_old - beta_update, type = "2") < tolerance) {
        break
      }
    }
  }else if(algorithm == "parallel Jacobi" | algorithm == 3) {
    cl = parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    for (i in 1: max_it){
      beta_update_old = beta_update
      beta_update = unlist(foreach::foreach(j = 1:n, .multicombine = TRUE) %dopar% {
(Y[j] - (X[j,-j] %*%beta_update[-j]))/D[j,j]})
if (norm(beta_update_old - beta_update, type = "2") < tolerance) {
  break
}
      #print(i)
    }
    parallel::stopCluster(cl)
  }else{
    print("Unknown algorithm.")
  }

  #print(beta_update)
  return(beta_update)
}

#' Algorithmic Leveraging
#'
#' This function implements algorithmic leveraging for
#' linear regression using uniform and leverage score, and the user can
#' specific subsampling of rows.

#'
#' @param X The design matrix
#' @param Y The response variable vector
#' @param algorithm Including "uniform" and "leverage", user can choose to
#' use 1 and 2, respectively，for these two choices. Default setting is "leverage".
#' @param sub_rows Number of subsamples to generate. Default setting is
#' 100.
#'
#' @return The parameter estimators.
#' @export
#'
#' @examples n = 500
#' X = matrix(rt(n, 6),n)
#' eps = rnorm(n)#error term
#' Y = - X + eps
#' algo_leverage(X,Y)
algo_leverage = function(X,
                         Y,
                         algorithm = "leverage",
                         sub_rows = 100) {
  X = as.matrix(X)
  Y = as.matrix(Y)
  n = nrow(X)
  if (sub_rows <= 0){
    return(print("sub_row should be a positive number."))
  }
  if (nrow(X) != nrow(Y) && nrow(X) != ncol(Y)){
    return(print("The dimensions of the two inputs don't match."))
  }
  if (algorithm == "uniform" | algorithm == 1) {
    uni_sample = sample(n, size = sub_rows, replace = TRUE)
    X_unif = X[uni_sample,]
    Y_unif = Y[uni_sample,]
    beta_unif = solve(t(X_unif) %*% X_unif) %*% t(X_unif) %*% Y_unif
    return(beta_unif)
  }else if(algorithm == "leverage" | algorithm == 2) {
    pi_lev_full_unnorm = apply(X, 1, function(x)
      x %*% solve(t(X) %*% X) %*% x)
    pi_lev_full = pi_lev_full_unnorm / as.vector(sqrt(pi_lev_full_unnorm %*% pi_lev_full_unnorm))
    lev_sample = sample(n,
                        size = sub_rows,
                        replace = TRUE,
                        prob = pi_lev_full)
    X_lev = X[lev_sample, , drop = F]
    Y_lev = Y[lev_sample,]
    cnt_element = sapply(1:n, function(x)
      sum(lev_sample == x))
    w_part = cnt_element / sub_rows * pi_lev_full
    w = diag(w_part)
    beta_blev = solve(t(X) %*% w %*% X) %*% t(X) %*% w %*% Y
    return(beta_blev)
  } else{
    print("Unknown algorithm.")
  }
}

#' Elastic Net by Coordinate Descent Algorithm
#'
#' This function fits elastic net to data using coordinate descent algorithm.
#'
#' @param X The design matrix
#' @param Y The response vector
#' @param lambda The hyperparameter for the restriction term(s).
#' By default it's 0.1.
#' @param alpha The hyperparameter of the elastic net. By default alpha = 0,
#' which is of the form of ridge regression.
#' @param max_it The maximum number of iteration. By default it's 10^6.
#' @param tolerance  Return when difference between two consecutive steps of beta are
#' smaller then this value.  By default it's 10^(-4)
#' @return The parameter estimates
#' @export
#'
#' @examples p = 20
#' n  = 50
#' beta_0 = rep(0, p)
#' alpha = c(0, 0.5, 1)
#' beta_partial = c(2, 0,-2, 0, 1, 0, -1, 0)
#' beta = c(beta_partial, rep(0, 12))
#' corr_mat = diag(20)
#' corr_mat[1, 2] = 0.8
#' corr_mat[2, 1] = 0.8
#' corr_mat[5, 6] = 0.8
#' corr_mat[6, 5] = 0.8
#' eps = matrix(rnorm(n), n)
#' X = MASS::mvrnorm(n, rep(0, p), Sigma = corr_mat)
#' Y = X[,] %*% beta + eps
#' elnet_coord(X,Y)
#' print(elnet_coord(X,Y))
elnet_coord = function(X, Y, lambda = 0.1, alpha = 0.5, max_it = 10^6, tolerance = 10^(-4)) {
  X = as.matrix(X)
  Y = as.matrix(Y)
  n = nrow(X)
  p = ncol(X)
  if (nrow(X) != nrow(Y) && nrow(X) != ncol(Y)){
    return(print("The dimensions of the two inputs don't match."))
  }
  if (alpha > 1  | alpha < 0){
    return(print("alpha should between 0 and 1 (inclusive)."))
  }
  r = rep(0, n)
  beta = rep(0,p)
  k = 1
  while (k <= max_it) {
    beta_prev = beta
    for (j in 1:p) {

      r = Y - X[, -j] %*% beta[-j]
      z = mean(X[, j] * r)
      beta[j] = soft_thres_helper(z, lambda, alpha, X[, j])
    }
    if (norm(beta-beta_prev,'2') < tolerance) {
      break
    }
    k = k + 1
  }

  return(beta)
}

soft_thres_helper = function(z, lambda, alpha, x_col) {
  if ((abs(z) - lambda * alpha) < 0) {
    return(0)
  } else{
    return(sign(z) * (abs(z) - lambda * alpha) / (lambda * (1 - alpha) + mean(x_col^2)))
  }
}

