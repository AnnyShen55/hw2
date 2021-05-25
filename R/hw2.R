

#' This function solves a linear system using Gauss-Seidel or Jacobi method,
#' if choose to do parallel, it allows user to specify how many
#' cores to use for parallel implementation.
#'
#' @param X The design matrix
#' @param Y The response variable vector
#' @param algorithm includes "Gauss-Seidel", "Jocobi", "parallel Gauss-Seidel",
#' and "parallel Jacobi". Can use integers 1, 2, 3, and 4 respectively.
#' If not specified, it would do "Gauss-Seidel".
#' @importFrom foreach %dopar%
#' @param ncores #the number of the cores that the user would like to use for
#' the parallel. If not specified, it would use the number of the system processors.
#' @param max_it #max number of iteration. Default setting is 10^4.
#'
#' @return the estimator vector of the parameters.
#' It would pop up error messages when requirements for either convergence,
#'shape of the design matrix, and algorithm name doesn't meet.
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
#' solve_ols(X,b)
solve_ols = function(X, Y,
                     algorithm = "Gauss-Seidel",
                     ncores = as.numeric(Sys.getenv("NUMBER_OF_PROCESSORS", "2")), max_it = 10^4) {
  X = as.matrix(X)
  Y = as.matrix(Y)
  if(nrow(X) != ncol(X)){

    return(print("Design matrix should be a square matrix."))
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
  x_update = rep(0, n)

  if (algorithm == "Gauss-Seidel" | algorithm == 1) {

    for (j in 1:max_it) {
      x_update_old = x_update
      x_update = solve(L + D) %*% (Y - U %*% x_update)

      if (norm(x_update_old - x_update, type = "2") < 10 ^ (-4)) {
        break
      }

    }
  }else if(algorithm == "Jocobi" | algorithm == 2) {
    for (j in 1:max_it) {
      x_update_old = x_update
      x_update = solve(D) %*% (Y - (L + U) %*% x_update)

      if (norm(x_update_old - x_update, type = "2") < 10 ^ (-4)) {
        break
      }
    }
  }else if(algorithm == "parallel Gauss-Seidel" | algorithm == 3) {
    cl = parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    output = foreach::foreach(j = 1:max_it, .multicombine = TRUE) %dopar% {
      x_update = solve(L + D) %*% (Y - U %*% x_update)

    }
    x_update = unlist(output)
    parallel::stopCluster(cl)
  }else if(algorithm == "parallel Jacobi" | algorithm == 4) {
    cl = parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    output = foreach::foreach(j = 1:max_it, .multicombine = TRUE) %dopar% {
      x_update = solve(D) %*% (Y - (L + U) %*% x_update)
    }
    x_update = unlist(output)
    parallel::stopCluster(cl)
  }else{
    print("Unknown algorithm.")
  }

  #print(x_update)
  return(x_update)
}

#' This function implements algorithmic leveraging for
#' linear regression using uniform and leverage score, and the user can
#' specific subsampling of rows.

#'
#' @param X The design matrix
#' @param Y The response variable vector
#' @param algorithm including "uniform" and "leverage", user can choose to
#' use 1 and 2, respectively. Default setting is "leverage".
#' @param subsampling_rows number of subsamples to generate. Default setting is
#' 100.
#'
#' @return the beta estimators.
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
                         subsampling_rows = 100) {
  X = as.matrix(X)
  Y = as.matrix(Y)
  n = nrow(X)
  if (algorithm == "uniform" | algorithm == 1) {
    uni_sample = sample(n, size = subsampling_rows, replace = TRUE)
    X_unif = X[uni_sample,]
    Y_unif = Y[uni_sample,]
    beta_unif = solve(t(X_unif) %*% X_unif) %*% t(X_unif) %*% Y_unif
    return(beta_unif)
  }else if(algorithm == "leverage" | algorithm == 2) {
    pi_lev_full_unnorm = apply(X, 1, function(x)
      x %*% solve(t(X) %*% X) %*% x)
    pi_lev_full = pi_lev_full_unnorm / sqrt(pi_lev_full_unnorm %*% pi_lev_full_unnorm)
    lev_sample = sample(n,
                        size = subsampling_rows,
                        replace = TRUE,
                        prob = pi_lev_full)
    X_lev = X[lev_sample, , drop = F]
    Y_lev = Y[lev_sample,]
    cnt_element = sapply(1:n, function(x)
      sum(lev_sample == x))
    w_part = cnt_element / subsampling_rows * pi_lev_full
    w = diag(w_part)
    beta_blev = solve(t(X) %*% w %*% X) %*% t(X) %*% w %*% Y
    return(beta_blev)
  } else{
    print("Unknown algorithm.")
  }
}

#' This function fits elastic net to data using coordinate descent algorithm.
#'
#' @param X The design matrix
#' @param Y The response vector
#' @param lambda the hyperparameter for the restriction term(s).
#' By default it's 0.5.
#' @param alpha #the hyperparameter of the elastic net. By default alpha = 0,
#' which is of the form of ridge regression.
#' @param max_it # the maximum number of iteration. By default it's 10^4
#' @param tolerance  # return when difference between two step of beta are
#' smaller then this value.
#' @return # the estimates of the betas
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

