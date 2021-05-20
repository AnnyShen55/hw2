

#' This function solves a linear system using Gauss-Seidel or Jacobi method,
#' if choose to do parallel, it allows user to specify how many
#' cores to use for parallel implementation.
#'
#' @param X The design matrix
#' @param Y The response variable vector
#' @param algorithm includes "Gauss-Seidel", "Jocobi", "parallel Gauss-Seidel",
#' and "parallel Jacobi". Can use integers 1, 2, 3, and 4 respectively.
#' If not specified, it would do "Gauss-Seidel".
#' @param ncores #the number of the cores that the user would like to use for
#' the parallel. If not specified, it would use the number of the system processors.
#' @param max_it #max number of iteration. Default setting is 10^4.
#'
#' @return the estimator vector of the parameters.
#' @export
#'
#' @examples n = 100
#' D = diag(rep(3, n))
#' U = rbind(cbind(rep(0, n - 1), diag(rep(-1, n - 1))), rep(0, n))
#' L = t(U)
#' v = rep(c(1, 0), as.integer(n / 2))
#' b = (L + D + U) %*% v
#' solve_ols(D+U+L,b,3)
solve_ols = function(X, Y,
                     algorithm = "Gauss-Seidel",
                     ncores = as.numeric(Sys.getenv("NUMBER_OF_PROCESSORS")), max_it = 10^4) {
  if(nrow(X) != ncol(X)){
    stop("design matrix should be a square matrix.")
  }

  #get D, U and L
  D = diag(X)
  U = upper.tri(X, diag = FALSE)
  L = lower.tri(X, diag = FALSE)

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
    output = foreach::foreach(j = 1:max_it, .multicombine = TRUE) %foreach::dopar% {
      x_update = solve(L + D) %*% (Y - U %*% x_update)
      x_update
    }
    x_update = unlist(output)
  }else if(algorithm == "parallel Jacobi" | algorithm == 4) {
    cl = parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    output = foreach::foreach(j = 1:max_it, .multicombine = TRUE) %foreach::dopar% {
      x_update = solve(D) %*% (Y - (L + U) %*% x_update)
      x_update
    }
    x_update = unlist(output)
  }else{
    stop("unknown algorithm.")
  }

  return(x_update)
}

algo_leverage = function(X,
                         Y,
                         algorithm = "leverage",
                         subsampling_rows = 100) {
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
    stop("unknown algorithm.")
  }
}

elnet_coord = function(X, Y, lambda = 0.5, alpha = 0, max_it = 10^4) {
  n = nrow(X)
  p = ncol(X)

  r = rep(0, n)
  beta = rep(0,p)
  k = 1
  while (k <= max_it) {
    beta_prev = beta
    for (j in 1:p) {
      for (i in 1:n) {
        r[j] = Y[i] - X[i, -j] %*% beta[-j]
      }
      beta[j] = soft_thres_helper(1 / n * X[, j] %*% r, lambda, alpha)
    }
    if (norm(beta_prev - beta, type = "2") < 10 ^ (-4)) {
      break
    }
    k = k + 1
  }
  return(beta)
}

soft_thres_helper = function(z, lambda, alpha) {
  if ((abs(z) - lambda * alpha) <= 0) {
    return(0)
  } else{
    return(sign(z) * (abs(z) - lambda * alpha) / (2 * lambda * (1 - alpha) + 1))
  }
}

