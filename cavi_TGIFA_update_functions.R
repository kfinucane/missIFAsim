source("V-TGIFA_utils.R")

update_eta <- function(Y, n, p, k_star, mean_sigma_inv,
                       mean_lambda, var_lambda, mean_mu) {
  var_eta <- lapply(seq_len(n), function(i) {
    armaInv(
     Reduce("+", lapply(1:p, function(j) {
        mean_sigma_inv[j, 1] * (tcrossprod(matrix(mean_lambda[j, ])) + var_lambda[[j]])
      })) + diag(k_star)
    )
  })

  mean_eta <- lapply(seq_len(n), function(i) {
    var_eta[[i]] %*% (Rfast::transpose(mean_lambda) %*%
                        diag(mean_sigma_inv[, 1]) %*%
                        (matrix(Y[i, ]) - mean_mu)
    )
  })

  return(list(
    "var_eta" = var_eta,
    "mean_eta" = matrix(unlist(mean_eta), byrow = TRUE, nrow = n)
  ))
}


update_lambda <- function(Y, n, p, k_star, mean_sigma_inv,
                          mean_eta, var_eta, mean_mu, mean_phi, mean_delta) {
  var_lambda <- lapply(seq_len(p), function(j) {
    armaInv(
      mean_sigma_inv[j, 1] * Reduce("+", lapply(seq_len(n), function(i) {
        (var_eta[[i]] + tcrossprod(matrix(mean_eta[i, ])))
      })) + diag(mean_phi[j, ] * cumprod(mean_delta))
    )
  })

  mean_lambda <- lapply(seq_len(p), function(j) {
    (var_lambda[[j]] * mean_sigma_inv[j, 1]) %*%
      Reduce("+", lapply(seq_len(n), function(i) {
        matrix(mean_eta[i, ]) * (Y[i, j] - mean_mu[j, 1])
      }))
  })

  return(list(
    "var_lambda" = var_lambda,
    "mean_lambda" = matrix(unlist(mean_lambda), byrow = TRUE, nrow = p)
  ))
}


update_sigma_inv <- function(Y, n, p, b_sigma, mean_mu, mean_lambda,
                             mean_eta, var_lambda, var_eta, var_mu) {
  rate_sigma_inv <- sapply(seq_len(p), function(j) {
    b_sigma +
      0.5 * (
        Reduce("+", lapply(seq_len(n), function(i) {
          (Y[i, j] - (mean_mu[j, 1] + t(matrix(mean_lambda[j, ])) %*% matrix(mean_eta[i, ])))^2
        })) +
          sum(diag(
            var_lambda[[j]] %*% Reduce("+", lapply(seq_len(n), function(i) {
              (var_eta[[i]] + tcrossprod(matrix(mean_eta[i, ])))
            }))
          )) +
          t(matrix(mean_lambda[j, ])) %*% Reduce("+", lapply(seq_len(n), function(i) {
            var_eta[[i]]
          })) %*% matrix(mean_lambda[j, ]) +
          var_mu[j, j] * n
      )
  })

  return(matrix(rate_sigma_inv))
}


update_mu <- function(Y, n, p, mean_sigma_inv, mu_varphi, var_eta, mean_lambda,
                      mean_eta, mu_tilde) {
  var_mu <- armaInv(
    n * diag(mean_sigma_inv[, 1]) + diag(mu_varphi)
  )

  mean_mu <- var_mu %*% (
    Reduce("+", lapply(seq_len(n), function(i) {
      diag(mean_sigma_inv[, 1]) %*% (matrix(Y[i, ]) - mean_lambda %*% matrix(mean_eta[i, ]))
    })) +
      diag(mu_varphi) %*% mu_tilde
  )

  return(list(
    "var_mu" = var_mu,
    "mean_mu" = mean_mu
  ))
}

update_alpha <- function(Y, R, trunc.point) {

  a_shape_alpha <- sum((1 - R) * (Y >= trunc.point)) + 1  # scalar

  return(a_shape_alpha)

}


update_phi <- function(p, k_star, kappa_2, mean_delta, var_lambda, mean_lambda){

  rate_phi <- sapply(seq_len(p), function(j){
    sapply(seq_len(k_star), function(h){
      kappa_2 + 0.5 * (
        cumprod(mean_delta)[h] *
          (var_lambda[[j]][h, h] + mean_lambda[j, h]^2)
      )
    })
  })
  return(t(rate_phi))
}



update_delta <- function(p, k_star, mean_phi, shape_delta, mean_delta,
                         var_lambda, mean_lambda) {

  new_mean_delta <- mean_delta

  rate_delta <- matrix(NA, nrow = k_star, ncol = 1)

  for (l in seq_len(k_star)) {

    new_rate <- 1 + 0.5 * (

      Reduce("+", lapply(seq_len(p), function(j){

        Reduce("+", lapply(l:k_star, function(f){
          mean_phi[j, f] * ((cumprod(new_mean_delta[ , 1])[f] / new_mean_delta[l , 1])) *
            (var_lambda[[j]][f, f] + mean_lambda[j, f]^2)
        }))

      }))

    )

    rate_delta[l, 1] <- new_rate

    new_mean_delta[l, 1] <- shape_delta[l, 1] / new_rate

  }

  return(list("rate_delta" = rate_delta,
              "mean_delta" = new_mean_delta))

}


# imputation function
imputer <- function(theta_current, miss_ind, miss_vec, n, p, k.star, trunc.point){

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # generate parameters
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  eta_gen <- matrix(NA, nrow = n, ncol = k.star)
  lambda_gen <- matrix(NA, nrow = p, ncol = k.star)
  sigma_inv_gen <- matrix(NA, nrow = p, ncol = 1)

  for (j in seq_len(p)) {

    lambda_gen[j, ] <- Rfast::rmvnorm(n = 1,
                                      mu = theta_current$mean_lambda[j, ],
                                      sigma = theta_current$var_lambda[[j]],
                                      seed = get_seed()
    )

    sigma_inv_gen[j, ] <- rgamma(n = 1,
                                 shape = theta_current$shape_sigma_inv[j, ],
                                 rate = theta_current$rate_sigma_inv[j, ])

  }

  for (i in seq_len(n)) {

    eta_gen[i, ] <- Rfast::rmvnorm(n = 1,
                                   mu = theta_current$mean_eta[i, ],
                                   sigma = theta_current$var_eta[[i]],
                                   seed = get_seed()
    )

  }

  mu_gen <- matrix(Rfast::rmvnorm(n = 1,
                                  mu = as.vector(theta_current$mean_mu),
                                  sigma = theta_current$var_mu,
                                  seed = get_seed()),
                   ncol = 1)

  alpha_gen <- rbeta(n = 1,
                     shape1 = theta_current$a_shape_alpha,
                     shape2 = theta_current$b_shape_alpha)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # impute
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Z <- rep(NA, length(miss_vec))
  imputed_values <- rep(NA, length(miss_vec))

  if (dim(miss_ind)[1] != 0) {

    impute_mean <- sweep(Rfast::mat.mult(lambda_gen, Rfast::transpose(eta_gen)), 1, mu_gen, "+")  # p x n

    for (point_ind in 1:dim(miss_ind)[1]) {

      impute_var <- (1 / sigma_inv_gen)  # 1 x p

      # ptruncnorm calcs lower.tail
      lower_prob <- truncnorm::ptruncnorm(trunc.point,
                                          a = 0,
                                          mean = impute_mean[miss_ind[point_ind, ][2], miss_ind[point_ind, ][1]],
                                          sd = sqrt(impute_var[miss_ind[point_ind, ][2], 1]))

      # as truncation of dist is 0, Inf
      # and lower_prob is 0, trunc.point
      higher_prob <- 1 - lower_prob

      zij_prob <- (lower_prob) / (lower_prob + alpha_gen * higher_prob)

      zij <- rbinom(1, 1, zij_prob)

      Z[point_ind] <- zij

      if (zij == 0) {

        # impute as MAR
        imputed_values[point_ind] <-
          truncnorm::rtruncnorm(1, a = trunc.point,
                                b = Inf,
                                mean = impute_mean[miss_ind[point_ind, ][2], miss_ind[point_ind, ][1]],
                                sd = sqrt(impute_var[miss_ind[point_ind, ][2], 1]))

      } else {

        # impute as MNAR
        imputed_values[point_ind] <-
          truncnorm::rtruncnorm(1, a = 0, b = trunc.point,
                                mean = impute_mean[miss_ind[point_ind, ][2], miss_ind[point_ind, ][1]],
                                sd = sqrt(impute_var[miss_ind[point_ind, ][2], 1]))

      }

    }

  }

  return(list(
    "Z" = Z,
    "imputed_values" = imputed_values))

}

cavi_updater <- function(input_theta, data_working, n, p, k.star, R, trunc.point,
                         mu_tilde,
                         mu_varphi,
                         kappa_1 = 3L,
                         kappa_2 = 2L,
                         a_sigma = 1L,
                         b_sigma = 0.25,
                         a_1 = 2.1,
                         a_2 = 3.1) {

  # extract parameters
  var_eta <- input_theta$var_eta
  mean_eta <- input_theta$mean_eta
  var_lambda <- input_theta$var_lambda
  mean_lambda <- input_theta$mean_lambda
  var_mu <- input_theta$var_mu
  mean_mu <- input_theta$mean_mu
  shape_sigma_inv <- input_theta$shape_sigma_inv
  rate_sigma_inv <- input_theta$rate_sigma_inv
  a_shape_alpha <- input_theta$a_shape_alpha
  b_shape_alpha <- input_theta$b_shape_alpha
  shape_phi <- input_theta$shape_phi
  rate_phi <- input_theta$rate_phi
  shape_delta <- input_theta$shape_delta
  rate_delta <- input_theta$rate_delta

  # create objects
  mean_sigma_inv <- shape_sigma_inv / rate_sigma_inv

  mean_delta <- shape_delta / rate_delta  # k x 1

  mean_phi <- shape_phi / rate_phi  # p x k


  # update latent factor score
  new_eta_params <- update_eta(data_working, n, p, k.star, mean_sigma_inv,
                               mean_lambda, var_lambda, mean_mu)

  var_eta <- new_eta_params$var_eta

  mean_eta <- new_eta_params$mean_eta


  # update factor loadings
  new_lambda_params <- update_lambda(data_working, n, p, k.star,
                                     mean_sigma_inv, mean_eta, var_eta,
                                     mean_mu, mean_phi, mean_delta)

  var_lambda <- new_lambda_params$var_lambda

  mean_lambda <- new_lambda_params$mean_lambda

  # update variance of idiosyncratic errors
  rate_sigma_inv <- update_sigma_inv(data_working, n, p, b_sigma,
                                     mean_mu, mean_lambda,
                                     mean_eta, var_lambda, var_eta, var_mu)

  mean_sigma_inv <- shape_sigma_inv / rate_sigma_inv

  # update mean
  new_mu_params <- update_mu(data_working, n, p, mean_sigma_inv,
                             mu_varphi, var_eta, mean_lambda, mean_eta,
                             mu_tilde)

  var_mu <- new_mu_params$var_mu

  mean_mu <- new_mu_params$mean_mu

  # update MAR probability
  a_shape_alpha <- update_alpha(data_working, R, trunc.point)

  # update local shrinkage params
  rate_phi <- update_phi(p, k.star, kappa_2, mean_delta, var_lambda, mean_lambda)

  mean_phi <- shape_phi / rate_phi

  # update global shrinkage params
  new_delta_params <- update_delta(p, k.star, mean_phi, shape_delta,
                                   mean_delta, var_lambda, mean_lambda)

  rate_delta <- new_delta_params$rate_delta

  mean_delta <- new_delta_params$mean_delta

  theta <- list("shape_sigma_inv" = shape_sigma_inv,
                "rate_sigma_inv" = rate_sigma_inv,
                "shape_delta" = shape_delta,
                "rate_delta" = rate_delta,
                "shape_phi" = shape_phi,
                "rate_phi" = rate_phi,
                "mean_eta" = mean_eta,
                "var_eta" = var_eta,
                "mean_mu" = mean_mu,
                "var_mu" = var_mu,
                "mean_lambda" = mean_lambda,
                "var_lambda" = var_lambda,
                "a_shape_alpha" = a_shape_alpha,
                "b_shape_alpha" = b_shape_alpha)

  return(theta)

}

