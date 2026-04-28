col_rifa_gibbs <- function(input_data, coding = NA, log.data = FALSE,
                           n.iters = 10000, k.star = 5,
                           verbose = TRUE, burn = 5000, thin = 5,
                           kappa_1 = 3L, kappa_2 = 2L,
                           a_sigma = 1L, b_sigma = 0.25, a_1 = 2.1, a_2 = 3.1,
                           nu = 3L,
                           return_chain = FALSE,
                           cred_level = 0.95,
                           var_mult = 1.0,
                           use_true_params = FALSE, generated_data = NULL,
                           init_imp = "half_min") {

  library(Rfast)
  source("rIFA_utils.R")

  # record inputs to include in output
  inputs <- list("input_data" = input_data,
                 "coding" = coding,
                 "log.data" = log.data,
                 "n.iters" = n.iters,
                 "k.star" = k.star,
                 "verbose" = verbose,
                 "burn" = burn,
                 "thin" = thin,
                 "kappa_1" = kappa_1,
                 "kappa_2" = kappa_2,
                 "a_sigma" = a_sigma,
                 "b_sigma" = b_sigma,
                 "a_1" = a_1,
                 "a_2" = a_2,
                 "nu" = nu,
                 "cred_level" = cred_level,
                 "var_mult" = var_mult,
                 "init_imp" = init_imp)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # formatting data as needed
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  data <- rIFA_prep(input_data, coding = NA, log_trans = log.data)$data
  R <- rIFA_prep(input_data, coding = NA, log_trans = log.data)$missingness

  # data dimensions
  n <- nrow(data)
  p <- ncol(data)

  # truncation point
  trunc.point <- min(data, na.rm = TRUE)

  # indices of missing data
  # col 1 contains row indices in data
  # col 2 contains col indices in data
  # so miss_ind[1, ] are the indices of the first missing point in the data
  miss_ind <- which(R == 0, arr.ind = TRUE)
  miss_vec <- which(R == 0)

  # number of missing points
  n_missing <- length(miss_vec)

  # ~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~
  # initialisation
  # ~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~

  if (verbose) {
    print("Setting initial values")
  }

  if (n_missing != 0) {

    if (use_true_params) {

      # have initial imputation be true values
      data_working <- data
      data_working[is.na(data_working)] <- generated_data$data_gen$data[is.na(data_working)]

    } else {

      if (init_imp == "SVD") {

        # initialise imputation with svd imputation
        start_svd <- softImpute::softImpute(data, rank.max = k.star, type = "svd")

        start_imp <- start_svd$u %*% diag(start_svd$d) %*% t(start_svd$v)

        initial_dataset <- data
        initial_dataset[is.na(initial_dataset)] <- start_imp[is.na(initial_dataset)]

      } else if (init_imp == "half_min") {

        # initialise imputation with half-min
        initial_dataset <- half_min_imp(data, miss_ind)

      }

      data_working <- initial_dataset

    }

  } else {

    # use full data
    data_working <- data

  }

  # if (use_true_params) {
  #
  #   mu_tilde <- matrix(generated_data$data_gen$mu_tilde)
  #   mu_varphi <- generated_data$data_gen$mu_varphi
  #
  # } else {
  #
  #   if (n_missing == 0) {
  #
  #     mu_tilde <- matrix(apply(data, 2, mean, na.rm = TRUE))
  #     mu_varphi <- rep(1L, p)
  #
  #   }
  #
  #   # initialise mu hyperparams
  #   mu_tilde <- matrix(apply(data, 2, mean, na.rm = TRUE)) * 0.9  # multiplied downwards
  #   mu_varphi <- 1 / (0.05 * apply(data, 2, mean, na.rm = TRUE))
  #   mu_varphi[(1:p)[-unique(miss_ind[ , 2])]] <- 1L
  #
  # }

  # use pca to initialise factors and loadings
  pca_res <- prcomp(data_working, scale. = FALSE, center = FALSE)

  lambda <- pca_res$rotation[ , 1:k.star]  # p x k.star
  rownames(lambda) <- NULL
  colnames(lambda) <- NULL

  eta <- matrix(NA, nrow = n, ncol = k.star)

  for (i in 1:n) {

    eta[i, ] <- Rfast::rmvnorm(n = 1,
                               mu = rep(0, k.star),
                               sigma = diag(k.star),
                               seed = get_seed())

  }  # n x k

  mu <- matrix(apply(data, 2, mean, na.rm = TRUE)) - apply(Rfast::mat.mult(lambda, t(eta)), 1, mean)
  mu_tilde <- mu
  mu_varphi <- 1 / (0.03 * apply(data, 2, function(x) abs(range(x, na.rm = TRUE)[2] - range(x, na.rm = TRUE)[1])))
  # mu_varphi <- 1 / (0.05 * apply(data, 2, mean, na.rm = TRUE))
  mu_varphi[(1:p)[-unique(miss_ind[,2])]] <- 1L

  alpha <- runif(1, 0, length(miss_vec) / (n * p))

  # eps <- matrix(NA, nrow = n, ncol = p)
  #
  # for (i in 1:n) {
  #
  #   eps[i, ] <- data_working[i, ] - mu - (lambda %*% matrix(eta[i ,]))
  #
  # }
  #
  # sigma_inv <- matrix(1 / diag(cov(eps, use = "complete")))
  #
  # Sigma_inv <- diag(sigma_inv[ , 1])

  sigma_inv <- matrix(1 / (apply(data, 2, var, na.rm = TRUE) * var_mult))

  Sigma_inv <- diag(sigma_inv[ , 1])

  phi <- matrix(rgamma(n = p * k.star, shape = kappa_1, rate = kappa_2),  # p x k.star
                nrow = p)

  delta <- matrix(rep(0, k.star), nrow = k.star)  # k.star x 1
  delta[1, 1] <- rgamma(n = 1, shape = a_1, rate = 1)
  delta[-1, 1] <- truncdist::rtrunc(spec = "gamma",
                                    n = k.star - 1,
                                    shape = a_2,
                                    rate = 1,
                                    a = 1)

  tau <- matrix(cumprod(delta))  # k.star x 1

  gamma <- matrix(rgamma(n, shape = nu / 2, rate = nu / 2),
                  nrow = n, ncol = 1)  # n x 1

  # ~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~
  # set up storage
  # ~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~

  # imputed dataset
  store_data <- list()

  # alpha - prob of MAR
  store_alpha <- list()

  # sigma
  store_sigma_inv <- list()

  # lambda - factor loadings
  store_lambda <- list()

  # eta
  store_eta <- list()

  # mu - mean
  store_mu <- list()

  # effective factors
  store_k_t <- list()

  # phi
  store_phi <- list()

  # delta
  store_delta <- list()

  # gamma
  store_gamma <- list()

  # Z
  store_Z <- list()

  for (m in seq_len(n.iters)) {


    if (verbose) {
      # create a readout
      statement <- paste("rIFA process running. Now on iteration ", m, " of ", n.iters)
      print(statement)
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # lambda update
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    gamma_eta <- sweep(transpose(eta), 2, c(gamma), FUN = "*")  # k.star x n

    eta_t_eta <- mat.mult(gamma_eta , eta)  # k.star x k.star

    for (j in 1:p) {

      if (k.star == 1) {

        lambda_var <- armaInv(sigma_inv[j] * eta_t_eta + diag(c(phi[j, ] * tau), nrow = 1))  # k.star x k.star

      } else {

        lambda_var <- armaInv(sigma_inv[j] * eta_t_eta + diag(c(phi[j, ] * tau)))  # k.star x k.star

      }

      lambda_mean <- mat.mult(lambda_var, mat.mult(
        gamma_eta, as.matrix(data_working[ , j] - mu[j])) * sigma_inv[j])  # k.star x 1

      lambda[j, ] <- rmvnorm(n = 1, mu = as.vector(lambda_mean),
                             sigma = lambda_var, seed = get_seed())

    }

    # lambda is a p x k.star matrix

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # eta update
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # collapsed gibbs sampler
    lambda_t_Sigma <- mat.mult(Rfast::transpose(lambda), Sigma_inv)

    eta_var <- armaInv(diag(k.star) + mat.mult(lambda_t_Sigma, lambda))  # k.star x k.star

    for (i in 1:n) {

      yi_minus_mu <- as.matrix(data_working[i, ] - mu)

      eta_mean <- mat.mult(eta_var, mat.mult(lambda_t_Sigma, yi_minus_mu))  # k.star x 1

      eta_var_mult <- ( 1 / (nu + p) ) * ( nu + mat.mult(mat.mult(Rfast::transpose(yi_minus_mu), Sigma_inv), yi_minus_mu) )

      eta[i , ] <-  mvnfast::rmvt(n = 1,
                                  mu = as.vector(eta_mean),
                                  sigma = as.numeric(eta_var_mult) * eta_var,
                                  df = nu + p)

    }

    # eta is a n x k.star matrix

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # sigma_inv update
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    lambda_eta_t <- mat.mult(lambda, Rfast::transpose(eta))  # p x n

    sigma_inv <- matrix(rgamma(p, shape = a_sigma + n/2,
                               rate = b_sigma + 0.5 * colsums(sweep((sweep(data_working, 2, mu, FUN = "-") -
                                                                       Rfast::transpose(lambda_eta_t))^2,
                                                                    1, gamma, "*"))))

    # rate correct here as per relationship between Ga with shape and rate and IG with shape and scale

    Sigma_inv <- diag(sigma_inv[ , 1])

    # sigma_inv is a p x 1 matrix
    # Sigma_inv is a p x p diag matrix

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # mu update
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    mu_temp1 <- mu_varphi * diag(p)

    mu_temp2 <- sum(gamma) * Sigma_inv

    mu_var <- armaInv(mu_temp1 + mu_temp2)

    # element-wise multiplication is correct here
    mu_temp3 <- colsums( sweep(matrix(rep(sigma_inv, n), nrow = n, byrow = TRUE), 1, gamma, "*") * (
      data_working - mat.mult(eta, transpose(lambda)))
    )

    mu_temp3 <- matrix(mu_temp3)

    mu_mean <- mat.mult(mu_var, mu_temp3 + mat.mult(mu_temp1, mu_tilde))

    mu <- rmvnorm(n = 1, mu = as.vector(mu_mean),
                  sigma = mu_var, seed = get_seed())


    mu <- matrix(mu, nrow = p) # p x 1

    rm(mu_temp1)
    rm(mu_temp2)
    rm(mu_temp3)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # phi update
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    for (h in 1:k.star) {

      phi[ , h] <- rgamma(p, shape = rep(0.5 + kappa_1, p), rate = kappa_2 + 0.5 * tau[h] * lambda[ , h]^2)

    }

    # phi is a p x k.star matrix

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # delta and tau update
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    delta[1, 1] <- rgamma(1,
                          shape = a_1 + 0.5 * k.star * p,
                          rate = 1 + 0.5 * sum((tau / delta[1, 1]) * colsums(lambda^2 * phi)))

    tau <- matrix(cumprod(delta))

    if (k.star > 1) {

      for (h in 2:k.star) {

        delta_shape <- a_2 + 0.5 * p * (k.star - h + 1)
        delta_rate <- 1 + 0.5 * sum(((tau / delta[h, 1]) * colsums(lambda^2 * phi))[h:k.star])

        delta[h, 1] <- truncdist::rtrunc(spec = "gamma",
                                         n = 1,
                                         shape = delta_shape,
                                         rate = delta_rate,
                                         a = 1)

        tau <- matrix(cumprod(delta))

      }

    }


    # delta is a k.star x 1 matrix
    # tau is a k.star x 1 matrix

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # gamma update
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    gamma_shape <- (p + k.star + nu) / 2

    for (i in seq_len(n)) {

      gamma_temp1 <- matrix(data_working[i, ]) - mu - lambda_eta_t[ , i]

      gamma_rate <- 0.5 * (mat.mult(transpose(gamma_temp1), mat.mult(Sigma_inv, gamma_temp1)) +
                             mat.mult(transpose(matrix(eta[i , ])), matrix(eta[i , ])) + nu)

      gamma[i, 1] <- rgamma(1, shape = gamma_shape, rate = as.numeric(gamma_rate))

    }

    # gamma is a n x 1 matrix

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # alpha update
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    alpha_shape1 <- sum((1 - R) * (data_working >= trunc.point)) + 1

    alpha_shape2 <- sum((R) * (data_working >= trunc.point)) + 1

    alpha <- rbeta(1, shape1 = alpha_shape1, shape2 = alpha_shape2)

    # alpha is a scalar

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # imputation
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Z <- rep(NA, length(miss_vec))

    if (dim(miss_ind)[1] != 0) {

      impute_mean <- sweep(mat.mult(lambda, Rfast::transpose(eta)), 1, mu, "+")  # p x n

      for (point_ind in 1:dim(miss_ind)[1]) {

        impute_var <- (1 / gamma[miss_ind[point_ind, ][1]]) * (1 / sigma_inv)  # p x 1

        lower_prob <- pnorm(trunc.point,
                            mean = impute_mean[miss_ind[point_ind, ][2], miss_ind[point_ind, ][1]],
                            sd = sqrt(impute_var[miss_ind[point_ind, ][2], 1]),
                            lower = TRUE)

        higher_prob <- 1 - lower_prob

        zij_prob <- (lower_prob) / (lower_prob + alpha * higher_prob)

        zij <- rbinom(1, 1, zij_prob)

        Z[point_ind] <- zij

        if (zij == 0) {

          # impute as MAR
          data_working[miss_ind[point_ind, ][1], miss_ind[point_ind, ][2]] <-
            truncnorm::rtruncnorm(1, a = trunc.point,
                                  b = Inf,
                                  mean = impute_mean[miss_ind[point_ind, ][2], miss_ind[point_ind, ][1]],
                                  sd = sqrt(impute_var[miss_ind[point_ind, ][2], 1]))

        } else {

          # impute as MNAR
          data_working[miss_ind[point_ind, ][1], miss_ind[point_ind, ][2]] <-
            truncnorm::rtruncnorm(1, a = -Inf, b = trunc.point,
                                  mean = impute_mean[miss_ind[point_ind, ][2], miss_ind[point_ind, ][1]],
                                  sd = sqrt(impute_var[miss_ind[point_ind, ][2], 1]))

        }

      }

    }

    # ~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~
    # update storage
    # ~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~

    if (m %% thin == 0) {

      # imputed dataset
      # store_data[[m / 5]] <- data_working
      store_data[[m / thin]] <- data_working[miss_vec]

      # alpha - prob of MAR
      store_alpha[[m / thin]] <- alpha

      # sigma
      store_sigma_inv[[m / thin]] <- sigma_inv

      # lambda - factor loadings
      store_lambda[[m / thin]] <- lambda

      # eta
      store_eta[[m / thin]] <- eta

      # mu - mean
      store_mu[[m / thin]] <- mu

      # effective factors
      store_k_t[[m / thin]] <- k.star

      # phi
      store_phi[[m / thin]] <- phi

      # delta
      store_delta[[m / thin]] <- delta

      # gamma
      store_gamma[[m / thin]] <- gamma

      # Z
      store_Z[[m / thin]] <- Z

    }

  }

  rIFA_res <- list("store_data" = store_data,
                   "store_alpha" = store_alpha,
                   "store_sigma_inv" = store_sigma_inv,
                   "store_lambda" = store_lambda,
                   "store_eta" = store_eta,
                   "store_mu" = store_mu,
                   "store_k_t" = store_k_t,
                   "store_phi" = store_phi,
                   "store_delta" = store_delta,
                   "store_gamma" = store_gamma,
                   "store_Z" = store_Z,
                   "input_data" = input_data,
                   "input_missingness" = R,
                   "var_mult" = var_mult,
                   "inputs" = inputs
  )

  formatted_res <- rIFA_res_format(rIFA_res, burn = burn, thin = thin,
                                   cred_level = cred_level)

  # exponentiate if data were logged
  if (log.data) {

    formatted_res$imputed_dataset <- exp(formatted_res$imputed_dataset)
    formatted_res$imputation_info$imputed_val <- exp(formatted_res$imputation_info$imputed_val)
    formatted_res$imputation_info$cred_int_upper <- exp(formatted_res$imputation_info$cred_int_upper)
    formatted_res$imputation_info$cred_int_lower <- exp(formatted_res$imputation_info$cred_int_lower)

  }

  output <- list("imputed_dataset" = formatted_res$imputed_dataset,
                 "imputation_info" = formatted_res$imputation_info)

  if (return_chain) {

    output$chain <- rIFA_res

  }

  return(output)

}


# # load testing data
# # from tIFA simulation study
# load("results/p_200_datasets/gaussian_dist_dataset_p_200_1.RData")
#
#
# test_data <- generated_data$data_gen$data
# test_data_w_missing <- generated_data$data_w_missing
#
# start_time <- Sys.time()
# test <- col_rifa_gibbs(input_data = test_data_w_missing,
#                        coding = NA,
#                        n.iters = 1000,
#                        log.data = FALSE,
#                        use_true_params = FALSE,
#                        # generated_data = generated_data,
#                        k.star = 3,
#                        verbose = TRUE,
#                        burn = 500, thin = 5,
#                        return_chain = TRUE,
#                        var_mult = 0.6)
# end_time <- Sys.time()
#
# table(test$imputation_info$miss_mech)
#
# true_miss_desig <- ifelse(which(is.na(generated_data$data_w_missing)) %in% generated_data$MNAR_inds,
#                           "MNAR",
#                           "MAR")
#
# mnar_conds <- which(true_miss_desig == "MNAR")
# mar_conds <- which(true_miss_desig == "MAR")
#
# sum(test$imputation_info$miss_mech == true_miss_desig) / length(true_miss_desig)
#
# sum(test$imputation_info$miss_mech[mnar_conds] == true_miss_desig[mnar_conds]) / length(true_miss_desig[mnar_conds])
#
# sum(test$imputation_info$miss_mech[mar_conds] == true_miss_desig[mar_conds]) / length(true_miss_desig[mar_conds])


