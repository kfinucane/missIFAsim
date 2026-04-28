

VI_TGIFA_model <- function(input_data, coding = NA, n.iters_min,
                        n.iters_max, k.star = 5, tolerance = 1,
                        imp.iters_min = 2, imp.iters_max = 10,
                        imp_pc_tol = 0.1, n_imp = 10,
                        verbose = TRUE,
                        kappa_1 = 3L, kappa_2 = 2L,
                        a_sigma = 1L, b_sigma = 0.25, a_1 = 2.1, a_2 = 3.1,
                        return_multiple = FALSE, return_params = TRUE,
                        use_true_params = FALSE, generated_data = NULL,
                        init_imp = "half_min") {

  source("VI-TGIFA_utils.R")
  source("cavi_TGIFA_update_functions.R")

  # record inputs to include in output
  inputs <- list("input_data" = input_data,
                 "coding" = coding,
                 "n.iters_min" = n.iters_min,
                 "n.iters_max" = n.iters_max,
                 "k.star" = k.star,
                 "tolerance" = tolerance,
                 "imp.iters_min" = imp.iters_min,
                 "imp.iters_max" = imp.iters_max,
                 "imp_pc_tol" = imp_pc_tol,
                 "n_imp" = n_imp,
                 "verbose" = verbose,
                 "kappa_1" = kappa_1,
                 "kappa_2" = kappa_2,
                 "a_sigma" = a_sigma,
                 "b_sigma" = b_sigma,
                 "a_1" = a_1,
                 "a_2" = a_2,
                 "return_multiple" = return_multiple,
                 "return_params" = return_params)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # formatting data as needed
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  data <- TGIFA_prep(input_data, coding = NA)$data
  R <- TGIFA_prep(input_data, coding = NA)$missingness

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

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # initialisation
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if (verbose) {
    print("Setting initial values")
  }

  if (n_missing != 0) {

    if (use_true_params) {

      # have initial imputation be true values
      data_working <- data
      data_working[is.na(data_working)] <- generated_data$data_gen$data[is.na(data_working)]

    } else {

      if (init_imp == "half_min") {

        # initialise imputation with half-min
        initial_dataset <- half_min_imp(data, miss_ind)
        data_working <- initial_dataset

      } else if (init_imp == "SVD") {

        # initialise imputation with svd imputation
        start_svd <- softImpute::softImpute(data, rank.max = k.star, type = "svd")

        start_imp <- start_svd$u %*% diag(start_svd$d) %*% t(start_svd$v)

        initial_dataset <- data
        initial_dataset[is.na(initial_dataset)] <- start_imp[is.na(initial_dataset)]
        data_working <- initial_dataset

      }

    }

  } else {

    # use full data
    data_working <- data

  }

  if (use_true_params) {

    mu_tilde <- matrix(generated_data$data_gen$mu_tilde)
    mu_varphi <- generated_data$data_gen$mu_varphi

  } else {

    if (n_missing == 0) {

      mu_tilde <- matrix(apply(data, 2, mean, na.rm = TRUE))
      mu_varphi <- rep(1L, p)

    }

    # initialise mu hyperparams
    mu_tilde <- matrix(apply(data, 2, mean, na.rm = TRUE)) * 0.9  # multiplied downwards
    mu_varphi <- 1 / (0.05 * apply(data, 2, mean, na.rm = TRUE))
    mu_varphi[(1:p)[-unique(miss_ind[ , 2])]] <- 1L

  }

  # fix params that do not vary per-iteration
  shape_delta <- matrix(NA, nrow = k.star, ncol = 1)  # k x 1
  shape_delta[1, 1] <- a_1 + ( p * k.star ) / 2
  shape_delta[-1, 1] <- a_2 + ( p * ( k.star - (2:k.star) + 1) ) / 2

  shape_phi <- matrix(kappa_1 + 0.5, nrow = p, ncol = k.star)  # p x k

  shape_sigma_inv <- matrix(a_sigma + 0.5 * n, nrow = p, ncol = 1)  # p x 1

  # initialise remaining parameters
  mean_mu <- mu_tilde  # p x 1

  var_mu <- (1 / mu_varphi) * diag(p)  # p x p

  rate_phi <- matrix(kappa_2, nrow = p, ncol = k.star)  # p x k

  rate_delta <- matrix(1, nrow = k.star, ncol = 1)  # k x 1

  a_shape_alpha <- sum((1 - R) * (data_working >= trunc.point)) + 1  # scalar

  b_shape_alpha <- sum((R) * (data_working >= trunc.point)) + 1  # scalar

  # start at CAVI res for Gaussian FA initialisation where possible
  # using capture.output to suppress internal messages that print
  capture.output(initial_params <- VIMSFA::cavi_fa(data_working, J = k.star,
                                                   scale = FALSE, min_iter = 5,
                                                   max_iter = 5, verbose = 0),
                file = nullfile())

  rate_sigma_inv <- matrix(initial_params$rate_psi, nrow = p, ncol = 1)  # p x 1

  mean_lambda <- initial_params$mean_lambda  # matrix p x k with k-length mean_lambda_p in row j of p

  var_lambda <- initial_params$var_lambda  # list of length p of k x k matrices

  mean_eta <- t(initial_params$mean_l)  # matrix n x k with k-length mean_eta_i in row i of n

  var_eta <- lapply(seq_len(n), function(i) initial_params$var_l)  # list of length n of k x k matrices

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # VI iterations
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # iteration counter for imputation
  m_imp <- 0

  starting_Z <- rep(0, n_missing)
  starting_Z[which(data_working[miss_vec] < trunc.point)] <- 1

  while (TRUE) {

    m_imp <- m_imp + 1

    if (verbose) {

      print(paste0("Imputation iteration: ", m_imp))

    }

    theta <- list("var_eta" = var_eta,
                  "mean_eta" = mean_eta,
                  "var_lambda" = var_lambda,
                  "mean_lambda" = mean_lambda,
                  "var_mu" = var_mu,
                  "mean_mu" = mean_mu,
                  "shape_sigma_inv" = shape_sigma_inv,
                  "rate_sigma_inv" = rate_sigma_inv,
                  "a_shape_alpha" = a_shape_alpha,
                  "b_shape_alpha" = b_shape_alpha,
                  "shape_phi" = shape_phi,
                  "rate_phi" = rate_phi,
                  "shape_delta" = shape_delta,
                  "rate_delta" = rate_delta)

    # iteration counter for CAVI
    m = 0

    # VI interations
    while (TRUE) {

      m <- m + 1

      if (verbose) {

        print(paste("CAVI iteration", m, "starting"))

      }


      # store initial values of this iteration for comparison later
      theta_init <- theta

      theta <- cavi_updater(theta_init, data_working, n, p, k.star, R, trunc.point,
                            mu_tilde = mu_tilde, mu_varphi = mu_varphi,
                            kappa_1 = kappa_1, kappa_2 = kappa_2,
                            a_sigma = a_sigma, b_sigma = b_sigma,
                            a_1 = a_1, a_2 = a_2)

      # assess convergence

      # get sum of squared differences from prev iteration
      # for lambda and sigma_inv
      mean_lambda_change <- sum_sq_dif_calc(theta$mean_lambda, theta_init$mean_lambda)
      var_lambda_change <- sum_sq_dif_calc(unlist(theta$var_lambda), unlist(theta_init$var_lambda))
      rate_sigma_inv_change <- sum_sq_dif_calc(theta$rate_sigma_inv, theta_init$rate_sigma_inv)

      # value to check against tolerance
      conv_check_val <- sqrt(mean_lambda_change + var_lambda_change + rate_sigma_inv_change) / (p * (k.star + 1))

      if (verbose) {

        print(paste0(
          "mean lambda change: ", round(mean_lambda_change),
          ", var_lambda_change: ", round(var_lambda_change),
          ", rate_sigma_inv change: ", round(rate_sigma_inv_change)
        ))

        print(paste0("convergence checking value: ", conv_check_val))

      }

      if (m >= n.iters_min) {

        # check if convergence or stopping criteria reached
        if (conv_check_val < tolerance) {
          break
        } else if (m > n.iters_max){
          break
        }

      }

    }  # end CAVI while

    if (n_missing > 0){

      # generate datasets
      indiv_imputer_res <- lapply(1:n_imp, function(x) imputer(theta, miss_ind, miss_vec, n, p, k.star, trunc.point))

      # extract Z per dataset
      # matrix n_missing x n_imp
      indiv_Z <- sapply(1:n_imp, function(x) indiv_imputer_res[[x]]$Z)

      # get modal Z designation
      # vector length n_missing
      modal_Z <- apply(indiv_Z, 1, Mode)

      # which values per imputation are modal Z
      # matrix n_missing x n_imp
      vals_to_use <- apply(indiv_Z, 2, function(x) x == modal_Z)

      # calculate uncertainty in the Z designation
      # vector length n_missing
      z_uncertainty <- 1 - (apply(vals_to_use, 1, sum) / n_imp)

      # extract individual imputed values
      # matrix of n_missing x n_imp size
      indiv_imputations <- sapply(1:n_imp, function(x) indiv_imputer_res[[x]]$imputed_values)

      # get mean of imputed values
      # using those that match the modal Z
      # vector of length n_missing
      modal_imp <- sapply(1:n_missing, function(x) mean(indiv_imputations[x , vals_to_use[x, ]]))

      # get quantile endpoints
      # matrix 2 x n_missing
      modal_quantiles <- sapply(1:n_missing, function(x) quantile(indiv_imputations[x , vals_to_use[x, ]], probs = c(0.025, 0.975)))

    }

    hamming_dist <- e1071::hamming.distance(starting_Z, modal_Z)

    if (verbose) {

      print(paste0("Hamming distance over number of missing entries is: ",
                   hamming_dist / n_missing))

    }


    # change this to solution with final decision on how many imps
    # ie which indiv_imputations to check

    if (m_imp >= imp.iters_min) {

      # check if convergence or stopping criteria reached
      if ( (hamming_dist / n_missing) < imp_pc_tol) {
        break
      } else if (m_imp == imp.iters_max){
        break
      }

    }

    # change indexing if needed
    starting_Z <- modal_Z

    # starting imputation for next iteration
    # modal draws

    data_working[miss_vec] <- modal_imp


  }  # end imp while

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # format for output
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # format data

  imputation_info <- data.frame(entry_row = miss_ind[ , 1],
                                entry_col = miss_ind[ , 2],
                                imputed_val = modal_imp,
                                cred_int_lower = modal_quantiles[1, ],
                                cred_int_upper = modal_quantiles[2, ],
                                miss_mech = ifelse(modal_Z == 1, "MNAR", "MAR"),
                                miss_mech_unc = z_uncertainty
                                )

  output <- list("imputed_dataset" = data_working,
                 "imputation_info" = imputation_info,
                 "inputs" = inputs)

  if (return_params) {

    output$params <- theta

  }

  # return full datasets for the n_imp imputations performed
  # ie return n_imp datasets with different imputated values
  if (return_multiple) {

    all_datasets <- lapply(1:n_imp, function(x) {
      data_working[miss_vec] <- indiv_imputations[ , x]
      return(data_working)
    })

    output$all_datasets <- all_datasets

  }

  return(output)

}

# directory <- "./results/testing_imputation"
#
# # load 18 x 400 example data
# load(paste0(directory, "/TGIFA_testing1.RData"))
#
# cavi_res <- V_TGIFA_model(input_data = generated_data$data_w_missing, k.star = 3,
#                         n.iters_min = 5, n.iters_max = 10, tolerance = 1,
#                         imp.iters_min = 2, imp.iters_max = 3, imp_pc_tol = 0.025,
#                         return_multiple = TRUE, return_params = TRUE,
#                         use_true_params = FALSE, generated_data = generated_data)
#
# cavi_res$imputation_info
#
# plot(cavi_res$final_imp, generated_data$data_gen$data[which(is.na(generated_data$data_w_missing))])
# abline(a = 0, b = 1, col = "red", lwd = 2)
#
# # check recovery
# check <- rep(NA, length(cavi_res$final_imp))
# truth <- generated_data$data_gen$data[which(is.na(generated_data$data_w_missing))]
#
# for (imp in 1:length(cavi_res$final_imp)){
#
#   check[imp] <- (truth[imp] > cavi_res$final_quantiles[1, imp]) & (truth[imp] < cavi_res$final_quantiles[2, imp])
#
# }
#
#
# sum(check) / 108
#
#
# save(cavi_res, file = paste0(directory, "/test_res.RData"))
