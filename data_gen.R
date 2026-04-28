# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# data generation function
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

data_sim <- function(dist, n, p, k = 5, sample_mean_in,
                     lambda, eta, broc_sized, var_mult){

  # function for seed for rcpp functions
  get_seed <- function() {

    sample.int(.Machine$integer.max, 1)

  }


  # define parameters for hyperpriors
  a_sigma <- 1L
  b_sigma <- 0.25

  if (dist == "t") {

    nu <- 3

    print("getting gamma")

    gamma <- matrix(rgamma(n = n, shape = nu/2, rate = nu/2))  # n x 1

  }

  print("getting eta")

  eta <- matrix(NA, nrow = n, ncol = k)

  if (dist == "gaussian") {

    for (i in 1:n) {

      eta[i, ] <- Rfast::rmvnorm(n = 1,
                                 mu = rep(0, k),
                                 sigma = diag(k),
                                 seed = get_seed())

    }  # n x k

  } else if (dist == "t") {

    for (i in 1:n) {

      eta[i, ] <- Rfast::rmvnorm(n = 1,
                                 mu = rep(0, k),
                                 sigma = (1 / gamma[i, 1]) * diag(k),
                                 seed = get_seed())

    }  # n x k

  }

  print("starting mu")

  mu_tilde <- as.vector(sample_mean_in) - apply(Rfast::mat.mult(lambda, t(eta)), 1, mean)  # p x 1
  mu_varphi <- 1 / (apply(broc_sized, 2, mean, na.rm = TRUE) * 0.05)  # vector of length p

  mu <- Rfast::rmvnorm(n = 1, mu = mu_tilde,
                       sigma = (1/mu_varphi) * diag(p),
                       seed = get_seed())

  mu <- matrix(t(mu))

  print("getting Sigma")

  Sigma <- diag(apply(broc_sized, 2, var, na.rm = TRUE) * var_mult)

  print("starting data")

  # data simulation
  y <- matrix(NA, nrow = n, ncol = p)

  if (dist == "gaussian") {

    for (i in 1:n) {

      y[i, ] <- Rfast::rmvnorm(n = 1,
                               mu = as.vector(mu) + (lambda %*% matrix(eta[i ,])),
                               sigma = Sigma,
                               seed = get_seed())

    }

  } else if (dist == "t") {

    for (i in 1:n) {

      y[i, ] <- Rfast::rmvnorm(n = 1,
                               mu = as.vector(mu) + (lambda %*% matrix(eta[i ,])),
                               sigma = (1 / gamma[i, 1]) * Sigma,
                               seed = get_seed())
    }

  }


  res <- list("data" = y,
              "eta" = eta,
              "lambda" = lambda,
              "Sigma" = Sigma,
              "mu_varphi" = mu_varphi,
              "mu_tilde" = mu_tilde,
              "mu" = mu)

  if (dist == "t") {

    res$nu = nu
    res$gamma = gamma

  }

  return(res)

}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# broc details function
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

faux_broc_sim <- function(dist_assumption, desired_n, desired_p, selected_factor_num, seeding,
                          prop_MNAR, prop_MAR, var_mult, LOD_quantile = 0) {

  inputs <- list("dist_assumption" = dist_assumption,
                 "desired_n" = desired_n,
                 "desired_p" = desired_p,
                 "selected_factor_num" = selected_factor_num,
                 "seeding" = seeding,
                 "prop_MNAR" = prop_MNAR,
                 "prop_MAR" = prop_MAR,
                 "var_mult" = var_mult,
                 "LOD_quantile" = LOD_quantile)

  # prop_MNAR is the proportion of MNAR columns desired
  # prop_MAR is the proportion of MAR VALUES desired
  # var_mult is the proportion of variance of the data
  # to be alloted to idiosyncratic errors

  set.seed(seeding)

  # load broc data at timepoint 0
  broc_t0 <- read.table("data/broc_example.csv", sep = ",")
  broc_t0 <- as.matrix(broc_t0)

  # missingness by variable
  missing_by_var <- apply(broc_t0, 2, function(x) sum(is.na(x)))
  prop_missing_by_var <- missing_by_var / nrow(broc_t0)

  # filter to only variables with no missingness
  false_full_broc_t0 <- as.matrix(broc_t0[ , which(prop_missing_by_var == 0)])

  # extract data mean for data gen
  sample_cols <- sort(sample(1:(dim(false_full_broc_t0)[2]), size = desired_p,
                             replace = FALSE))

  sample_mean_in <- apply(false_full_broc_t0, 2, mean, na.rm = TRUE)[sample_cols]

  broc_in <- false_full_broc_t0[ , sample_cols]

  sample_data_pca_res <- prcomp(false_full_broc_t0[ , sample_cols],
                                center = FALSE,
                                scale. = FALSE)

  # extract lambda from PCA for data gen
  lambda_in <- sample_data_pca_res$rotation[ , 1:selected_factor_num]  # p x k
  rownames(lambda_in) <- NULL
  colnames(lambda_in) <- NULL

  eta_in <- sample_data_pca_res$x[ , 1:selected_factor_num]

  # generate data based off 2 factors
  # judged from prop of variance plot
  data_gen <- data_sim(dist = dist_assumption,
                       n = desired_n,
                       p = desired_p,
                       k = selected_factor_num,
                       sample_mean_in = sample_mean_in,
                       lambda = lambda_in,
                       eta = eta_in,
                       broc_sized = broc_in,
                       var_mult = var_mult)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # for purposes of missingness allocation testing
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if (LOD_quantile > 0) {

    LOD <- quantile(data_gen$data, probs = LOD_quantile)

  } else {

    LOD <- min(data_gen$data) - 1  # 0

  }


  if (prop_MNAR == 0 & prop_MAR == 0 & LOD_quantile == 0) {

    MNAR_inds <- NA
    ind_MAR <- NA
    missingness_applied <- NA

  } else {

    missingness_applied <- data_gen$data

    if (prop_MNAR != 0) {

      # set proportion of MNAR columns desired

      # get number of columns MNAR to be applied to
      n_MNAR_cols <- ceiling(prop_MNAR * ncol(data_gen$data))

      # generate number of MNAR entries in these columns
      n_MNAR_per_col <- ceiling(runif(n = n_MNAR_cols,
                                      min = 0.1, max = 0.2) * dim(data_gen$data)[1])

      # order cols by column mins ascending
      # so essentially choosing the columns with the minimum values
      # and applying MNAR to these
      MNAR_cols <- order(apply(data_gen$data, 2, min, na.rm = TRUE), decreasing = FALSE)[1:n_MNAR_cols]

      # remove the allocated % of points from each col when ordered
      for (col in 1:length(MNAR_cols)) {

        ind <- order(missingness_applied[ , MNAR_cols[col]])[1:(n_MNAR_per_col[col])]

        missingness_applied[ind , MNAR_cols[col]] <- NA

      }

      # get indices of MNAR points for reference
      MNAR_inds <- which(is.na(missingness_applied))

    } else {

      MNAR_inds <- NA

    }

    if (LOD_quantile > 0) {

      # apply desired LOD
      missingness_applied[missingness_applied < LOD] <- NA

      # get indices of MNAR points for reference
      MNAR_inds <- which(is.na(missingness_applied))

    }


    if (prop_MAR != 0) {

      # set proportion of MAR ENTRIES desired
      # note this is not columns, separate to MNAR prop selection

      # convert prop to number of MAR points
      n_MAR <- prop_MAR * prod(dim(missingness_applied))

      # apply MAR to remaining points
      if (!is.na(MNAR_inds[1])) {

        ind_MAR <- sample((1:length(missingness_applied))[-(MNAR_inds)], n_MAR, replace = FALSE)

      } else {
        ind_MAR <- sample((1:length(missingness_applied)), n_MAR, replace = FALSE)

      }

      # remove allocated MAR points from dataset
      missingness_applied[ind_MAR] <- NA

    } else {

      ind_MAR <- NA

    }

  }

  generated_data <- list(data_gen = data_gen,
                         MNAR_inds = MNAR_inds,
                         MAR_inds = ind_MAR,
                         seeding = seeding,
                         data_w_missing = missingness_applied,
                         seeding = seeding,
                         prop_MNAR_cols = prop_MNAR,
                         prop_MAR_entries = prop_MAR,
                         prop_LOD = desired_LOD_quantile,
                         LOD = LOD,
                         inputs = inputs)

  return(generated_data)

}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# settings for simulation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dist_assumption <- "t"  # "gaussian" or "t"
n_samp <- 1  # how many datasets to generate
desired_n <- 18  # desired number of observations
desired_p <- 1391 # desired number of variables
directory <- paste0("./results/p_", desired_p, "_datasets")  # directory to save
desired_k_star <- 5  # desired k.star
desired_MNAR_prop <- 0.0  # VESTIGAL SO KEEP THIS AT 0
desired_MAR_prop <- 0.015 # prop of points above LOD to be MAR
desired_var_mult <- 0.6  # prop of var not explained by factors
desired_LOD_quantile <- 0.015  # if want no MNAR missingness set this to zero

# generate seeds
seeds <- ceiling(1000 * runif(n_samp))

for (i in 1:length(seeds)) {

  print(i)

  generated_data <- faux_broc_sim(dist_assumption = dist_assumption, desired_n = desired_n,
                                  desired_p = desired_p, selected_factor_num = desired_k_star,
                                  seeding = seeds[i], prop_MNAR = desired_MNAR_prop,
                                  prop_MAR = desired_MAR_prop,
                                  var_mult = desired_var_mult,
                                  LOD_quantile = desired_LOD_quantile)

  filename <- paste0(directory, "/", dist_assumption, "_dist_dataset_p_", desired_p, "_", i, ".RData")

  save(generated_data, file = filename)

}