# source necessary scripts
source("tGIFA_model.R")
source("tGIFA_utils.R")

# load necessary packages
library(softImpute)
library(Rfast)

sim_num <- 1 # which simulated dataset to run

p <- 1391  # number of variables
sim_data_dist <- "gaussian"  # distribution desired

# proposal distribution covariance multipliers
chosen_prop_mult <- list("mu_prop_mult" = 0.05,
                         "eta_prop_mult" = 3L,
                         "lambda_prop_mult" = 3L,
                         "sigma_inv_prop_mult" = 1.5)


# directory for loading and saving
directory <- paste0("results/p_", p, "_datasets")

# load simulated data
load(paste0(directory, "/", sim_data_dist, "_dist_dataset_p_", p, "_", sim_num, ".RData"))

# remove variables that have > 25% missingness
# missingness by variable
missing_by_var <- apply(generated_data$data_w_missing, 2, function(x) sum(is.na(x)))
prop_missing_by_var <- missing_by_var / nrow(generated_data$data_w_missing)

if (length(which(prop_missing_by_var > 0.25)) > 0) {
  generated_data$use_data$data <- generated_data$data_gen$data[ , -which(prop_missing_by_var > 0.25 )]
  generated_data$use_data$data_w_missing <- generated_data$data_w_missing[ , -which(prop_missing_by_var > 0.25 )]
  generated_data$use_data$eta <- generated_data$data_gen$eta
  generated_data$use_data$lambda <- generated_data$data_gen$lambda[-which(prop_missing_by_var > 0.25 ), ]
  generated_data$use_data$Sigma <- generated_data$data_gen$Sigma[-which(prop_missing_by_var > 0.25 ), -which(prop_missing_by_var > 0.25 )]
  generated_data$use_data$mu_varphi <- generated_data$data_gen$mu_varphi[-which(prop_missing_by_var > 0.25 )]
  generated_data$use_data$mu_tilde <- generated_data$data_gen$mu_tilde[-which(prop_missing_by_var > 0.25 )]
  generated_data$use_data$mu <- matrix(generated_data$data_gen$mu[-which(prop_missing_by_var > 0.25 ), ])
} else {
  generated_data$use_data$data <- generated_data$data_gen$data
  generated_data$use_data$data_w_missing <- generated_data$data_w_missing
  generated_data$use_data$eta <- generated_data$data_gen$eta
  generated_data$use_data$lambda <- generated_data$data_gen$lambda
  generated_data$use_data$Sigma <- generated_data$data_gen$Sigma
  generated_data$use_data$mu_varphi <- generated_data$data_gen$mu_varphi
  generated_data$use_data$mu_tilde <- generated_data$data_gen$mu_tilde
  generated_data$use_data$mu <- generated_data$data_gen$mu
}

# seed
set.seed(generated_data$seeding + sim_num)

# run TGIFA model
TGIFA_res <- TGIFA_model(input_data = generated_data$use_data$data_w_missing, n.iters = 10, k.star = 5,
                         burn = 0, thin = 5, return_chain = TRUE, prop_mult = chosen_prop_mult,
                         use_true_params = FALSE, generated_data = generated_data,
                         var_mult = 0.6,
                         init_imp = "half_min")

# file to save
filename <- paste0(directory, "/TGIFA_MCMC_", sim_data_dist, "_p_", p, "_", sim_num, "_1_res.RData")

# save
save(TGIFA_res, file = filename)

