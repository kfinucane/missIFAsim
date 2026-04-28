source("IFA_model.R")

p <- 1391
sim_data_dist <- "t"

# directory for loading and saving
directory <- paste0("results/p_", p, "_datasets")

sim_num <- 1

chain <- 1

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

set.seed(generated_data$seeding + sim_num + chain)

IFA_res <- ifa_gibbs(input_data = generated_data$use_data$data_w_missing,
                     coding = NA,
                     n.iters = 10000,
                     log.data = FALSE,
                     use_true_params = FALSE,
                     k.star = 5,
                     verbose = TRUE,
                     burn = 5000, thin = 5,
                     return_chain = TRUE,
                     var_mult = 0.6)

save(IFA_res, file = paste0(directory, "/IFA_MCMC_", sim_data_dist, "_p_", p, "_", sim_num, "_", chain, "_res.RData"))
