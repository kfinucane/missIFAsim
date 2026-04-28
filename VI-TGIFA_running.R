# source model script
source("VI-TGIFA_model.R")

p <- 1391  # number of variables
sim_data_dist <- "gaussian"  # distribution desired

# directory for loading and saving
directory <- paste0("results/p_", p, "_datasets")

# run over all simulations (generally fast)
for (sim_num in 1:10) {

  # load 18 x 1391 data
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

  # run method
  cavi_res <- VI_TGIFA_model(input_data = generated_data$use_data$data_w_missing, k.star = 5,
                          n.iters_min = 5, n.iters_max = 50, tolerance = 0.01,
                          imp.iters_min = 3, imp.iters_max = 10, n_imp = 100,
                          imp_pc_tol = 0.05,
                          return_multiple = TRUE, return_params = TRUE,
                          use_true_params = FALSE, generated_data = generated_data)

  # record end time
  end_time <- Sys.time()

  # save results and timing
  save(cavi_res, file = paste0(directory, "/VI-TGIFA_res_", sim_data_dist, "_p_", p, "_", sim_num, "_1_res.RData"))

}
