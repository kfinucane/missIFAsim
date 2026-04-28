#' Extract information from input data in preparation for imputation.
#'
#' @description
#' Extracts information from a dataset for use in a TGIFA model.
#'
#'
#' @param data The data to be processed. If not a matrix, will be transformed in to a matrix.
#' @param scaled `TRUE`/`FALSE.` Should the input data be scaled? Defaults to `FALSE.`
#' @param coding The coding used to indicate a data entry is missing. Defaults to `NA.`
#' @param pareto_scale `TRUE`/`FALSE.` Should the input data be Pareto scaled? Defaults to `FALSE.`
#' @param log_trans `TRUE`/`FALSE.` Should the input data be log transformed? Defaults to `FALSE.`
#'
#' @return A list containing three entries.
#'
#' * data: Contains the processed data in matrix format with missing entries coded as `NA`.
#'
#' * missingness_pattern: Contains the missingness pattern of the data for use in the TGIFA model.
#'
#' * divisors: Contains a vector of the divisors used to scale the input data, if applicable. If unscaled, the value will be `NULL`.
#'
#' @importFrom stats sd
#'
#' @noRd
TGIFA_prep <- function(data, scaled = FALSE, coding = NA, pareto_scale = FALSE,
                       log_trans = FALSE) {

  # convert to matrix as required
  if (!(is.matrix(data))){
    data <- as.matrix(data)
  }

  # how missing data is encoded
  data[data == coding] <- NA

  # centre or scale as required
  if (scaled){
    data <- scale(data, center = FALSE, scale = scaled)
  }

  divisors <- NULL

  if (pareto_scale) {

    divisors <- sqrt(apply(data, 2, sd, na.rm = TRUE))

    data <- sweep(data, MARGIN = 2, STATS = divisors, FUN = "/")
  }

  if (log_trans) {
    data <- log(data)
  }

  # generate missingness pattern
  missingness_pattern <- 1 - is.na(data)

  return(list("data" = data, "missingness_pattern" = missingness_pattern, "divisors" = divisors))
}


#' Seeding function for RCPP functions
#'
#' @description
#' Seeding function for RCPP functions when `set.seed()` is not appropriate.
#'
#'
#' @noRd
get_seed <- function() {

  sample.int(.Machine$integer.max, 1)

}


# matrix inverse in Cpp
Rcpp::cppFunction("arma::mat armaInv(const arma::mat & x) { return arma::inv(x); }", depends = "RcppArmadillo")

# function to get squared difference
# between two objects of the same size
# vectors or matrices
sum_sq_dif_calc <- function(object1, object2) {

  sq_dif <- (object1 - object2)^2

  sum_sq_dif <- sum(sq_dif)

  return(sum_sq_dif)

}


#' Mode
#'
#' @param x A vector.
#'
#' @return The mode of x. If x is multi-modal, the first mode (by index in x) is returned.
#'
#' @noRd
Mode <- function(x) {

  options <- unique(x)

  return(options[which.max(tabulate(match(x, options)))])

}

# function to perform half minimum imputation
half_min_imp <- function(data, miss_ind) {

  half_mins_per_col <- apply(data, 2, function(x) min(x, na.rm = TRUE) / 2)

  half_min_imp_data <- data

  for (ind in 1:dim(miss_ind)[1]) {

    half_min_imp_data[miss_ind[ind, 1], miss_ind[ind, 2]] <- half_mins_per_col[miss_ind[ind, 2]]

  }

  return(half_min_imp_data)

}

# calculate mean absolute error
MAE_calc <- function(true_vals, est_vals) {

  abs_dif <- abs(true_vals - est_vals)

  MAE <- sum(abs_dif) / length(true_vals)

  return(MAE)

}
