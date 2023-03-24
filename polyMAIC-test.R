###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###
###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###
#
# Preamble ----------------------------------------------------------------
#
# Code Purpose:   Test application of polyMAIC
#  Based on https://becarispublishing.com/doi/10.2217/cer-2021-0266
#
#  Alsop, Jonathan C., and Lawrence O. Pont. "Matching-adjusted indirect
#  comparison via a polynomial-based non-linear optimization method." Journal of
#  Comparative Effectiveness Research 11.8 (2022): 551-561.
###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###
###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###
###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###
###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###
###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###

library(dplyr)
library(nloptr)
objective_function_gen <- function(x, IPD1, continuous_vars = NULL, binary_vars = NULL) {
  # Prepare the matrix with the columns for the intercept, continuous, and binary variables
  X <- matrix(1, nrow(IPD1), 1)


  # Add columns for continuous variables if not NULL
  if (!is.null(continuous_vars)) {
    for (i in 1:length(continuous_vars)) {
      x_i <- IPD1[[continuous_vars[i]]]
      x_i_normalized <- (x_i - min(x_i)) / (max(x_i) - min(x_i))
      X <- cbind(X, x_i_normalized, x_i_normalized^2, x_i_normalized^3, x_i_normalized^4)
    }
  }

  # Add columns for binary variables if not NULL
  if (!is.null(binary_vars)) {
    for (i in 1:length(binary_vars)) {
      x_i <- IPD1[[binary_vars[i]]]
      X <- cbind(X, x_i)
    }
  }

  # Calculate the objective function value
  w <- exp(X %*% x)
  sw <- sum(w)
  ess <- sw^2 / sum(w^2)
  res <- -ess

  return(dplyr::lst(res, w))
}


constraint_r2_gen <- function(x, IPD1, continuous_vars = NULL, continuous_tolerances = NULL, binary_vars = NULL, binary_tolerances = NULL) {
  # Calculate weights using the objective_function_gen
  weights <- objective_function_gen(x, IPD1, continuous_vars, binary_vars)$w

  N <- nrow(IPD1)
  sw <- sum(weights)

  # Initialize constraints list
  constraints <- list()

  # Add constraints for continuous variables
  if (!is.null(continuous_vars) && !is.null(continuous_tolerances)) {
    for (i in 1:length(continuous_vars)) {
      x_i <- IPD1[[continuous_vars[i]]]
      wm <- sum(x_i * weights) / sw
      wsd <- sqrt(sum(weights * (x_i - wm)^2) / (sw - 1))

      r_mean <- abs(wm - continuous_tolerances[[i]]$mean) - continuous_tolerances[[i]]$mto
      r_sd <- abs(wsd - continuous_tolerances[[i]]$sd) - continuous_tolerances[[i]]$sto

      constraints <- append(constraints, list(r_mean, r_sd))
    }
  }

  # Add constraints for binary variables
  if (!is.null(binary_vars) && !is.null(binary_tolerances)) {
    for (i in 1:length(binary_vars)) {
      x_i <- IPD1[[binary_vars[i]]]
      wm <- sum(x_i * weights) / sw

      r_binary <- abs(wm - binary_tolerances[[i]]$p) - binary_tolerances[[i]]$mto
      constraints <- append(constraints, list(r_binary))
    }
  }

  # Add constraint for the sum of weights
  constraints <- append(constraints, list(abs(sw - N)))

  return(unlist(constraints))
}




############################################################################## #
############################################################################## #
#
# 2. Test generic version of function----
#
#     Section Notes
#
############################################################################## #
############################################################################## #

set.seed(1)
IPD1 <- data.frame(
  AGE = rnorm(1000, mean = 60, sd = 7),
  VAR2 = rnorm(1000, mean = 20, sd = 5),
  SEXF = rbinom(1000, 1, prob = 0.5)
)

continuous_vars <- c("AGE", "VAR2")
binary_vars <- c("SEXF")

# Define the continuous variable tolerances
continuous_tolerances <- list(
  list(mean = 55, sd = 6, mto = 0.005, sto = 0.005),
  list(mean = 18, sd = 5, mto = 0.01, sto = 0.01)
)

# Define the binary variable tolerances
binary_tolerances <- list(
  list(p = 0.6, mto = 0.01))




# Prepare initial_values
n_continuous <- length(continuous_vars)
n_binary <- length(binary_vars)
n_params <- 1 + n_continuous * 4 + n_binary
initial_values <- rep(0, n_params)

# Call the optimizer
opt_result <- nloptr(
  x0 = initial_values,
  eval_f = function(x) objective_function_gen(x, IPD1, continuous_vars, binary_vars)$res,
  lb = rep(-Inf, n_params),
  ub = rep(Inf, n_params),
  eval_g_ineq = function(x) constraint_r2_gen(x, IPD1, continuous_vars, continuous_tolerances, binary_vars, binary_tolerances),
  opts = list(
    maxeval = 70000,
    algorithm = "NLOPT_LN_COBYLA",
    "xtol_rel" = 1.0e-6
  )
)

# Check results
print(opt_result)

# Compute w using the optimized x
objective_result <- objective_function_gen(opt_result$solution, IPD1, continuous_vars = c("AGE", "VAR2"), binary_vars = c("SEXF"))
w <- objective_result$w

# Add computed weights to IPD1
IPD1$w <- w
sum(w)
# Print IPD1 dataset
print(IPD1)

# Calculate the summary statistics
summary_stats <- IPD1 %>%
  dplyr::select(AGE, VAR2, SEXF, w) %>%
  tidyr::pivot_longer(-w) %>%
  dplyr::group_by(name) %>%
  dplyr::summarise(
    point = matrixStats::weightedMean(value, w),
    sd = matrixStats::weightedSd(value, w)
  )

summary_stats
