library(medoutcon)
library(stringr)

### g-computation
estimate_manual <- function(data, ymodel, mmodel, l_a_model) {
  tempdat <- data

  g_m_model <- glm(formula = mmodel, family = "binomial", data = tempdat)

  # Exposition à A1
  data_a1 <- tempdat
  data_a1$a <- 1
  g_m_1 <- predict(g_m_model, newdata = data_a1, type = "response")

  # Exposition à A0
  data_a0 <- tempdat
  data_a0$a <- 0
  g_m_0 <- predict(g_m_model, newdata = data_a0, type = "response")

  # Q function
  q_model <- glm(ymodel, family = "binomial", data = tempdat)

  # do(M = 0) ou do(M = 1)
  data_m0 <- data_m1 <- tempdat
  data_m1$m <- 1
  data_m0$m <- 0

  q_pred_m1 <- predict(q_model, newdata = data_m1, type = "response")
  q_pred_m0 <- predict(q_model, newdata = data_m0, type = "response")

  tempdat$q_gamma_a1 <- q_pred_m1 * g_m_1 + q_pred_m0 * (1 - g_m_1)
  tempdat$q_gamma_a0 <- q_pred_m1 * g_m_0 + q_pred_m0 * (1 - g_m_0)

  q_model_a1 <- glm(
    formula = paste("q_gamma_a1", l_a_model, sep = "~"),
    family = "quasibinomial", data = tempdat
  )
  q_model_a0 <- glm(
    formula = paste("q_gamma_a0", l_a_model, sep = "~"),
    family = "quasibinomial", data = tempdat
  )

  q_pred_a1_gamma_a1 <- predict(q_model_a1, newdata = data_a1, type = "response")
  q_pred_a1_gamma_a0 <- predict(q_model_a0, newdata = data_a1, type = "response")
  q_pred_a0_gamma_a0 <- predict(q_model_a0, newdata = data_a0, type = "response")

  psi_mrnde <- mean(q_pred_a1_gamma_a0) - mean(q_pred_a0_gamma_a0)
  psi_mrnie <- mean(q_pred_a1_gamma_a1) - mean(q_pred_a1_gamma_a0)

  return(c(
    psi_mrnde,
    psi_mrnie
  ))
}


### IPTW
iptw_direct_indirect <- function(data) {
  g_a <- glm(a ~ 1, family = "binomial", data = data)
  g_a_l0 <- glm(a ~ w1 + w2 + w3, family = "binomial", data = data)
  g_m_a <- glm(m ~ a, family = "binomial", data = data)
  g_m_l <- glm(m ~ w1 + w2 + w3 + a + z, family = "binomial", data = data)

  pred_g1_a <- predict(g_a, type = "response")
  pred_g0_a <- 1 - pred_g1_a

  pred_g1_a_l0 <- predict(g_a_l0, type = "response")
  pred_g0_a_l0 <- 1 - pred_g1_a_l0

  pred_g1_m_a <- predict(g_m_a, type = "response")
  pred_g0_m_a <- 1 - pred_g1_m_a

  pred_g1_m_l <- predict(g_m_l, type = "response")
  pred_g0_m_l <- 1 - pred_g1_m_l

  g_a <- gm_a <- g_a_l <- gm_al <- rep(NA, nrow(data))
  g_a[data$a == 1] <- pred_g1_a[data$a == 1]
  g_a[data$a == 0] <- pred_g0_a[data$a == 0]
  g_a_l[data$a == 1] <- pred_g1_a_l0[data$a == 1]
  g_a_l[data$a == 0] <- pred_g0_a_l0[data$a == 0]
  gm_a[data$m == 1] <- pred_g1_m_a[data$m == 1]
  gm_a[data$m == 0] <- pred_g0_m_a[data$m == 0]
  gm_al[data$m == 1] <- pred_g1_m_l[data$m == 1]
  gm_al[data$m == 0] <- pred_g0_m_l[data$m == 0]

  sw <- (g_a * gm_a) / (g_a_l * gm_al)

  msm_y <- glm(y ~ a + m + a:m,
    family = "gaussian",
    data = data, weights = sw
  )
  msm_m <- glm(m ~ a + w1 + w2 + w3, family = "binomial", data = data)

  iptw_edn_l0 <- msm_y$coefficients["a"] +
    (msm_y$coefficients["a:m"] *
      plogis(rep(msm_m$coefficients["(Intercept)"], nrow(data)) +
        msm_m$coefficients["w1"] * data$w1 +
        msm_m$coefficients["w2"] * data$w2 +
        msm_m$coefficients["w3"] * data$w3))

  iptw_ein_l0 <- (msm_y$coefficients["m"] + msm_y$coefficients["a:m"]) *
    (plogis(rep(msm_m$coefficients["(Intercept)"], nrow(data)) +
      msm_m$coefficients["a"] +
      msm_m$coefficients["w1"] * data$w1 +
      msm_m$coefficients["w2"] * data$w2 +
      msm_m$coefficients["w3"] * data$w3) -
      plogis(rep(msm_m$coefficients["(Intercept)"], nrow(data)) +
        msm_m$coefficients["w1"] * data$w1 +
        msm_m$coefficients["w2"] * data$w2 +
        msm_m$coefficients["w3"] * data$w3))

  iptw_edn <- mean(iptw_edn_l0)
  iptw_ein <- mean(iptw_ein_l0)

  return(c(
    iptw_edn,
    iptw_ein
  ))
}


#### TMLE
tmle_direct_indirect <- function(data, w_names, m_names) {
  w_names <- str_subset(colnames(data), "w")
  m_names <- str_subset(colnames(data), "m")

  dir_tmle <- medoutcon(
    W = data[, w_names],
    A = data$a,
    Z = data$z,
    M = data[, m_names],
    Y = data$y,
    effect = "direct",
    u_learners = sl3::Lrnr_glm_fast$new(),
    v_learners = sl3::Lrnr_glm_fast$new(),
    estimator = "tmle"
  )

  ind_tmle <- medoutcon(
    W = data[, w_names],
    A = data$a,
    Z = data$z,
    M = data[, m_names],
    Y = data$y,
    effect = "indirect",
    u_learners = sl3::Lrnr_glm_fast$new(),
    v_learners = sl3::Lrnr_glm_fast$new(),
    estimator = "tmle"
  )

  return(c(
    dir_tmle$theta,
    ind_tmle$theta
  ))
}


#### one-step
onestep_direct_indirect <- function(data) {
  w_names <- str_subset(colnames(data), "w")
  m_names <- str_subset(colnames(data), "m")

  dir_onestep <- medoutcon(
    W = data[, w_names],
    A = data$a,
    Z = data$z,
    M = data[, m_names],
    Y = data$y,
    effect = "direct",
    u_learners = sl3::Lrnr_glm_fast$new(),
    v_learners = sl3::Lrnr_glm_fast$new(),
    estimator = "onestep"
  )

  ind_onestep <- medoutcon(
    W = data[, w_names],
    A = data$a,
    Z = data$z,
    M = data[, m_names],
    Y = data$y,
    effect = "indirect",
    u_learners = sl3::Lrnr_glm_fast$new(),
    v_learners = sl3::Lrnr_glm_fast$new(),
    estimator = "onestep"
  )

  return(c(
    dir_onestep$theta,
    ind_onestep$theta
  ))
}


# Constants
true_sde <- 0.064
true_sie <- 0.0112

ymodel <- "y ~ w1 + w2 + w3 + a + z + m"
mmodel <- "m ~ w1 + w2 + w3 + a"
l_a_model <- "w1 + w2 + w3 + a"

file_path <- "../Data/"
results_path <- "./Results/"
sim_path <- "quantitative_simulations/"

set.seed(42)
n_sim <- 500
idx <- sample(1:1000, n_sim)


# Get bias estimates from files
get_bias_estimates <- function(n_sim, path, FUN = NULL, ...) { # nolint
  bias_estimates <- matrix(nrow = n_sim, ncol = 2)
  colnames(bias_estimates) <- c("SDE", "SIE")

  for (i in 1:n_sim) {
    if (i %% 10 == 0 | i == 1) {
      print(paste("Simulation ", i, " of ", n_sim))
    }

    data <- read.csv(paste0(path, "data_", idx[i], ".csv", sep = ""))

    if (!str_detect(path, "rudolph")) {
      data <- subset(data, select = -y_qol)
      colnames(data) <- c("w1", "w2", "a", "z", "m", "y")

      results <- FUN(data, ...)
      bias_estimates[i, "SDE"] <- results[1] - true_sde
      bias_estimates[i, "SIE"] <- results[2] - true_sie
    } else {
      results <- FUN(data, ...)
      bias_estimates[i, "SDE"] <- results[1]
      bias_estimates[i, "SIE"] <- results[2]
    }
  }

  return(as.data.frame(bias_estimates))
}

# Plot bias estimates as boxplot
plot_bias_estimates <- function(bias_estimates) {
  par(mfrow = c(1, 2))
  # SDE
  boxplot(bias_estimates$SDE,
    main = "Bias estimates of SDE",
    ylab = "Bias",
    col = "steelblue2",
    border = "black"
  )
  # abline(h = 0, col = "black", lty = "dashed")

  # SIE
  boxplot(bias_estimates$SIE,
    main = "Bias estimates of SIE",
    ylab = "Bias",
    col = "steelblue2",
    border = "black"
  )
  # abline(h = 0, col = "black", lty = "dashed")
}


## g-comp
bias_estimates_gcomp <- get_bias_estimates(
  n_sim, paste(file_path, sim_path, sep = ""),
  FUN = estimate_manual,
  ymodel, mmodel, l_a_model
)
plot_bias_estimates(bias_estimates_gcomp)

write.csv(
  bias_estimates_gcomp,
  paste0(results_path, "comparisons/", "biases_gcomp_rud_", n_sim, ".csv", sep = ""),
  row.names = FALSE
)

## IPTW
bias_estimates_iptw <- get_bias_estimates(
  n_sim, paste(file_path, sim_path, sep = ""),
  FUN = iptw_direct_indirect
)
plot_bias_estimates(bias_estimates_iptw)

write.csv(
  bias_estimates_iptw,
  paste(results_path, "comparisons/", "biases_iptw_rud_", n_sim, ".csv", sep = ""),
  row.names = FALSE
)

## TMLE
bias_estimates_tmle <- get_bias_estimates(
  n_sim, paste(file_path, sim_path, sep = ""),
  FUN = tmle_direct_indirect
)
plot_bias_estimates(bias_estimates_tmle)

write.csv(
  bias_estimates_tmle,
  paste(results_path, "comparisons/", "biases_tmle_rud_", n_sim, ".csv", sep = ""),
  row.names = FALSE
)

## one-step
bias_estimates_onestep <- get_bias_estimates(
  n_sim, paste(file_path, sim_path, sep = ""),
  FUN = onestep_direct_indirect
)
plot_bias_estimates(bias_estimates_onestep)

write.csv(
  bias_estimates_onestep,
  paste(results_path, "comparisons/", "biases_onestep_rud_", n_sim, ".csv", sep = ""),
  row.names = FALSE
)


### Final comparison
# Plot final results
plot_comparison <- function(biases) {
  boxplot(biases,
    main = "Biases",
    ylab = "Bias",
    col = "steelblue2",
    border = "black"
  )
  abline(h = 0, col = "black", lty = "dashed")
}

results_path <- "./Results/"
n_sim <- 750

bias_estimates_gcomp <- read.csv(
  paste(results_path, "comparisons/", "biases_gcomp_quant_", n_sim, ".csv", sep = "")
)
bias_estimates_iptw <- read.csv(
  paste(results_path, "comparisons/", "biases_iptw_quant_", n_sim, ".csv", sep = "")
)
bias_estimates_tmle <- read.csv(
  paste(results_path, "comparisons/", "biases_tmle_quant_", n_sim, ".csv", sep = "")
)
bias_estimates_onestep <- read.csv(
  paste(results_path, "comparisons/", "biases_onestep_quant_", n_sim, ".csv", sep = "")
)

# SDE
sde_biases <- matrix(nrow = n_sim, ncol = 4)
colnames(sde_biases) <- c("SDE_gcomp", "SDE_IPTW", "SDE_TMLE", "SDE_onestep")

sde_biases[, 1] <- bias_estimates_gcomp$SDE
sde_biases[, 2] <- bias_estimates_iptw$SDE
sde_biases[, 3] <- bias_estimates_tmle$SDE
sde_biases[, 4] <- bias_estimates_onestep$SDE

write.csv(
  sde_biases,
  paste(results_path, "comparisons/", "biases_sde_quant_", n_sim, ".csv", sep = ""),
  row.names = FALSE
)

# SIE
sie_biases <- matrix(nrow = n_sim, ncol = 4)
colnames(sie_biases) <- c("SIE_gcomp", "SIE_IPTW", "SIE_TMLE", "SIE_onestep")

sie_biases[, 1] <- bias_estimates_gcomp$SIE
sie_biases[, 2] <- bias_estimates_iptw$SIE
sie_biases[, 3] <- bias_estimates_tmle$SIE
sie_biases[, 4] <- bias_estimates_onestep$SIE

write.csv(
  sie_biases,
  paste(results_path, "comparisons/", "biases_sie_quant_", n_sim, ".csv", sep = ""),
  row.names = FALSE
)

plot_comparison(sde_biases)
plot_comparison(sie_biases)