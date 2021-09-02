## Estimation manuelle (g-comp)
estimate_manual <- function(data, sl_library) {
  tempdat <- data

  y <- tempdat$m
  x <- subset(tempdat, select = c("w1", "w2", "a"))
  g_m_model <- SuperLearner(y, x,
    family = binomial(),
    SL.library = sl_library
  )

  # Exposition à A1
  data_a1 <- tempdat
  data_a1$a <- 1
  g_m_1 <- predict(g_m_model, newdata = data_a1, type = "response")

  # Exposition à A0
  data_a0 <- tempdat
  data_a0$a <- 0
  g_m_0 <- predict(g_m_model, newdata = data_a0, type = "response")

  # Q function
  y <- tempdat$y
  x <- subset(tempdat, select = -c(y))
  q_model <- SuperLearner(y, x,
    family = binomial(),
    SL.library = sl_library
  )

  # do(M = 0) ou do(M = 1)
  data_m0 <- data_m1 <- tempdat
  data_m1$m <- 1
  data_m0$m <- 0

  q_pred_m1 <- predict(q_model, newdata = data_m1, type = "response")
  q_pred_m0 <- predict(q_model, newdata = data_m0, type = "response")

  q_preds_m1 <- q_pred_m1$pred
  q_preds_m0 <- q_pred_m0$pred

  tempdat$q_gamma_a1 <- q_preds_m1 * g_m_1$pred + q_preds_m0 * (1 - g_m_1$pred)
  tempdat$q_gamma_a0 <- q_preds_m1 * g_m_0$pred + q_preds_m0 * (1 - g_m_0$pred)

  y <- tempdat$q_gamma_a1
  x <- subset(tempdat, select = c("w1", "w2", "a"))
  q_model_a1 <- SuperLearner(y, x,
    family = "quasibinomial",
    SL.library = sl_library
  )

  y <- tempdat$q_gamma_a0
  x <- subset(tempdat, select = c("w1", "w2", "a"))
  q_model_a0 <- SuperLearner(y, x,
    family = "quasibinomial",
    SL.library = sl_library
  )

  q_pred_a1_gamma_a1 <- predict(q_model_a1, data_a1, type = "response")
  q_pred_a1_gamma_a0 <- predict(q_model_a0, data_a1, type = "response")
  q_pred_a0_gamma_a0 <- predict(q_model_a0, data_a0, type = "response")

  psi_m_rnde <- mean(q_pred_a1_gamma_a0$pred) - mean(q_pred_a0_gamma_a0$pred)
  psi_m_rnie <- mean(q_pred_a1_gamma_a1$pred) - mean(q_pred_a1_gamma_a0$pred)

  return(c(
    psi_m_rnde,
    psi_m_rnie
  ))
}


## Article Kara E. Rudolph (TMLE)
estimate_rudolph <- function(data, sl_library) {
  obsdat <- data
  tmpdat <- obsdat
  dfa1 <- dfa0 <- dfm1 <- dfm0 <- obsdat

  dfa1$a <- dfm1$m <- 1
  dfa0$a <- dfm0$m <- 0

  y <- obsdat$z
  x <- obsdat[, c("w1", "w2", "a")]
  sl_glm_zfit <- SuperLearner(y, x,
    family = binomial(),
    SL.library = sl_library
  )
  pred_z_a0 <- predict(sl_glm_zfit, dfa0, type = "response")
  pred_z_a1 <- predict(sl_glm_zfit, dfa1, type = "response")

  preds_z_a0 <- pred_z_a0$pred
  preds_z_a1 <- pred_z_a1$pred

  y <- obsdat$m
  x <- subset(tmpdat, select = -c(m, y))
  sl_glm_mfit <- SuperLearner(y, x,
    family = binomial(),
    SL.library = sl_library
  )
  pred_m_a0 <- predict(sl_glm_mfit, dfa0, type = "response")
  pred_m_a1 <- predict(sl_glm_mfit, dfa1, type = "response")

  preds_m_a0 <- pred_m_a0$pred
  preds_m_a1 <- pred_m_a1$pred

  gm_a0 <- preds_m_a0 * preds_z_a0 + preds_m_a0 * (1 - preds_z_a0)
  gm_a1 <- preds_m_a1 * preds_z_a1 + preds_m_a1 * (1 - preds_z_a1)

  y <- obsdat$a
  x <- obsdat[, c("w1", "w2")]
  sl_glm_afit <- SuperLearner(y, x,
    family = binomial(),
    SL.library = sl_library
  )
  pred_sl <- predict(sl_glm_afit, type = "response")
  pred_ps_a1 <- I(tmpdat$a == 1) / pred_sl$pred
  pred_ps_a0 <- I(tmpdat$a == 0) / (1 - pred_sl$pred)

  y <- obsdat$m
  x <- subset(tmpdat, select = -c(m, y))
  sl_glm_mfit <- SuperLearner(y, x,
    family = binomial(),
    SL.library = sl_library
  )
  sl_mazw <- predict(sl_glm_mfit, type = "response")
  sl_psm <- I(tmpdat$m == 1) * sl_mazw$pred + I(tmpdat$m == 0) * (1 - sl_mazw$pred)

  y <- obsdat$y
  x <- subset(tmpdat, select = -c(y))
  sl_glm_yfit <- SuperLearner(y, x,
    family = binomial(),
    SL.library = sl_library
  )
  tmpdat$qyinit <- cbind(
    predict(sl_glm_yfit, newdata = tmpdat, type = "response")$pred,
    predict(sl_glm_yfit, newdata = dfm0, type = "response")$pred,
    predict(sl_glm_yfit, newdata = dfm1, type = "response")$pred
  )

  # ----------------------------
  # Get E(Y_{1, g_0})
  # ----------------------------
  b <- I(tmpdat$m) == 1
  c <- I(tmpdat$m) == 0
  tmpdat$h_a1_gm_a0 <- ((b * gm_a0 + c * (1 - gm_a0)) / sl_psm) * pred_ps_a1

  epsilon_a1_gm_a0 <- coef(glm(y ~ 1,
    weights = h_a1_gm_a0, offset = qlogis(qyinit[, 1]),
    family = "quasibinomial", data = tmpdat
  ))

  tmpdat$qyup_m0_a1_gm_a0 <- plogis(qlogis(tmpdat$qyinit[, 2]) + epsilon_a1_gm_a0)
  tmpdat$qyup_m1_a1_gm_a0 <- plogis(qlogis(tmpdat$qyinit[, 3]) + epsilon_a1_gm_a0)

  tmpdat$q_a1_gm_a0 <- tmpdat$qyup_m0_a1_gm_a0 * (1 - gm_a0) + tmpdat$qyup_m1_a1_gm_a0 * gm_a0

  df <- tmpdat
  y <- df[df$a == 1, ]$q_a1_gm_a0
  x <- subset(tmpdat, select = -c(q_a1_gm_a0))
  x <- subset(x, a == 1, select = c("w1", "w2"))

  q_a1_g0_fit <- SuperLearner(y, x,
    family = "quasibinomial",
    SL.library = sl_library
  )
  q_a1_g0 <- predict(q_a1_g0_fit, newdata = tmpdat, type = "response")

  # Necessary if A is non random
  epsilon_z_a1_gm_a0 <- coef(glm(q_a1_gm_a0 ~ 1,
    weights = pred_ps_a1, offset = qlogis(q_a1_g0$pred),
    family = "quasibinomial", data = tmpdat
  ))

  q_zup_a1_gm_a0 <- plogis(qlogis(q_a1_g0$pred) + epsilon_z_a1_gm_a0)

  # Calculate TMLE
  tmle_a1_m0 <- mean(q_zup_a1_gm_a0)

  # ----------------------------
  # Get E(Y_{0, g_0})
  # ----------------------------
  b <- tmpdat$m
  tmpdat$h_a0_gm_a0 <- ((b * gm_a0 + (1 - b) * (1 - gm_a0)) / sl_psm) * pred_ps_a0

  epsilon_a0_gm_a0 <- coef(glm(y ~ 1,
    weights = h_a0_gm_a0, offset = qlogis(qyinit[, 1]),
    family = "quasibinomial", data = tmpdat
  ))

  tmpdat$qyup_m0_a0_gm_a0 <- plogis(qlogis(tmpdat$qyinit[, 2]) + epsilon_a0_gm_a0)
  tmpdat$qyup_m1_a0_gm_a0 <- plogis(qlogis(tmpdat$qyinit[, 3]) + epsilon_a0_gm_a0)

  tmpdat$q_a0_gm_a0 <- tmpdat$qyup_m0_a0_gm_a0 * (1 - gm_a0) + tmpdat$qyup_m1_a0_gm_a0 * gm_a0

  df <- tmpdat
  y <- df[df$a == 0, ]$q_a0_gm_a0
  x <- subset(tmpdat, select = -c(q_a0_gm_a0))
  x <- subset(x, a == 0, select = c("w1", "w2"))

  q_a0_g0_fit <- SuperLearner(y, x,
    family = "quasibinomial",
    SL.library = sl_library
  )
  q_a0_g0 <- predict(q_a0_g0_fit, newdata = tmpdat, type = "response")

  # Necessary if A is non random
  epsilon_z_a0_gm_a0 <- coef(glm(q_a0_gm_a0 ~ 1,
    weights = pred_ps_a0, offset = qlogis(q_a0_g0$pred),
    family = "quasibinomial", data = tmpdat
  ))

  q_zup_a0_gm_a0 <- plogis(qlogis(q_a0_g0$pred) + epsilon_z_a0_gm_a0)

  # Calculate TMLE
  tmle_a0_m0 <- mean(q_zup_a0_gm_a0)

  # ----------------------------
  # Get E(Y_{1, g_1})
  # ----------------------------
  tmpdat$h_a1_gm_a1 <- ((b * gm_a1 + (1 - b) * (1 - gm_a1)) / sl_psm) * pred_ps_a1

  epsilon_a1_gm_a1 <- coef(glm(y ~ 1,
    weights = h_a1_gm_a1, offset = qlogis(qyinit[, 1]),
    family = "quasibinomial", data = tmpdat
  ))

  tmpdat$qyup_m0_a1_gm_a1 <- plogis(qlogis(tmpdat$qyinit[, 2]) + epsilon_a1_gm_a1)
  tmpdat$qyup_m1_a1_gm_a1 <- plogis(qlogis(tmpdat$qyinit[, 3]) + epsilon_a1_gm_a1)

  tmpdat$q_a1_gm_a1 <- tmpdat$qyup_m0_a1_gm_a1 * (1 - gm_a1) + tmpdat$qyup_m1_a1_gm_a1 * gm_a1

  df <- tmpdat
  y <- df[df$a == 1, ]$q_a1_gm_a1
  x <- subset(tmpdat, select = -c(q_a1_gm_a1))
  x <- subset(x, a == 1, select = c("w1", "w2"))

  q_a1_g1_fit <- SuperLearner(y, x,
    family = "quasibinomial",
    SL.library = sl_library
  )
  q_a1_g1 <- predict(q_a1_g1_fit, newdata = tmpdat, type = "response")

  # Necessary if A is non random
  epsilon_z_a1_gm_a1 <- coef(glm(q_a1_gm_a1 ~ 1,
    weights = pred_ps_a1, offset = qlogis(q_a1_g1$pred),
    family = "quasibinomial", data = tmpdat
  ))

  q_zup_a1_gm_a1 <- plogis(qlogis(q_a1_g1$pred) + epsilon_z_a1_gm_a1)

  # Calculate TMLE
  tmle_a1_m1 <- mean(q_zup_a1_gm_a1)

  # ----------------------------
  # Estimate variances using EIC
  # ----------------------------
  tmpdat$qyup_a1_g0 <- plogis(qlogis(tmpdat$qyinit[, 1]) + epsilon_a1_gm_a0)
  tmpdat$qyup_a1_g1 <- plogis(qlogis(tmpdat$qyinit[, 1]) + epsilon_a1_gm_a1)
  tmpdat$qyup_a0_g0 <- plogis(qlogis(tmpdat$qyinit[, 1]) + epsilon_a0_gm_a0)

  # EIC for E(Y_{1, g_0})
  eic1_a1_g0 <- tmpdat$h_a1_gm_a0 * (tmpdat$y - tmpdat$qyup_a1_g0)
  eic2_a1_g0 <- pred_ps_a1 * (tmpdat$q_a1_gm_a0 - q_zup_a1_gm_a0)
  eic3_a1_g0 <- q_zup_a1_gm_a0 - tmle_a1_m0

  eic_a1_g0 <- eic1_a1_g0 + eic2_a1_g0 + eic3_a1_g0

  # EIC for E(Y_{1, g_1})
  eic1_a1_g1 <- tmpdat$h_a1_gm_a1 * (tmpdat$y - tmpdat$qyup_a1_g1)
  eic2_a1_g1 <- pred_ps_a1 * (tmpdat$q_a1_gm_a1 - q_zup_a1_gm_a1)
  eic3_a1_g1 <- q_zup_a1_gm_a1 - tmle_a1_m1

  eic_a1_g1 <- eic1_a1_g1 + eic2_a1_g1 + eic3_a1_g1

  # EIC for E(Y_{0, g_0})
  eic1_a0_g0 <- tmpdat$h_a0_gm_a0 * (tmpdat$y - tmpdat$qyup_a0_g0)
  eic2_a0_g0 <- pred_ps_a0 * (tmpdat$q_a0_gm_a0 - q_zup_a0_gm_a0)
  eic3_a0_g0 <- q_zup_a0_gm_a0 - tmle_a0_m0

  eic_a0_g0 <- eic1_a0_g0 + eic2_a0_g0 + eic3_a0_g0

  # Mediation params., variances and 95% IC
  sde_tmle <- tmle_a1_m0 - tmle_a0_m0
  sde_eic <- eic_a1_g0 - eic_a0_g0
  var_sde_eic <- var(sde_eic) / nrow(tmpdat)

  sie_tmle <- tmle_a1_m1 - tmle_a1_m0
  sie_eic <- eic_a1_g1 - eic_a1_g0
  var_sie_eic <- var(sie_eic) / nrow(tmpdat)

  return(c(
    sde_tmle,
    sie_tmle
  ))
}


### Medoutcon (TMLE)
tmle_direct_indirect <- function(data, learner) {
  w_names <- str_subset(colnames(data), "w")
  m_names <- str_subset(colnames(data), "m")

  dir_tmle <- medoutcon(
    W = data[, w_names],
    A = data$a,
    Z = data$z,
    M = data[, m_names],
    Y = data$y,
    effect = "direct",
    u_learners = learner,
    v_learners = learner,
    estimator = "tmle"
  )

  ind_tmle <- medoutcon(
    W = data[, w_names],
    A = data$a,
    Z = data$z,
    M = data[, m_names],
    Y = data$y,
    effect = "indirect",
    u_learners = learner,
    v_learners = learner,
    estimator = "tmle"
  )

  return(c(
    dir_tmle$theta,
    ind_tmle$theta
  ))
}


### Medoutcon (one-step)
onestep_direct_indirect <- function(data, learner) {
  w_names <- str_subset(colnames(data), "w")
  m_names <- str_subset(colnames(data), "m")

  dir_os <- medoutcon(
    W = data[, w_names],
    A = data$a,
    Z = data$z,
    M = data[, m_names],
    Y = data$y,
    effect = "direct",
    u_learners = learner,
    v_learners = learner,
    estimator = "onestep"
  )

  ind_os <- medoutcon(
    W = data[, w_names],
    A = data$a,
    Z = data$z,
    M = data[, m_names],
    Y = data$y,
    effect = "indirect",
    u_learners = learner,
    v_learners = learner,
    estimator = "onestep"
  )

  return(c(
    dir_os$theta,
    ind_os$theta
  ))
}


### IPTW
iptw_direct_indirect <- function(data, sl_library) {
  g_a <- glm(a ~ 1, family = "binomial", data = data)

  y <- data$a
  x <- subset(data, select = c("w1", "w2"))
  g_a_l0 <- SuperLearner(y, x,
    family = binomial(),
    SL.library = sl_library
  )

  y <- data$m
  x <- subset(data, select = a)
  g_m_a <- SuperLearner(y, x,
    family = binomial(),
    SL.library = sl_library
  )

  x <- subset(data, select = c("w1", "w2", "a", "z"))
  g_m_l <- SuperLearner(y, x,
    family = binomial(),
    SL.library = sl_library
  )

  pred_g1_a <- predict(g_a, type = "response")
  pred_g0_a <- 1 - pred_g1_a

  pred_g1_a_l0 <- predict(g_a_l0, type = "response")$pred
  pred_g0_a_l0 <- 1 - pred_g1_a_l0

  pred_g1_m_a <- predict(g_m_a, type = "response")$pred
  pred_g0_m_a <- 1 - pred_g1_m_a

  pred_g1_m_l <- predict(g_m_l, type = "response")$pred
  pred_g0_m_l <- 1 - pred_g1_m_l

  ga <- gm_a <- ga_l <- gm_al <- rep(NA, nrow(data))

  ga[data$a == 1] <- pred_g1_a[data$a == 1]
  ga[data$a == 0] <- pred_g0_a[data$a == 0]
  ga_l[data$a == 1] <- pred_g1_a_l0[data$a == 1]
  ga_l[data$a == 0] <- pred_g0_a_l0[data$a == 0]
  gm_a[data$m == 1] <- pred_g1_m_a[data$m == 1]
  gm_a[data$m == 0] <- pred_g0_m_a[data$m == 0]
  gm_al[data$m == 1] <- pred_g1_m_l[data$m == 1]
  gm_al[data$m == 0] <- pred_g0_m_l[data$m == 0]

  sw <- (ga * gm_a) / (ga_l * gm_al)

  msm_y <- glm(y ~ a + m + a:m, family = "gaussian", data = data, weights = sw)
  msm_m <- glm(m ~ a + w1 + w2, family = "binomial", data = data, weights = (ga / ga_l))

  iptw_edn_l0 <- msm_y$coefficients["a"] +
    (msm_y$coefficients["a:m"] *
      plogis(rep(msm_m$coefficients["(Intercept)"], nrow(data)) +
        msm_m$coefficients["w1"] * data$w1 +
        msm_m$coefficients["w2"] * data$w2))

  iptw_ein_l0 <- (msm_y$coefficients["m"] + msm_y$coefficients["a:m"]) *
    (plogis(rep(msm_m$coefficients["(Intercept)"], nrow(data)) +
      msm_m$coefficients["a"] +
      msm_m$coefficients["w1"] * data$w1 +
      msm_m$coefficients["w2"] * data$w2) -
      plogis(rep(msm_m$coefficients["(Intercept)"], nrow(data)) +
        msm_m$coefficients["w1"] * data$w1 +
        msm_m$coefficients["w2"] * data$w2))

  iptw_edn <- mean(iptw_edn_l0)
  iptw_ein <- mean(iptw_ein_l0)

  return(c(
    iptw_edn,
    iptw_ein
  ))
}



libs <- c(
  "SL.glm", "SL.glm.interaction", "SL.ranger",
  "SL.bayesglm", "SL.speedlm"
)

learners <- c(
  sl3::Lrnr_bayesglm$new(), sl3::Lrnr_glm_fast$new(),
  sl3::Lrnr_ranger$new(), sl3::Lrnr_nnet$new()
)
learner_names <- c("bayes_glm", "glm_fast", "ranger", "nnet")

# Constants
true_sde <- 0.064
true_sie <- 0.0112

file_path <- "../Data/"
results_path <- "./Results/"
sim_path <- "simulations_rudolph/"
n_sim <- 200


# Get bias estimates from files
get_bias_estimates <- function(n_sim, path, FUN, ...) { # nolint
  suppressPackageStartupMessages({
    library(medoutcon)
    library(SuperLearner)
    library(ranger)
    library(arm)
    library(speedglm)
    library(stringr)
    library(sl3)
  })

  set.seed(42)
  idx <- sample(1:1000, n_sim)

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

      results <- try(FUN(data, ...))
      bias_estimates[i, "SDE"] <- results[1] - true_sde
      bias_estimates[i, "SIE"] <- results[2] - true_sie
    } else {
      results <- try(FUN(data, ...))
      bias_estimates[i, "SDE"] <- results[1]
      bias_estimates[i, "SIE"] <- results[2]
    }
  }

  return(as.data.frame(bias_estimates))
}

# Plot bias estimates as boxplot
plot_bias_estimates <- function(sde_bias, sie_bias, n_sim, title, sub) {
  par(mfrow = c(1, 2))
  # SDE
  boxplot(sde_bias,
    main = paste("Bias estimates of SDE with", title, "\nand n_sim = ", n_sim),
    ylab = "Bias",
    col = "steelblue2",
    border = "black",
    sub = sub
  )
  abline(h = 0, col = "black", lty = "dashed")

  # SIE
  boxplot(sie_bias,
    main = paste("Bias estimates of SIE with", title, "\nand n_sim = ", n_sim),
    ylab = "Bias",
    col = "steelblue2",
    border = "black",
    sub = sub
  )
  abline(h = 0, col = "black", lty = "dashed")
}


## g-comp
gcomp_sl_sde_estimates <- matrix(nrow = n_sim, ncol = length(libs))
gcomp_sl_sie_estimates <- matrix(nrow = n_sim, ncol = length(libs))
colnames(gcomp_sl_sde_estimates) <- libs
colnames(gcomp_sl_sie_estimates) <- libs

write.csv(gcomp_sl_sde_estimates,
  paste(results_path, "sl_comparisons/", "gcomp_sde_conf_", n_sim, ".csv", sep = ""),
  row.names = FALSE
)
write.csv(gcomp_sl_sie_estimates,
  paste(results_path, "sl_comparisons/", "gcomp_sie_conf_", n_sim, ".csv", sep = ""),
  row.names = FALSE
)

for (i in seq_len(length(libs))) {
  cat(paste("\n######## Library ", i, " of ", length(libs), " ########\n"))

  bias_estimates_gcomp <- get_bias_estimates(
    n_sim, paste(file_path, sim_path, sep = ""),
    FUN = estimate_manual,
    libs[i]
  )

  gcomp_sl_sde_estimates <- read.csv(paste(
    results_path, "sl_comparisons/", "gcomp_sde_conf_", n_sim, ".csv",
    sep = ""
  ))
  gcomp_sl_sie_estimates <- read.csv(paste(
    results_path, "sl_comparisons/", "gcomp_sie_conf_", n_sim, ".csv",
    sep = ""
  ))

  gcomp_sl_sde_estimates[, i] <- bias_estimates_gcomp[, 1]
  gcomp_sl_sie_estimates[, i] <- bias_estimates_gcomp[, 2]

  plot_bias_estimates(
    gcomp_sl_sde_estimates, gcomp_sl_sie_estimates,
    n_sim, "g-computation", "quantitative confounder"
  )

  write.csv(
    gcomp_sl_sde_estimates,
    paste(results_path, "sl_comparisons/", "gcomp_sde_conf_", n_sim, ".csv", sep = ""),
    row.names = FALSE
  )
  write.csv(
    gcomp_sl_sie_estimates,
    paste(results_path, "sl_comparisons/", "gcomp_sie_conf_", n_sim, ".csv", sep = ""),
    row.names = FALSE
  )
}


## IPTW
iptw_sl_sde_estimates <- matrix(nrow = n_sim, ncol = length(libs))
iptw_sl_sie_estimates <- matrix(nrow = n_sim, ncol = length(libs))
colnames(iptw_sl_sde_estimates) <- libs
colnames(iptw_sl_sie_estimates) <- libs

write.csv(iptw_sl_sde_estimates,
  paste(results_path, "sl_comparisons/", "iptw_sde_conf_", n_sim, ".csv", sep = ""),
  row.names = FALSE
)
write.csv(iptw_sl_sie_estimates,
  paste(results_path, "sl_comparisons/", "iptw_sie_conf_", n_sim, ".csv", sep = ""),
  row.names = FALSE
)

for (i in seq_len(length(libs))) {
  cat(paste("\n######## Library ", i, " of ", length(libs), " ########\n"))

  bias_estimates_iptw <- get_bias_estimates(
    n_sim, paste(file_path, sim_path, sep = ""),
    FUN = iptw_direct_indirect,
    libs[i]
  )

  iptw_sl_sde_estimates <- read.csv(paste(
    results_path, "sl_comparisons/", "iptw_sde_conf_", n_sim, ".csv",
    sep = ""
  ))
  iptw_sl_sie_estimates <- read.csv(paste(
    results_path, "sl_comparisons/", "iptw_sie_conf_", n_sim, ".csv",
    sep = ""
  ))

  iptw_sl_sde_estimates[, i] <- bias_estimates_iptw[, 1]
  iptw_sl_sie_estimates[, i] <- bias_estimates_iptw[, 2]

  plot_bias_estimates(
    iptw_sl_sde_estimates, iptw_sl_sie_estimates,
    n_sim, "IPTW", "quantitative confounder"
  )

  write.csv(
    iptw_sl_sde_estimates,
    paste(results_path, "sl_comparisons/", "iptw_sde_conf_", n_sim, ".csv", sep = ""),
    row.names = FALSE
  )
  write.csv(
    iptw_sl_sie_estimates,
    paste(results_path, "sl_comparisons/", "iptw_sie_conf_", n_sim, ".csv", sep = ""),
    row.names = FALSE
  )
}


## TMLE
tmle_sl_sde_estimates <- matrix(nrow = n_sim, ncol = length(learners))
tmle_sl_sie_estimates <- matrix(nrow = n_sim, ncol = length(learners))
colnames(tmle_sl_sde_estimates) <- learner_names
colnames(tmle_sl_sie_estimates) <- learner_names

write.csv(tmle_sl_sde_estimates,
  paste(results_path, "sl_comparisons/", "tmle_sde_rud_", n_sim, ".csv", sep = ""),
  row.names = FALSE
)
write.csv(tmle_sl_sie_estimates,
  paste(results_path, "sl_comparisons/", "tmle_sie_rud_", n_sim, ".csv", sep = ""),
  row.names = FALSE
)

ind <- 1
for (learner in learners) {
  cat(paste("\n######## Library ", ind, " of ", length(learners), " ########\n"))

  bias_estimates_tmle <- get_bias_estimates(
    n_sim, paste(file_path, sim_path, sep = ""),
    FUN = tmle_direct_indirect,
    # FUN = estimate_rudolph,
    learner
  )

  tmle_sl_sde_estimates <- read.csv(paste(
    results_path, "sl_comparisons/", "tmle_sde_rud_", n_sim, ".csv",
    sep = ""
  ))
  tmle_sl_sie_estimates <- read.csv(paste(
    results_path, "sl_comparisons/", "tmle_sie_rud_", n_sim, ".csv",
    sep = ""
  ))

  tmle_sl_sde_estimates[, ind] <- bias_estimates_tmle[, 1]
  tmle_sl_sie_estimates[, ind] <- bias_estimates_tmle[, 2]

  plot_bias_estimates(
    tmle_sl_sde_estimates, tmle_sl_sie_estimates,
    n_sim, "TMLE", "modèle Rudolph"
  )

  write.csv(
    tmle_sl_sde_estimates,
    paste(results_path, "sl_comparisons/", "tmle_sde_rud_", n_sim, ".csv", sep = ""),
    row.names = FALSE
  )
  write.csv(
    tmle_sl_sie_estimates,
    paste(results_path, "sl_comparisons/", "tmle_sie_rud_", n_sim, ".csv", sep = ""),
    row.names = FALSE
  )
  ind <- ind + 1
}


## one-step
onestep_sl_sde_estimates <- matrix(nrow = n_sim, ncol = length(learners))
onestep_sl_sie_estimates <- matrix(nrow = n_sim, ncol = length(learners))
colnames(onestep_sl_sde_estimates) <- learner_names
colnames(onestep_sl_sie_estimates) <- learner_names

write.csv(onestep_sl_sde_estimates,
  paste(results_path, "sl_comparisons/", "onestep_sde_rud_", n_sim, ".csv", sep = ""),
  row.names = FALSE
)
write.csv(onestep_sl_sie_estimates,
  paste(results_path, "sl_comparisons/", "onestep_sie_rud_", n_sim, ".csv", sep = ""),
  row.names = FALSE
)

ind <- 1
for (learner in learners) {
  cat(paste("\n######## Library ", ind, " of ", length(learners), " ########\n"))

  bias_estimates_onestep <- get_bias_estimates(
    n_sim, paste(file_path, sim_path, sep = ""),
    FUN = onestep_direct_indirect,
    learner
  )

  onestep_sl_sde_estimates <- read.csv(paste(
    results_path, "sl_comparisons/", "onestep_sde_rud_", n_sim, ".csv",
    sep = ""
  ))
  onestep_sl_sie_estimates <- read.csv(paste(
    results_path, "sl_comparisons/", "onestep_sie_rud_", n_sim, ".csv",
    sep = ""
  ))

  onestep_sl_sde_estimates[, ind] <- bias_estimates_onestep[, 1]
  onestep_sl_sie_estimates[, ind] <- bias_estimates_onestep[, 2]

  plot_bias_estimates(
    onestep_sl_sde_estimates, onestep_sl_sie_estimates,
    n_sim, "one-step", "modèle Rudolph"
  )

  write.csv(
    onestep_sl_sde_estimates,
    paste(results_path, "sl_comparisons/", "onestep_sde_rud_", n_sim, ".csv", sep = ""),
    row.names = FALSE
  )
  write.csv(
    onestep_sl_sie_estimates,
    paste(results_path, "sl_comparisons/", "onestep_sie_rud_", n_sim, ".csv", sep = ""),
    row.names = FALSE
  )
  ind <- ind + 1
}


### Compile results
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
n_sim <- 100

# SDE
sde_biases <- matrix(nrow = n_sim, ncol = 4)
colnames(sde_biases) <- c("SDE_gcomp", "SDE_IPTW", "SDE_TMLE", "SDE_onestep")

sde_bias_gcomp <- read.csv(
  paste(results_path, "sl_comparisons/", "gcomp_sde_", n_sim, ".csv", sep = "")
)
sde_bias_iptw <- read.csv(
  paste(results_path, "sl_comparisons/", "iptw_sde_", n_sim, ".csv", sep = "")
)
sde_bias_tmle <- read.csv(
  paste(results_path, "sl_comparisons/", "tmle_sde_", n_sim, ".csv", sep = "")
)
sde_bias_onestep <- read.csv(
  paste(results_path, "sl_comparisons/", "onestep_sde_", n_sim, ".csv", sep = "")
)

sde_biases[, 1] <- sde_bias_gcomp$SDE
sde_biases[, 2] <- sde_bias_iptw$SDE
sde_biases[, 3] <- sde_bias_tmle$SDE
sde_biases[, 4] <- sde_bias_onestep$SDE

# SIE
sie_biases <- matrix(nrow = n_sim, ncol = 4)
colnames(sie_biases) <- c("SIE_gcomp", "SIE_IPTW", "SIE_TMLE", "SIE_onestep")

sie_biases[, 1] <- bias_estimates_gcomp$SIE
sie_biases[, 2] <- bias_estimates_iptw$SIE
sie_biases[, 3] <- bias_estimates_tmle$SIE
sie_biases[, 4] <- bias_estimates_onestep$SIE

plot_comparison(sde_biases)
plot_comparison(sie_biases)