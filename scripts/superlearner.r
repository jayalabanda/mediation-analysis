library(SuperLearner)
library(ranger)
library(arm)
library(speedglm)
library(medoutcon)
library(tidyverse)

libs <- c(
  "SL.glm", "SL.glm.interaction", "SL.ranger",
  "SL.bayesglm", "SL.lm", "SL.speedlm"
)
learners <- c(sl3::Lrnr_bayesglm$new(), sl3::Lrnr_glm_fast$new(), sl3::Lrnr_ranger$new())

set.seed(42)
n_sim <- 50
idx <- sample(1:1000, n_sim)

file_path <- "../Data/"
data <- data.frame(read.csv(paste(file_path, "data_sim.csv", sep = "")))
data <- subset(data, select = -c(Y_qol))
head(data)


## Estimation manuelle (g-comp)
estimate_manual <- function(data, sl_library) {
  tempdat <- data
  colnames(tempdat) <- c("w1", "w2", "w3", "a", "z", "m", "y")

  y <- tempdat$m
  x <- subset(tempdat, select = c("w1", "w2", "w3", "a"))
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
  x <- subset(tempdat, select = c("w1", "w2", "w3", "a"))
  q_model_a1 <- SuperLearner(y, x,
    family = "quasibinomial",
    SL.library = sl_library
  )

  y <- tempdat$q_gamma_a0
  x <- subset(tempdat, select = c("w1", "w2", "w3", "a"))
  q_model_a0 <- SuperLearner(y, x,
    family = "quasibinomial",
    SL.library = sl_library
  )

  q_pred_a1_gamma_a1 <- predict(q_model_a1, data_a1, type = "response")
  q_pred_a1_gamma_a0 <- predict(q_model_a0, data_a1, type = "response")
  q_pred_a0_gamma_a0 <- predict(q_model_a0, data_a0, type = "response")

  psi_m_rnde <- mean(q_pred_a1_gamma_a0$pred) - mean(q_pred_a0_gamma_a0$pred)
  psi_m_rnie <- mean(q_pred_a1_gamma_a1$pred) - mean(q_pred_a1_gamma_a0$pred)

  results <- list(
    psi_m_rnde = psi_m_rnde,
    psi_m_rnie = psi_m_rnie
  )
  return(results)
}

results <- estimate_manual(data, "SL.speedlm")
results

### SL.glm.interaction
# $psi_m_rnde
# [1] 0.06221807

# $psi_m_rnie
# [1] 0.009689907

### SL.ranger
# [1] 0.03603389

# $psi_m_rnie
# [1] 0.005867094

### SL.bayesglm
# $psi_m_rnde
# [1] 0.062521

# $psi_m_rnie
# [1] 0.009826797

### SL.glm
# $psi_m_rnde
# [1] 0.0625841

# $psi_m_rnie
# [1] 0.009845864

### SL.lm
# $psi_m_rnde
# [1] 0.06600732

# $psi_m_rnie
# [1] 0.008698918

### SL.speedlm
# $psi_m_rnde
# [1] 0.06600732

# $psi_m_rnie
# [1] 0.008698918


## Article Kara E. Rudolph (TMLE)
estimate_rudolph <- function(data, sl_library) {
  obsdat <- data
  colnames(obsdat) <- c("w1", "w2", "a", "z", "m", "y")

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

  results <- data.frame(
    sde = sde_tmle,
    sde_var = var_sde_eic,
    sie = sie_tmle,
    sie_var = var_sie_eic,
    sde_lb = sde_tmle - 1.96 * sqrt(var_sde_eic), # sde lower bound
    sde_ub = sde_tmle + 1.96 * sqrt(var_sde_eic), # sde upper bound
    sie_lb = sie_tmle - 1.96 * sqrt(var_sie_eic), # sie lower bound
    sie_ub = sie_tmle + 1.96 * sqrt(var_sie_eic) # sie upper bound
  )
  return(results)
}

results <- estimate_rudolph(data, "SL.speedglm")
results$sde
results$sie

### SL.glm
# > results$sde
# [1] 0.08026604

# > results$sie
# [1] 0.008459422

### SL.glm.interaction
# > results$sde
# [1] 0.07943351

# > results$sie
# [1] 0.008547797

### SL.ranger
# > results$sde
# [1] 0.07039337

# > results$sie
# [1] 0.006884582

### SL.bayesglm
# > results$sde
# [1] 0.08018934

# > results$sie
# [1] 0.008446274

### SL.lm
# > results$sde
# [1] 0.08432383

# > results$sie
# [1] 0.007735236

### SL.speedlm
# > results$sde
# [1] 0.08432383

# > results$sie
# [1] 0.007735236


### Medoutcon (one-step)
ind_dir_effects_medoutcon <- function(data, w_names, m_names, learner) {
  dir_os <- medoutcon(
    W = data[, w_names],
    A = data$A,
    Z = data$Z,
    M = data[, m_names],
    Y = data$Y,
    effect = "direct",
    estimator = "onestep",
    u_learners = learner,
    v_learners = learner
  )

  ind_os <- medoutcon(
    W = data[, w_names],
    A = data$A,
    Z = data$Z,
    M = data[, m_names],
    Y = data$Y,
    effect = "indirect",
    estimator = "onestep",
    u_learners = learner,
    v_learners = learner
  )

  return(list(
    dir_result_os = dir_os,
    ind_result_os = ind_os
  ))
}


### IPTW
iptw_direct_indirect <- function(data, sl_library) {
  g_a <- glm(a ~ 1, family = "binomial", data = data)

  y <- data$a
  x <- subset(data, select = c("w1", "w2"))
  g_a_l0 <- SuperLearner(y, x, family = binomial(), SL.library = sl_library)

  y <- data$m
  x <- subset(data, select = c(a))
  g_m_a <- SuperLearner(y, x, family = binomial(), SL.library = sl_library)

  x <- subset(data, select = c("w1", "w2", "a", "z"))
  g_m_l <- SuperLearner(y, x, family = binomial(), SL.library = sl_library)

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

  return(list(iptw_edn = iptw_edn, iptw_ein = iptw_ein))
}

## Save results
save_path <- "./Results/estimates/"

#### Rudolph article: TMLE
file_path <- "../Data/simulations_rudolph/"

results_sde <- matrix(nrow = n_sim, ncol = length(libs))
results_sie <- matrix(nrow = n_sim, ncol = length(libs))
colnames(results_sde) <- libs
colnames(results_sie) <- libs

for (i in 1:n_sim) {
  print(paste0("Simulation ", i))

  data <- read.csv(paste0(file_path, "data_", idx[i], ".csv"))
  for (j in seq_len(length(libs))) {
    results <- estimate_rudolph(data, libs[j])
    results_sde[i, j] <- results$sde
    results_sie[i, j] <- results$sie
  }
}

results_sde
results_sie

write.csv(results_sde, paste(save_path, "estimates_sde_SL_rud.csv"), row.names = FALSE)
write.csv(results_sie, paste(save_path, "estimates_sie_SL_rud.csv"), row.names = FALSE)

boxplot(results_sde)
boxplot(results_sie)


##### G-comp
file_path <- "../Data/simulations_rudolph/"

results_sde <- matrix(nrow = n_sim, ncol = length(libs))
results_sie <- matrix(nrow = n_sim, ncol = length(libs))
colnames(results_sde) <- libs
colnames(results_sie) <- libs

for (i in 1:n_sim) {
  print(paste0("Simulation ", i))

  data <- read.csv(paste0(file_path, "data_", idx[i], ".csv"))
  for (j in seq_len(length(libs))) {
    results <- estimate_manual(data, libs[j])
    results_sde[i, j] <- results$psi_m_rnde
    results_sie[i, j] <- results$psi_m_rnie
  }
}

results_sde
results_sie

write.csv(results_sde, paste(save_path, "estimates_sde_SL.csv"), row.names = FALSE)
write.csv(results_sie, paste(save_path, "estimates_sie_SL.csv"), row.names = FALSE)

boxplot(results_sde)
boxplot(results_sie)


#### Positivity
file_path <- "../Data/new_simulations/"

results_sde <- matrix(nrow = n_sim, ncol = length(libs))
results_sie <- matrix(nrow = n_sim, ncol = length(libs))
colnames(results_sde) <- libs
colnames(results_sie) <- libs

for (i in 1:n_sim) {
  print(paste0("Simulation ", i))

  data <- read.csv(paste0(file_path, "data_", idx[i], ".csv"))
  data <- subset(data, select = -c(y_qol))

  for (j in seq_len(length(libs))) {
    results <- estimate_rudolph(data, libs[j])
    results_sde[i, j] <- results$sde
    results_sie[i, j] <- results$sie
  }
}

results_sde
results_sie

write.csv(results_sde, paste(save_path, "estimates_sde_posit.csv"), row.names = FALSE)
write.csv(results_sie, paste(save_path, "estimates_sie_posit.csv"), row.names = FALSE)

boxplot(results_sde)
boxplot(results_sie)


#### Medoutcon (One-step)
file_path <- "../Data/new_simulations/"

results_sde <- matrix(nrow = n_sim, ncol = length(learners))
results_sie <- matrix(nrow = n_sim, ncol = length(learners))
colnames(results_sde) <- c("Bayesglm", "Glm_fast", "Ranger")
colnames(results_sie) <- c("Bayesglm", "Glm_fast", "Ranger")

for (i in 1:n_sim) {
  print(paste0("Simulation ", i))

  data <- read.csv(paste0(file_path, "data_", idx[i], ".csv"))
  data <- subset(data, select = -c(y_qol))
  colnames(data) <- c("W_1", "W_2", "A", "Z", "M_1", "Y")

  w_names <- str_subset(colnames(data), "W")
  m_names <- str_subset(colnames(data), "M")

  for (j in seq_len(length(learners))) {
    results <- ind_dir_effects_medoutcon(data, w_names, m_names, learners[j])
    results_sde[i, j] <- results$dir_result_os$theta
    results_sie[i, j] <- results$ind_result_os$theta
  }
}

results_sde
results_sie

write.csv(results_sde, paste(save_path, "estimates_sde_posi_moc.csv"), row.names = FALSE)
write.csv(results_sie, paste(save_path, "estimates_sie_posi_moc.csv"), row.names = FALSE)

boxplot(results_sde)
boxplot(results_sie)


##### IPTW positivity
file_path <- "../Data/new_simulations/"

results_sde <- matrix(nrow = n_sim, ncol = length(libs))
results_sie <- matrix(nrow = n_sim, ncol = length(libs))
colnames(results_sde) <- libs
colnames(results_sie) <- libs

for (i in 1:n_sim) {
  print(paste0("Simulation ", i))

  data <- read.csv(paste0(file_path, "data_", idx[i], ".csv"))
  data <- subset(data, select = -c(w3))
  colnames(data) <- c("w1", "w2", "a", "z", "m", "y")

  for (j in seq_len(length(libs))) {
    try(results <- iptw_direct_indirect(data, libs[j]))
    results_sde[i, j] <- results$iptw_edn
    results_sie[i, j] <- results$iptw_ein
  }
}

results_sde
results_sie

write.csv(results_sde, paste(save_path, "estimates_sde_posit_iptw.csv"), row.names = FALSE)
write.csv(results_sie, paste(save_path, "estimates_sie_posit_iptw.csv"), row.names = FALSE)

boxplot(results_sde)
boxplot(results_sie)


##### IPTW positivity with quantitative variables
file_path <- "../Data/quantitative_simulations/"

results_sde <- matrix(nrow = n_sim, ncol = length(libs))
results_sie <- matrix(nrow = n_sim, ncol = length(libs))
colnames(results_sde) <- libs
colnames(results_sie) <- libs

for (i in 1:n_sim) {
  print(paste0("Simulation ", i))

  data <- read.csv(paste0(file_path, "data_", idx[i], ".csv"))
  data <- subset(data, select = -c(y_qol))
  colnames(data) <- c("w1", "w2", "a", "z", "m", "y")

  for (j in seq_len(length(libs))) {
    try(results <- iptw_direct_indirect(data, libs[j]))
    results_sde[i, j] <- results$iptw_edn
    results_sie[i, j] <- results$iptw_ein
  }
}

results_sde
results_sie

write.csv(results_sde, paste(save_path, "estimates_sde_posit_quant_iptw.csv"), row.names = FALSE)
write.csv(results_sie, paste(save_path, "estimates_sie_posit_quant_iptw.csv"), row.names = FALSE)

boxplot(results_sde)
boxplot(results_sie)