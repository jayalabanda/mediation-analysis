file_path <- "../Data/"
data <- data.frame(read.csv(paste(file_path, "data_sim.csv", sep = "")))
data <- subset(data, select = -c(Y_qol))
head(data)

estimate_effects <- function(data, treatment, covariates, control, mediator,
                             outcome, amodel, ymodel, mmodel, zmodel, qmodel) {
  obsdat <- data
  colnames(obsdat)[grep(treatment, colnames(obsdat))] <- "a"
  colnames(obsdat)[grep(control, colnames(obsdat))] <- "z"
  colnames(obsdat)[grep(mediator, colnames(obsdat))] <- "m"
  colnames(obsdat)[grep(outcome, colnames(obsdat))] <- "y"

  for (i in seq_len(length(covariates))) {
    colnames(obsdat)[grep(covariates[i], colnames(obsdat))] <- paste0("w", i)
  }

  tmpdat <- obsdat

  mfit <- glm(mmodel, data = obsdat, family = "binomial")
  yfit <- glm(ymodel, data = obsdat, family = "binomial")

  dfa1 <- dfa0 <- dfm1 <- dfm0 <- obsdat

  dfa1$a <- dfm1$m <- 1
  dfa0$a <- dfm0$m <- 0

  zfit <- glm(formula = zmodel, family = "binomial", data = obsdat)
  z_a0 <- predict(zfit, newdata = dfa0, type = "response")
  z_a1 <- predict(zfit, newdata = dfa1, type = "response")

  mfit <- glm(formula = mmodel, family = "binomial", data = obsdat)
  m_a0 <- predict(mfit, newdata = dfa0, type = "response")
  m_a1 <- predict(mfit, newdata = dfa1, type = "response")

  gm_a0 <- m_a0 * z_a0 + m_a0 * (1 - z_a0)
  gm_a1 <- m_a1 * z_a1 + m_a1 * (1 - z_a1)

  afit <- glm(formula = amodel, family = "binomial", data = tmpdat)
  pred <- predict(afit, type = "response")
  ps_a1 <- I(tmpdat$a == 1) / pred
  ps_a0 <- I(tmpdat$a == 0) / (1 - pred)

  mfit <- glm(formula = mmodel, family = "binomial", data = tmpdat)
  mazw <- predict(mfit, type = "response")
  psm <- I(tmpdat$m == 1) * mazw + I(tmpdat$m == 0) * (1 - mazw)

  yfit <- glm(formula = ymodel, family = "binomial", data = tmpdat)
  tmpdat$qyinit <- cbind(
    predict(yfit, newdata = tmpdat, type = "response"),
    predict(yfit, newdata = dfm0, type = "response"),
    predict(yfit, newdata = dfm1, type = "response")
  )

  # ----------------------------
  # Get E(Y_{1, g_0})
  # ----------------------------
  b <- I(tmpdat$m) == 1
  c <- I(tmpdat$m) == 0
  tmpdat$h_a1_gm_a0 <- ((b * gm_a0 + c * (1 - gm_a0)) / psm) * ps_a1

  epsilon_a1_gm_a0 <- coef(glm(y ~ 1,
    weights = h_a1_gm_a0, offset = qlogis(qyinit[, 1]),
    family = "quasibinomial", data = tmpdat
  ))

  tmpdat$qyup_m0_a1_gm_a0 <-
    plogis(qlogis(tmpdat$qyinit[, 2]) + epsilon_a1_gm_a0)
  tmpdat$qyup_m1_a1_gm_a0 <-
    plogis(qlogis(tmpdat$qyinit[, 3]) + epsilon_a1_gm_a0)

  tmpdat$q_a1_gm_a0 <-
    tmpdat$qyup_m0_a1_gm_a0 * (1 - gm_a0) + tmpdat$qyup_m1_a1_gm_a0 * gm_a0

  q_a1_g0_fit <- glm(
    formula = paste("q_a1_gm_a0", qmodel, sep = "~"), # zmodel
    family = "quasibinomial",
    data = tmpdat[tmpdat$a == 1, ]
  )
  q_a1_g0 <- predict(q_a1_g0_fit, newdata = tmpdat, type = "response")

  # Necessary if A is non random
  epsilon_z_a1_gm_a0 <- coef(glm(q_a1_gm_a0 ~ 1,
    weights = ps_a1, offset = qlogis(q_a1_g0),
    family = "quasibinomial", data = tmpdat
  ))
  q_zup_a1_gm_a0 <- plogis(qlogis(q_a1_g0) + epsilon_z_a1_gm_a0)

  # Calculate TMLE
  tmle_a1_m0 <- mean(q_zup_a1_gm_a0)

  # ----------------------------
  # Get E(Y_{0, g_0})
  # ----------------------------
  b <- tmpdat$m
  tmpdat$h_a0_gm_a0 <- ((b * gm_a0 + (1 - b) * (1 - gm_a0)) / psm) * ps_a0

  epsilon_a0_gm_a0 <- coef(glm(y ~ 1,
    weights = h_a0_gm_a0, offset = qlogis(qyinit[, 1]),
    family = "quasibinomial", data = tmpdat
  ))

  tmpdat$qyup_m0_a0_gm_a0 <-
    plogis(qlogis(tmpdat$qyinit[, 2]) + epsilon_a0_gm_a0)
  tmpdat$qyup_m1_a0_gm_a0 <-
    plogis(qlogis(tmpdat$qyinit[, 3]) + epsilon_a0_gm_a0)

  tmpdat$q_a0_gm_a0 <-
    tmpdat$qyup_m0_a0_gm_a0 * (1 - gm_a0) + tmpdat$qyup_m1_a0_gm_a0 * gm_a0

  q_a0_g0_fit <- glm(
    formula = paste("q_a0_gm_a0", qmodel, sep = "~"), # zmodel
    family = "quasibinomial",
    data = tmpdat[tmpdat$a == 0, ]
  )
  q_a0_g0 <- predict(q_a0_g0_fit, newdata = tmpdat, type = "response")

  # Necessary if A is non random
  epsilon_z_a0_gm_a0 <- coef(glm(q_a0_gm_a0 ~ 1,
    weights = ps_a0, offset = qlogis(q_a0_g0),
    family = "quasibinomial", data = tmpdat
  ))
  q_zup_a0_gm_a0 <- plogis(qlogis(q_a0_g0) + epsilon_z_a0_gm_a0)

  # Calculate TMLE
  tmle_a0_m0 <- mean(q_zup_a0_gm_a0)

  # ----------------------------
  # Get E(Y_{1, g_1})
  # ----------------------------
  tmpdat$h_a1_gm_a1 <- ((b * gm_a1 + (1 - b) * (1 - gm_a1)) / psm) * ps_a1

  epsilon_a1_gm_a1 <- coef(glm(y ~ 1,
    weights = h_a1_gm_a1, offset = qlogis(qyinit[, 1]),
    family = "quasibinomial", data = tmpdat
  ))

  tmpdat$qyup_m0_a1_gm_a1 <-
    plogis(qlogis(tmpdat$qyinit[, 2]) + epsilon_a1_gm_a1)
  tmpdat$qyup_m1_a1_gm_a1 <-
    plogis(qlogis(tmpdat$qyinit[, 3]) + epsilon_a1_gm_a1)

  tmpdat$q_a1_gm_a1 <-
    tmpdat$qyup_m0_a1_gm_a1 * (1 - gm_a1) + tmpdat$qyup_m1_a1_gm_a1 * gm_a1

  q_a1_g1_fit <- glm(
    formula = paste("q_a1_gm_a1", qmodel, sep = "~"), # zmodel
    family = "quasibinomial",
    data = tmpdat[tmpdat$a == 1, ]
  )
  q_a1_g1 <- predict(q_a1_g1_fit, newdata = tmpdat, type = "response")

  # Necessary if A is non random
  epsilon_z_a1_gm_a1 <- coef(glm(q_a1_gm_a1 ~ 1,
    weights = ps_a1, offset = qlogis(q_a1_g1),
    family = "quasibinomial", data = tmpdat
  ))
  q_zup_a1_gm_a1 <- plogis(qlogis(q_a1_g1) + epsilon_z_a1_gm_a1)

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
  eic2_a1_g0 <- ps_a1 * (tmpdat$q_a1_gm_a0 - q_zup_a1_gm_a0)
  eic3_a1_g0 <- q_zup_a1_gm_a0 - tmle_a1_m0

  eic_a1_g0 <- eic1_a1_g0 + eic2_a1_g0 + eic3_a1_g0

  # EIC for E(Y_{1, g_1})
  eic1_a1_g1 <- tmpdat$h_a1_gm_a1 * (tmpdat$y - tmpdat$qyup_a1_g1)
  eic2_a1_g1 <- ps_a1 * (tmpdat$q_a1_gm_a1 - q_zup_a1_gm_a1)
  eic3_a1_g1 <- q_zup_a1_gm_a1 - tmle_a1_m1

  eic_a1_g1 <- eic1_a1_g1 + eic2_a1_g1 + eic3_a1_g1

  # EIC for E(Y_{0, g_0})
  eic1_a0_g0 <- tmpdat$h_a0_gm_a0 * (tmpdat$y - tmpdat$qyup_a0_g0)
  eic2_a0_g0 <- ps_a0 * (tmpdat$q_a0_gm_a0 - q_zup_a0_gm_a0)
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
    cbind(
      sde = sde_tmle,
      sde_var = var_sde_eic,
      sie = sie_tmle,
      sie_var = var_sie_eic,
      sde_lb = sde_tmle - 1.96 * sqrt(var_sde_eic), # sde lower bound
      sde_ub = sde_tmle + 1.96 * sqrt(var_sde_eic), # sde upper bound
      sie_lb = sie_tmle - 1.96 * sqrt(var_sie_eic), # sie lower bound
      sie_ub = sie_tmle + 1.96 * sqrt(var_sie_eic) # sie upper bound
    )
  )

  return(results)
}

covariates <- c("L0_male", "L0_parent_low_educ_lv")

amodel <- "a ~ w1 + w2"
zmodel <- "z ~ w1 + w2 + a"
mmodel <- "m ~ w1 + w2 + a + z"
ymodel <- "y ~ w1 + w2 + a + z + m"
qmodel <- "w1 + w2"

results <- estimate_effects(
  data, "A0_ace", covariates,
  "L1", "M_smoking", "Y_death",
  amodel, ymodel, mmodel, zmodel, qmodel
)

results
# sde     0.06180492
# sde_var 0.0001999548
# sie     0.009527786
# sie_var 1.011022e-05
# sde_lb  0.03408947
# sde_ub  0.08952037
# sie_lb  0.00329566
# sie_ub  0.01575991

## Bootstrap
n_sim <- 200
n_boot <- 200
true_sde <- 0.0625841
true_sie <- 0.009845864
sim_data_path <- "../Data/simulations/"
covariates <- c("l0_male", "l0_parent_low_educ_lv")

# Direct effect
estimates_sde_rud <- matrix(NA, ncol = 3, nrow = n_sim)
colnames(estimates_sde_rud) <- c("sde_rud", "sd_rud", "cov_rud")
write.csv(
  estimates_sde_rud,
  file = paste(file_path, "estimates_sde_rud.csv", sep = ""),
  row.names = FALSE
)

# Indirect effect
estimates_sie_rud <- matrix(NA, ncol = 3, nrow = n_sim)
colnames(estimates_sie_rud) <- c("sie_rud", "sd_rud", "cov_rud")
write.csv(
  estimates_sie_rud,
  file = paste(file_path, "estimates_sie_rud.csv", sep = ""),
  row.names = FALSE
)

start_time <- Sys.time()
for (i in 1:n_sim) {
  print(paste0("Simulation n°: ", i))
  data_sim <- read.csv(paste0(sim_data_path, "data_", i, ".csv", sep = ""))
  data_sim <- subset(data_sim, select = -c(y_qol))

  rud_estimates <- estimate_effects(
    data = data_sim,
    treatment = "a0_ace",
    covariates = covariates,
    control = "l1",
    mediator = "m_smoking",
    outcome = "y_death",
    amodel = amodel,
    ymodel = ymodel,
    mmodel = mmodel,
    zmodel = zmodel,
    qmodel = qmodel
  )

  estimates_sde_rud <- read.csv(
    paste(file_path, "estimates_sde_rud.csv", sep = "")
  )
  estimates_sde_rud[i, "sde_rud"] <- rud_estimates$sde

  estimates_sie_rud <- read.csv(
    paste(file_path, "estimates_sie_rud.csv", sep = "")
  )
  estimates_sie_rud[i, "sie_rud"] <- rud_estimates$sie

  boot_sde_estimates <- data.frame(matrix(NA, nrow = n_boot, ncol = 1))
  colnames(boot_sde_estimates) <- "boot_sde_estimate"

  boot_sie_estimates <- data.frame(matrix(NA, nrow = n_boot, ncol = 1))
  colnames(boot_sie_estimates) <- "boot_sie_estimate"

  for (b in 1:n_boot) {
    idx <- sample(seq(1, nrow(data_sim)), replace = TRUE)
    boot_data <- data_sim[idx, ]

    if (b %% 100 == 0 & b != 0) {
      print(paste0("Bootstrap number: ", b))
    }

    rud_boot <- estimate_effects(
      data = boot_data,
      treatment = "a0_ace",
      covariates = covariates,
      control = "l1",
      mediator = "m_smoking",
      outcome = "y_death",
      amodel = amodel,
      ymodel = ymodel,
      mmodel = mmodel,
      zmodel = zmodel,
      qmodel = qmodel
    )

    boot_sde_estimates[b, ] <- rud_boot$sde
    boot_sie_estimates[b, ] <- rud_boot$sie
  }

  boot_sde_est <- boot_sde_estimates$boot_sde_estimate
  estimates_sde_rud[i, "sd_rud"] <- sd(boot_sde_est)
  estimates_sde_rud[i, "cov_rud"] <- as.numeric(
    true_sde >= estimates_sde_rud[i, 1] - 1.96 * sd(boot_sde_est) &
      true_sde <= estimates_sde_rud[i, 1] + 1.96 * sd(boot_sde_est)
  )

  boot_sie_est <- boot_sie_estimates$boot_sie_estimate
  estimates_sie_rud[i, "sd_rud"] <- sd(boot_sie_est)
  estimates_sie_rud[i, "cov_rud"] <- as.numeric(
    true_sie >= estimates_sie_rud[i, 1] - 1.96 * sd(boot_sie_est) &
      true_sie <= estimates_sie_rud[i, 1] + 1.96 * sd(boot_sie_est)
  )

  write.csv(estimates_sde_rud,
    file = paste(file_path, "estimates_sde_rud.csv", sep = ""),
    row.names = FALSE
  )
  write.csv(estimates_sie_rud,
    file = paste(file_path, "estimates_sie_rud.csv", sep = ""),
    row.names = FALSE
  )
}
end_time <- Sys.time()
diff <- end_time - start_time
diff

estimates_sde_rud <- read.csv(
  paste(file_path, "estimates_sde_rud.csv", sep = "")
)
estimates_sde_rud <- data.frame(estimates_sde_rud)
head(estimates_sde_rud)

estimates_sie_rud <- read.csv(
  paste(file_path, "estimates_sie_rud.csv", sep = "")
)
estimates_sie_rud <- data.frame(estimates_sie_rud)
head(estimates_sie_rud)

## Calculate results

# bias
bias_sde <- mean(estimates_sde_rud$sde_rud) - true_sde
bias_sie <- mean(estimates_sie_rud$sie_rud) - true_sie

# variance & standard error
n_rows <- nrow(data_sim)
var_sde <- mean(
  (estimates_sde_rud$sde_rud - mean(estimates_sde_rud$sde_rud))^2
) * n_rows / (n_rows - 1)
var_sie <- mean(
  (estimates_sie_rud$sie_rud - mean(estimates_sie_rud$sie_rud))^2
) * n_rows / (n_rows - 1)

se_sde <- sqrt(var_sde)
se_sie <- sqrt(var_sie)

# Standardized bias
sd_bias_sde <- bias_sde / se_sde
sd_bias_sie <- bias_sie / se_sie

# MSE
mse_sde <- var_sde + bias_sde^2
mse_sie <- var_sie + bias_sie^2

# Average estimated standard error
av_estimated_se_sde <- sqrt(mean(estimates_sde_rud$sd_rud^2))
av_estimated_se_sie <- sqrt(mean(estimates_sie_rud$sd_rud^2))

# Coverage
cov_sde <- mean(estimates_sde_rud$cov_rud)
cov_sie <- mean(estimates_sie_rud$cov_rud)

# Summarize results
results_sde <- data.frame(
  bias = bias_sde, # -0.0003190485
  variance = var_sde, # 0.0002560941
  STD = se_sde, # 0.01600294
  std.bias = sd_bias_sde, # -0.01993687
  MSE = mse_sde, # 0.0002561959
  av.est.std = av_estimated_se_sde, # 0.01501659
  coverage = cov_sde # 0.93
)

results_sde
write.csv(
  results_sde,
  paste(file_path, "results_sde_rud.csv", sep = ""),
  row.names = FALSE
)

results_sie <- data.frame(
  bias = bias_sie, # 0.000745668
  variance = var_sie, # 1.637447e-05
  STD = se_sie, # 0.004046538
  std.bias = sd_bias_sie, # 0.1842731
  MSE = mse_sie, # 1.693049e-05
  av.est.std = av_estimated_se_sie, # 0.004210313
  coverage = cov_sie # 0.965
)

results_sie
write.csv(
  results_sie,
  paste(file_path, "results_sie_rud.csv", sep = ""),
  row.names = FALSE
)

# n_sim = 100, n_boot = 500
# 4h

# Avec M binaire
#         bias     variance        STD  std.bias          MSE
# 1 0.00396171 0.0002467648 0.01570875 0.2521977 0.0002624599
# av.est.std coverage
# 0.01472793     0.93

#           bias     variance         STD   std.bias          MSE
# 1 -0.001949152 9.766798e-06 0.003125188 -0.6236911 1.356599e-05
# av.est.std coverage
# 0.003218264    0.88


# n_sim = 200, n_boot = 200
# 3h14

# Avec M binaire
#         bias     variance        STD  std.bias          MSE
# 1 0.00260686 0.0002538018 0.01593116 0.1636328 0.0002605975
# av.est.std coverage
#  0.0147856    0.925

#           bias     variance         STD   std.bias          MSE
# 1 -0.002078757 9.433828e-06 0.003071454 -0.6767991 1.375506e-05
#    av.est.std  coverage
# 1 0.003223341   0.855


################################################
################################################
## Medoutcon
library(medoutcon)
library(tidyverse)
library(data.table)
library(hal9001)
library(sl3)

file_path <- "../Data/"
data <- data.frame(read.csv(paste(file_path, "data_sim.csv", sep = "")))
data <- subset(data, select = -c(Y_qol))
head(data)

colnames(data)[1] <- "W_1"
colnames(data)[2] <- "W_2"
colnames(data)[3] <- "A"
colnames(data)[4] <- "Z"
colnames(data)[5] <- "M"
colnames(data)[6] <- "Y"
head(data)

w_names <- str_subset(colnames(data), "W")
m_names <- str_subset(colnames(data), "M")

os_de <- medoutcon(
  W = data[, w_names],
  A = data$A,
  Z = data$Z,
  M = data[, m_names],
  Y = data$Y,
  effect = "direct",
  estimator = "onestep",
  u_learners = sl3::Lrnr_glm_fast$new(),
  v_learners = sl3::Lrnr_glm_fast$new()
)

os_de
# Interventional Direct Effect
# Estimator: onestep
# Estimate: 0.061
# Std. Error: 0.014
# 95% CI: [0.032, 0.088]

tmle_de <- medoutcon(
  W = data[, w_names],
  A = data$A,
  Z = data$Z,
  M = data[, m_names],
  Y = data$Y,
  effect = "direct",
  estimator = "tmle",
  u_learners = sl3::Lrnr_glm_fast$new(),
  v_learners = sl3::Lrnr_glm_fast$new()
)

tmle_de
# Interventional Direct Effect
# Estimator: tmle
# Estimate: 0.034
# Std. Error: 0.014
# 95% CI: [0.005, 0.062]

os_ie <- medoutcon(
  W = data[, w_names],
  A = data$A,
  Z = data$Z,
  M = data[, m_names],
  Y = data$Y,
  effect = "indirect",
  estimator = "onestep",
  u_learners = sl3::Lrnr_glm_fast$new(),
  v_learners = sl3::Lrnr_glm_fast$new()
)

os_ie
# Interventional Indirect Effect
# Estimator: onestep
# Estimate: 0.011
# Std. Error: 0.004
# 95% CI: [0.004, 0.018]

tmle_ie <- medoutcon(
  W = data[, w_names],
  A = data$A,
  Z = data$Z,
  M = data[, m_names],
  Y = data$Y,
  effect = "indirect",
  estimator = "tmle",
  u_learners = sl3::Lrnr_glm_fast$new(),
  v_learners = sl3::Lrnr_glm_fast$new()
)

tmle_ie
# Interventional Indirect Effect
# Estimator: tmle
# Estimate: 0.009
# Std. Error: 0.004
# 95% CI: [0.001, 0.016]

# Create practical function
ind_dir_effects_medoutcon <- function(data, w_names, m_names) {
  dir_os <- medoutcon(
    W = data[, w_names],
    A = data$A,
    Z = data$Z,
    M = data[, m_names],
    Y = data$Y,
    effect = "direct",
    estimator = "onestep"
    # u_learners = sl3::Lrnr_glm_fast$new(),
    # v_learners = sl3::Lrnr_glm_fast$new()
  )

  dir_tmle <- medoutcon(
    W = data[, w_names],
    A = data$A,
    Z = data$Z,
    M = data[, m_names],
    Y = data$Y,
    effect = "direct",
    estimator = "tmle"
    # u_learners = sl3::Lrnr_glm_fast$new(),
    # v_learners = sl3::Lrnr_glm_fast$new()
  )

  ind_os <- medoutcon(
    W = data[, w_names],
    A = data$A,
    Z = data$Z,
    M = data[, m_names],
    Y = data$Y,
    effect = "indirect",
    estimator = "onestep"
    # u_learners = sl3::Lrnr_glm_fast$new(),
    # v_learners = sl3::Lrnr_glm_fast$new()
  )

  ind_tmle <- medoutcon(
    W = data[, w_names],
    A = data$A,
    Z = data$Z,
    M = data[, m_names],
    Y = data$Y,
    effect = "indirect",
    estimator = "tmle"
    # u_learners = sl3::Lrnr_glm_fast$new(),
    # v_learners = sl3::Lrnr_glm_fast$new()
  )

  return(list(
    dir_result_os = dir_os,
    dir_result_tmle = dir_tmle,
    ind_result_os = ind_os,
    ind_result_tmle = ind_tmle
  ))
}

# Define expit function
expit <- function(x) {
  return(1 / (1 + exp(-x)))
}

## Multiple simulaitions
n_sim <- 200
true_sde <- 0.0625841
true_sie <- 0.009845864
sim_data_path <- "../Data/simulations/"

# Direct effect
estimates_sde_moc <- matrix(NA, ncol = 6, nrow = n_sim)
colnames(estimates_sde_moc) <- c(
  "sde_os", "sd_os", "cov_os",
  "sde_tmle", "sd_tmle", "cov_tmle"
)
write.csv(
  estimates_sde_moc,
  file = paste(file_path, "estimates_sde_moc.csv", sep = ""),
  row.names = FALSE
)

# Indirect effect
estimates_sie_moc <- matrix(NA, ncol = 6, nrow = n_sim)
colnames(estimates_sie_moc) <- c(
  "sie_os", "sd_os", "cov_os",
  "sie_tmle", "sd_tmle", "cov_tmle"
)
write.csv(
  estimates_sie_moc,
  file = paste(file_path, "estimates_sie_moc.csv", sep = ""),
  row.names = FALSE
)

start_time <- Sys.time()
for (i in 1:n_sim) {
  print(paste0("Simulation n°: ", i))
  data_sim <- read.csv(paste0(sim_data_path, "data_", i, ".csv", sep = ""))
  data_sim <- subset(data_sim, select = -c(y_qol))

  colnames(data_sim)[1] <- "W_1"
  colnames(data_sim)[2] <- "W_2"
  colnames(data_sim)[3] <- "A"
  colnames(data_sim)[4] <- "Z"
  colnames(data_sim)[5] <- "M_1"
  colnames(data_sim)[6] <- "Y"

  # data_sim[, "M_2"] <- rbinom(
  #   nrow(data_sim), 1,
  #   expit(data_sim$W_1 + data_sim$W_2 + data_sim$A - data_sim$Z +
  #     data_sim$A * data_sim$Z - 0.3 * data_sim$A * data_sim$W_2)
  # )

  # data_sim[, "M_3"] <- rbinom(
  #   nrow(data_sim), 1,
  #   expit(data_sim$W_1 - data_sim$W_2 + 0.2 * data_sim$A + data_sim$Z +
  #     0.1 * data_sim$A * data_sim$Z - 0.4 * data_sim$A * data_sim$W_2)
  # )

  w_names <- str_subset(colnames(data_sim), "W")
  m_names <- str_subset(colnames(data_sim), "M")

  # Calculate effects from function
  results_dir_ind <- ind_dir_effects_medoutcon(data_sim, w_names, m_names)
  # dir_result_os
  # dir_result_tmle
  # ind_result_os
  # ind_result_tmle

  estimates_sde_moc <- read.csv(
    paste(file_path, "estimates_sde_moc.csv", sep = "")
  )
  estimates_sde_moc[i, "sde_os"] <- results_dir_ind$dir_result_os$theta
  estimates_sde_moc[i, "sde_tmle"] <- results_dir_ind$dir_result_tmle$theta

  estimates_sie_moc <- read.csv(
    paste(file_path, "estimates_sie_moc.csv", sep = "")
  )
  estimates_sie_moc[i, "sie_os"] <- results_dir_ind$ind_result_os$theta
  estimates_sie_moc[i, "sie_tmle"] <- results_dir_ind$ind_result_tmle$theta

  stand_dev_dir <- sqrt(results_dir_ind$dir_result_os$var)
  estimates_sde_moc[i, "sd_os"] <- stand_dev_dir
  estimates_sde_moc[i, "cov_os"] <- as.numeric(
    true_sde >= estimates_sde_moc[i, 1] - 1.96 * stand_dev_dir &
      true_sde <= estimates_sde_moc[i, 1] + 1.96 * stand_dev_dir
  )

  stand_dev_dir <- sqrt(results_dir_ind$dir_result_tmle$var)
  estimates_sde_moc[i, "sd_tmle"] <- stand_dev_dir
  estimates_sde_moc[i, "cov_tmle"] <- as.numeric(
    true_sde >= estimates_sde_moc[i, 4] - 1.96 * stand_dev_dir &
      true_sde <= estimates_sde_moc[i, 4] + 1.96 * stand_dev_dir
  )

  stand_dev_ind <- sqrt(results_dir_ind$ind_result_os$var)
  estimates_sie_moc[i, "sd_os"] <- stand_dev_ind
  estimates_sie_moc[i, "cov_os"] <- as.numeric(
    true_sie >= estimates_sie_moc[i, 1] - 1.96 * stand_dev_ind &
      true_sie <= estimates_sie_moc[i, 1] + 1.96 * stand_dev_ind
  )

  stand_dev_ind <- sqrt(results_dir_ind$ind_result_tmle$var)
  estimates_sie_moc[i, "sd_tmle"] <- stand_dev_ind
  estimates_sie_moc[i, "cov_tmle"] <- as.numeric(
    true_sie >= estimates_sie_moc[i, 4] - 1.96 * stand_dev_ind &
      true_sie <= estimates_sie_moc[i, 4] + 1.96 * stand_dev_ind
  )

  write.csv(
    estimates_sde_moc,
    file = paste(file_path, "estimates_sde_moc.csv", sep = ""),
    row.names = FALSE
  )
  write.csv(
    estimates_sie_moc,
    file = paste(file_path, "estimates_sie_moc.csv", sep = ""),
    row.names = FALSE
  )
}
end_time <- Sys.time()
diff <- end_time - start_time

estimates_sde_moc <- read.csv(
  paste(file_path, "estimates_sde_moc.csv", sep = "")
)
estimates_sde_moc <- data.frame(estimates_sde_moc)
head(estimates_sde_moc)

estimates_sie_moc <- read.csv(
  paste(file_path, "estimates_sie_moc.csv", sep = "")
)
estimates_sie_moc <- data.frame(estimates_sie_moc)
head(estimates_sie_moc)

## Calculate results

# estimate
sde_estimate_os <- mean(estimates_sde_moc$sde_os)
sde_estimate_tmle <- mean(estimates_sde_moc$sde_tmle)
sie_estimate_os <- mean(estimates_sie_moc$sie_os)
sie_estimate_tmle <- mean(estimates_sie_moc$sie_tmle)

# bias
bias_sde_os <- mean(estimates_sde_moc$sde_os) - true_sde
bias_sde_tmle <- mean(estimates_sde_moc$sde_tmle) - true_sde
bias_sie_os <- mean(estimates_sie_moc$sie_os) - true_sie
bias_sie_tmle <- mean(estimates_sie_moc$sie_tmle) - true_sie

# variance & standard error
n_rows <- nrow(data_sim)
var_sde_os <- mean(
  (estimates_sde_moc$sde_os - mean(estimates_sde_moc$sde_os))^2
) * n_rows / (n_rows - 1)
var_sie_os <- mean(
  (estimates_sie_moc$sie_os - mean(estimates_sie_moc$sie_os))^2
) * n_rows / (n_rows - 1)
var_sde_tmle <- mean(
  (estimates_sde_moc$sde_tmle - mean(estimates_sde_moc$sde_tmle))^2
) * n_rows / (n_rows - 1)
var_sie_tmle <- mean(
  (estimates_sie_moc$sie_tmle - mean(estimates_sie_moc$sie_tmle))^2
) * n_rows / (n_rows - 1)

se_sde_os <- sqrt(var_sde_os)
se_sde_tmle <- sqrt(var_sde_tmle)
se_sie_os <- sqrt(var_sie_os)
se_sie_tmle <- sqrt(var_sie_tmle)

# Standardized bias
sd_bias_sde_os <- bias_sde_os / se_sde_os
sd_bias_sde_tmle <- bias_sde_tmle / se_sde_tmle
sd_bias_sie_os <- bias_sie_os / se_sie_os
sd_bias_sie_tmle <- bias_sie_tmle / se_sie_tmle

# MSE
mse_sde_os <- var_sde_os + bias_sde_os^2
mse_sde_tmle <- var_sde_tmle + bias_sde_tmle^2
mse_sie_os <- var_sie_os + bias_sie_os^2
mse_sie_tmle <- var_sie_tmle + bias_sie_tmle^2

# Average estimated standard error
av_estimated_se_sde_os <- sqrt(mean(estimates_sde_moc$sd_os^2))
av_estimated_se_sde_tmle <- sqrt(mean(estimates_sde_moc$sd_tmle^2))
av_estimated_se_sie_os <- sqrt(mean(estimates_sie_moc$sd_os^2))
av_estimated_se_sie_tmle <- sqrt(mean(estimates_sie_moc$sd_tmle^2))

# Coverage
cov_sde_os <- mean(estimates_sde_moc$cov_os)
cov_sde_tmle <- mean(estimates_sde_moc$cov_tmle)
cov_sie_os <- mean(estimates_sie_moc$cov_os)
cov_sie_tmle <- mean(estimates_sie_moc$cov_tmle)

# Summarize results
results_sde <- data.frame(
  ## One-step
  estimate_os = sde_estimate_os, # 0.06209866
  bias_os = bias_sde_os, # -0.000485435
  variance_os = var_sde_os, # 0.0002692855
  STD_os = se_sde_os, # 0.01640992
  std_bias_os = sd_bias_sde_os, # -0.0295818
  MSE_os = mse_sde_os, # 0.0002695211
  av_est_std_os = av_estimated_se_sde_os, # 0.0153701
  coverage_os = cov_sde_os, # 0.935
  ## TMLE
  estimate_tmle = sde_estimate_tmle, # 0.0525999
  bias_tmle = bias_sde_tmle, # -0.009984198
  variance_tmle = var_sde_tmle, # 0.0003756425
  STD_tmle = se_sde_tmle, # 0.0193815
  std_bias_tmle = sd_bias_sde_tmle, # -0.5151406
  MSE_tmle = mse_sde_tmle, # 0.0004753267
  av_est_std_tmle = av_estimated_se_sde_tmle, # 0.01535946
  coverage_tmle = cov_sde_tmle, # 0.85
)

results_sde
write.csv(
  results_sde,
  paste(file_path, "results_sde_moc.csv", sep = ""),
  row.names = FALSE
)

results_sie <- data.frame(
  ## One-step
  estimate_os = sie_estimate_os, # 0.01086758
  bias_os = bias_sie_os, # 0.001021716
  variance_os = var_sie_os, # 1.868794e-05
  STD_os = se_sie_os, # 0.004322954
  std_bias_os = sd_bias_sie_os, # 0.2363468
  MSE_os = mse_sie_os, # 1.973184e-05
  av_est_std.os = av_estimated_se_sie_os, # 0.004543823
  coverage_os = cov_sie_os, # 0.96
  ## TMLE
  estimate_tmle = sie_estimate_tmle, # 0.01091418
  bias_tmle = bias_sie_tmle, # 0.001068317
  variance_tmle = var_sie_tmle, # 1.637269e-05
  STD_tmle = se_sie_tmle, # 0.004046318
  std_bias_tmle = sd_bias_sie_tmle, # 0.2640221
  MSE_tmle = mse_sie_tmle, # 1.751399e-05
  av_est_std_tmle = av_estimated_se_sie_tmle, # 0.004375105
  coverage_tmle = cov_sie_tmle, # 0.965
)

results_sie
write.csv(
  results_sie,
  paste(file_path, "results_sie_moc.csv", sep = ""),
  row.names = FALSE
)

# Utiliser Lrnr glm fast dans u_learners et v_learners accélère le calcul
# Par contre utiliser Lrnr glmnet est beaucoup plus lent (comparable avec hal)

# Plus on ajoute de médiateurs, plus les résultats sont variables
# Biais, MSE et variance plus élevés.
# Pas possible de faire plus de variables de confusion avec medoutcon.
# Pas de différence dans le temps d'exécution.


# avec n_sim = 100
# 28 mins

# Avec M binaire
# sde_estimate_os 0.06327047
# bias_sde_os 0.0006863745
# var_sde_os 0.0002774893
# se_sde_os 0.01665801
# sd_bias_sde_os 0.04120387
# mse_sde_os 0.0002779604
# av_estimated_se_sde_os 0.01539581
# cov_sde_os 0.94

# sde_estimate_tmle 0.05303517
# bias_sde_tmle -0.009548932
# var_sde_tmle 0.0003930228
# se_sde_tmle 0.0198248
# sd_bias_sde_tmle -0.4816659
# mse_sde_tmle 0.0004842049
# av_estimated_se_sde_tmle 0.01537876
# cov_sde_tmle 0.81

# sie_estimate_os 0.01103531
# bias_sie_os 0.001189445
# var_sie_os 1.97582e-05
# se_sie_os 0.00444502
# sd_bias_sie_os 0.2675906
# mse_sie_os 2.117298e-05
# av_estimated_se_sie_os 0.004523985
# cov_sie_os 0.99

# sie_estimate_tmle 0.01092913
# bias_sie_tmle 0.001083267
# var_sie_tmle 1.649212e-05
# se_sie_tmle 0.00406105
# sd_bias_sie_tmle 0.2667455
# mse_sie_tmle 1.766559e-05
# av_estimated_se_sie_tmle 0.004350968
# cov_sie_tmle 0.97


# avec n_sim = 200
# 66 minutes

# Avec M binaire
# sde_estimate_os 0.06213466
# bias_sde_os -0.0004494406
# var_sde_os 0.0002702406
# se_sde_os 0.016439
# sd_bias_sde_os -0.0273399
# mse_sde_os 0.0002704426
# av_estimated_se_sde_os 0.01540242
# cov_sde_os 0.945

# sde_estimate_tmle 0.05272032
# bias_sde_tmle -0.00986378
# var_sde_tmle 0.0003822515
# se_sde_tmle 0.01955125
# sd_bias_sde_tmle -0.5045088
# mse_sde_tmle 0.0004795457
# av_estimated_se_sde_tmle 0.01539767
# cov_sde_tmle 0.835

# sie_estimate_os 0.01057567
# bias_sie_os 0.0007298047
# var_sie_os 1.804866e-05
# se_sie_os 0.004248371
# sd_bias_sie_os 0.1717846
# mse_sie_os 1.858127e-05
# av_estimated_se_sie_os 0.004537201
# cov_sie_os 0.965

# sie_estimate_tmle 0.01065447
# bias_sie_tmle 0.000808602
# var_sie_tmle 1.664688e-05
# se_sie_tmle 0.004080059
# sd_bias_sie_tmle 0.1981839
# mse_sie_tmle 1.730072e-05
# av_estimated_se_sie_tmle 0.004373933
# cov_sie_tmle 0.965


# avec lrnr hal, n_sim = 100
# 60 minutes

# Avec M binaire
# sde_estimate_os 0.06326683
# bias_sde_os 0.0006827346
# var_sde_os 0.0002709446
# se_sde_os 0.01646039
# sd_bias_sde_os 0.04147741
# mse_sde_os 0.0002714107
# av_estimated_se_sde_os 0.01535658
# cov_sde_os 0.94

# sde_estimate_tmle 0.05332462
# bias_sde_tmle -0.00925948
# var_sde_tmle 0.0003779163
# se_sde_tmle 0.01944007
# sd_bias_sde_tmle -0.476309
# mse_sde_tmle 0.0004636542
# av_estimated_se_sde_tmle 0.01534995
# cov_sde_tmle 0.86

# sie_estimate_os 0.01106628
# bias_sie_os 0.001220411
# var_sie_os 1.935189e-05
# se_sie_os 0.004399078
# sd_bias_sie_os 0.2774243
# mse_sie_os 2.084129e-05
# av_estimated_se_sie_os 0.004528392
# cov_sie_os 0.91

# sie_estimate_tmle 0.01099387
# bias_sie_tmle 0.001148005
# var_sie_tmle 1.739804e-05
# se_sie_tmle 0.004171096
# sd_bias_sie_tmle 0.2752287
# mse_sie_tmle 1.871596e-05
# av_estimated_se_sie_tmle 0.004365288
# cov_sie_tmle 0.96


# avec lrnr hal, n_sim = 200
# 127 minutes

# Avec M binaire
# sde_estimate_os 0.06211439
# bias_sde_os -0.000469712
# var_sde_os 0.0002727083
# se_sde_os 0.01651388
# sd_bias_sde_os -0.02844346
# mse_sde_os 0.000272929
# av_estimated_se_sde_os 0.01537274
# cov_sde_os 0.935

# sde_estimate_tmle 0.05266297
# bias_sde_tmle -0.009921131
# var_sde_tmle 0.0003772731
# se_sde_tmle 0.01942352
# sd_bias_sde_tmle -0.5107792
# mse_sde_tmle 0.000475702
# av_estimated_se_sde_tmle 0.01536134
# cov_sde_tmle 0.825

# sie_estimate_os 0.01077314
# bias_sie_os 0.0009272753
# var_sie_os 1.741075e-05
# se_sie_os 0.004172619
# sd_bias_sie_os 0.2222286
# mse_sie_os 1.827059e-05
# av_estimated_se_sie_os 0.004543031
# cov_sie_os 0.98

# sie_estimate_tmle 0.01095812
# bias_sie_tmle 0.00111226
# var_sie_tmle 1.616048e-05
# se_sie_tmle 0.00402001
# sd_bias_sie_tmle 0.2766809
# mse_sie_tmle 1.73976e-05
# av_estimated_se_sie_tmle 0.004373798
# cov_sie_tmle 0.97