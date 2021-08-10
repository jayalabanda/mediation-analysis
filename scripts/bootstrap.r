sim_param_time_varying_l <- function(a_m_interaction = NULL) {
  # L0
  p_l0_male <- 0.5
  p_l0_parent_low_educ_lv <- 0.65

  b_a <- 0.05 # reference prevalence is 5%
  b_male_a <- 0.04 # + 0.04 for the effect of l0_male -> a0_ace
  b_parent_educ_a <- 0.06 # +0.06 for effect of l0_parent_low_educ_lv -> a0_ace

  b_l1 <- 0.30 # reference prevalence is 30%
  b_male_l1 <- -0.05 # - 0.05 for the effect of l0_male -> l1
  b_parent_l1 <- +0.08 # + 0.08 for the effect of l0_parent_low_educ_lv -> l1
  b_a_l1 <- +0.2 # +0.2 for the effect of a0_ace -> l1

  b_m <- 0.2 # reference prevalence is 20%
  b_male_m <- 0.05 # +0.05 for the effect of l0_male -> m_smoking
  b_parent_educ_m <- 0.06 # +0.06 effect of l0_parent_low_educ_lv -> m_smoking
  b_a_m <- 0.1 # +0.10 for the effect of a0_ace -> m_smoking
  b_l1_m <- 0.2 # +0.2 for the effect of l1 -> m_smoking

  b_y <- 0.1 # reference prevalence is 10%
  b_male_y <- 0.06 # +0.06 for the effect of l0_male -> Y
  b_parent_educ_y <- 0.04 # +0.04 for the effect of l0_parent_low_educ_lv -> Y
  b_a_y <- 0.05 # 0.05 for the effect of a0_ace -> Y
  b_l1_y <- 0.07 # +0.07 for the effect of l1 -> Y
  b_m_y <- 0.08 # 0.08 for the effect of m_smoking -> Y
  b_am_y <- 0.03 # 0.03 for the interaction effect a0_ace * m_smoking -> Y

  mu_y <- 75 # reference mean for QoL
  c_male_y <- -1 # -1 for the effect of l0_male -> Y
  c_parent_educ_y <- -3 # -3 for the effect of l0_parent_low_educ_lv -> Y
  c_a_y <- -4 # -4 for the effect of a0_ace -> Y
  c_l1_y <- -5 # -5 for the effect of l1 -> Y
  c_m_y <- -9 # -9 for the effect of m_smoking -> Y
  c_am_y <- -5 # - 5 for the interaction effect a0_ace * m_smoking  -> Y
  sd_y <- 10 # standard deviation of the residuals

  # A*M interaction ?
  a_m_inter <- a_m_interaction

  coef <- c(
    p_l0_male = p_l0_male, p_l0_parent_low_educ_lv = p_l0_parent_low_educ_lv,
    b_a = b_a, b_male_a = b_male_a, b_parent_educ_a = b_parent_educ_a,
    b_l1 = b_l1, b_male_l1 = b_male_l1, b_parent_l1 = b_parent_l1,
    b_a_l1 = b_a_l1, b_m = b_m, b_male_m = b_male_m,
    b_parent_educ_m = b_parent_educ_m, b_l1_m = b_l1_m, b_a_m = b_a_m,
    b_y = b_y, b_male_y = b_male_y, b_parent_educ_y = b_parent_educ_y,
    b_a_y = b_a_y, b_l1_y = b_l1_y, b_m_y = b_m_y, b_am_y = b_am_y,
    mu_y = mu_y, c_male_y = c_male_y, c_parent_educ_y = c_parent_educ_y,
    c_a_y = c_a_y, c_l1_y = c_l1_y, c_m_y = c_m_y, c_am_y = c_am_y,
    sd_y = sd_y, a_m_inter = a_m_inter
  )

  return(coef)
}

gen_data_time_varying_l <- function(N, a_m_inter) {
  # input parameters: sample size N and presence of A*M interaction

  b <- sim_param_time_varying_l(a_m_interaction = a_m_inter)

  # baseline confounders: l0_parent_low_educ_lv & l0_male
  l0_male <- rbinom(N, size = 1, prob = b["p_l0_male"])
  l0_parent_low_educ_lv <- rbinom(N,
    size = 1,
    prob = b["p_l0_parent_low_educ_lv"]
  )

  # exposure: a0_ace
  a0_ace <- rbinom(N, size = 1, prob = b["b_a"] +
    b["b_male_a"] * l0_male +
    b["b_parent_educ_a"] * l0_parent_low_educ_lv)

  # intermediate confounder between m_smoking and Y, not affected by A0 l1
  l1 <- rbinom(N, size = 1, prob = b["b_l1"] +
    b["b_male_l1"] * l0_male +
    b["b_parent_l1"] * l0_parent_low_educ_lv +
    b["b_a_l1"] * a0_ace)

  # mediator: m_smoking
  m_smoking <- rbinom(N, size = 1, prob = b["b_m"] +
    b["b_male_m"] * l0_male +
    b["b_parent_educ_m"] * l0_parent_low_educ_lv +
    b["b_a_m"] * a0_ace +
    b["b_l1_m"] * l1)

  # y_death
  y_death <- rbinom(N, size = 1, prob = b["b_y"] +
    b["b_male_y"] * l0_male +
    b["b_parent_educ_y"] * l0_parent_low_educ_lv +
    b["b_a_y"] * a0_ace +
    b["b_l1_y"] * l1 +
    b["b_m_y"] * m_smoking +
    b["b_am_y"] * a0_ace * m_smoking * a_m_inter)

  # y_qol
  y_qol <- (b["mu_y"] +
    b["c_male_y"] * l0_male +
    b["c_parent_educ_y"] * l0_parent_low_educ_lv +
    b["c_a_y"] * a0_ace +
    b["c_l1_y"] * l1 +
    b["c_m_y"] * m_smoking +
    b["c_am_y"] * a0_ace * m_smoking * a_m_inter) +
    rnorm(N, mean = 0, sd = b["sd_y"])

  # data.frame
  data_sim <- data.frame(
    l0_male, l0_parent_low_educ_lv, a0_ace,
    l1, m_smoking, y_death, y_qol
  )

  return(data_sim)
}

estimate_manual <- function(data, ymodel, mmodel, qmodel, A, M) {
  tempdat <- data
  colnames(tempdat)[grep(A, colnames(tempdat))] <- "a"
  colnames(tempdat)[grep(M, colnames(tempdat))] <- "m"

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
    formula = paste("q_gamma_a1", qmodel, sep = "~"),
    family = "quasibinomial", data = tempdat
  )
  q_model_a0 <- glm(
    formula = paste("q_gamma_a0", qmodel, sep = "~"),
    family = "quasibinomial", data = tempdat
  )

  q_pred_a1_gamma_a1 <- predict(q_model_a1,
    newdata = data_a1, type = "response"
  )
  q_pred_a1_gamma_a0 <- predict(q_model_a0,
    newdata = data_a1, type = "response"
  )
  q_pred_a0_gamma_a0 <- predict(q_model_a0,
    newdata = data_a0, type = "response"
  )

  psi_mr_nde <- mean(q_pred_a1_gamma_a0) - mean(q_pred_a0_gamma_a0)
  psi_mr_nie <- mean(q_pred_a1_gamma_a1) - mean(q_pred_a1_gamma_a0)

  return(list(
    psi_mr_nde = psi_mr_nde,
    psi_mr_nie = psi_mr_nie
  ))
}

create_data_table <- function(data, covars, l1, A, M, outcome) {
  # binary outcome
  require(data.table)
  len <- nrow(data)
  tempdat <- data
  colnames(tempdat)[grep(A, colnames(tempdat))] <- "a"
  colnames(tempdat)[grep(M, colnames(tempdat))] <- "m"

  data_stremr <- data.frame(tempdat)
  data_stremr$ID <- 1:len
  data_stremr$t0 <- rep(0, len)
  data_stremr$t1 <- rep(1, len)
  data_stremr$alive <- rep(0, len)

  new_df <- data.frame(data[, covars])

  data_dt_0 <- data.frame(
    ID = data_stremr$ID,
    t = data_stremr$t0,
    l1 = rep(NA, len),
    A.t = data_stremr$a,
    Y.tplus1 = data_stremr$alive,
    A.tminus1 = rep(NA, len)
  )

  data_dt_1 <- data.frame(
    ID = data_stremr$ID,
    t = data_stremr$t1,
    l1 = data_stremr[, l1],
    A.t = data_stremr$m,
    Y.tplus1 = data_stremr[, outcome],
    A.tminus1 = data_stremr$a
  )

  data_dt_0 <- cbind(data_dt_0, new_df)
  data_dt_1 <- cbind(data_dt_1, new_df)

  # To data.table
  data_dt <- rbind(data_dt_0, data_dt_1)
  data_dt <- data.table::data.table(data_dt)
  data_dt <- data_dt[order(ID, t)]

  data_dt[, ("p_a1_gamma_a0") := NA]
  data_dt[, ("p_a0_gamma_a0") := NA]
  data_dt[, ("p_a1_gamma_a1") := NA]

  return(data_dt)
}

estimate_stremr <- function(data, mmodel, covars, l1, A, M, outcome, q_forms) {
  require(stremr)
  tempdat <- create_data_table(data, covars, l1, A, M, outcome)

  colnames(data)[grep(A, colnames(data))] <- "a"
  colnames(data)[grep(M, colnames(data))] <- "m"

  g_m_model <- glm(formula = mmodel, family = "binomial", data = data)

  # Exposition à A = 1
  data_a1 <- data
  data_a1$a <- 1
  g_m_1 <- predict(g_m_model, newdata = data_a1, type = "response")

  # Exposition à A = 0
  data_a0 <- data
  data_a0$a <- 0
  g_m_0 <- predict(g_m_model, newdata = data_a0, type = "response")

  tempdat$p_a1_gamma_a0[tempdat$t == 0] <- rep(1, nrow(data))
  tempdat$p_a1_gamma_a0[tempdat$t == 1] <- g_m_0

  tempdat$p_a0_gamma_a0[tempdat$t == 0] <- rep(0, nrow(data))
  tempdat$p_a0_gamma_a0[tempdat$t == 1] <- g_m_0

  tempdat$p_a1_gamma_a1[tempdat$t == 0] <- rep(1, nrow(data))
  tempdat$p_a1_gamma_a1[tempdat$t == 1] <- g_m_1

  # Import data for stremr
  o_data <- stremr::importData(
    tempdat,
    ID = "ID",
    t_name = "t",
    covars = append(covars, l1),
    TRT = "A.t",
    OUTCOME = "Y.tplus1",
    CENS = NULL, MONITOR = NULL
  )

  # g-computation
  # E(Y_{a=1, \Gamma_a=0})
  gcomp_y_a1_g0 <- stremr::fit_GCOMP(
    o_data,
    tvals = c(0:1),
    intervened_TRT = "p_a1_gamma_a0",
    Qforms = q_forms
  )

  # E(Y_{a=0, \Gamma_a=0})
  gcomp_y_a0_g0 <- stremr::fit_GCOMP(
    o_data,
    tvals = c(0:1),
    intervened_TRT = "p_a0_gamma_a0",
    Qforms = q_forms
  )

  # E(Y_{a=1, \Gamma_a=1})
  gcomp_y_a1_g1 <- stremr::fit_GCOMP(
    o_data,
    tvals = c(0:1),
    intervened_TRT = "p_a1_gamma_a1",
    Qforms = q_forms
  )

  psi_mr_nde_gcomp <-
    gcomp_y_a1_g0$estimates$cum.inc[2] - gcomp_y_a0_g0$estimates$cum.inc[2]
  psi_mr_nie_gcomp <-
    gcomp_y_a1_g1$estimates$cum.inc[2] - gcomp_y_a1_g0$estimates$cum.inc[2]

  return(list(
    psi_mr_nde_gcomp = psi_mr_nde_gcomp,
    psi_mr_nie_gcomp = psi_mr_nie_gcomp
  ))
}

###### Bootstrap

file_path <- "../Data/simulations/"

n_sim <- 100
n_boot <- 200
true_sde <- 0.0625841
true_sie <- 0.009845864

q_forms <- c(
  "Qkplus1 ~ l0_male + l0_parent_low_educ_lv + A.t",
  "Qkplus1 ~ l0_male + l0_parent_low_educ_lv + A.t + l1 + A.tminus1"
)

ymodel <- "y_death ~ l0_male + l0_parent_low_educ_lv + a + l1 + m"
mmodel <- "m ~ l0_male + l0_parent_low_educ_lv + a"
qmodel <- "l0_male + l0_parent_low_educ_lv + a"

## Direct effect
estimates_sde_gcomp <- matrix(NA, ncol = 3, nrow = n_sim)
colnames(estimates_sde_gcomp) <- c("sde_gcomp", "sd_gcomp", "cov_gcomp")
write.csv(
  estimates_sde_gcomp,
  file = paste(file_path, "estimates_sde_gcomp.csv", sep = ""),
  row.names = FALSE
)

## Indirect effect
estimates_sie_gcomp <- matrix(NA, ncol = 3, nrow = n_sim)
colnames(estimates_sie_gcomp) <- c("sie_gcomp", "sd_gcomp", "cov_gcomp")
write.csv(
  estimates_sie_gcomp,
  file = paste(file_path, "estimates_sie_gcomp.csv", sep = ""),
  row.names = FALSE
)

## Both using stremr
estimates_stremr <- matrix(NA, ncol = 6, nrow = n_sim)
colnames(estimates_stremr) <- c(
  "sde_gcomp", "sd_dir_gcomp", "cov_dir_gcomp",
  "sie_gcomp", "sd_ind_gcomp", "cov_ind_gcomp"
)
write.csv(
  estimates_stremr,
  file = paste(file_path, "estimates_stremr.csv", sep = ""),
  row.names = FALSE
)

set.seed(42)
idx <- sample(1:1000, size = n_sim)
sim <- 1

start_time <- Sys.time()
for (i in idx) {
  print(paste0("Simulation n°: ", sim))
  sim <- sim + 1
  data_sim <- read.csv(paste0(file_path, "data_", i, ".csv", sep = ""))
  data_sim <- subset(data_sim, select = -c(y_qol))

  gcomp_estimates <- estimate_manual(
    data = data_sim,
    ymodel = ymodel,
    mmodel = mmodel,
    qmodel = qmodel,
    A = "a0_ace",
    M = "m_smoking"
  )

  stremr_estimates <- estimate_stremr(
    data = data_sim, mmodel = mmodel,
    covars = c("l0_male", "l0_parent_low_educ_lv"),
    l1 = "l1", A = "a0_ace", M = "m_smoking",
    outcome = "y_death", q_forms = q_forms
  )

  estimates_sde_gcomp <- read.csv(
    paste(file_path, "estimates_sde_gcomp.csv", sep = "")
  )
  estimates_sde_gcomp[i, "sde_gcomp"] <- gcomp_estimates$psi_mr_nde

  estimates_sie_gcomp <- read.csv(
    paste(file_path, "estimates_sie_gcomp.csv", sep = "")
  )
  estimates_sie_gcomp[i, "sie_gcomp"] <- gcomp_estimates$psi_mr_nie

  estimates_stremr <- read.csv(
    paste(file_path, "estimates_stremr.csv", sep = "")
  )
  estimates_stremr[i, "sde_gcomp"] <- stremr_estimates$psi_mr_nde_gcomp
  estimates_stremr[i, "sie_gcomp"] <- stremr_estimates$psi_mr_nie_gcomp

  boot_sde_estimates <- data.frame(matrix(NA, nrow = n_boot, ncol = 1))
  colnames(boot_sde_estimates) <- "boot_sde_estimate"

  boot_sie_estimates <- data.frame(matrix(NA, nrow = n_boot, ncol = 1))
  colnames(boot_sie_estimates) <- "boot_sie_estimate"

  boot_stremr_estimates <- data.frame(matrix(NA, nrow = n_boot, ncol = 2))
  colnames(boot_stremr_estimates) <- c("boot_sde_estimate", "boot_sie_estimate")

  for (b in 1:n_boot) {
    idx <- sample(seq(1, nrow(data_sim)), replace = TRUE)
    boot_data <- data_sim[idx, ]

    if (b %% 100 == 0 & b != 0) {
      print(paste0("Bootstrap number: ", b))
    }

    gcomp_boot <- estimate_manual(
      data = boot_data,
      ymodel = ymodel,
      mmodel = mmodel,
      qmodel = qmodel,
      A = "a0_ace",
      M = "m_smoking"
    )

    stremr_boot <- estimate_stremr(
      data = boot_data, mmodel = mmodel,
      covars = c("l0_male", "l0_parent_low_educ_lv"),
      l1 = "l1", A = "a0_ace", M = "m_smoking",
      outcome = "y_death", q_forms = q_forms
    )

    boot_sde_estimates[b, ] <- gcomp_boot$psi_mr_nde
    boot_sie_estimates[b, ] <- gcomp_boot$psi_mr_nie
    boot_stremr_estimates[b, 1] <- stremr_boot$psi_mr_nde_gcomp
    boot_stremr_estimates[b, 2] <- stremr_boot$psi_mr_nie_gcomp
  }

  boot_sde_est <- boot_sde_estimates$boot_sde_estimate
  estimates_sde_gcomp[i, "sd_gcomp"] <- sd(boot_sde_est)
  estimates_sde_gcomp[i, "cov_gcomp"] <- as.numeric(
    true_sde >= estimates_sde_gcomp[i, 1] - 1.96 * sd(boot_sde_est) &
      true_sde <= estimates_sde_gcomp[i, 1] + 1.96 * sd(boot_sde_est)
  )

  boot_sie_est <- boot_sie_estimates$boot_sie_estimate
  estimates_sie_gcomp[i, "sd_gcomp"] <- sd(boot_sie_est)
  estimates_sie_gcomp[i, "cov_gcomp"] <- as.numeric(
    true_sie >= estimates_sie_gcomp[i, 1] - 1.96 * estimates_sie_gcomp[i, 2] &
      true_sie <= estimates_sie_gcomp[i, 1] + 1.96 * estimates_sie_gcomp[i, 2]
  )

  boot_sde_est <- boot_stremr_estimates$boot_sde_estimate
  boot_sie_est <- boot_stremr_estimates$boot_sie_estimate
  estimates_stremr[i, "sd_dir_gcomp"] <- sd(boot_sde_est)
  estimates_stremr[i, "sd_ind_gcomp"] <- sd(boot_sie_est)

  estimates_stremr[i, "cov_dir_gcomp"] <- as.numeric(
    true_sde >= estimates_stremr[i, 1] - 1.96 * estimates_stremr[i, 2] &
      true_sde <= estimates_stremr[i, 1] + 1.96 * estimates_stremr[i, 2]
  )
  estimates_stremr[i, "cov_ind_gcomp"] <- as.numeric(
    true_sie >= estimates_stremr[i, 4] - 1.96 * estimates_stremr[i, 5] &
      true_sie <= estimates_stremr[i, 4] + 1.96 * estimates_stremr[i, 5]
  )

  write.csv(estimates_sde_gcomp,
    file = paste(file_path, "estimates_sde_gcomp.csv", sep = ""),
    row.names = FALSE
  )
  write.csv(estimates_sie_gcomp,
    file = paste(file_path, "estimates_sie_gcomp.csv", sep = ""),
    row.names = FALSE
  )
  write.csv(estimates_stremr,
    file = paste(file_path, "estimates_stremr.csv", sep = ""),
    row.names = FALSE
  )
}
end_time <- Sys.time()
diff <- end_time - start_time
diff

estimates_sde_gcomp <- read.csv(
  paste(file_path, "estimates_sde_gcomp.csv", sep = "")
)
estimates_sde_gcomp <- data.frame(estimates_sde_gcomp)
head(estimates_sde_gcomp)

estimates_sie_gcomp <- read.csv(
  paste(file_path, "estimates_sie_gcomp.csv", sep = "")
)
estimates_sie_gcomp <- data.frame(estimates_sie_gcomp)
head(estimates_sie_gcomp)

estimates_stremr <- read.csv(
  paste(file_path, "estimates_stremr.csv", sep = "")
)
estimates_stremr <- data.frame(estimates_stremr)
head(estimates_stremr)

## Calculate results

# bias
bias_sde <- mean(estimates_sde_gcomp$sde_gcomp) - true_sde
bias_sie <- mean(estimates_sie_gcomp$sie_gcomp) - true_sie
bias_sde_stremr <- mean(estimates_stremr$sde_gcomp) - true_sde
bias_sie_stremr <- mean(estimates_stremr$sie_gcomp) - true_sie

# variance & standard error
n_rows <- nrow(data_sim)
var_sde <- mean(
  (estimates_sde_gcomp$sde_gcomp - mean(estimates_sde_gcomp$sde_gcomp))^2
) * n_rows / (n_rows - 1)
var_sie <- mean(
  (estimates_sie_gcomp$sie_gcomp - mean(estimates_sie_gcomp$sie_gcomp))^2
) * n_rows / (n_rows - 1)
var_sde_stremr <- mean(
  (estimates_stremr$sde_gcomp - mean(estimates_stremr$sde_gcomp))^2
) * n_rows / (n_rows - 1)
var_sie_stremr <- mean(
  (estimates_stremr$sie_gcomp - mean(estimates_stremr$sie_gcomp))^2
) * n_rows / (n_rows - 1)

se_sde <- sqrt(var_sde)
se_sie <- sqrt(var_sie)
se_sde_stremr <- sqrt(var_sde_stremr)
se_sie_stremr <- sqrt(var_sie_stremr)

# Standardized bias
sd_bias_sde <- bias_sde / se_sde
sd_bias_sie <- bias_sie / se_sie
sd_bias_sde_stremr <- bias_sde_stremr / se_sde_stremr
sd_bias_sie_stremr <- bias_sie_stremr / se_sie_stremr

# MSE
mse_sde <- var_sde + bias_sde^2
mse_sie <- var_sie + bias_sie^2
mse_sde_stremr <- var_sde_stremr + bias_sde_stremr^2
mse_sie_stremr <- var_sie_stremr + bias_sie_stremr^2

# Average estimated standard error
av_estimated_se_sde <- sqrt(mean(estimates_sde_gcomp$sd_gcomp^2))
av_estimated_se_sie <- sqrt(mean(estimates_sie_gcomp$sd_gcomp^2))
av_estimated_se_sde_stremr <- sqrt(mean(estimates_stremr$sd_dir_gcomp^2))
av_estimated_se_sie_stremr <- sqrt(mean(estimates_stremr$sd_ind_gcomp^2))

# Coverage
cov_sde <- mean(estimates_sde_gcomp$cov_gcomp)
cov_sie <- mean(estimates_sie_gcomp$cov_gcomp)
cov_sde_stremr <- mean(estimates_stremr$cov_dir_gcomp)
cov_sie_stremr <- mean(estimates_stremr$cov_ind_gcomp)

# Summarize results
results_sde <- data.frame(
  bias = bias_sde, # -0.002334099
  variance = var_sde, # 0.0001804081
  STD = se_sde, # 0.01343161
  std.bias = sd_bias_sde, # -0.1737766
  MSE = mse_sde, # 0.0001858561
  av.est.std = av_estimated_se_sde, # 0.01383005
  coverage = cov_sde # 0.96
)

results_sde
write.csv(results_sde, paste(file_path, "results_sde.csv", sep = ""))

results_sie <- data.frame(
  bias = bias_sie, # 0.002782764
  variance = var_sie, # 3.847359e-06
  STD = se_sie, # 0.001961469
  std.bias = sd_bias_sie, # 1.418714
  MSE = mse_sie, # 1.159113e-05
  av.est.std = av_estimated_se_sie, # 0.002075798
  coverage = cov_sie # 0.77
)

results_sie
write.csv(results_sie, paste(file_path, "results_sie.csv", sep = ""))

results_stremr_sde <- data.frame(
  bias = bias_sde_stremr, # -0.002334099
  variance = var_sde_stremr, # 0.0001804081
  STD = se_sde_stremr, # 0.01343161
  std.bias = sd_bias_sde_stremr, # -0.1737766
  MSE = mse_sde_stremr, # 0.0001858561
  av.est.std = av_estimated_se_sde_stremr, # 0.01383005
  coverage = cov_sde_stremr # 0.96
)

results_stremr_sde
write.csv(
  results_stremr_sde,
  paste(file_path, "results_stremr_sde.csv", sep = "")
)

results_stremr_sie <- data.frame(
  bias = bias_sie_stremr, # 0.002782764
  variance = var_sie_stremr, # 3.847359e-06
  STD = se_sie_stremr, # 0.001961469
  std.bias = sd_bias_sie_stremr, # 1.418714
  MSE = mse_sie_stremr, # 1.159113e-05
  av.est.std = av_estimated_se_sie_stremr, # 0.002075798
  coverage = cov_sie_stremr # 0.77
)

results_stremr_sie
write.csv(
  results_stremr_sie,
  paste(file_path, "results_stremr_sie.csv", sep = "")
)