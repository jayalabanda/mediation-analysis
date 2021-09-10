### IPTW
iptw_direct_indirect <- function(data, inter) {
  g_a <- glm(a ~ 1, family = "binomial", data = data)
  g_a_l0 <- glm(a ~ w1 + w2, family = "binomial", data = data)
  g_m_a <- glm(m ~ a, family = "binomial", data = data)
  g_m_l <- glm(m ~ w1 + w2 + a + z, family = "binomial", data = data)

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
  msm_m <- glm(m ~ a + w1 + w2, family = "binomial", data = data)

  if (inter) {
    msm_y <- glm(y ~ a + m + a:m,
      family = "gaussian",
      data = data, weights = sw
    )

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
  } else {
    msm_y <- glm(y ~ a + m,
      family = "gaussian",
      data = data, weights = sw
    )

    iptw_edn_l0 <- msm_y$coefficients["a"]

    iptw_ein_l0 <- msm_y$coefficients["m"] *
      (plogis(rep(msm_m$coefficients["(Intercept)"], nrow(data)) +
        msm_m$coefficients["a"] +
        msm_m$coefficients["w1"] * data$w1 +
        msm_m$coefficients["w2"] * data$w2) -
        plogis(rep(msm_m$coefficients["(Intercept)"], nrow(data)) +
          msm_m$coefficients["w1"] * data$w1 +
          msm_m$coefficients["w2"] * data$w2))
  }

  iptw_edn <- mean(iptw_edn_l0)
  iptw_ein <- mean(iptw_ein_l0)

  return(list(
    iptw_edn = iptw_edn,
    iptw_ein = iptw_ein
  ))
}

file_path <- "../Data/"
results_path <- "./Results/"
sim_path <- "simulations/"
ests_path <- "estimates/"

true_sde <- 0.064
true_sie <- 0.0112

n_sim <- 750
n_boot <- 200


## Direct effect
estimates_sde_iptw <- matrix(NA, ncol = 3, nrow = n_sim)
colnames(estimates_sde_iptw) <- c("sie_iptw", "sd_iptw", "cov_iptw")
write.csv(
  estimates_sde_iptw,
  file = paste(results_path, ests_path, "estimates_sde_iptw.csv", sep = ""),
  row.names = FALSE
)

## Indirect effect
estimates_sie_iptw <- matrix(NA, ncol = 3, nrow = n_sim)
colnames(estimates_sie_iptw) <- c("sie_iptw", "sd_iptw", "cov_iptw")
write.csv(
  estimates_sie_iptw,
  file = paste(results_path, ests_path, "estimates_sie_iptw.csv", sep = ""),
  row.names = FALSE
)


set.seed(42)
idx <- sample(1:1000, size = n_sim)
sim <- 1

for (i in idx) {
  if (sim %% 10 == 0 | sim == 1) {
    print(paste0("Simulation n°: ", sim))
  }
  sim <- sim + 1
  data_sim <- read.csv(paste0(file_path, sim_path, "data_", i, ".csv"))
  data_sim <- subset(data_sim, select = -c(y_qol))
  colnames(data_sim) <- c("w1", "w2", "a", "z", "m", "y")

  iptw_estimates <- iptw_direct_indirect(
    data = data_sim,
    inter = FALSE
  )

  estimates_sde_iptw <- read.csv(
    paste(results_path, ests_path, "estimates_sde_iptw.csv", sep = "")
  )
  estimates_sde_iptw[i, "sie_iptw"] <- iptw_estimates$iptw_edn

  estimates_sie_iptw <- read.csv(
    paste(results_path, ests_path, "estimates_sie_iptw.csv", sep = "")
  )
  estimates_sie_iptw[i, "sie_iptw"] <- iptw_estimates$iptw_ein

  boot_sde_estimates <- data.frame(matrix(NA, nrow = n_boot, ncol = 1))
  colnames(boot_sde_estimates) <- "boot_sde_estimate"

  boot_sie_estimates <- data.frame(matrix(NA, nrow = n_boot, ncol = 1))
  colnames(boot_sie_estimates) <- "boot_sie_estimate"

  for (b in 1:n_boot) {
    ind <- sample(seq(1, nrow(data_sim)), replace = TRUE)
    boot_data <- data_sim[ind, ]

    if (b %% 100 == 0 & b != 0) {
      print(paste0("Bootstrap number: ", b))
    }

    gcomp_boot <- iptw_direct_indirect(
      data = boot_data,
      inter = FALSE
    )

    boot_sde_estimates[b, ] <- gcomp_boot$iptw_edn
    boot_sie_estimates[b, ] <- gcomp_boot$iptw_ein
  }

  boot_sde_est <- boot_sde_estimates$boot_sde_estimate
  estimates_sde_iptw[i, "sd_iptw"] <- sd(boot_sde_est)
  estimates_sde_iptw[i, "cov_iptw"] <- as.numeric(
    true_sde >= estimates_sde_iptw[i, 1] - 1.96 * sd(boot_sde_est) &
      true_sde <= estimates_sde_iptw[i, 1] + 1.96 * sd(boot_sde_est)
  )

  boot_sie_est <- boot_sie_estimates$boot_sie_estimate
  estimates_sie_iptw[i, "sd_iptw"] <- sd(boot_sie_est)
  estimates_sie_iptw[i, "cov_iptw"] <- as.numeric(
    true_sie >= estimates_sie_iptw[i, 1] - 1.96 * estimates_sie_iptw[i, 2] &
      true_sie <= estimates_sie_iptw[i, 1] + 1.96 * estimates_sie_iptw[i, 2]
  )

  write.csv(estimates_sde_iptw,
    file = paste(results_path, ests_path, "estimates_sde_iptw.csv", sep = ""),
    row.names = FALSE
  )
  write.csv(estimates_sie_iptw,
    file = paste(results_path, ests_path, "estimates_sie_iptw.csv", sep = ""),
    row.names = FALSE
  )
}


estimates_sde_iptw <- read.csv(
  paste(results_path, ests_path, "estimates_sde_iptw.csv", sep = "")
)
estimates_sde_iptw <- data.frame(estimates_sde_iptw)
head(estimates_sde_iptw)

library(tidyr)
library(magrittr)

estimates_sde_iptw %<>%
  drop_na()

estimates_sie_iptw <- read.csv(
  paste(results_path, ests_path, "estimates_sie_iptw.csv", sep = "")
)
estimates_sie_iptw <- data.frame(estimates_sie_iptw)
head(estimates_sie_iptw)

estimates_sie_iptw %<>%
  drop_na()

## Calculate results

# bias
bias_sde <- mean(estimates_sde_iptw$sie_iptw) - true_sde
bias_sie <- mean(estimates_sie_iptw$sie_iptw) - true_sie

# variance & standard error
n_rows <- nrow(data_sim)
var_sde <- mean(
  (estimates_sde_iptw$sie_iptw - mean(estimates_sde_iptw$sie_iptw))^2
) * n_rows / (n_rows - 1)
var_sie <- mean(
  (estimates_sie_iptw$sie_iptw - mean(estimates_sie_iptw$sie_iptw))^2
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
av_estimated_se_sde <- sqrt(mean(estimates_sde_iptw$sd_iptw^2))
av_estimated_se_sie <- sqrt(mean(estimates_sie_iptw$sd_iptw^2))

# Coverage
cov_sde <- mean(estimates_sde_iptw$cov_iptw)
cov_sie <- mean(estimates_sie_iptw$cov_iptw)

# Summarize results
results_sde <- data.frame(
  bias = bias_sde,
  variance = var_sde,
  STD = se_sde,
  std_bias = sd_bias_sde,
  MSE = mse_sde,
  av_est_std = av_estimated_se_sde,
  coverage = cov_sde
)

results_sde
# write.csv(results_sde, paste(results_path, "results/", "results_sde_iptw.csv", sep = ""))

results_sie <- data.frame(
  bias = bias_sie,
  variance = var_sie,
  STD = se_sie,
  std_bias = sd_bias_sie,
  MSE = mse_sie,
  av_est_std = av_estimated_se_sie,
  coverage = cov_sie
)

results_sie
# write.csv(results_sie, paste(results_path, "results/", "results_sie_iptw.csv", sep = ""))

# n_sim = 100
# sde
#          bias     variance        STD  std_bias          MSE av_est_std
# 1 0.002079584 0.0002917651 0.01708113 0.1217474 0.0002960898 0.01552183
#   coverage
# 1     0.92

# sie
#            bias     variance         STD   std_bias          MSE  av_est_std
# 1 -0.0002560235 2.891264e-06 0.001700372 -0.1505691 2.956812e-06 0.001803405
#   coverage
# 1     0.91


# n_sim = 250
# sde
#           bias     variance        STD     std_bias          MSE av_est_std
# 1 -3.99406e-05 0.0002605039 0.01614013 -0.002474614 0.0002605054 0.01539173
#   coverage
# 1    0.936

# sie
#            bias     variance         STD    std_bias          MSE  av_est_std
# 1 -0.0001880437 3.573513e-06 0.001890374 -0.09947432 3.608874e-06 0.001814843
#   coverage
# 1    0.908


# n_sim = 500
# sde
#           bias     variance        STD    std_bias          MSE av_est_std
# 1 -0.000475612 0.0002500413 0.01581269 -0.03007786 0.0002502675 0.01535048
#   coverage
# 1    0.936

# sie
#            bias     variance         STD    std_bias          MSE  av_est_std
# 1 -0.0001618821 3.322718e-06 0.001822832 -0.08880802 3.348923e-06 0.001822781
#   coverage
# 1     0.92


# n_sim = 750
# sde
#            bias     variance        STD    std_bias          MSE av_est_std
# 1 -0.0003881096 0.0002465893 0.01570316 -0.02471538 0.0002467399 0.01536024
#    coverage
# 1 0.9426667

# sie
#            bias     variance         STD    std_bias          MSE  av_est_std
# 1 -0.0001004024 3.310287e-06 0.001819419 -0.05518375 3.320368e-06 0.001829653
#    coverage
# 1 0.9293333


# n_sim = 1000
# sde
#           bias     variance        STD    std_bias          MSE av_est_std
# 1 0.0001060799 0.0002424025 0.01556928 0.006813412 0.0002424137 0.01536758
#   coverage
# 1    0.949

# sie
#            bias    variance         STD   std_bias         MSE  av_est_std
# 1 -0.0001302327 3.25807e-06 0.001805012 -0.0721506 3.27503e-06 0.001828228
#   coverage
# 1    0.936