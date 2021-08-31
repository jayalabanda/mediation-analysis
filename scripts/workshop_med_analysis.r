n <- 1e6
w <- rnorm(n)
a <- rbinom(n, 1, 0.5)
m <- rnorm(n, w + a)
y <- rnorm(n, w + a + m)

lm_y <- lm(y ~ m + a + w)

pred_y1 <- predict(lm_y, newdata = data.frame(a = 1, m = m, w = w))
pred_y2 <- predict(lm_y, newdata = data.frame(a = 0, m = m, w = w))

pseudo <- pred_y1 - pred_y2
lm_pseudo <- lm(pseudo ~ a + w)

pred_pseudo <- predict(lm_pseudo, newdata = data.frame(a = 0, w = w))
mean(pred_pseudo)

n <- 1e6
w <- rnorm(n)
a <- rbinom(n, 1, 0.5)
z <- rbinom(n, 1, 0.5 + 0.2 * a)
m <- rnorm(n, w + a - z)
y <- rnorm(n, w + a + z + m)

lm_y <- lm(y ~ m + a + z + w)
pred_a1z0 <- predict(lm_y, newdata = data.frame(m = m, a = 1, z = 0, w = w))
pred_a1z1 <- predict(lm_y, newdata = data.frame(m = m, a = 1, z = 1, w = w))

prob_z <- lm(z ~ a)
pred_z <- predict(prob_z, newdata = data.frame(a = 1))

pseudo_out <- pred_a1z0 * (1 - pred_z) + pred_a1z1 * pred_z

fit_pseudo <- lm(pseudo_out ~ a + w)
pred_pseudo <- predict(fit_pseudo, newdata = data.frame(a = 0, w = w))

mean(pred_pseudo)

# Stochastic effects
n <- 1e6
w <- rnorm(n)
a <- rbinom(n, 1, plogis(1 + w))
m <- rnorm(n, w + a)
y <- rnorm(n, w + a + m)

fit_y1 <- lm(y ~ m + a + w)
fit_y2 <- lm(y ~ a + w)

pred_y1_a1 <- predict(fit_y1, newdata = data.frame(a = 1, m, w))
pred_y1_a0 <- predict(fit_y1, newdata = data.frame(a = 0, m, w))
pred_y2_a1 <- predict(fit_y2, newdata = data.frame(a = 1, w))
pred_y2_a0 <- predict(fit_y2, newdata = data.frame(a = 0, w))

pseudo_a1 <- pred_y2_a1 - pred_y1_a1
pseudo_a0 <- pred_y2_a0 - pred_y1_a0

pscore_fit <- glm(a ~ w, family = "binomial")
pscore <- predict(pscore_fit, type = "response")
pscore_delta <- 2 * pscore / (2 * pscore + 1 - pscore)

plot(pscore, pscore_delta,
  xlab = "Observed prop. score",
  ylab = "Prop. score under intervention"
)
abline(0, 1)

odds <- (pscore_delta / (1 - pscore_delta)) / (pscore / (1 - pscore))
summary(odds)

indirect <- pseudo_a1 * pscore_delta + pseudo_a0 * (1 - pscore_delta)

mean(indirect)

direct <- (pred_y1_a1 - y) * pscore_delta +
  (pred_y1_a0 - y) * (1 - pscore_delta)

mean(direct)

###############################################
###############################################
# Our data
file_path <- "../Data/"
data <- data.frame(read.csv(paste(file_path, "data_sim.csv", sep = "")))
data <- subset(data, select = -c(Y_qol)) # remove Y_qol
head(data)

# Y_{1, G_0} : A = 1, M = G_0
# Y_{0, G_0} : A = 0, M = G_0
# Y_{1, G_1} : A = 1, M = G_1

# Int. indirect effect : E[Y_{1, G_1} - Y_{1, G_0}]
# Int. direct effect : E[Y_{1, G_0} - Y_{0, G_0}]

library(stringr)

workshop_estimates <- function(data) {
  w_names <- str_subset(colnames(data), "w")
  w <- as.matrix(data[, w_names])
  a <- data$a
  z <- data$z
  m <- data$m
  y <- data$y

  # (a, a') = (1, 0)
  lm_y <- lm(y ~ m + a + z + w)
  pred_a1z0 <- predict(lm_y,
    newdata = data.frame(m = m, a = 1, z = 0, w = w)
  )
  pred_a1z1 <- predict(lm_y,
    newdata = data.frame(m = m, a = 1, z = 1, w = w)
  )

  prob_z <- lm(z ~ a)
  pred_z <- predict(prob_z, newdata = data.frame(a = 1))

  pseudo_out_1_0 <- pred_a1z0 * (1 - pred_z) + pred_a1z1 * pred_z

  fit_pseudo <- lm(pseudo_out_1_0 ~ a + w)
  pred_pseudo_1_0 <- predict(fit_pseudo,
    newdata = data.frame(a = 0, w = w)
  )

  res_1_0 <- mean(pred_pseudo_1_0)

  # (a, a') = (1, 1)
  lm_y <- lm(y ~ m + a + z + w)
  pred_a1z0 <- predict(lm_y,
    newdata = data.frame(m = m, a = 1, z = 0, w = w)
  )
  pred_a1z1 <- predict(lm_y,
    newdata = data.frame(m = m, a = 1, z = 1, w = w)
  )

  prob_z <- lm(z ~ a)
  pred_z <- predict(prob_z, newdata = data.frame(a = 1))

  pseudo_out_1_1 <- pred_a1z0 * (1 - pred_z) + pred_a1z1 * pred_z

  fit_pseudo <- lm(pseudo_out_1_1 ~ a + w)
  pred_pseudo_1_1 <- predict(fit_pseudo,
    newdata = data.frame(a = 1, w = w)
  )

  res_1_1 <- mean(pred_pseudo_1_1)

  # (a, a') = (0, 0)
  lm_y <- lm(y ~ m + a + z + w)
  pred_a0z0 <- predict(lm_y,
    newdata = data.frame(m = m, a = 0, z = 0, w = w)
  )
  pred_a0z1 <- predict(lm_y,
    newdata = data.frame(m = m, a = 0, z = 1, w = w)
  )

  prob_z <- lm(z ~ a)
  pred_z <- predict(prob_z, newdata = data.frame(a = 0))

  pseudo_out_0_0 <- pred_a0z0 * (1 - pred_z) + pred_a0z1 * pred_z

  fit_pseudo <- lm(pseudo_out_0_0 ~ a + w)
  pred_pseudo_0_0 <- predict(fit_pseudo,
    newdata = data.frame(a = 0, w = w)
  )

  res_0_0 <- mean(pred_pseudo_0_0)

  dir_effect <- res_1_0 - res_0_0
  ind_effect <- res_1_1 - res_1_0

  return(c(
    dir_effect,
    ind_effect
  ))
}

colnames(data) <- c("w1", "w2", "a", "z", "m", "y")

# Results
res <- workshop_estimates(data)
res
# $ind_effect
# [1] 0.008698918

# $dir_effect
# [1] 0.06637155

## Bootstrap
n_sim <- 500
n_boot <- 500
true_sde <- 0.0624
true_sie <- 0.0112
sim_data_path <- "../Data/simulations/"

# Direct effect
estimates_sde_ws <- matrix(NA, ncol = 3, nrow = n_sim)
colnames(estimates_sde_ws) <- c("sde", "sd", "cov")
write.csv(
  estimates_sde_ws,
  file = paste(file_path, "estimates_sde_ws.csv", sep = ""),
  row.names = FALSE
)

# Indirect effect
estimates_sie_ws <- matrix(NA, ncol = 3, nrow = n_sim)
colnames(estimates_sie_ws) <- c("sie", "sd", "cov")
write.csv(
  estimates_sie_ws,
  file = paste(file_path, "estimates_sie_ws.csv", sep = ""),
  row.names = FALSE
)

start_time <- Sys.time()
for (i in 1:n_sim) {
  print(paste0("Simulation n°: ", i))
  data_sim <- read.csv(paste0(sim_data_path, "data_", i, ".csv", sep = ""))
  data_sim <- subset(data_sim, select = -c(y_qol))

  colnames(data_sim) <- c("w1", "w2", "a", "z", "m", "y")

  # data_sim$m <- # sapply(data_sim$m, bin_to_quant, simplify = "array")

  results <- workshop_estimates(data_sim)

  estimates_sde_ws <- read.csv(
    paste(file_path, "estimates_sde_ws.csv", sep = "")
  )
  estimates_sde_ws[i, "sde"] <- results$dir_effect

  estimates_sie_ws <- read.csv(
    paste(file_path, "estimates_sie_ws.csv", sep = "")
  )
  estimates_sie_ws[i, "sie"] <- results$ind_effect

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

    ws_boot <- workshop_estimates(boot_data)

    boot_sde_estimates[b, ] <- ws_boot$dir_effect
    boot_sie_estimates[b, ] <- ws_boot$ind_effect
  }

  boot_sde_est <- boot_sde_estimates$boot_sde_estimate
  estimates_sde_ws[i, "sd"] <- sd(boot_sde_est)
  estimates_sde_ws[i, "cov"] <- as.numeric(
    true_sde >= estimates_sde_ws[i, 1] - 1.96 * sd(boot_sde_est) &
      true_sde <= estimates_sde_ws[i, 1] + 1.96 * sd(boot_sde_est)
  )

  boot_sie_est <- boot_sie_estimates$boot_sie_estimate
  estimates_sie_ws[i, "sd"] <- sd(boot_sie_est)
  estimates_sie_ws[i, "cov"] <- as.numeric(
    true_sie >= estimates_sie_ws[i, 1] - 1.96 * sd(boot_sie_est) &
      true_sie <= estimates_sie_ws[i, 1] + 1.96 * sd(boot_sie_est)
  )

  write.csv(estimates_sde_ws,
    file = paste(file_path, "estimates_sde_ws.csv", sep = ""),
    row.names = FALSE
  )
  write.csv(estimates_sie_ws,
    file = paste(file_path, "estimates_sie_ws.csv", sep = ""),
    row.names = FALSE
  )
}
end_time <- Sys.time()
diff <- (end_time - start_time) / 60

estimates_sde_ws <- read.csv(
  paste(file_path, "estimates_sde_ws.csv", sep = "")
)
estimates_sde_ws <- data.frame(estimates_sde_ws)
head(estimates_sde_ws)

estimates_sie_ws <- read.csv(
  paste(file_path, "estimates_sie_ws.csv", sep = "")
)
estimates_sie_ws <- data.frame(estimates_sie_ws)
head(estimates_sie_ws)

## Calculate results

# estimate
sde_estimate <- mean(estimates_sde_ws$sde)
sie_estimate <- mean(estimates_sie_ws$sie)

# bias
bias_sde <- mean(estimates_sde_ws$sde) - true_sde
bias_sie <- mean(estimates_sie_ws$sie) - true_sie

# variance & standard error
n_rows <- nrow(data_sim)
var_sde <- mean(
  (estimates_sde_ws$sde - mean(estimates_sde_ws$sde))^2
) * n_rows / (n_rows - 1)
var_sie <- mean(
  (estimates_sie_ws$sie - mean(estimates_sie_ws$sie))^2
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
av_estimated_se_sde <- sqrt(mean(estimates_sde_ws$sd^2))
av_estimated_se_sie <- sqrt(mean(estimates_sie_ws$sd^2))

# Coverage
cov_sde <- mean(estimates_sde_ws$cov)
cov_sie <- mean(estimates_sie_ws$cov)

# Summarize results
results_sde <- data.frame(
  sde_estimate = sde_estimate, # 0.06379728
  bias = bias_sde, # 0.001213182
  variance = var_sde, # 0.0002157029
  STD = se_sde, # 0.01468683
  std.bias = sd_bias_sde, # 0.0826034
  MSE = mse_sde, # 0.0002171747
  av.est.std = av_estimated_se_sde, # 0.01447822
  coverage = cov_sde # 0.944
)

results_sde
write.csv(
  results_sde,
  paste(file_path, "results_sde_ws.csv", sep = ""),
  row.names = FALSE
)

results_sie <- data.frame(
  sie_estimate = sie_estimate, # 0.01119133
  bias = bias_sie, # 0.001345465
  variance = var_sie, # 3.298703e-06
  STD = se_sie, # 0.001816233
  std.bias = sd_bias_sie, # 0.7407997
  MSE = mse_sie, # 5.108979e-06
  av.est.std = av_estimated_se_sie, # 0.001838103
  coverage = cov_sie # 0.914
)

results_sie
write.csv(
  results_sie,
  paste(file_path, "results_sie_ws.csv", sep = ""),
  row.names = FALSE
)

# Pas de différence remarquable entre les deux méthodes
# que ce soit avec M binaire ou M continu sauf pour coverage

# n_sim = 100, n_boot = 500
# 38 minutes

# Avec M continu
#   sde_estimate        bias     variance       STD  std.bias         MSE
# 1    0.0672141 0.004630001 0.0002363751 0.0153745 0.3011481 0.000257812
#   av.est.std coverage
# 1 0.01447884     0.91

#   sie_estimate         bias     variance         STD  std.bias          MSE
# 1  0.008088085 -0.001757779 2.249273e-06 0.001499758 -1.172042 5.339059e-06
#    av.est.std coverage
# 1 0.001528076     0.75


# n_sim = 200, n_boot = 500
# 72 minutes

# Avec M continu
#   sde_estimate        bias     variance        STD  std.bias          MSE
# 1   0.06575144 0.003167339 0.0002426105 0.01557596 0.2033479 0.0002526426
#   av.est.std coverage
# 1 0.01439175    0.905

#   sie_estimate         bias     variance         STD  std.bias          MSE
# 1  0.008128275 -0.001717589 2.271654e-06 0.001507201 -1.139589 5.221766e-06
#    av.est.std coverage
# 1 0.001545938    0.735


# n_sim = 500, n_boot = 500
# 2h43 minutes

# Avec M continu
#   sde_estimate        bias     variance        STD  std.bias          MSE
# 1   0.06663224 0.004048135 0.0002292207 0.01514004 0.2673795 0.0002456081
#   av.est.std coverage
# 1  0.0144567    0.924

#   sie_estimate         bias     variance         STD std.bias          MSE
# 1  0.008059886 -0.001785978 2.142665e-06 0.001463785 -1.22011 5.332382e-06
#    av.est.std coverage
# 1 0.001531339    0.758


# n_sim = 100, n_boot = 500
# 40 minutes

# Avec M binaire
#   sde_estimate        bias     variance        STD   std.bias          MSE
# 1   0.06410054 0.001516443 0.0002363297 0.01537302 0.09864314 0.0002386293
#   av.est.std coverage
# 1 0.01453757     0.93

#   sie_estimate        bias     variance         STD  std.bias          MSE
# 1   0.01117657 0.001330706 3.260057e-06 0.001805563 0.7370037 5.030837e-06
#    av.est.std coverage
# 1 0.001836964     0.91


# n_sim = 200, n_boot = 500
# 66 minutes

# Avec M binaire
#   sde_estimate         bias     variance       STD    std.bias          MSE
# 1   0.06262966 4.555886e-05 0.0002472692 0.0157248 0.002897262 0.0002472713
#   av.est.std coverage
# 1 0.01453101    0.945

#   sie_estimate        bias     variance         STD  std.bias          MSE
# 1   0.01122455 0.001378684 3.327952e-06 0.001824268 0.7557467 5.228722e-06
#    av.est.std coverage
# 1 0.001837334     0.91


# n_sim = 500, n_boot = 500
# 2h53 minutes

# Avec M binaire
#   sde_estimate         bias     variance        STD   std.bias          MSE
# 1   0.06347579 0.0008916896 0.0002292661 0.01514153 0.05889031 0.0002300612
#   av.est.std coverage
# 1  0.0144751    0.938

#   sie_estimate        bias     variance         STD  std.bias         MSE
# 1   0.01118948 0.001343612 3.206816e-06 0.001790758 0.7503035 5.01211e-06
#    av.est.std coverage
# 1 0.001837824     0.91



#### Biases
file_path <- "../Data/"
sim_path <- "new_simulations/"

true_sde <- 0.064
true_sie <- 0.0112

set.seed(42)
n_sim <- 500
idx <- sample(1:1000, size = n_sim)

bias_estimates <- matrix(nrow = n_sim, ncol = 2)
colnames(bias_estimates) <- c("SDE", "SIE")

for (i in 1:n_sim) {
  if (i %% 10 == 0 | i == 1) {
    print(paste0("Simulation ", i, " of ", n_sim))
  }

  data <- read.csv(paste0(file_path, sim_path, "data_", idx[i], ".csv", sep = ""))
  if (!str_detect(sim_path, "rudolph")) {
    data <- subset(data, select = -y_qol)
    colnames(data) <- c("w1", "w2", "a", "z", "m", "y")

    res <- workshop_estimates(data)
    bias_estimates[i, 1] <- res[1] - true_sde
    bias_estimates[i, 2] <- res[2] - true_sie
  } else {
    res <- workshop_estimates(data)
    bias_estimates[i, 1] <- res[1]
    bias_estimates[i, 2] <- res[2]
  }
}

boxplot(bias_estimates)
abline(h = 0, col = ifelse(!str_detect(sim_path, "rudolph"), "black", NA))

mean(bias_estimates[, 1])
mean(bias_estimates[, 2])

# [1] -0.0002119083
# [1] 6.550189e-05

# Rudolph
# [1] 0.09879082
# [1] 0.03257756