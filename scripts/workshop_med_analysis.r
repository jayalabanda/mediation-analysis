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

# Our data
file_path <- "../Data/"
data <- data.frame(read.csv(paste(file_path, "data_sim.csv", sep = "")))
data <- subset(data, select = -c(Y_qol)) # remove Y_qol
head(data)

colnames(data)[1] <- "w1"
colnames(data)[2] <- "w2"
colnames(data)[3] <- "a"
colnames(data)[4] <- "z"
colnames(data)[5] <- "m"
colnames(data)[6] <- "y"
head(data)

# function to convert
# binary vector to quantitative
bin_to_quant <- function(x) {
  x[x == 0] <- runif(1, min = 0, max = 0.5)
  x[x == 1] <- runif(1, min = 0.5, max = 1)
  return(x)
}

w1 <- data$w1
w2 <- data$w2
a <- data$a
z <- data$z
m <- data$m
y <- data$y

# Y_{1, G_0} : A = 1, M = G_0
# Y_{0, G_0} : A = 0, M = G_0
# Y_{1, G_1} : A = 1, M = G_1

# Int. indirect effect : E[Y_{1, G_1} - Y_{1, G_0}]
# Int. direct effect : E[Y_{1, G_0} - Y_{0, G_0}]

workshop_estimates <- function(data) {
  w1 <- data$w1
  w2 <- data$w2
  a <- data$a
  z <- data$z
  m <- data$m
  # m <- # lapply(data$m, bin_to_quant)
  y <- data$y

  # (a, a') = (1, 0)
  lm_y <- lm(y ~ m + a + z + w1 + w2)
  pred_a1z0 <- predict(lm_y,
    newdata = data.frame(m = m, a = 1, z = 0, w1 = w1, w2 = w2)
  )
  pred_a1z1 <- predict(lm_y,
    newdata = data.frame(m = m, a = 1, z = 1, w1 = w1, w2 = w2)
  )

  prob_z <- lm(z ~ a)
  pred_z <- predict(prob_z, newdata = data.frame(a = 1))

  pseudo_out <- pred_a1z0 * (1 - pred_z) + pred_a1z1 * pred_z

  fit_pseudo <- lm(pseudo_out ~ a + w1 + w2)
  pred_pseudo_1_0 <- predict(fit_pseudo,
    newdata = data.frame(a = 0, w1 = w1, w2 = w2)
  )

  res_1_0 <- mean(pred_pseudo_1_0)

  # (a, a') = (1, 1)
  lm_y <- lm(y ~ m + a + z + w1 + w2)
  pred_a1z0 <- predict(lm_y,
    newdata = data.frame(m = m, a = 1, z = 0, w1 = w1, w2 = w2)
  )
  pred_a1z1 <- predict(lm_y,
    newdata = data.frame(m = m, a = 1, z = 1, w1 = w1, w2 = w2)
  )

  prob_z <- lm(z ~ a)
  pred_z <- predict(prob_z, newdata = data.frame(a = 1))

  pseudo_out <- pred_a1z0 * (1 - pred_z) + pred_a1z1 * pred_z

  fit_pseudo <- lm(pseudo_out ~ a + w1 + w2)
  pred_pseudo_1_1 <- predict(fit_pseudo,
    newdata = data.frame(a = 1, w1 = w1, w2 = w2)
  )

  res_1_1 <- mean(pred_pseudo_1_1)

  # (a, a') = (0, 0)
  lm_y <- lm(y ~ m + a + z + w1 + w2)
  pred_a0z0 <- predict(lm_y,
    newdata = data.frame(m = m, a = 0, z = 0, w1 = w1, w2 = w2)
  )
  pred_a0z1 <- predict(lm_y,
    newdata = data.frame(m = m, a = 0, z = 1, w1 = w1, w2 = w2)
  )

  prob_z <- lm(z ~ a)
  pred_z <- predict(prob_z, newdata = data.frame(a = 0))

  pseudo_out <- pred_a0z0 * (1 - pred_z) + pred_a0z1 * pred_z

  fit_pseudo <- lm(pseudo_out ~ a + w1 + w2)
  pred_pseudo_0_0 <- predict(fit_pseudo,
    newdata = data.frame(a = 0, w1 = w1, w2 = w2)
  )

  res_0_0 <- mean(pred_pseudo_0_0)

  return(list(
    ind_effect = res_1_1 - res_1_0,
    dir_effect = res_1_0 - res_0_0
  ))
}

# Results
res <- workshop_estimates(data)
res

## Bootstrap
n_sim <- 1000
n_boot <- 500
true_sde <- 0.0625841
true_sie <- 0.009845864
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

for (i in 1:n_sim) {
  print(paste0("Simulation n°: ", i))
  data_sim <- read.csv(paste0(sim_data_path, "data_", i, ".csv", sep = ""))
  data_sim <- subset(data_sim, select = -c(y_qol))
  colnames(data_sim)[1] <- "w1"
  colnames(data_sim)[2] <- "w2"
  colnames(data_sim)[3] <- "a"
  colnames(data_sim)[4] <- "z"
  colnames(data_sim)[5] <- "m"
  colnames(data_sim)[6] <- "y"

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
    boot_data <- data_sim[idx,]

    if (b %% 100 == 0 & b != 0) {
      print(paste0("Bootstrap number: ", b))
    }

    ws_boot <- workshop_estimates(boot_data)

    boot_sde_estimates[b,] <- ws_boot$dir_effect
    boot_sie_estimates[b,] <- ws_boot$ind_effect
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
  (estimates_sde_ws$sde - mean(estimates_sde_ws$sde)) ^ 2
) * n_rows / (n_rows - 1)
var_sie <- mean(
  (estimates_sie_ws$sie - mean(estimates_sie_ws$sie)) ^ 2
) * n_rows / (n_rows - 1)

se_sde <- sqrt(var_sde)
se_sie <- sqrt(var_sie)

# Standardized bias
sd_bias_sde <- bias_sde / se_sde
sd_bias_sie <- bias_sie / se_sie

# MSE
mse_sde <- var_sde + bias_sde ^ 2
mse_sie <- var_sie + bias_sie ^ 2

# Average estimated standard error
av_estimated_se_sde <- sqrt(mean(estimates_sde_ws$sd ^ 2))
av_estimated_se_sie <- sqrt(mean(estimates_sie_ws$sd ^ 2))

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
# que ce soit M binaire ou M continu
