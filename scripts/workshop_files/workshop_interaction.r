workshop_estimates <- function(data) {
  require(stringr)

  w_names <- str_subset(colnames(data), "w")
  w <- as.matrix(data[, w_names])
  a <- data$a
  z <- data$z
  m <- data$m
  y <- data$y

  # (a, a') = (1, 0)
  # lm_y <- lm(y ~ m + a + z + w) # a:m
  lm_y <- lm(y ~ z + w + a * m) # note : a * m <=> a + m + a:m
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
  # lm_y <- lm(y ~ m + a + z + w) # a:m
  lm_y <- lm(y ~ z + w + a * m)
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
  # lm_y <- lm(y ~ m + a + z + w)
  lm_y <- lm(y ~ z + w + a * m)
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

data_path <- "../Data/simulations/"
data_path_inter <- "../Data/simulations_inter/"

data_1 <- read.csv(paste0(data_path, "data_100.csv"))
data_2 <- read.csv(paste0(data_path_inter, "data_100.csv"))

colnames(data_1) <- colnames(data_2) <- c("w1", "w2", "a", "z", "m", "y")

res <- workshop_estimates(data_1)
res
# [1] 0.05930487 0.01354218

res <- workshop_estimates(data_2)
res
# [1] 0.04104866 0.02594974

true_sde <- 0.0624
true_sie <- 0.0112
true_sde_inter <- 0.073882
true_sie_inter <- 0.0154

set.seed(42)
n_sim <- 1000
idx <- sample(1:1000, size = n_sim)

bias_estimates <- matrix(nrow = n_sim, ncol = 4)
colnames(bias_estimates) <- c("SDE_no_int", "SIE_no_int", "SDE_int", "SIE_int")

for (i in 1:n_sim) {
  data_1 <- read.csv(paste0(data_path, "data_", idx[i], ".csv"))
  data_2 <- read.csv(paste0(data_path_inter, "data_", idx[i], ".csv"))
  colnames(data_1) <- colnames(data_2) <- c("w1", "w2", "a", "z", "m", "y")

  res_1 <- workshop_estimates(data_1)
  res_2 <- workshop_estimates(data_2)

  bias_estimates[i, 1] <- res_1[1] - true_sde
  bias_estimates[i, 2] <- res_1[2] - true_sie
  bias_estimates[i, 3] <- res_2[1] - true_sde_inter
  bias_estimates[i, 4] <- res_2[2] - true_sie_inter
}

boxplot(bias_estimates)