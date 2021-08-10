# define expit function
expit <- function(x) {
  return(1 / (1 + exp(-x)))
}

set.seed(1000)
n <- 10000
w_1 <- rbinom(n, 1, 0.6)
w_2 <- rbinom(n, 1, 0.3)
w_3 <- rbinom(n, 1, 0.2 + (w_1 + w_2) / 3)
a <- rbinom(n, 1, expit(0.25 * (w_1 + w_2 + w_3) + 3 * w_1 * w_2 - 2))
z <- rbinom(n, 1, expit((w_1 + w_2 + w_3) / 3 - a - a * w_3 - 0.25))
m <- rbinom(n, 1, expit(w_1 + w_2 + a - z + a * z - 0.3 * a * w_2))
y <- rbinom(n, 1, expit((a - z + m - a * z) / (w_1 + w_2 + w_3 + 1)))

data <- data.frame(w_1, w_2, w_3, a, z, m, y)
head(data)

iptw_direct_indirect <- function(data, covariates, outcome) {
  g_a <- glm(a ~ 1, family = "binomial", data = data)
  g_a_l0 <- glm(a ~ w_1 + w_2 + w_3, family = "binomial", data = data)
  g_m_a <- glm(m ~ a, family = "binomial", data = data)
  g_m_l <- glm(m ~ w_1 + w_2 + w_3 + a + z, family = "binomial", data = data)

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
  msm_m <- glm(m ~ a + w_1 + w_2 + w_3, family = "binomial", data = data)

  iptw_edn_l0 <- msm_y$coefficients["a"] +
    (msm_y$coefficients["a:m"] *
      plogis(rep(msm_m$coefficients["(Intercept)"], nrow(data)) +
        msm_m$coefficients["w_1"] * data$w_1 +
        msm_m$coefficients["w_2"] * data$w_2 +
        msm_m$coefficients["w_3"] * data$w_3))

  iptw_ein_l0 <- (msm_y$coefficients["m"] + msm_y$coefficients["a:m"]) *
    (plogis(rep(msm_m$coefficients["(Intercept)"], nrow(data)) +
      msm_m$coefficients["a"] +
      msm_m$coefficients["w_1"] * data$w_1 +
      msm_m$coefficients["w_2"] * data$w_2 +
      msm_m$coefficients["w_3"] * data$w_3) -
      plogis(rep(msm_m$coefficients["(Intercept)"], nrow(data)) +
        msm_m$coefficients["w_1"] * data$w_1 +
        msm_m$coefficients["w_2"] * data$w_2 +
        msm_m$coefficients["w_3"] * data$w_3))

  iptw_edn <- mean(iptw_edn_l0)
  iptw_ein <- mean(iptw_ein_l0)

  return(list(iptw_edn = iptw_edn, iptw_ein = iptw_ein))
}

res <- iptw_direct_indirect(data)
res
# $iptw_edn
# [1] 0.1354538

# $iptw_ein
# [1] 0.02732996

iptw_direct_indirect_bis <- function(data) {
  g_a <- glm(a ~ 1, family = "binomial", data = data)
  g_a_l0 <- glm(a ~ w_1 + w_2, family = "binomial", data = data)
  g_m_a <- glm(m ~ a, family = "binomial", data = data)
  g_m_l <- glm(m ~ w_1 + w_2 + a + z, family = "binomial", data = data)

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
  msm_m <- glm(m ~ a + w_1 + w_2, family = "binomial", data = data)

  iptw_edn_l0 <- msm_y$coefficients["a"] +
    (msm_y$coefficients["a:m"] *
      plogis(rep(msm_m$coefficients["(Intercept)"], nrow(data)) +
        msm_m$coefficients["w_1"] * data$w_1 +
        msm_m$coefficients["w_2"] * data$w_2))

  iptw_ein_l0 <- (msm_y$coefficients["m"] + msm_y$coefficients["a:m"]) *
    (plogis(rep(msm_m$coefficients["(Intercept)"], nrow(data)) +
      msm_m$coefficients["a"] +
      msm_m$coefficients["w_1"] * data$w_1 +
      msm_m$coefficients["w_2"] * data$w_2) -
      plogis(rep(msm_m$coefficients["(Intercept)"], nrow(data)) +
        msm_m$coefficients["w_1"] * data$w_1 +
        msm_m$coefficients["w_2"] * data$w_2))

  iptw_edn <- mean(iptw_edn_l0)
  iptw_ein <- mean(iptw_ein_l0)

  return(list(iptw_edn = iptw_edn, iptw_ein = iptw_ein))
}

# Speed
file_path <- "../Data/simulations/"
n_sim <- 1000
set.seed(42)
idx <- sample(1:1000, size = n_sim)
sim <- 1

start_time <- Sys.time()
for (i in idx) {
  print(paste("Simulation ", sim))
  data <- read.csv(paste0(file_path, "data_", i, ".csv", sep = ""))
  data <- subset(data, select = -c(y_qol))
  colnames(data) <- c("w_1", "w_2", "a", "z", "m", "y")

  res <- iptw_direct_indirect_bis(data)
  sim <- sim + 1
}
end_time <- Sys.time()
diff <- end_time - start_time
diff

# avec n_sim = 100
# 13.81 s

# avec n_sim = 100
# 26.91 s

# avec n_sim = 500
# 68.42 s

# avec n_sim = 1000
# 140.54 s

###############
library(medoutcon)

start_time_1 <- Sys.time()
ests_dir_os <- medoutcon(
  W = data.frame(w_1, w_2, w_3),
  A = a,
  Z = z,
  M = m,
  Y = y,
  effect = "direct",
  u_learners = sl3::Lrnr_glm_fast$new(),
  v_learners = sl3::Lrnr_glm_fast$new(),
  estimator = "onestep"
)
end_time_1 <- Sys.time()
diff_1 <- end_time_1 - start_time_1
diff_1

ests_dir_os
# Interventional Direct Effect
# Estimator: onestep
# Estimate: 0.126
# Std. Error: 0.021
# 95% CI: [0.085, 0.167]

start_time_2 <- Sys.time()
ests_dir_tmle <- medoutcon(
  W = data.frame(w_1, w_2, w_3),
  A = a,
  Z = z,
  M = m,
  Y = y,
  effect = "direct",
  u_learners = sl3::Lrnr_glm_fast$new(),
  v_learners = sl3::Lrnr_glm_fast$new(),
  estimator = "tmle"
)
end_time_2 <- Sys.time()
diff_2 <- end_time_2 - start_time_2
diff_2

ests_dir_tmle
# Interventional Direct Effect
# Estimator: tmle
# Estimate: 0.059
# Std. Error: 0.021
# 95% CI: [0.019, 0.1]

start_time_3 <- Sys.time()
ests_ind_os <- medoutcon(
  W = data.frame(w_1, w_2, w_3),
  A = a,
  Z = z,
  M = m,
  Y = y,
  effect = "indirect",
  u_learners = sl3::Lrnr_glm_fast$new(),
  v_learners = sl3::Lrnr_glm_fast$new(),
  estimator = "onestep"
)
end_time_3 <- Sys.time()
diff_3 <- end_time_3 - start_time_3
diff_3

ests_ind_os
# Interventional Indirect Effect
# Estimator: onestep
# Estimate: 0.033
# Std. Error: 0.012
# 95% CI: [0.009, 0.056]

start_time_4 <- Sys.time()
ests_ind_tmle <- medoutcon(
  W = data.frame(w_1, w_2, w_3),
  A = a,
  Z = z,
  M = m,
  Y = y,
  effect = "indirect",
  u_learners = sl3::Lrnr_glm_fast$new(),
  v_learners = sl3::Lrnr_glm_fast$new(),
  estimator = "tmle"
)
end_time_4 <- Sys.time()
diff_4 <- end_time_4 - start_time_4
diff_4

ests_ind_tmle
# Interventional Indirect Effect
# Estimator: tmle
# Estimate: 0.028
# Std. Error: 0.012
# 95% CI: [0.005, 0.051]

# les estimations pour l'effet direct sont bonnes avec onestep, pas avec tmle
# en revanche, les estimations de l'effet indirect sont bonnes avec tmle
# et moins bonnes avec onestep

file_path <- "../Data/"
data <- data.frame(read.csv(paste(file_path, "data_sim.csv", sep = "")))
data <- subset(data, select = -c(Y_qol))
head(data)

res <- iptw_direct_indirect(data)
res

file_path <- "../Data/simulations/"
set.seed(42)
idx <- sample(1:1000, size = 500)

start_time <- Sys.time()
for (i in 1:1000) {
  print(paste0("Simulation ", i))
  data <- data.frame(read.csv(paste(file_path, "data_", i, ".csv", sep = "")))
  data <- subset(data, select = -c(y_qol))

  ests_dir_os <- medoutcon(
    W = data.frame(data[, 1], data[, 2]),
    A = data[, 3],
    Z = data[, 4],
    M = data[, 5],
    Y = data[, 6],
    effect = "indirect",
    u_learners = sl3::Lrnr_glm_fast$new(),
    v_learners = sl3::Lrnr_glm_fast$new(),
    estimator = "tmle"
  )
}
end_time <- Sys.time()
diff <- end_time - start_time
diff

# 100 samples, direct effect
# one-step: 6.01 mins
# tmle: 6.78 mins

# 200 samples, direct effect
# one-step: 16.83 mins
# tmle: 15.05 mins

# 500 samples, direct effect
# one-step: 53.67 mins
# tmle: 62.46 mins

# 1000 samples, direct effect
# one-step: 167.16 mins
# tmle: 105.70 mins


# 100 samples, indirect effect
# one-step: 8.63 mins
# tmle: 8.70 mins

# 200 samples, indirect effect
# one-step: 16.40 mins
# tmle: 15.73 mins

# 500 samples, indirect effect
# one-step: 49.75 mins
# tmle: 49.16 mins

# 1000 samples, indirect effect
# one-step: 100.31 mins
# tmle: 107.15 mins