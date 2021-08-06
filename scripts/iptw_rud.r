# define expit function
expit <- function(x) {
  return(1 / (1 + exp(-x)))
}

set.seed(420)
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

iptw_direct_indirect <- function(data) {
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

library(medoutcon)

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

ests_dir_os
# Interventional Direct Effect
# Estimator: onestep
# Estimate: 0.126
# Std. Error: 0.021
# 95% CI: [0.085, 0.167]

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

ests_dir_tmle
# Interventional Direct Effect
# Estimator: tmle
# Estimate: 0.059
# Std. Error: 0.021
# 95% CI: [0.019, 0.1]

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

ests_ind_os
# Interventional Indirect Effect
# Estimator: onestep
# Estimate: 0.033
# Std. Error: 0.012
# 95% CI: [0.009, 0.056]

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

ests_ind_tmle
# Interventional Indirect Effect
# Estimator: tmle
# Estimate: 0.028
# Std. Error: 0.012
# 95% CI: [0.005, 0.051]

# les timations pour l'effet direct sont bonnes avec onestep, mais pas avec tmle
# en revanche, les estimations de l'effet indirect sont bonnes avec tmle
# et moins bonnes avec onestep

file_path <- "../Data/"
data <- data.frame(read.csv(paste(file_path, "data_sim.csv", sep = "")))
data <- subset(data, select = -c(Y_qol))
head(data)

res <- iptw_direct_indirect(data)
res

# TODO: corriger la fonction pour qu'elle s'adapte à ces données
# TODO: comparer avec les résultats de medoutcon