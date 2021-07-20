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

file_path <- "../../Data/simulations/"

for (i in 1:1000) {
  print(paste0("Simulation n°:", i))
  data <- gen_data_time_varying_l(N = 10000, a_m_inter = 0)
  write.csv(
    data,
    file = paste0(file_path, "data_", i, ".csv"), row.names = FALSE
  )
  rm(data)
}