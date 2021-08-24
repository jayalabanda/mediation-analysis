library(medoutcon)
library(sl3)
library(tidyverse)

# Create practical function
ind_dir_effects_medoutcon <- function(data, w_names, m_names) {
  dir_os <- medoutcon(
    W = data[, w_names],
    A = data$A,
    Z = data$Z,
    M = data[, m_names],
    Y = data$Y,
    effect = "direct",
    u_learners = sl3::Lrnr_glm_fast$new(),
    v_learners = sl3::Lrnr_glm_fast$new(),
    estimator = "onestep"
  )

  dir_tmle <- medoutcon(
    W = data[, w_names],
    A = data$A,
    Z = data$Z,
    M = data[, m_names],
    Y = data$Y,
    effect = "direct",
    u_learners = sl3::Lrnr_glm_fast$new(),
    v_learners = sl3::Lrnr_glm_fast$new(),
    estimator = "tmle"
  )

  ind_os <- medoutcon(
    W = data[, w_names],
    A = data$A,
    Z = data$Z,
    M = data[, m_names],
    Y = data$Y,
    effect = "indirect",
    u_learners = sl3::Lrnr_glm_fast$new(),
    v_learners = sl3::Lrnr_glm_fast$new(),
    estimator = "onestep"
  )

  ind_tmle <- medoutcon(
    W = data[, w_names],
    A = data$A,
    Z = data$Z,
    M = data[, m_names],
    Y = data$Y,
    effect = "indirect",
    u_learners = sl3::Lrnr_glm_fast$new(),
    v_learners = sl3::Lrnr_glm_fast$new(),
    estimator = "tmle"
  )

  return(list(
    dir_result_os = dir_os,
    dir_result_tmle = dir_tmle,
    ind_result_os = ind_os,
    ind_result_tmle = ind_tmle
  ))
}

data_path <- "../Data/new_simulations/"
true_sde <- 0.067885
true_sie <- 0.0154

n_sim <- 20
idx <- sample(1:1000, size = n_sim)
med_results <- matrix(nrow = n_sim, ncol = 4)
colnames(med_results) <- c("sde_os", "sde_tmle", "sie_os", "sie_tmle")

for (i in 1:n_sim) {
  print(paste0("Simulation n°: ", i))

  data_bis <- read.csv(paste0(data_path, "data_", idx[i], ".csv", sep = ""))
  data_bis <- subset(data_bis, select = -c(y_qol))
  colnames(data_bis) <- c("W_1", "W_2", "A", "Z", "M_1", "Y")

  w_names <- str_subset(colnames(data_bis), "W")
  m_names <- str_subset(colnames(data_bis), "M")

  # Calculate effects from function
  results_dir_ind_bis <- ind_dir_effects_medoutcon(data_bis, w_names, m_names)
  med_results[i, 1] <- results_dir_ind_bis$dir_result_os$theta - true_sde
  med_results[i, 2] <- results_dir_ind_bis$dir_result_tmle$theta - true_sde
  med_results[i, 3] <- results_dir_ind_bis$ind_result_os$theta - true_sie
  med_results[i, 4] <- results_dir_ind_bis$ind_result_tmle$theta - true_sie
}

med_results <- data.frame(med_results)
boxplot(med_results)
abline(h = 0, col = "red")

write.csv(med_results,
  file = "./Results/medoutcon_positivity.csv",
  row.names = FALSE
)