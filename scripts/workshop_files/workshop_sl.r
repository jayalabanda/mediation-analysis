workshop_estimates_sl <- function(data, sl_library) {
  require(stringr)
  require(SuperLearner)

  w_names <- str_subset(colnames(data), "w")
  w <- as.matrix(data[, w_names])
  a <- data$a
  z <- data$z
  m <- data$m
  y <- data$y

  # (a, a') = (1, 0)
  # lm_y <- lm(y ~ z + w + a * m)
  x <- data.frame(a = a, z = z, w = w, m = m)
  lm_y <- SuperLearner(y, x,
    family = "binomial",
    SL.library = sl_library
  )

  newdat <- data.frame(a = a, z = z, w = w, m = m)
  newdat$a <- 1
  newdat$z <- 0
  pred_a1z0 <- predict(lm_y,
    newdata = newdat
  )$pred

  newdat$z <- 1
  pred_a1z1 <- predict(lm_y,
    newdata = newdat
  )$pred

  y <- z
  x <- as.data.frame(a)
  # prob_z <- lm(z ~ a)
  prob_z <- SuperLearner(y, x,
    family = "binomial",
    SL.library = sl_library
  )

  newdat <- data.frame(a = a)
  newdat$a <- 1
  pred_z <- predict(prob_z, newdata = newdat)$pred

  pseudo_out_1_0 <- pred_a1z0 * (1 - pred_z) + pred_a1z1 * pred_z

  y <- pseudo_out_1_0
  x <- data.frame(a = a, w = w)
  # fit_pseudo <- lm(pseudo_out_1_0 ~ a + w)
  fit_pseudo <- SuperLearner(y, x,
    family = "binomial",
    SL.library = sl_library
  )

  newdat <- data.frame(a = a, w = w)
  newdat$a <- 0
  pred_pseudo_1_0 <- predict(fit_pseudo,
    newdata = newdat
  )$pred

  res_1_0 <- mean(pred_pseudo_1_0)

  # (a, a') = (1, 1)
  # lm_y <- lm(y ~ z + w + a * m)
  y <- data$y
  x <- data.frame(a = a, z = z, w = w, m = m)
  lm_y <- SuperLearner(y, x,
    family = "binomial",
    SL.library = sl_library
  )

  newdat <- data.frame(a = a, z = z, w = w, m = m)
  newdat$a <- 1
  newdat$z <- 0
  pred_a1z0 <- predict(lm_y,
    newdata = newdat
  )$pred

  newdat$z <- 1
  pred_a1z1 <- predict(lm_y,
    newdata = newdat
  )$pred

  y <- z
  x <- as.data.frame(a)
  # prob_z <- lm(z ~ a)
  prob_z <- SuperLearner(y, x,
    family = "binomial",
    SL.library = sl_library
  )

  newdat <- data.frame(a = a)
  newdat$a <- 1
  pred_z <- predict(prob_z, newdata = newdat)$pred

  pseudo_out_1_1 <- pred_a1z0 * (1 - pred_z) + pred_a1z1 * pred_z

  y <- pseudo_out_1_1
  x <- data.frame(a = a, w = w)
  # fit_pseudo <- lm(pseudo_out_1_1 ~ a + w)
  fit_pseudo <- SuperLearner(y, x,
    family = "binomial",
    SL.library = sl_library
  )

  newdat <- data.frame(a = a, w = w)
  newdat$a <- 1
  pred_pseudo_1_1 <- predict(fit_pseudo,
    newdata = newdat
  )$pred

  res_1_1 <- mean(pred_pseudo_1_1)

  # (a, a') = (0, 0)
  # lm_y <- lm(y ~ z + w + a * m)
  y <- data$y
  x <- data.frame(a = a, z = z, w = w, m = m)
  lm_y <- SuperLearner(y, x,
    family = "binomial",
    SL.library = sl_library
  )

  newdat <- data.frame(a = a, z = z, w = w, m = m)
  newdat$a <- 0
  newdat$z <- 0
  pred_a0z0 <- predict(lm_y,
    newdata = newdat
  )$pred

  newdat$z <- 1
  pred_a0z1 <- predict(lm_y,
    newdata = newdat
  )$pred

  y <- z
  x <- as.data.frame(a)
  # prob_z <- lm(z ~ a)
  prob_z <- SuperLearner(y, x,
    family = "binomial",
    SL.library = sl_library
  )

  newdat <- data.frame(a = a)
  newdat$a <- 0
  pred_z <- predict(prob_z, newdata = newdat)$pred

  pseudo_out_0_0 <- pred_a0z0 * (1 - pred_z) + pred_a0z1 * pred_z

  y <- pseudo_out_0_0
  x <- data.frame(a = a, w = w)
  # fit_pseudo <- lm(pseudo_out_0_0 ~ a + w)
  fit_pseudo <- SuperLearner(y, x,
    family = "binomial",
    SL.library = sl_library
  )

  newdat <- data.frame(a = a, w = w)
  newdat$a <- 0
  pred_pseudo_0_0 <- predict(fit_pseudo,
    newdata = newdat
  )$pred

  res_0_0 <- mean(pred_pseudo_0_0)

  dir_effect <- res_1_0 - res_0_0
  ind_effect <- res_1_1 - res_1_0

  return(c(
    dir_effect,
    ind_effect
  ))
}

data_path <- "../Data/simulations/"
libs <- c(
  "SL.glm", "SL.bayesglm", "SL.glm.interaction",
  "SL.speedlm", "SL.mean", "SL.earth"
)

true_sde <- 0.0624
true_sie <- 0.0112
true_sde_rud <- 0.124793
true_sie_rud <- 0.03026875

n_sims <- c(100, 200, 300, 400, 500)
for (n_sim in n_sims) {
  set.seed(42)
  idx <- sample(1:1000, size = n_sim)

  bias_estimates <- matrix(nrow = n_sim, ncol = 2)
  colnames(bias_estimates) <- c("SDE", "SIE")

  ws_sl_sde_estimates <- matrix(nrow = n_sim, ncol = length(libs))
  ws_sl_sie_estimates <- matrix(nrow = n_sim, ncol = length(libs))
  colnames(ws_sl_sde_estimates) <- libs
  colnames(ws_sl_sie_estimates) <- libs

  write.csv(ws_sl_sde_estimates,
    paste0("scripts/workshop_files/sde_ws_", n_sim, ".csv"),
    row.names = FALSE
  )
  write.csv(ws_sl_sie_estimates,
    paste0("scripts/workshop_files/sie_ws_", n_sim, ".csv"),
    row.names = FALSE
  )

  for (j in seq_len(length(libs))) {
    cat(paste0("\n######### Library ", j, " of ", length(libs), " #########\n"))
    cat(paste0(libs[j], "\n"))

    for (i in 1:n_sim) {
      if (i %% 10 == 0 | i == 1) {
        print(paste0("Simulation ", i, " of ", n_sim))
      }

      data <- read.csv(paste0(data_path, "data_", idx[i], ".csv"))

      if (str_detect(data_path, "rudolph")) {
        res_1 <- try(workshop_estimates_sl(data, libs[j]))
        bias_estimates[i, 1] <- res_1[1] - true_sde_rud
        bias_estimates[i, 2] <- res_1[2] - true_sie_rud
      } else {
        data <- subset(data, select = -y_qol)
        colnames(data) <- c("w1", "w2", "a", "z", "m", "y")

        res_1 <- try(workshop_estimates_sl(data, libs[j]))
        bias_estimates[i, 1] <- res_1[1] - true_sde
        bias_estimates[i, 2] <- res_1[2] - true_sie
      }
    }

    ws_sl_sde_estimates <- read.csv(paste0(
      "scripts/workshop_files/sde_ws_", n_sim, ".csv"
    ))
    ws_sl_sie_estimates <- read.csv(paste0(
      "scripts/workshop_files/sie_ws_", n_sim, ".csv"
    ))

    ws_sl_sde_estimates[, j] <- bias_estimates[, 1]
    ws_sl_sie_estimates[, j] <- bias_estimates[, 2]

    write.csv(
      ws_sl_sde_estimates,
      paste0("scripts/workshop_files/sde_ws_", n_sim, ".csv"),
      row.names = FALSE
    )
    write.csv(
      ws_sl_sie_estimates,
      paste0("scripts/workshop_files/sie_ws_", n_sim, ".csv"),
      row.names = FALSE
    )

    par(mfrow = c(1, 2))
    boxplot(ws_sl_sde_estimates)
    boxplot(ws_sl_sie_estimates)
  }
}
