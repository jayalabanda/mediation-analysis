# Create visual representation
plot.dag <- function(exposures, outcomes) {
  require(dagitty)
  require(ggdag)

  dag <- dagitty("dag {
    L_0 -> A
    L_0 -> L_1
    L_0 -> M
    L_0 -> Y
    L_1 -> Y
    L_1 -> M
    A -> M
    A -> L_1
    A -> Y
    M -> Y
    A -> L_1 -> Y
    A -> L_1 -> M -> Y
    A -> M -> Y
  }")

  exposures(dag) <- exposures
  outcomes(dag) <- outcomes

  tidy_dag <- tidy_dagitty(dag)

  p <- ggdag(tidy_dag) + theme_dag()
  print(p)
}

# Plot the DAG
plot.dag(exposures = c("A"), outcomes = c("Y"))

#####################################################
#####################################################

# Manual computation
estimate.manual <- function(data, ymodel, mmodel, l_a_model, A, M) {
  tempdat <- data
  colnames(tempdat)[grep(A, colnames(tempdat))] <- "a"
  colnames(tempdat)[grep(M, colnames(tempdat))] <- "m"

  g_M_model <- glm(formula = mmodel, family = "binomial", data = tempdat)

  # Exposition à A1
  data_A1 <- tempdat
  data_A1$a <- 1
  g_M_1 <- predict(g_M_model, newdata = data_A1, type = "response")

  # Exposition à A0
  data_A0 <- tempdat
  data_A0$a <- 0
  g_M_0 <- predict(g_M_model, newdata = data_A0, type = "response")

  # Q function
  Q_model <- glm(ymodel, family = "binomial", data = tempdat)

  # do(M = 0) ou do(M = 1)
  data_M0 <- data_M1 <- tempdat
  data_M1$m <- 1
  data_M0$m <- 0

  Q_pred_M1 <- predict(Q_model, newdata = data_M1, type = "response")
  Q_pred_M0 <- predict(Q_model, newdata = data_M0, type = "response")

  tempdat$Q_gamma_A1 <- Q_pred_M1 * g_M_1 + Q_pred_M0 * (1 - g_M_1)
  tempdat$Q_gamma_A0 <- Q_pred_M1 * g_M_0 + Q_pred_M0 * (1 - g_M_0)

  Q_model_A1 <- glm(
    formula = paste("Q_gamma_A1", l_a_model, sep = "~"),
    family = "quasibinomial", data = tempdat
  )
  Q_model_A0 <- glm(
    formula = paste("Q_gamma_A0", l_a_model, sep = "~"),
    family = "quasibinomial", data = tempdat
  )

  Q_pred_A1_gamma_A1 <- predict(Q_model_A1, newdata = data_A1, type = "response")
  Q_pred_A1_gamma_A0 <- predict(Q_model_A0, newdata = data_A1, type = "response")
  Q_pred_A0_gamma_A0 <- predict(Q_model_A0, newdata = data_A0, type = "response")

  Psi_mrNDE <- mean(Q_pred_A1_gamma_A0) - mean(Q_pred_A0_gamma_A0)
  Psi_mrNIE <- mean(Q_pred_A1_gamma_A1) - mean(Q_pred_A1_gamma_A0)

  return(list(
    Psi_mrNDE = Psi_mrNDE,
    Psi_mrNIE = Psi_mrNIE
  ))
}

# Retrieve simulated data
file_path <- "../Data/"
data <- read.csv(paste(file_path, "data_sim.csv", sep = ""))
head(data)

bin_to_quant <- function(x) {
  x[x == 0] <- runif(1, min = 0, max = 0.5)
  x[x == 1] <- runif(1, min = 0.5, max = 1)
  return(x)
}

# data$L1 <- # sapply(data$L1, bin_to_quant, simplify = "array")

# Create models
ymodel <- "Y_death ~ L0_male + L0_parent_low_educ_lv + a + L1 + m"
mmodel <- "m ~ L0_male + L0_parent_low_educ_lv + a"
l_a_model <- "L0_male + L0_parent_low_educ_lv + a"

# Results
res <- estimate.manual(
  data = data, ymodel = ymodel, mmodel = mmodel,
  l_a_model = l_a_model, A = "A0_ace", M = "M_smoking"
)

Psi_mrNDE <- res$Psi_mrNDE
Psi_mrNDE

Psi_mrNIE <- res$Psi_mrNIE
Psi_mrNIE


## Bootstrap
file_path <- "../Data/simulations/"
n_sim <- 1000
set.seed(42)
idx <- sample(1:1000, n_sim)
sim <- 1

ymodel <- "y_death ~ l0_male + l0_parent_low_educ_lv + a + l1 + m"
mmodel <- "m ~ l0_male + l0_parent_low_educ_lv + a"
l_a_model <- "l0_male + l0_parent_low_educ_lv + a"

start_time <- Sys.time()
for (i in idx) {
  print(paste0("Simulation ", sim))
  data <- read.csv(paste(file_path, "data_", i, ".csv", sep = ""))
  res <- estimate.manual(
    data = data, ymodel = ymodel, mmodel = mmodel,
    l_a_model = l_a_model, A = "a0_ace", M = "m_smoking"
  )
  sim <- sim + 1
}
end_time <- Sys.time()
diff <- end_time - start_time
diff

# avec n_sim = 100
# 12.28 s

# avec n_sim = 200
# 24.58 s

# avec n_sim = 500
# 61.42 s

# avec n_sim = 1000
# 121.49 s

#####################################################
#####################################################

# Results with stremr package
create.data.table <- function(data, covars, L1, TRT, M, outcome) {
  # binary outcome
  require(data.table)

  # data_stremr <- subset(data, select = -c(Y_qol))
  data_stremr <- data.frame(data)
  data_stremr$ID <- 1:nrow(data_stremr)
  data_stremr$t0 <- rep(0, nrow(data_stremr))
  data_stremr$t1 <- rep(1, nrow(data_stremr))
  data_stremr$alive <- rep(0, nrow(data_stremr))

  new.df <- data.frame(data[, covars])

  dataDT0 <- data.frame(
    ID = data_stremr$ID,
    t = data_stremr$t0,
    L1 = rep(NA, nrow(data_stremr)),
    A.t = data_stremr[, TRT],
    Y.tplus1 = data_stremr$alive,
    A.tminus1 = rep(NA, nrow(data_stremr))
  )

  dataDT1 <- data.frame(
    ID = data_stremr$ID,
    t = data_stremr$t1,
    L1 = data_stremr[, L1],
    A.t = data_stremr[, M],
    Y.tplus1 = data_stremr[, outcome],
    A.tminus1 = data_stremr[, TRT]
  )

  dataDT0 <- cbind(dataDT0, new.df)
  dataDT1 <- cbind(dataDT1, new.df)

  # To data.table
  dataDT <- rbind(dataDT0, dataDT1)
  dataDT <- data.table(dataDT)
  dataDT <- dataDT[order(ID, t)]

  dataDT[, ("p_A1_Gamma_A0") := NA]
  dataDT[, ("p_A0_Gamma_A0") := NA]
  dataDT[, ("p_A1_Gamma_A1") := NA]

  return(dataDT)
}

estimate.stremr <- function(data, mmodel, covars, L1, TRT, M, outcome, Qforms) {
  require(stremr)
  tempdat <- create.data.table(data, covars, L1, TRT, M, outcome)

  g_M_model <- glm(formula = mmodel, family = "binomial", data = data)

  # Exposition à A = 1
  data.A1 <- data
  data.A1[, TRT] <- 1
  g_M_1 <- predict(g_M_model, newdata = data.A1, type = "response")

  # Exposition à A = 0
  data.A0 <- data
  data.A0[, TRT] <- 0
  g_M_0 <- predict(g_M_model, newdata = data.A0, type = "response")

  tempdat$p_A1_gamma_A0[tempdat$t == 0] <- rep(1, nrow(data))
  tempdat$p_A1_gamma_A0[tempdat$t == 1] <- g_M_0

  tempdat$p_A0_gamma_A0[tempdat$t == 0] <- rep(0, nrow(data))
  tempdat$p_A0_gamma_A0[tempdat$t == 1] <- g_M_0

  tempdat$p_A1_gamma_A1[tempdat$t == 0] <- rep(1, nrow(data))
  tempdat$p_A1_gamma_A1[tempdat$t == 1] <- g_M_1

  # Import data for stremr
  OData <- importData(
    tempdat,
    ID = "ID",
    t_name = "t",
    covars = append(covars, L1),
    TRT = "A.t",
    OUTCOME = "Y.tplus1",
    CENS = NULL, MONITOR = NULL
  )

  # g-computation
  # E(Y_{a=1, \Gamma_a=0})
  gcomp_Y_A1_G0 <- fit_GCOMP(OData,
    tvals = c(0:1),
    intervened_TRT = "p_A1_gamma_A0",
    Qforms = Qforms
  )

  # E(Y_{a=0, \Gamma_a=0})
  gcomp_Y_A0_G0 <- fit_GCOMP(OData,
    tvals = c(0:1),
    intervened_TRT = "p_A0_gamma_A0",
    Qforms = Qforms
  )

  # E(Y_{a=1, \Gamma_a=1})
  gcomp_Y_A1_G1 <- fit_GCOMP(OData,
    tvals = c(0:1),
    intervened_TRT = "p_A1_gamma_A1",
    Qforms = Qforms
  )

  Psi_mrNDE_GCOMP <- gcomp_Y_A1_G0$estimates$cum.inc[2] - gcomp_Y_A0_G0$estimates$cum.inc[2]
  Psi_mrNIE_GCOMP <- gcomp_Y_A1_G1$estimates$cum.inc[2] - gcomp_Y_A1_G0$estimates$cum.inc[2]

  return(list(
    Psi_mrNDE_GCOMP = Psi_mrNDE_GCOMP,
    Psi_mrNIE_GCOMP = Psi_mrNIE_GCOMP
  ))
}

Qforms <- c(
  "Qkplus1 ~ L0_male + L0_parent_low_educ_lv + A.t",
  "Qkplus1 ~ L0_male + L0_parent_low_educ_lv + A.t + L1 + A.tminus1"
)

# Results
res <- estimate.stremr(
  data = data, mmodel = mmodel,
  covars = c("L0_male", "L0_parent_low_educ_lv"),
  L1 = "L1", TRT = "A0_ace", M = "M_smoking",
  outcome = "Y_death", Qforms = Qforms
)

Psi_mrNDE_GCOMP <- res$Psi_mrNDE_GCOMP
Psi_mrNDE_GCOMP

Psi_mrNIE_GCOMP <- res$Psi_mrNIE_GCOMP
Psi_mrNIE_GCOMP

#####################################################
#             Article Kara E. Rudolph               #
#####################################################

# Calculer les effets stochastiques directs et indirects selon le modèle :
#     - W = f(Uw)
#     - A = f(W, Ua)
#     - Z = f(W, A, Uz)
#     - M = f(W, A, Z, Um)
#     - Y = f(W, A, Z, M, Uy)

# Ici,  W = (L0_male, L0_parent_low_educ_lv)
#       A = A0_ace
#       Z = L1
#       M = M_smoking
#       Y = Y_death

medtmle <- function(data, covars, A, Z, M, outcome, amodel, zmodel, mmodel, ymodel, qmodel) {
  tmpdat <- obsdat <- data
  # Replace column names
  colnames(obsdat)[grep(A, colnames(obsdat))] <- "a"
  colnames(obsdat)[grep(Z, colnames(obsdat))] <- "z"
  colnames(obsdat)[grep(M, colnames(obsdat))] <- "m"
  colnames(obsdat)[grep(outcome, colnames(obsdat))] <- "y"

  for (i in 1:length(covars)) {
    colnames(obsdat)[grep(covars[i], colnames(obsdat))] <- paste0("w", i)
  }

  # Create data frames
  dfa1 <- dfa0 <- dfa1z1 <- dfa1z0 <- dfa0z1 <- dfa0z0 <- dfm1 <- dfm0 <- tmpdat

  dfa1$a <- dfa1z1$a <- dfa1z1$z <- dfa1z0$a <- dfa0z1$z <- dfm1$m <- 1
  dfa0$a <- dfa1z0$z <- dfa0z1$a <- dfa0z0$a <- dfa0z0$z <- dfm0$m <- 0

  zfit <- glm(formula = zmodel, family = "binomial", data = obsdat)
  zA0 <- predict(zfit, newdata = dfa0, type = "response")
  zA1 <- predict(zfit, newdata = dfa1, type = "response")

  # pas besoin de zA0 et zA1 ?

  mfit <- glm(formula = mmodel, family = "binomial", data = obsdat)
  mZ1A0 <- predict(mfit, newdata = dfa0z1, type = "response")
  mZ0A0 <- predict(mfit, newdata = dfa0z0, type = "response")
  mZ1A1 <- predict(mfit, newdata = dfa1z1, type = "response") # enlever
  mZ0A1 <- predict(mfit, newdata = dfa1z0, type = "response") # enlever

  # TODO
  # juste a = 1 et a = 0

  # Empirical dist. of M_{A=0} | Z
  gmA0 <- mZ1A0 * zA0 + mZ0A0 * (1 - zA0)
  # Empirical dist. of M_{A=1} | Z
  gmA1 <- mZ1A1 * zA1 + mZ0A1 * (1 - zA1)

  # m sans ajuster sur z
  # dfa0, dfa1

  # Weights
  afit <- glm(formula = amodel, family = "binomial", data = tmpdat)
  pred <- predict(afit, type = "response")
  psA1 <- I(tmpdat$a == 1) / pred
  psA0 <- I(tmpdat$a == 0) / (1 - pred)

  mfit <- glm(formula = mmodel, family = "binomial", data = tmpdat)
  mazw <- predict(mfit, type = "response")
  psm <- I(tmpdat$m == 1) * mazw + I(tmpdat$m == 0) * (1 - mazw)

  # Initial estimate of outcome/Y
  yfit <- glm(formula = ymodel, family = "binomial", data = tmpdat)
  tmpdat$qyinit <- cbind(
    predict(yfit, newdata = tmpdat, type = "response"),
    predict(yfit, newdata = dfm0, type = "response"),
    predict(yfit, newdata = dfm1, type = "response")
  )
  # ----------------------------
  # Get E(Y_{1, g_0})
  # ----------------------------

  tmpdat$hA1gmA0 <- ((I(tmpdat$m == 1) * gmA0 + I(tmpdat$m == 0) * (1 - gmA0)) / psm) * psA1

  epsilon.A1gmA0 <- coef(glm(y ~ 1,
    weights = hA1gmA0, offset = qlogis(qyinit[, 1]),
    family = "quasibinomial", data = tmpdat
  ))
  tmpdat$qyupM0A1gmA0 <- plogis(qlogis(tmpdat$qyinit[, 2]) + epsilon.A1gmA0)
  tmpdat$qyupM1A1gmA0 <- plogis(qlogis(tmpdat$qyinit[, 3]) + epsilon.A1gmA0)

  tmpdat$QA1gmA0 <- tmpdat$qyupM0A1gmA0 * (1 - gmA0) + tmpdat$qyupM1A1gmA0 * gmA0

  QA1g0.fit <- glm(
    formula = paste("QA1gmA0", qmodel, sep = "~"), # zmodel ?
    family = "quasibinomial",
    data = tmpdat[tmpdat$a == 1, ]
  )
  QA1g0 <- predict(QA1g0.fit, newdata = tmpdat, type = "response")

  # Necessary if A is non random
  epsilon.zA1gmA0 <- coef(glm(QA1gmA0 ~ 1,
    weights = psA1, offset = qlogis(QA1g0),
    family = "quasibinomial", data = tmpdat
  ))
  QzupA1gmA0 <- plogis(qlogis(QA1g0) + epsilon.zA1gmA0)

  # Calculate TMLE
  tmleA1M0 <- mean(QzupA1gmA0)

  # ----------------------------
  # Get E(Y_{0, g_0})
  # ----------------------------
  tmpdat$hA0gmA0 <- ((tmpdat$m * gmA0 + (1 - tmpdat$m) * (1 - gmA0)) / psm) * psA0

  epsilon.A0gmA0 <- coef(glm(y ~ 1,
    weights = hA0gmA0, offset = qlogis(qyinit[, 1]),
    family = "quasibinomial", data = tmpdat
  ))
  tmpdat$qyupM0A0gmA0 <- plogis(qlogis(tmpdat$qyinit[, 2]) + epsilon.A0gmA0)
  tmpdat$qyupM1A0gmA0 <- plogis(qlogis(tmpdat$qyinit[, 3]) + epsilon.A0gmA0)

  tmpdat$QA0gmA0 <- tmpdat$qyupM0A0gmA0 * (1 - gmA0) + tmpdat$qyupM1A0gmA0 * gmA0

  QA0g0.fit <- glm(
    formula = paste("QA0gmA0", qmodel, sep = "~"),
    family = "quasibinomial",
    data = tmpdat[tmpdat$a == 0, ]
  )
  QA0g0 <- predict(QA0g0.fit, newdata = tmpdat, type = "response")

  # Necessary if A is non random
  epsilon.zA0gmA0 <- coef(glm(QA0gmA0 ~ 1,
    weights = psA0, offset = qlogis(QA0g0),
    family = "quasibinomial", data = tmpdat
  ))
  QzupA0gmA0 <- plogis(qlogis(QA0g0) + epsilon.zA0gmA0)

  # Calculate TMLE
  tmleA0M0 <- mean(QzupA0gmA0)

  # ----------------------------
  # Get E(Y_{1, g_1})
  # ----------------------------
  tmpdat$hA1gmA1 <- ((tmpdat$m * gmA1 + (1 - tmpdat$m) * (1 - gmA1)) / psm) * psA1

  epsilon.A1gmA1 <- coef(glm(y ~ 1,
    weights = hA1gmA1, offset = qlogis(qyinit[, 1]),
    family = "quasibinomial", data = tmpdat
  ))
  tmpdat$qyupM0A1gmA1 <- plogis(qlogis(tmpdat$qyinit[, 2]) + epsilon.A1gmA1)
  tmpdat$qyupM1A1gmA1 <- plogis(qlogis(tmpdat$qyinit[, 3]) + epsilon.A1gmA1)

  tmpdat$QA1gmA1 <- tmpdat$qyupM0A1gmA1 * (1 - gmA1) + tmpdat$qyupM1A1gmA1 * gmA1

  QA1g1.fit <- glm(
    formula = paste("QA1gmA1", qmodel, sep = "~"),
    family = "quasibinomial",
    data = tmpdat[tmpdat$a == 1, ]
  )
  QA1g1 <- predict(QA1g1.fit, newdata = tmpdat, type = "response")

  # Necessary if A is non random
  epsilon.zA1gmA1 <- coef(glm(QA1gmA1 ~ 1,
    weights = psA1, offset = qlogis(QA1g1),
    family = "quasibinomial", data = tmpdat
  ))
  QzupA1gmA1 <- plogis(qlogis(QA1g1) + epsilon.zA1gmA1)

  # Calculate TMLE
  tmleA1M1 <- mean(QzupA1gmA1)

  # ----------------------------
  # Estimate variances using EIC
  # ----------------------------
  tmpdat$qyupA1g0 <- plogis(qlogis(tmpdat$qyinit[, 1]) + epsilon.A1gmA0)
  tmpdat$qyupA1g1 <- plogis(qlogis(tmpdat$qyinit[, 1]) + epsilon.A1gmA1)
  tmpdat$qyupA0g0 <- plogis(qlogis(tmpdat$qyinit[, 1]) + epsilon.A0gmA0)

  # EIC for E(Y_{1, g_0})
  eic1A1g0 <- tmpdat$hA1gmA0 * (tmpdat$y - tmpdat$qyupA1g0)
  eic2A1g0 <- psA1 * (tmpdat$QA1gmA0 - QzupA1gmA0)
  eic3A1g0 <- QzupA1gmA0 - tmleA1M0

  eicA1g0 <- eic1A1g0 + eic2A1g0 + eic3A1g0

  # EIC for E(Y_{1, g_1})
  eic1A1g1 <- tmpdat$hA1gmA1 * (tmpdat$y - tmpdat$qyupA1g1)
  eic2A1g1 <- psA1 * (tmpdat$QA1gmA1 - QzupA1gmA1)
  eic3A1g1 <- QzupA1gmA1 - tmleA1M1

  eicA1g1 <- eic1A1g1 + eic2A1g1 + eic3A1g1

  # EIC for E(Y_{0, g_0})
  eic1A0g0 <- tmpdat$hA0gmA0 * (tmpdat$y - tmpdat$qyupA0g0)
  eic2A0g0 <- psA0 * (tmpdat$QA0gmA0 - QzupA0gmA0)
  eic3A0g0 <- QzupA0gmA0 - tmleA0M0

  eicA0g0 <- eic1A0g0 + eic2A0g0 + eic3A0g0

  # Mediation params., variances and 95% IC
  sde_tmle <- tmleA1M0 - tmleA0M0
  sde_eic <- eicA1g0 - eicA0g0
  var_sde_eic <- var(sde_eic) / nrow(tmpdat)

  sie_tmle <- tmleA1M1 - tmleA1M0
  sie_eic <- eicA1g1 - eicA1g0
  var_sie_eic <- var(sie_eic) / nrow(tmpdat)

  results <- data.frame(
    cbind(
      sde = sde_tmle,
      sde_var = var_sde_eic,
      sie = sie_tmle,
      sie_var = var_sie_eic,
      sde_lb = sde_tmle - 1.96 * sqrt(var_sde_eic), # sde lower bound
      sde_ub = sde_tmle + 1.96 * sqrt(var_sde_eic), # sde upper bound
      sie_lb = sie_tmle - 1.96 * sqrt(var_sie_eic), # sie lower bound
      sie_ub = sie_tmle + 1.96 * sqrt(var_sie_eic) # sie upper bound
    )
  )
}

file_path <- "../../Data/"
data <- data.frame(read.csv(paste(file_path, "data_sim.csv", sep = "")))
data <- subset(data, select = -c(Y_qol)) # remove Y_qol
head(data)

amodel <- "a ~ w1 + w2"
zmodel <- "z ~ w1 + w2 + a"
mmodel <- "m ~ w1 + w2 + a"
ymodel <- "y ~ w1 + w2 + a + z + m"
qmodel <- "w1 + w2"

res <- medtmle(
  data = data,
  covars <- c("L0_male", "L0_parent_low_educ_lv"),
  A = "A0_ace",
  Z = "L1",
  M = "M_smoking",
  outcome = "Y_death",
  amodel = amodel,
  zmodel = zmodel,
  mmodel = mmodel,
  ymodel = ymodel,
  qmodel = qmodel
)

res