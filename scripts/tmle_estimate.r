library(data.table)
library(stremr)

file_path <- "../../Data/"
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

qmodel <- "w1 + w2 + a"
mmodel <- "m ~ w1 + w2 + a"
ymodel <- "y ~ w1 + w2 + a + z + m"

g_M_model <- glm(mmodel, family = "binomial", data = data)

data.A1 <- data
data.A1$a <- 1
g_M_1 <- predict(g_M_model, newdata = data.A1, type = "response")

data.A0 <- data
data.A0$a <- 0
g_M_0 <- predict(g_M_model, newdata = data.A0, type = "response")
Q2_model <- glm(ymodel, family = "binomial", data = data)

data.M0 <- data.M1 <- data
data.M1$m <- 1
data.M0$m <- 0

Q2_pred_M1 <- predict(Q2_model, newdata = data.M1, type = "response")
Q2_pred_M0 <- predict(Q2_model, newdata = data.M0, type = "response")

Q2_Gamma_A1 <- Q2_pred_M1 * g_M_1 + Q2_pred_M0 * (1 - g_M_1)

Q2_Gamma_A0 <- Q2_pred_M1 * g_M_0 + Q2_pred_M0 * (1 - g_M_0)

Q1_model_A1 <- glm(
  paste("Q2_Gamma_A1", qmodel, sep = "~"),
  family = "quasibinomial", data = data
)
Q1_model_A0 <- glm(
  paste("Q2_Gamma_A0", qmodel, sep = "~"),
  family = "quasibinomial", data = data
)

Q1_pred_A1_Gamma_A1 <- predict(Q1_model_A1,
  newdata = data.A1, type = "response"
)

Q1_pred_A1_Gamma_A0 <- predict(Q1_model_A0,
  newdata = data.A1, type = "response"
)
Q1_pred_A0_Gamma_A0 <- predict(Q1_model_A0,
  newdata = data.A0, type = "response"
)

Psi_mrNDE <- mean(Q1_pred_A1_Gamma_A0) - mean(Q1_pred_A0_Gamma_A0)
Psi_mrNDE
# [1] 0.0625841

Psi_mrNIE <- mean(Q1_pred_A1_Gamma_A1) - mean(Q1_pred_A1_Gamma_A0)
Psi_mrNIE
# [1] 0.009845864

ymodel <- "Y_death ~ L0_male + L0_parent_low_educ_lv + A0_ace + L1 + M_smoking"
mmodel <- "M_smoking ~ L0_male + L0_parent_low_educ_lv + A0_ace"
amodel <- "A.t ~ L0_male + L0_parent_low_educ_lv"

create.data.table <- function(data, covars, L1, TRT, M, outcome) {
  # binary outcome
  require(data.table)

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

estimate.stremr <- function(data, amodel, mmodel, covars, L1, TRT, M, outcome, Qforms) {
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

  OData <- fitPropensity(
    OData,
    gform_TRT = amodel
  )

  # TMLE
  # E(Y_{a=1, \Gamma_a=0})
  gcomp_Y_A1_G0 <- fit_TMLE(OData,
    tvals = c(0:1),
    intervened_TRT = "p_A1_gamma_A0",
    Qforms = Qforms
  )

  # E(Y_{a=0, \Gamma_a=0})
  gcomp_Y_A0_G0 <- fit_TMLE(OData,
    tvals = c(0:1),
    intervened_TRT = "p_A0_gamma_A0",
    Qforms = Qforms
  )

  # E(Y_{a=1, \Gamma_a=1})
  gcomp_Y_A1_G1 <- fit_TMLE(OData,
    tvals = c(0:1),
    intervened_TRT = "p_A1_gamma_A1",
    Qforms = Qforms
  )

  Psi_mrNDE_GCOMP <-
    gcomp_Y_A1_G0$estimates$cum.inc[2] - gcomp_Y_A0_G0$estimates$cum.inc[2]
  Psi_mrNIE_GCOMP <-
    gcomp_Y_A1_G1$estimates$cum.inc[2] - gcomp_Y_A1_G0$estimates$cum.inc[2]

  return(list(
    Psi_mrNDE_GCOMP = Psi_mrNDE_GCOMP,
    Psi_mrNIE_GCOMP = Psi_mrNIE_GCOMP
  ))
}

data <- data.frame(read.csv(paste(file_path, "data_sim.csv", sep = "")))
data <- subset(data, select = -c(Y_qol)) # remove Y_qol
head(data)

Qforms <- c(
  "Qkplus1 ~ L0_male + L0_parent_low_educ_lv + A.t",
  "Qkplus1 ~ L0_male + L0_parent_low_educ_lv + A.t + L1 + A.tminus1"
)

# Results
res <- estimate.stremr(
  data = data, amodel = amodel, mmodel = mmodel,
  covars = c("L0_male", "L0_parent_low_educ_lv"),
  L1 = "L1", TRT = "A0_ace", M = "M_smoking",
  outcome = "Y_death", Qforms = Qforms
)

Psi_mrNDE_GCOMP <- res$Psi_mrNDE_GCOMP
Psi_mrNDE_GCOMP
# [1] 0.06167823

Psi_mrNIE_GCOMP <- res$Psi_mrNIE_GCOMP
Psi_mrNIE_GCOMP
# [1] 0.00948735