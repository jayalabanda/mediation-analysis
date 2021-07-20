#### simulation de Benkeser et al. A nonparametric super-efficient estimator of the ATE.

rm(list = ls())

dataGenerator <- function(N, gamma, a) {
  W1 <- runif(n = N, min = -1.5, max = 1.5)
  W2 <- runif(n = N, min = -1.5, max = 1.5)
  W3 <- runif(n = N, min = -1.5, max = 1.5)
  W4 <- runif(n = N, min = -1.5, max = 1.5)
  W5 <- runif(n = N, min = -1.5, max = 1.5)
  W6 <- runif(n = N, min = -1.5, max = 1.5)
  W7 <- runif(n = N, min = -1.5, max = 1.5)

  W8 <- rbinom(n = N, size = 1, prob = 0.5)

  if (is.null(a)) {
    A <- rbinom(n = N, size = 1, prob = plogis(0.5 * gamma - gamma * W8 +
      2^(1 - 1) * W1 + # = 1 * W1
      2^(1 - 2) * W2 + # = 0.5 * W2
      2^(1 - 3) * W3 + # = 0.25 * W3
      2^(1 - 4) * W4 + # = 0.125 * W4
      2^(1 - 5) * W5 + # = 0.0625 * W5
      2^(1 - 6) * W6 + # = 0.03125 * W6
      2^(1 - 7) * W7)) # = 0.015625 * W7
  } else {
    A <- rep(a, N)
  }

  Y <- rnorm(n = N, mean = A - (2^(1 - 1) * W1 +
    2^(1 - 2) * W2 +
    2^(1 - 3) * W3 +
    2^(1 - 4) * W4 +
    2^(1 - 5) * W5 +
    2^(1 - 6) * W6 +
    2^(1 - 7) * W7), sd = 1)
  data <- data.frame(W1, W2, W3, W4, W5, W6, W7, W8, A, Y)
  if (is.null(a)) {
    return(data)
  } else {
    return(Psi = mean(data$Y))
  }
}

set.seed(1234)
data.gamma0 <- dataGenerator(N = 100000, gamma = 0, a = NULL)
data.gamma3 <- dataGenerator(N = 100000, gamma = 3, a = NULL)
data.gamma6 <- dataGenerator(N = 100000, gamma = 6, a = NULL)


g.gamma0.W <- glm(A ~ W1 + W2 + W3 + W4 + W5 + W6 + W7, family = "binomial", data = data.gamma0)
g.gamma3.W <- glm(A ~ W1 + W2 + W3 + W4 + W5 + W6 + W7 + W8, family = "binomial", data = data.gamma3)
g.gamma6.W <- glm(A ~ W1 + W2 + W3 + W4 + W5 + W6 + W7 + W8, family = "binomial", data = data.gamma6)

pred0.g1.W <- predict(g.gamma0.W, data = data.gamma0, type = "response")
pred3.g1.W <- predict(g.gamma3.W, data = data.gamma3, type = "response")
pred6.g1.W <- predict(g.gamma6.W, data = data.gamma6, type = "response")

boxplot(pred0.g1.W, pred3.g1.W, pred6.g1.W)
summary(data.frame(pred0.g1.W, pred3.g1.W, pred6.g1.W))
#      pred0.g1.W        pred3.g1.W        pred6.g1.W
# Min.   :0.05715   Min.   :0.01298   Min.   :0.003361
# 1st Qu.:0.31905   1st Qu.:0.18472   1st Qu.:0.046552
# Median :0.50105   Median :0.50406   Median :0.311928
# Mean   :0.50091   Mean   :0.50162   Mean   :0.497710
# 3rd Qu.:0.68294   3rd Qu.:0.81767   3rd Qu.:0.952308
# Max.   :0.94732   Max.   :0.98477   Max.   :0.996964

# ça donne effectivement des scores de propension compris entre (0.05;0.95) , (0.01;0.99) , (0.003, 0.997)

# valeurs théoriques : est égale à 1 dans les 3 cas +++
# set.seed(1234)
# Psi.gamma0 <- dataGenerator(N = 1e6, gamma = 0, a = 1)
# Psi.gamma3 <- dataGenerator(N = 1e6, gamma = 3, a = 1)
# Psi.gamma6 <- dataGenerator(N = 1e6, gamma = 6, a = 1)


# exemple avec N=1000 et gamma = 3
# obsdata <- dataGenerator(N = 1000, gamma = 3, a = NULL)

### fonction d'estimation avec CTMLE et TMLE
TMLE_CTMLE.estimator <- function(gamma, obsdata) {

  # 0. transform Y : Y* = (Y-a)/(b-a)
  minY <- min(obsdata$Y) # on peut aussi proposer une valeur théorique
  maxY <- max(obsdata$Y) # on peut aussi proposer une valeur théorique
  obsdata$Y <- (obsdata$Y - minY) / (maxY - minY) # normalise Y à une valeur entre 0 et 1

  # 1. estimate OR (outcome regression)
  # here we use correctly specified glm
  # OR = E_{P_0}(Y | A=1, W=w) => on ne s'intéresse à la régression de Y | W que pour A=1 +++
  est_or <- glm(Y ~ W1 + W2 + W3 + W4 + W5 + W6 + W7, data = obsdata[obsdata$A == 1, ], family = "quasibinomial")

  # 2. predict outcome
  # on logistic scale
  obsdata$logit.Qbar_n.W <- predict(est_or, newdata = obsdata)

  # and on probability scale
  obsdata$Qbar_n.W <- plogis(obsdata$logit.Qbar_n.W)

  # 3. estimate adaptive PS
  ifelse(gamma == 0,
    est_ps_tmle <- glm(A ~ W1 + W2 + W3 + W4 + W5 + W6 + W7, family = "binomial", data = obsdata),
    est_ps_tmle <- glm(A ~ W1 + W2 + W3 + W4 + W5 + W6 + W7 + W8, family = "binomial", data = obsdata)
  )

  # here we use natural splines
  library(splines)
  est_ps_ctmle <- glm(A ~ ns(Qbar_n.W, df = 2), family = "binomial", data = obsdata)

  # here we use hal9001
  # library(hal9001)
  # est_ps_ctmle <- fit_hal(X = obsdata$Qbar_n.W, Y = obsdata$A, fit_type = "glmnet", family = "binomial", yolo = FALSE)
  # ??? la fonction HAL ne fonctionne pas vraiment bien j'ai l'impression...

  # 4. predict PS
  # on probability scale
  obsdata$Gbar_n.W.Qbar_n.tmle <- predict(est_ps_tmle, type = "response")
  obsdata$Gbar_n.W.Qbar_n.ctmle <- predict(est_ps_ctmle, new_data = obsdata, type = "response")

  # 5. fit OR working model
  # define covariate at observed values of A
  obsdata$H_n.AW.tmle <- obsdata$A / obsdata$Gbar_n.W.Qbar_n.tmle
  obsdata$H_n.AW.ctmle <- obsdata$A / obsdata$Gbar_n.W.Qbar_n.ctmle

  # define covariate at A = 1
  obsdata$H_n.1W.tmle <- 1 / obsdata$Gbar_n.W.Qbar_n.tmle
  obsdata$H_n.1W.ctmle <- 1 / obsdata$Gbar_n.W.Qbar_n.ctmle

  # fit working model
  est_wm.tmle <- glm(Y ~ -1 + offset(logit.Qbar_n.W) + H_n.AW.tmle, family = "quasibinomial", data = obsdata) # note: -1 removes the intercept term
  # get coefficient
  eps_hash_n.tmle <- est_wm.tmle$coefficients

  est_wm.ctmle <- glm(Y ~ -1 + offset(logit.Qbar_n.W) + H_n.AW.ctmle, family = "quasibinomial", data = obsdata)
  eps_hash_n.ctmle <- est_wm.ctmle$coefficients

  # 6. target OR estimate
  # generate prediction from working model with
  # covariate H_n.1W
  Qbar_hash_n.tmle <- plogis(obsdata$logit.Qbar_n.W + eps_hash_n.tmle * obsdata$H_n.1W.tmle)
  Qbar_hash_n.ctmle <- plogis(obsdata$logit.Qbar_n.W + eps_hash_n.ctmle * obsdata$H_n.1W.ctmle)

  # 7. compute plug -in
  psi_hash_n.tmle <- mean(Qbar_hash_n.tmle)
  # reverse transformation (continuous outcome)
  psi_hash_n.tmle <- (maxY - minY) * psi_hash_n.tmle + minY
  # à noter que pour la variance il faudra aussi la transformer en sens inverse : \sigma^2 = (b-a)^2.\sigma*^2

  psi_hash_n.ctmle <- mean(Qbar_hash_n.ctmle)
  psi_hash_n.ctmle <- (maxY - minY) * psi_hash_n.ctmle + minY


  ### To compute the one-step estimator, we replace steps 5-7 as follows.
  # 5. evaluate influence function
  Dstar.O.tmle <- obsdata$H_n.AW.tmle * (obsdata$Y - obsdata$Qbar_n.W) + obsdata$Qbar_n.W - mean(obsdata$Qbar_n.W)
  Dstar.O.ctmle <- obsdata$H_n.AW.ctmle * (obsdata$Y - obsdata$Qbar_n.W) + obsdata$Qbar_n.W - mean(obsdata$Qbar_n.W)

  # 6. compute one - step estimator
  psi_hashos_n.tmle <- mean(obsdata$Qbar_n.W) + mean(Dstar.O.tmle)
  psi_hashos_n.tmle <- (maxY - minY) * psi_hashos_n.tmle + minY

  psi_hashos_n.ctmle <- mean(obsdata$Qbar_n.W) + mean(Dstar.O.ctmle)
  psi_hashos_n.ctmle <- (maxY - minY) * psi_hashos_n.ctmle + minY

  return(list(
    psi_hash_n.tmle = psi_hash_n.tmle,
    psi_hash_n.ctmle = psi_hash_n.ctmle,
    psi_hashos_n.tmle = psi_hashos_n.tmle,
    psi_hashos_n.ctmle = psi_hashos_n.ctmle
  ))
}

set.seed(1234)
simdata <- dataGenerator(N = 1000, gamma = 3, a = NULL)
psi.gamma3 <- TMLE_CTMLE.estimator(gamma = 3, obsdata = simdata)
psi.gamma3
# $psi_hash_n.tmle
# [1] 0.9940709

# $psi_hash_n.ctmle
# [1] 1.00563

# $psi_hashos_n.tmle
# [1] 0.9933039

# $psi_hashos_n.ctmle
# [1] 1.005546