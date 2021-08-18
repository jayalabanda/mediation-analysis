

rm(list=ls())

#############################################################################################################################################
#############################################################################################################################################
##### 0) Code pour définir les paramètres des équations structurelles
##### en présence de confusion intermédiaire influencé par l'exposition initiale
#############################################################################################################################################
#############################################################################################################################################

sim.param.time.varying.L <- function(A.M.interaction = NULL) {
  # L0
  p_L0_male <- 0.5
  p_L0_parent_low_educ_lv <- 0.65
  
  # A: A0_ace <- rbinom( 0.05 + 0.04 * L0_male + 0.06 * L0_parent_low_educ_lv ) 
  b_A <- 0.05   # reference prevalence is 5%
  b_male_A <- 0.04  # + 0.04 for the effect of L0_male -> A0_ace
  b_parent_educ_A <- 0.06  # +0.06 for the effect of L0_parent_low_educ_lv -> A0_ace
  
  # L1: L1 <- rbinom( 0.30 - 0.05 * L0_male + 0.08 * L0_parent_low_educ_lv + 0.2 * A0_ace ) 
  b_L1 <- 0.30   # reference prevalence is 30%
  b_male_L1 <- -0.05  # - 0.05 for the effect of L0_male -> L1
  b_parent_L1 <- +0.08 # + 0.08 for the effect of L0_parent_low_educ_lv -> L1
  b_A_L1 <- +0.2 # +0.2 for the effect of A0_ace -> L1
  
  # M: M_smoking <- rbinom( 0.2 + 0.05 * L0_male + 0.06 * L0_parent_low_educ_lv + 0.2 * L1 + 0.1 * A0_ace ) 
  b_M <- 0.2 # reference prevalence is 20%
  b_male_M <- 0.05 # +0.05 for the effect of L0_male -> M_smoking
  b_parent_educ_M <- 0.06 # +0.06 for the effect of L0_parent_low_educ_lv -> M_smoking
  b_A_M <- 0.1 # +0.10 for the effect of A0_ace -> M_smoking
  b_L1_M <- 0.2 # +0.2 for the effect of L1 -> M_smoking
  
  # Y binary: rbinom( 0.10 + 0.06 * L0_male + 0.04 * L0_parent_low_educ_lv + 0.05 * A0_ace + 0.07 * L1 + 0.08 * M_smoking +
  #                   0.03 * A0_ace * M_smoking * A.M.inter ) 
  b_Y <- 0.1 # reference prevalence is 10%
  b_male_Y <- 0.06 # +0.06 for the effect of L0_male -> Y
  b_parent_educ_Y <- 0.04 # +0.04 for the effect of L0_parent_low_educ_lv -> Y
  b_A_Y <- 0.05 # 0.05 for the effect of A0_ace -> Y
  b_L1_Y <- 0.07 # +0.07 for the effect of L1 -> Y
  b_M_Y <- 0.08 # 0.08 for the effect of M_smoking -> Y
  b_AM_Y <- 0.03 # 0.03 for the interaction effect A0_ace * M_smoking -> Y
  
  # Y continuous: (75 - 1 * L0_male - 3 * L0_parent_low_educ_lv - 4 * A0_ace -3.5 * L1 - 9 * M_smoking + 
  #             -5 * A0_ace * M_smoking * A.M.inter ) + rnorm(N, mean = 0, sd = 10)
  mu_Y <- 75 # reference mean for QoL
  c_male_Y <- -1 # -1 for the effect of L0_male -> Y
  c_parent_educ_Y <- -3 # -3 for the effect of L0_parent_low_educ_lv -> Y
  c_A_Y <- -4 # -4 for the effect of A0_ace -> Y
  c_L1_Y <- -5 # -5 for the effect of L1 -> Y
  c_M_Y <- -9 # -9 for the effect of M_smoking -> Y
  c_AM_Y <- -5  # - 5 for the interaction effect A0_ace * M_smoking  -> Y
  sd_Y <- 10 # standard deviation of the residuals
  
  # A*M interaction ?
  A.M.inter <- A.M.interaction
  
  coef <- c( p_L0_male = p_L0_male, p_L0_parent_low_educ_lv = p_L0_parent_low_educ_lv, 
             b_A = b_A, b_male_A = b_male_A, b_parent_educ_A = b_parent_educ_A, 
             b_L1 = b_L1, b_male_L1 = b_male_L1, b_parent_L1 = b_parent_L1, b_A_L1 = b_A_L1,
             b_M = b_M, b_male_M = b_male_M, b_parent_educ_M = b_parent_educ_M, b_L1_M = b_L1_M, b_A_M = b_A_M,
             b_Y = b_Y, b_male_Y = b_male_Y, b_parent_educ_Y = b_parent_educ_Y, b_A_Y = b_A_Y, b_L1_Y = b_L1_Y, b_M_Y = b_M_Y, b_AM_Y = b_AM_Y,
             mu_Y = mu_Y, c_male_Y = c_male_Y, c_parent_educ_Y = c_parent_educ_Y, c_A_Y = c_A_Y, c_L1_Y = c_L1_Y, c_M_Y = c_M_Y, c_AM_Y = c_AM_Y, 
             sd_Y = sd_Y, A.M.inter = A.M.inter)
  
  return(coef)
}


#############################################################################################################################################
#############################################################################################################################################
### 1) Calcul de l'effet total théorique  : ATE = E(Y_{a=1}) - E(Y_{a=0})
#############################################################################################################################################
#############################################################################################################################################

true.ATE.time.var.conf <- function(interaction = NULL) {
  b <- sim.param.time.varying.L(A.M.interaction = interaction)
  
  # binary outcome (death)
  S <- cbind(expand.grid(c(0,1),c(0,1),c(0,1), c(0,1)), rep(NA,n=2^4))
  colnames(S) <- list("male","parent_educ","L1","M","sum")
  for (n in 1:16) {
    S[n,"sum"] <- ( ( ( b["b_Y"] + 
                          b["b_male_Y"] * S[n,"male"] + 
                          b["b_parent_educ_Y"] * S[n,"parent_educ"] + 
                          b["b_A_Y"] * 1 + 
                          b["b_L1_Y"] * S[n,"L1"] +
                          b["b_M_Y"] * S[n,"M"] +
                          b["b_AM_Y"] * 1 * S[n,"M"] * b["A.M.inter"] ) *
                        (( b["b_M"] + 
                             b["b_male_M"] * S[n,"male"] + 
                             b["b_parent_educ_M"] * S[n,"parent_educ"] + 
                             b["b_L1_M"] * S[n,"L1"] +
                             b["b_A_M"] * 1 )^( S[n,"M"] )) *
                        (( 1 - (b["b_M"] + 
                                  b["b_male_M"] * S[n,"male"] + 
                                  b["b_parent_educ_M"] * S[n,"parent_educ"] + 
                                  b["b_L1_M"] * S[n,"L1"] +
                                  b["b_A_M"] * 1) )^( 1 - S[n,"M"] )) *
                        (( b["b_L1"] +
                             b["b_male_L1"] * S[n,"male"] +  
                             b["b_parent_L1"] * S[n,"parent_educ"] +
                             b["b_A_L1"] * 1)^( S[n,"L1"] )) *
                        (( 1 - ( b["b_L1"] +
                                   b["b_male_L1"] * S[n,"male"] +  
                                   b["b_parent_L1"] * S[n,"parent_educ"] +
                                   b["b_A_L1"] * 1))^( 1 - S[n,"L1"] )) ) - 
                      ( ( b["b_Y"] + 
                            b["b_male_Y"] * S[n,"male"] + 
                            b["b_parent_educ_Y"] * S[n,"parent_educ"] + 
                            b["b_A_Y"] * 0 + 
                            b["b_L1_Y"] * S[n,"L1"] +
                            b["b_M_Y"] * S[n,"M"] +
                            b["b_AM_Y"] * 0 * S[n,"M"] * b["A.M.inter"] ) *
                          (( b["b_M"] + 
                               b["b_male_M"] * S[n,"male"] + 
                               b["b_parent_educ_M"] * S[n,"parent_educ"] + 
                               b["b_L1_M"] * S[n,"L1"] +
                               b["b_A_M"] * 0 )^( S[n,"M"] )) *
                          (( 1 - (b["b_M"] + 
                                    b["b_male_M"] * S[n,"male"] + 
                                    b["b_parent_educ_M"] * S[n,"parent_educ"] + 
                                    b["b_L1_M"] * S[n,"L1"] +
                                    b["b_A_M"] * 0) )^( 1 - S[n,"M"] )) *
                          (( b["b_L1"] +
                               b["b_male_L1"] * S[n,"male"] +  
                               b["b_parent_L1"] * S[n,"parent_educ"] +
                               b["b_A_L1"] * 0)^( S[n,"L1"] )) *
                          (( 1 - ( b["b_L1"] +
                                     b["b_male_L1"] * S[n,"male"] +  
                                     b["b_parent_L1"] * S[n,"parent_educ"] +
                                     b["b_A_L1"] * 0))^( 1 - S[n,"L1"] )) ) ) *
      ((b["p_L0_male"])^(S[n,"male"])) * 
      ((1 - b["p_L0_male"])^(1 - S[n,"male"])) * 
      ((b["p_L0_parent_low_educ_lv"])^(S[n,"parent_educ"])) *
      ((1 - b["p_L0_parent_low_educ_lv"])^(1 - S[n,"parent_educ"])) 
  }
  
  ATE.death <- sum(S[,"sum"])
  
  # quantitative outcome (QoL)
  S <- cbind(expand.grid(c(0,1),c(0,1),c(0,1), c(0,1)), rep(NA,n=2^4))
  colnames(S) <- list("male","parent_educ","L1","M","sum")
  for (n in 1:16) {
    S[n,"sum"] <- ( ( ( b["mu_Y"] + 
                          b["c_male_Y"] * S[n,"male"] + 
                          b["c_parent_educ_Y"] * S[n,"parent_educ"] +
                          b["c_A_Y"] * 1 +
                          b["c_L1_Y"] * S[n,"L1"] +
                          b["c_M_Y"] * S[n,"M"] + 
                          b["c_AM_Y"] * 1 * S[n,"M"] * b["A.M.inter"] ) *
                        (( b["b_M"] + 
                             b["b_male_M"] * S[n,"male"] + 
                             b["b_parent_educ_M"] * S[n,"parent_educ"] + 
                             b["b_L1_M"] * S[n,"L1"] +
                             b["b_A_M"] * 1 )^( S[n,"M"] )) *
                        (( 1 - (b["b_M"] + 
                                  b["b_male_M"] * S[n,"male"] + 
                                  b["b_parent_educ_M"] * S[n,"parent_educ"] + 
                                  b["b_L1_M"] * S[n,"L1"] +
                                  b["b_A_M"] * 1) )^( 1 - S[n,"M"] )) *
                        (( b["b_L1"] +
                             b["b_male_L1"] * S[n,"male"] +  
                             b["b_parent_L1"] * S[n,"parent_educ"] +
                             b["b_A_L1"] * 1)^( S[n,"L1"] )) *
                        (( 1 - ( b["b_L1"] +
                                   b["b_male_L1"] * S[n,"male"] +  
                                   b["b_parent_L1"] * S[n,"parent_educ"] +
                                   b["b_A_L1"] * 1))^( 1 - S[n,"L1"] )) ) - 
                      ( ( b["mu_Y"] + 
                            b["c_male_Y"] * S[n,"male"] + 
                            b["c_parent_educ_Y"] * S[n,"parent_educ"] +
                            b["c_A_Y"] * 0 +
                            b["c_L1_Y"] * S[n,"L1"] +
                            b["c_M_Y"] * S[n,"M"] + 
                            b["c_AM_Y"] * 0 * S[n,"M"] * b["A.M.inter"] ) *
                          (( b["b_M"] + 
                               b["b_male_M"] * S[n,"male"] + 
                               b["b_parent_educ_M"] * S[n,"parent_educ"] + 
                               b["b_L1_M"] * S[n,"L1"] +
                               b["b_A_M"] * 0 )^( S[n,"M"] )) *
                          (( 1 - (b["b_M"] + 
                                    b["b_male_M"] * S[n,"male"] + 
                                    b["b_parent_educ_M"] * S[n,"parent_educ"] + 
                                    b["b_L1_M"] * S[n,"L1"] +
                                    b["b_A_M"] * 0) )^( 1 - S[n,"M"] )) ) *
                      (( b["b_L1"] +
                           b["b_male_L1"] * S[n,"male"] +  
                           b["b_parent_L1"] * S[n,"parent_educ"] +
                           b["b_A_L1"] * 0)^( S[n,"L1"] )) *
                      (( 1 - ( b["b_L1"] +
                                 b["b_male_L1"] * S[n,"male"] +  
                                 b["b_parent_L1"] * S[n,"parent_educ"] +
                                 b["b_A_L1"] * 0))^( 1 - S[n,"L1"] )) ) *
      ((b["p_L0_male"])^(S[n,"male"])) * 
      ((1 - b["p_L0_male"])^(1 - S[n,"male"])) * 
      ((b["p_L0_parent_low_educ_lv"])^(S[n,"parent_educ"])) *
      ((1 - b["p_L0_parent_low_educ_lv"])^(1 - S[n,"parent_educ"])) 
  }
  
  ATE.qol <- sum(S[,"sum"])
  
  return(list(ATE.death = ATE.death, ATE.qol = ATE.qol))
}

######################## résultats obtenus : 
# sans terme d'interaction
true.ATE2.no.inter <- true.ATE.time.var.conf(interaction = 0)
true.ATE2.no.inter
# $ATE.death
# [1] 0.0752
# 
# $ATE.qol
# [1] -6.26

# avec un terme d'interaction
true.ATE2.with.inter <- true.ATE.time.var.conf(interaction = 1)
true.ATE2.with.inter
# $ATE.death
# [1] 0.089282
# 
# $ATE.qol
# [1] -8.607



#############################################################################################################################################
#############################################################################################################################################
### 2) Calcul de l'effet direct et indirect attendus théoriquement 
###    dans les données simulées (marginal randomized direct and indirect effets)
#############################################################################################################################################
#############################################################################################################################################

true.marg.random.time.var <- function(interaction = NULL) {
  b <- sim.param.time.varying.L(A.M.interaction = interaction)
  
  # marginal distribution of M (conditionnellement à L0)
  M.S <- cbind(expand.grid(c(0,1),c(0,1),c(0,1),c(0,1),c(0,1)), rep(NA,n=2^5))
  colnames(M.S) <- list("male","parent_educ","L1","M","A","sum")
  
  for (n in 1:32) {
    M.S[n,"sum"] <- (( b["b_M"] +                                                              
                         b["b_male_M"] * M.S[n,"male"] + 
                         b["b_parent_educ_M"] * M.S[n,"parent_educ"] + 
                         b["b_L1_M"] * M.S[n,"L1"] +
                         b["b_A_M"] * M.S[n,"A"])^( M.S[n,"M"] )) *
      (( 1 - (b["b_M"] + 
                b["b_male_M"] * M.S[n,"male"] + 
                b["b_parent_educ_M"] * M.S[n,"parent_educ"] + 
                b["b_L1_M"] * M.S[n,"L1"] +
                b["b_A_M"] * M.S[n,"A"]) )^( 1 - M.S[n,"M"] )) *
      (( b["b_L1"] +                                                           
           b["b_male_L1"] * M.S[n,"male"] +  
           b["b_parent_L1"] * M.S[n,"parent_educ"] +
           b["b_A_L1"] * M.S[n,"A"])^( M.S[n,"L1"] )) *
      (( 1 - ( b["b_L1"] +
                 b["b_male_L1"] * M.S[n,"male"] +  
                 b["b_parent_L1"] * M.S[n,"parent_educ"] +
                 b["b_A_L1"] * M.S[n,"A"]))^( 1 - M.S[n,"L1"] )) 
  }
  
  M0.A0.L00 <- sum(M.S[M.S[,"M"]==0 & M.S[,"A"]==0 & M.S[,"male"]==0 & M.S[,"parent_educ"]==0,"sum"])
  M0.A0.L01 <- sum(M.S[M.S[,"M"]==0 & M.S[,"A"]==0 & M.S[,"male"]==0 & M.S[,"parent_educ"]==1,"sum"])
  M0.A0.L10 <- sum(M.S[M.S[,"M"]==0 & M.S[,"A"]==0 & M.S[,"male"]==1 & M.S[,"parent_educ"]==0,"sum"])
  M0.A0.L11 <- sum(M.S[M.S[,"M"]==0 & M.S[,"A"]==0 & M.S[,"male"]==1 & M.S[,"parent_educ"]==1,"sum"])
  
  M1.A0.L00 <- sum(M.S[M.S[,"M"]==1 & M.S[,"A"]==0 & M.S[,"male"]==0 & M.S[,"parent_educ"]==0,"sum"])
  M1.A0.L01 <- sum(M.S[M.S[,"M"]==1 & M.S[,"A"]==0 & M.S[,"male"]==0 & M.S[,"parent_educ"]==1,"sum"])
  M1.A0.L10 <- sum(M.S[M.S[,"M"]==1 & M.S[,"A"]==0 & M.S[,"male"]==1 & M.S[,"parent_educ"]==0,"sum"])
  M1.A0.L11 <- sum(M.S[M.S[,"M"]==1 & M.S[,"A"]==0 & M.S[,"male"]==1 & M.S[,"parent_educ"]==1,"sum"])
  
  M0.A1.L00 <- sum(M.S[M.S[,"M"]==0 & M.S[,"A"]==1 & M.S[,"male"]==0 & M.S[,"parent_educ"]==0,"sum"])
  M0.A1.L01 <- sum(M.S[M.S[,"M"]==0 & M.S[,"A"]==1 & M.S[,"male"]==0 & M.S[,"parent_educ"]==1,"sum"])
  M0.A1.L10 <- sum(M.S[M.S[,"M"]==0 & M.S[,"A"]==1 & M.S[,"male"]==1 & M.S[,"parent_educ"]==0,"sum"])
  M0.A1.L11 <- sum(M.S[M.S[,"M"]==0 & M.S[,"A"]==1 & M.S[,"male"]==1 & M.S[,"parent_educ"]==1,"sum"])
  
  M1.A1.L00 <- sum(M.S[M.S[,"M"]==1 & M.S[,"A"]==1 & M.S[,"male"]==0 & M.S[,"parent_educ"]==0,"sum"])
  M1.A1.L01 <- sum(M.S[M.S[,"M"]==1 & M.S[,"A"]==1 & M.S[,"male"]==0 & M.S[,"parent_educ"]==1,"sum"])
  M1.A1.L10 <- sum(M.S[M.S[,"M"]==1 & M.S[,"A"]==1 & M.S[,"male"]==1 & M.S[,"parent_educ"]==0,"sum"])
  M1.A1.L11 <- sum(M.S[M.S[,"M"]==1 & M.S[,"A"]==1 & M.S[,"male"]==1 & M.S[,"parent_educ"]==1,"sum"])
  
  # binary outcome (death)
  S <- cbind(expand.grid(c(0,1),c(0,1),c(0,1),c(0,1)), rep(NA,n=2^4), rep(NA,n=2^4), rep(NA,n=2^4))
  colnames(S) <- list("male","parent_educ","L1","M","sum.psi11", "sum.psi10", "sum.psi00")
  for (n in 1:16) {
    S[n,"sum.psi11"] <-  ( b["b_Y"] +                                            # A=1
                             b["b_male_Y"] * S[n,"male"] + 
                             b["b_parent_educ_Y"] * S[n,"parent_educ"] + 
                             b["b_A_Y"] * 1 + 
                             b["b_L1_Y"] * S[n,"L1"] +
                             b["b_M_Y"] * S[n,"M"] +
                             b["b_AM_Y"] * 1 * S[n,"M"] * b["A.M.inter"] ) *
      ((M1.A1.L00*(S[n,"male"]==0)*(S[n,"parent_educ"]==0) +                # A'=1
          M1.A1.L01*(S[n,"male"]==0)*(S[n,"parent_educ"]==1) +
          M1.A1.L10*(S[n,"male"]==1)*(S[n,"parent_educ"]==0) + 
          M1.A1.L11*(S[n,"male"]==1)*(S[n,"parent_educ"]==1) )^( S[n,"M"] )) *
      ((M0.A1.L00*(S[n,"male"]==0)*(S[n,"parent_educ"]==0) +                
          M0.A1.L01*(S[n,"male"]==0)*(S[n,"parent_educ"]==1) +
          M0.A1.L10*(S[n,"male"]==1)*(S[n,"parent_educ"]==0) + 
          M0.A1.L11*(S[n,"male"]==1)*(S[n,"parent_educ"]==1) )^( 1 - S[n,"M"] )) *
      (( b["b_L1"] +                                                           # A=1
           b["b_male_L1"] * S[n,"male"] +  
           b["b_parent_L1"] * S[n,"parent_educ"] +
           b["b_A_L1"] * 1)^( S[n,"L1"] )) *
      (( 1 - ( b["b_L1"] +
                 b["b_male_L1"] * S[n,"male"] +  
                 b["b_parent_L1"] * S[n,"parent_educ"] +
                 b["b_A_L1"] * 1))^( 1 - S[n,"L1"] )) *
      ((b["p_L0_male"])^(S[n,"male"])) * 
      ((1 - b["p_L0_male"])^(1 - S[n,"male"])) * 
      ((b["p_L0_parent_low_educ_lv"])^(S[n,"parent_educ"])) *
      ((1 - b["p_L0_parent_low_educ_lv"])^(1 - S[n,"parent_educ"])) 
    
    S[n,"sum.psi10"] <-  ( b["b_Y"] +                                            # A=1
                             b["b_male_Y"] * S[n,"male"] + 
                             b["b_parent_educ_Y"] * S[n,"parent_educ"] + 
                             b["b_A_Y"] * 1 + 
                             b["b_L1_Y"] * S[n,"L1"] +
                             b["b_M_Y"] * S[n,"M"] +
                             b["b_AM_Y"] * 1 * S[n,"M"] * b["A.M.inter"] ) *
      ((M1.A0.L00*(S[n,"male"]==0)*(S[n,"parent_educ"]==0) +                # A'=0
          M1.A0.L01*(S[n,"male"]==0)*(S[n,"parent_educ"]==1) +
          M1.A0.L10*(S[n,"male"]==1)*(S[n,"parent_educ"]==0) + 
          M1.A0.L11*(S[n,"male"]==1)*(S[n,"parent_educ"]==1) )^( S[n,"M"] )) *
      ((M0.A0.L00*(S[n,"male"]==0)*(S[n,"parent_educ"]==0) +                
          M0.A0.L01*(S[n,"male"]==0)*(S[n,"parent_educ"]==1) +
          M0.A0.L10*(S[n,"male"]==1)*(S[n,"parent_educ"]==0) + 
          M0.A0.L11*(S[n,"male"]==1)*(S[n,"parent_educ"]==1) )^( 1 - S[n,"M"] )) *
      (( b["b_L1"] +                                                           # A=1
           b["b_male_L1"] * S[n,"male"] +  
           b["b_parent_L1"] * S[n,"parent_educ"] +
           b["b_A_L1"] * 1)^( S[n,"L1"] )) *
      (( 1 - ( b["b_L1"] +
                 b["b_male_L1"] * S[n,"male"] +  
                 b["b_parent_L1"] * S[n,"parent_educ"] +
                 b["b_A_L1"] * 1))^( 1 - S[n,"L1"] )) *
      ((b["p_L0_male"])^(S[n,"male"])) * 
      ((1 - b["p_L0_male"])^(1 - S[n,"male"])) * 
      ((b["p_L0_parent_low_educ_lv"])^(S[n,"parent_educ"])) *
      ((1 - b["p_L0_parent_low_educ_lv"])^(1 - S[n,"parent_educ"])) 
    
    S[n,"sum.psi00"] <-  ( b["b_Y"] +                                            # A=0
                             b["b_male_Y"] * S[n,"male"] + 
                             b["b_parent_educ_Y"] * S[n,"parent_educ"] + 
                             b["b_A_Y"] * 0 + 
                             b["b_L1_Y"] * S[n,"L1"] +
                             b["b_M_Y"] * S[n,"M"] +
                             b["b_AM_Y"] * 0 * S[n,"M"] * b["A.M.inter"] ) *
      ((M1.A0.L00*(S[n,"male"]==0)*(S[n,"parent_educ"]==0) +                # A'=0
          M1.A0.L01*(S[n,"male"]==0)*(S[n,"parent_educ"]==1) +
          M1.A0.L10*(S[n,"male"]==1)*(S[n,"parent_educ"]==0) + 
          M1.A0.L11*(S[n,"male"]==1)*(S[n,"parent_educ"]==1) )^( S[n,"M"] )) *
      ((M0.A0.L00*(S[n,"male"]==0)*(S[n,"parent_educ"]==0) +                
          M0.A0.L01*(S[n,"male"]==0)*(S[n,"parent_educ"]==1) +
          M0.A0.L10*(S[n,"male"]==1)*(S[n,"parent_educ"]==0) + 
          M0.A0.L11*(S[n,"male"]==1)*(S[n,"parent_educ"]==1) )^( 1 - S[n,"M"] )) *
      (( b["b_L1"] +                                                           # A=0
           b["b_male_L1"] * S[n,"male"] +  
           b["b_parent_L1"] * S[n,"parent_educ"] +
           b["b_A_L1"] * 0)^( S[n,"L1"] )) *
      (( 1 - ( b["b_L1"] +
                 b["b_male_L1"] * S[n,"male"] +  
                 b["b_parent_L1"] * S[n,"parent_educ"] +
                 b["b_A_L1"] * 0))^( 1 - S[n,"L1"] )) *
      ((b["p_L0_male"])^(S[n,"male"])) * 
      ((1 - b["p_L0_male"])^(1 - S[n,"male"])) * 
      ((b["p_L0_parent_low_educ_lv"])^(S[n,"parent_educ"])) *
      ((1 - b["p_L0_parent_low_educ_lv"])^(1 - S[n,"parent_educ"])) 
  }
  
  mrNDE.death <- sum(S[,"sum.psi10"]) - sum(S[,"sum.psi00"])
  mrNIE.death <- sum(S[,"sum.psi11"]) - sum(S[,"sum.psi10"])
  
  # quantitative outcome (QoL)
  S <- cbind(expand.grid(c(0,1),c(0,1),c(0,1),c(0,1)), rep(NA,n=2^4), rep(NA,n=2^4), rep(NA,n=2^4))
  colnames(S) <- list("male","parent_educ","L1","M","sum.psi11", "sum.psi10", "sum.psi00")
  for (n in 1:16) {
    S[n,"sum.psi11"] <-  ( b["mu_Y"] +                                            # A=1
                             b["c_male_Y"] * S[n,"male"] + 
                             b["c_parent_educ_Y"] * S[n,"parent_educ"] + 
                             b["c_A_Y"] * 1 + 
                             b["c_L1_Y"] * S[n,"L1"] +
                             b["c_M_Y"] * S[n,"M"] +
                             b["c_AM_Y"] * 1 * S[n,"M"] * b["A.M.inter"] ) *
      ((M1.A1.L00*(S[n,"male"]==0)*(S[n,"parent_educ"]==0) +                # A'=1
          M1.A1.L01*(S[n,"male"]==0)*(S[n,"parent_educ"]==1) +
          M1.A1.L10*(S[n,"male"]==1)*(S[n,"parent_educ"]==0) + 
          M1.A1.L11*(S[n,"male"]==1)*(S[n,"parent_educ"]==1) )^( S[n,"M"] )) *
      ((M0.A1.L00*(S[n,"male"]==0)*(S[n,"parent_educ"]==0) +                
          M0.A1.L01*(S[n,"male"]==0)*(S[n,"parent_educ"]==1) +
          M0.A1.L10*(S[n,"male"]==1)*(S[n,"parent_educ"]==0) + 
          M0.A1.L11*(S[n,"male"]==1)*(S[n,"parent_educ"]==1) )^( 1 - S[n,"M"] )) *
      (( b["b_L1"] +                                                           # A=1
           b["b_male_L1"] * S[n,"male"] +  
           b["b_parent_L1"] * S[n,"parent_educ"] +
           b["b_A_L1"] * 1)^( S[n,"L1"] )) *
      (( 1 - ( b["b_L1"] +
                 b["b_male_L1"] * S[n,"male"] +  
                 b["b_parent_L1"] * S[n,"parent_educ"] +
                 b["b_A_L1"] * 1))^( 1 - S[n,"L1"] )) *
      ((b["p_L0_male"])^(S[n,"male"])) * 
      ((1 - b["p_L0_male"])^(1 - S[n,"male"])) * 
      ((b["p_L0_parent_low_educ_lv"])^(S[n,"parent_educ"])) *
      ((1 - b["p_L0_parent_low_educ_lv"])^(1 - S[n,"parent_educ"])) 
    
    S[n,"sum.psi10"] <-  ( b["mu_Y"] +                                            # A=1
                             b["c_male_Y"] * S[n,"male"] + 
                             b["c_parent_educ_Y"] * S[n,"parent_educ"] + 
                             b["c_A_Y"] * 1 +
                             b["c_L1_Y"] * S[n,"L1"] +
                             b["c_M_Y"] * S[n,"M"] +
                             b["c_AM_Y"] * 1 * S[n,"M"] * b["A.M.inter"] ) *
      ((M1.A0.L00*(S[n,"male"]==0)*(S[n,"parent_educ"]==0) +                # A'=0
          M1.A0.L01*(S[n,"male"]==0)*(S[n,"parent_educ"]==1) +
          M1.A0.L10*(S[n,"male"]==1)*(S[n,"parent_educ"]==0) + 
          M1.A0.L11*(S[n,"male"]==1)*(S[n,"parent_educ"]==1) )^( S[n,"M"] )) *
      ((M0.A0.L00*(S[n,"male"]==0)*(S[n,"parent_educ"]==0) +                
          M0.A0.L01*(S[n,"male"]==0)*(S[n,"parent_educ"]==1) +
          M0.A0.L10*(S[n,"male"]==1)*(S[n,"parent_educ"]==0) + 
          M0.A0.L11*(S[n,"male"]==1)*(S[n,"parent_educ"]==1) )^( 1 - S[n,"M"] )) *
      (( b["b_L1"] +                                                           # A=1
           b["b_male_L1"] * S[n,"male"] +  
           b["b_parent_L1"] * S[n,"parent_educ"] +
           b["b_A_L1"] * 1)^( S[n,"L1"] )) *
      (( 1 - ( b["b_L1"] +
                 b["b_male_L1"] * S[n,"male"] +  
                 b["b_parent_L1"] * S[n,"parent_educ"] +
                 b["b_A_L1"] * 1))^( 1 - S[n,"L1"] )) *
      ((b["p_L0_male"])^(S[n,"male"])) * 
      ((1 - b["p_L0_male"])^(1 - S[n,"male"])) * 
      ((b["p_L0_parent_low_educ_lv"])^(S[n,"parent_educ"])) *
      ((1 - b["p_L0_parent_low_educ_lv"])^(1 - S[n,"parent_educ"])) 
    
    S[n,"sum.psi00"] <-  ( b["mu_Y"] +                                            # A=0
                             b["c_male_Y"] * S[n,"male"] + 
                             b["c_parent_educ_Y"] * S[n,"parent_educ"] + 
                             b["c_A_Y"] * 0 + 
                             b["c_L1_Y"] * S[n,"L1"] +
                             b["c_M_Y"] * S[n,"M"] +
                             b["c_AM_Y"] * 0 * S[n,"M"] * b["A.M.inter"] ) *
      ((M1.A0.L00*(S[n,"male"]==0)*(S[n,"parent_educ"]==0) +                # A'=0
          M1.A0.L01*(S[n,"male"]==0)*(S[n,"parent_educ"]==1) +
          M1.A0.L10*(S[n,"male"]==1)*(S[n,"parent_educ"]==0) + 
          M1.A0.L11*(S[n,"male"]==1)*(S[n,"parent_educ"]==1) )^( S[n,"M"] )) *
      ((M0.A0.L00*(S[n,"male"]==0)*(S[n,"parent_educ"]==0) +                
          M0.A0.L01*(S[n,"male"]==0)*(S[n,"parent_educ"]==1) +
          M0.A0.L10*(S[n,"male"]==1)*(S[n,"parent_educ"]==0) + 
          M0.A0.L11*(S[n,"male"]==1)*(S[n,"parent_educ"]==1) )^( 1 - S[n,"M"] )) *
      (( b["b_L1"] +                                                           # A=0
           b["b_male_L1"] * S[n,"male"] +  
           b["b_parent_L1"] * S[n,"parent_educ"] +
           b["b_A_L1"] * 0)^( S[n,"L1"] )) *
      (( 1 - ( b["b_L1"] +
                 b["b_male_L1"] * S[n,"male"] +  
                 b["b_parent_L1"] * S[n,"parent_educ"] +
                 b["b_A_L1"] * 0))^( 1 - S[n,"L1"] )) *
      ((b["p_L0_male"])^(S[n,"male"])) * 
      ((1 - b["p_L0_male"])^(1 - S[n,"male"])) * 
      ((b["p_L0_parent_low_educ_lv"])^(S[n,"parent_educ"])) *
      ((1 - b["p_L0_parent_low_educ_lv"])^(1 - S[n,"parent_educ"])) 
  }
  
  mrNDE.qol <- sum(S[,"sum.psi10"]) - sum(S[,"sum.psi00"])
  mrNIE.qol <- sum(S[,"sum.psi11"]) - sum(S[,"sum.psi10"])
  
  return(list(mrNDE.death = mrNDE.death, mrNIE.death = mrNIE.death, 
              mrNDE.qol = mrNDE.qol, mrNIE.qol = mrNIE.qol))
}



######################## résultats obtenus : 
# sans terme d'interaction
true.marg.random2.no.inter <- true.marg.random.time.var(interaction = 0)
true.marg.random2.no.inter
# $mrNDE.death
# [1] 0.064
# 
# $mrNIE.death
# [1] 0.0112
# 
# $mrNDE.qol
# [1] -5
# 
# $mrNIE.qol
# [1] -1.26


# avec terme d'interaction
true.marg.random2.with.inter <- true.marg.random.time.var(interaction = 1)
true.marg.random2.with.inter
# $mrNDE.death
# [1] 0.073882
# 
# $mrNIE.death
# [1] 0.0154
# 
# $mrNDE.qol
# [1] -6.647
# 
# $mrNIE.qol
# [1] -1.96



# on a plus besoin des fonctions qui calculent les vrais paramètres : 
rm(true.ATE.time.var.conf, true.marg.random.time.var)



#############################################################################################################################################
#############################################################################################################################################
### 3) Fonction permettant de simuler une base de données
#############################################################################################################################################
#############################################################################################################################################

gen.data.time.varying.L <- function(N, A.M.inter) { # input parameters are the sample size N and the presence of A*M interaction with A.M.inter = 0 or 1
  
  b <- sim.param.time.varying.L(A.M.interaction = A.M.inter)
  
  # baseline confounders: parent's educational level=L0_parent_low_educ_lv & sex=L0_male
  L0_male <- rbinom(N, size = 1, prob = b["p_L0_male"]) 
  L0_parent_low_educ_lv <- rbinom(N, size = 1, prob = b["p_L0_parent_low_educ_lv"])  
  
  # exposure: A0_ace
  A0_ace <- rbinom(N, size = 1, prob =  b["b_A"] + 
                     b["b_male_A"] * L0_male + 
                     b["b_parent_educ_A"] * L0_parent_low_educ_lv ) 
  
  # intermediate confounder between M_smoking and Y, not affected by A0 L1
  L1 <- rbinom(N, size = 1, prob = b["b_L1"] + 
                 b["b_male_L1"] * L0_male + 
                 b["b_parent_L1"] * L0_parent_low_educ_lv + 
                 b["b_A_L1"] * A0_ace )
  
  # mediator: M_smoking
  M_smoking <- rbinom(N, size = 1, prob = b["b_M"] + 
                        b["b_male_M"] * L0_male + 
                        b["b_parent_educ_M"] * L0_parent_low_educ_lv + 
                        b["b_A_M"] * A0_ace +
                        b["b_L1_M"] * L1) 
  
  # Y_death 
  Y_death <- rbinom(N, size = 1, prob = b["b_Y"] + 
                      b["b_male_Y"] * L0_male + 
                      b["b_parent_educ_Y"] * L0_parent_low_educ_lv + 
                      b["b_A_Y"] * A0_ace + 
                      b["b_L1_Y"] * L1 +
                      b["b_M_Y"] * M_smoking +
                      b["b_AM_Y"] * A0_ace * M_smoking * A.M.inter ) 
  
  # Y_qol 
  Y_qol <- ( b["mu_Y"] + 
               b["c_male_Y"] * L0_male + 
               b["c_parent_educ_Y"] * L0_parent_low_educ_lv +
               b["c_A_Y"] * A0_ace +
               b["c_L1_Y"] * L1 +
               b["c_M_Y"] * M_smoking + 
               b["c_AM_Y"] * A0_ace * M_smoking * A.M.inter ) + 
    rnorm(N, mean = 0, sd = b["sd_Y"])
  
  # data.frame
  data.sim <- data.frame(L0_male, L0_parent_low_educ_lv, A0_ace, L1, M_smoking, Y_death, Y_qol)
  
  return( data.sim )
}


############ exemple de données simulées : 
set.seed(1234)
data.sim <- gen.data.time.varying.L(N=10000, A.M.inter=0)

summary(data.sim)





#############################################################################################################################################
#############################################################################################################################################
### 4) Estimation des rmNDE et rmNIE par g-computation iterative
#############################################################################################################################################
#############################################################################################################################################

### on veut calculer l'effet direct naturel (marginal randomized Natural direct effect)
### Psi_mrNDE = E(Y_{A=1,m=Gamma_0}) - E(Y_{A=0,m=Gamma_0})
### et l'effet indirect naturel (marginal randomized Natural indirect effect)
### Psi_mrNIE = E(Y_{A=1,m=Gamma_1}) - E(Y_{A=1,m=Gamma_0})

### a) Treatment mechanism (pour le médiateur)
###    il faut commencer par estimer la distribution contrefactuelle de M sous les scénarios fictifs do(A=1) et do(A=0)

#      a.i) Il s'agit de distributions marginales, conditionnellement à L0 (ce sont des scores de propension, on les notes g() par convention)
g_M_model <- glm(M_smoking ~ L0_male + L0_parent_low_educ_lv + A0_ace, family = "binomial", data = data.sim, )
# note : ici, on n'ajuste pas sur L1 car on veut l'effet marginal de A -> M, conditionnellement à L0

#      a.ii) Ensuite prédire la distribution attendue si toute la population est exposée à A=1
data.sim.A1 <- data.sim
data.sim.A1$A0_ace <- 1
g_M_1 <- predict(g_M_model, newdata = data.sim.A1, type="response")

#                      la distribution attendue si toute la population est exposée à A=0
data.sim.A0 <- data.sim
data.sim.A0$A0_ace <- 0
g_M_0 <- predict(g_M_model, newdata = data.sim.A0, type="response")


### b) Q function for the exposure to the mediator

###    b.i) On estime le modèle du critère de jugement en fonction de l'exposition au médiateur M
Q2_model <- glm(Y_death ~ L0_male + L0_parent_low_educ_lv + A0_ace + L1 + M_smoking, family = "binomial", data = data.sim)

###    b.ii) Pour chaque individu, prédire les résultats attendus si toute la population était exposée à do(M=0) ou à do(M=1)
data.sim.M0 <- data.sim.M1 <- data.sim
data.sim.M1$M_smoking  <- 1
data.sim.M0$M_smoking  <- 0

Q2_pred_M1 <- predict(Q2_model, newdata = data.sim.M1, type = "response")
Q2_pred_M0 <- predict(Q2_model, newdata = data.sim.M0, type = "response")

###    b.iii) Calculer une somme pondérée des Q2 prédits, avec les poids obtenus par g_M_1 ou par g_M_0 
# prédiction sous do(A=1)
Q2_Gamma_A1 <- Q2_pred_M1 * g_M_1 + Q2_pred_M0 * (1 - g_M_1)
Q2_Gamma_A1_bis <- Q2_pred_M1 * mean(g_M_1) + Q2_pred_M0 * (1 - mean(g_M_1))


# prédiction sous do(A=0)
Q2_Gamma_A0 <- Q2_pred_M1 * g_M_0 + Q2_pred_M0 * (1 - g_M_0)
Q2_Gamma_A0_bis <- Q2_pred_M1 * mean(g_M_0) + Q2_pred_M0 * (1 - mean(g_M_0))


### c) Q function for the exposure to the initial exposure 

###    c.i) Estimer Q1 à l'aide d'une régression quasi-logistique des valeurs Q2_gamma_Aa prédites à l'étape précédente
Q1_model_A1 <- glm(Q2_Gamma_A1 ~ L0_male + L0_parent_low_educ_lv + A0_ace, family = "quasibinomial", data = data.sim)
Q1_model_A0 <- glm(Q2_Gamma_A0 ~ L0_male + L0_parent_low_educ_lv + A0_ace, family = "quasibinomial", data = data.sim)

Q1_model_A1_bis <- glm(Q2_Gamma_A1_bis ~ L0_male + L0_parent_low_educ_lv + A0_ace, family = "quasibinomial", data = data.sim)
Q1_model_A0_bis <- glm(Q2_Gamma_A0_bis ~ L0_male + L0_parent_low_educ_lv + A0_ace, family = "quasibinomial", data = data.sim)

###    c.ii) Utiliser le modèle Q1_model_A1 pour prédire le résultat attendu si toute la population était exposée à do(A=1)
Q1_pred_A1_Gamma_A1 <- predict(Q1_model_A1, newdata = data.sim.A1, type = "response")
Q1_pred_A1_Gamma_A1_bis <- predict(Q1_model_A1_bis, newdata = data.sim.A1, type = "response")

###    c.iii) Utiliser le modèle Q1_model_A0 pour prédire le résultat attendu si toute la population était exposée à do(A=1) ou à do(A=0)
Q1_pred_A1_Gamma_A0 <- predict(Q1_model_A0, newdata = data.sim.A1, type = "response")
Q1_pred_A0_Gamma_A0 <- predict(Q1_model_A0, newdata = data.sim.A0, type = "response")

Q1_pred_A1_Gamma_A0_bis <- predict(Q1_model_A0_bis, newdata = data.sim.A1, type = "response")
Q1_pred_A0_Gamma_A0_bis <- predict(Q1_model_A0_bis, newdata = data.sim.A0, type = "response")

### d) Estimer les effets directs et indirects naturels, à partir de la moyenne des Q1 prédits à l'étape précédente : 

### Psi_mrNDE = E(Y_{A=1,m=Gamma_0}) - E(Y_{A=0,m=Gamma_0})
Psi_mrNDE <- mean(Q1_pred_A1_Gamma_A0) - mean(Q1_pred_A0_Gamma_A0)
Psi_mrNDE
# [1] 0.0625841

### Psi_mrNIE = E(Y_{A=1,m=Gamma_1}) - E(Y_{A=1,m=Gamma_0})
Psi_mrNIE <- mean(Q1_pred_A1_Gamma_A1) - mean(Q1_pred_A1_Gamma_A0)
Psi_mrNIE
# [1] 0.009845864

mean(Q1_pred_A1_Gamma_A0_bis) - mean(Q1_pred_A0_Gamma_A0_bis)
# [1] 0.06261304
mean(Q1_pred_A1_Gamma_A1_bis) - mean(Q1_pred_A1_Gamma_A0_bis)
# [1] 0.009793567
# donc on n'a pas exactement les mêmes résultats si on tire au sort dans G_{a*} ou dans G_{a*|L0} +++


######################### Question : 
######################### est-ce que l'estimation est plus précise si on tire au sort dans G_{a*|L0} plutôt que G_{a*} ?
g.comp.direct.indirect.estimation <- function(data) {
  ### a) Treatment mechanism (pour le médiateur)
  ###    il faut commencer par estimer la distribution contrefactuelle de M sous les scénarios fictifs do(A=1) et do(A=0)
  
  #      a.i) Il s'agit de distributions marginales, conditionnellement à L0 (ce sont des scores de propension, on les notes g() par convention)
  g_M_model <- glm(M_smoking ~ L0_male + L0_parent_low_educ_lv + A0_ace, family = "binomial", data = data, )
  # note : ici, on n'ajuste pas sur L1 car on veut l'effet marginal de A -> M, conditionnellement à L0
  
  #      a.ii) Ensuite prédire la distribution attendue si toute la population est exposée à A=1
  data.A1 <- data
  data.A1$A0_ace <- 1
  g_M_1 <- predict(g_M_model, newdata = data.A1, type="response")
  
  #                      la distribution attendue si toute la population est exposée à A=0
  data.A0 <- data
  data.A0$A0_ace <- 0
  g_M_0 <- predict(g_M_model, newdata = data.A0, type="response")
  
  
  ### b) Q function for the exposure to the mediator
  
  ###    b.i) On estime le modèle du critère de jugement en fonction de l'exposition au médiateur M
  Q2_model <- glm(Y_death ~ L0_male + L0_parent_low_educ_lv + A0_ace + L1 + M_smoking + A0_ace * M_smoking, family = "binomial", data = data)
  
  ###    b.ii) Pour chaque individu, prédire les résultats attendus si toute la population était exposée à do(M=0) ou à do(M=1)
  data.M0 <- data.M1 <- data
  data.M1$M_smoking  <- 1
  data.M0$M_smoking  <- 0
  
  Q2_pred_M1 <- predict(Q2_model, newdata = data.M1, type = "response")
  Q2_pred_M0 <- predict(Q2_model, newdata = data.M0, type = "response")
  
  ###    b.iii) Calculer une somme pondérée des Q2 prédits, avec les poids obtenus par g_M_1 ou par g_M_0 
  # prédiction sous do(A=1)
  Q2_Gamma_A1 <- Q2_pred_M1 * g_M_1 + Q2_pred_M0 * (1 - g_M_1)
  Q2_Gamma_A1_bis <- Q2_pred_M1 * mean(g_M_1) + Q2_pred_M0 * (1 - mean(g_M_1))
  
  
  # prédiction sous do(A=0)
  Q2_Gamma_A0 <- Q2_pred_M1 * g_M_0 + Q2_pred_M0 * (1 - g_M_0)
  Q2_Gamma_A0_bis <- Q2_pred_M1 * mean(g_M_0) + Q2_pred_M0 * (1 - mean(g_M_0))
  
  
  ### c) Q function for the exposure to the initial exposure 
  
  ###    c.i) Estimer Q1 à l'aide d'une régression quasi-logistique des valeurs Q2_gamma_Aa prédites à l'étape précédente
  Q1_model_A1 <- glm(Q2_Gamma_A1 ~ L0_male + L0_parent_low_educ_lv + A0_ace, family = "quasibinomial", data = data)
  Q1_model_A0 <- glm(Q2_Gamma_A0 ~ L0_male + L0_parent_low_educ_lv + A0_ace, family = "quasibinomial", data = data)
  
  Q1_model_A1_bis <- glm(Q2_Gamma_A1_bis ~ L0_male + L0_parent_low_educ_lv + A0_ace, family = "quasibinomial", data = data)
  Q1_model_A0_bis <- glm(Q2_Gamma_A0_bis ~ L0_male + L0_parent_low_educ_lv + A0_ace, family = "quasibinomial", data = data)
  
  ###    c.ii) Utiliser le modèle Q1_model_A1 pour prédire le résultat attendu si toute la population était exposée à do(A=1)
  Q1_pred_A1_Gamma_A1 <- predict(Q1_model_A1, newdata = data.A1, type = "response")
  Q1_pred_A1_Gamma_A1_bis <- predict(Q1_model_A1_bis, newdata = data.A1, type = "response")
  
  ###    c.iii) Utiliser le modèle Q1_model_A0 pour prédire le résultat attendu si toute la population était exposée à do(A=1) ou à do(A=0)
  Q1_pred_A1_Gamma_A0 <- predict(Q1_model_A0, newdata = data.A1, type = "response")
  Q1_pred_A0_Gamma_A0 <- predict(Q1_model_A0, newdata = data.A0, type = "response")
  
  Q1_pred_A1_Gamma_A0_bis <- predict(Q1_model_A0_bis, newdata = data.A1, type = "response")
  Q1_pred_A0_Gamma_A0_bis <- predict(Q1_model_A0_bis, newdata = data.A0, type = "response")
  
  ### d) Estimer les effets directs et indirects naturels, à partir de la moyenne des Q1 prédits à l'étape précédente : 
  
  ### Psi_mrNDE = E(Y_{A=1,m=Gamma_0}) - E(Y_{A=0,m=Gamma_0})
  Psi_mrNDE <- mean(Q1_pred_A1_Gamma_A0) - mean(Q1_pred_A0_Gamma_A0)
  
  ### Psi_mrNIE = E(Y_{A=1,m=Gamma_1}) - E(Y_{A=1,m=Gamma_0})
  Psi_mrNIE <- mean(Q1_pred_A1_Gamma_A1) - mean(Q1_pred_A1_Gamma_A0)
  
  Psi_mrNDE_bis <- mean(Q1_pred_A1_Gamma_A0_bis) - mean(Q1_pred_A0_Gamma_A0_bis)
  Psi_mrNIE_bis <- mean(Q1_pred_A1_Gamma_A1_bis) - mean(Q1_pred_A1_Gamma_A0_bis)
  
  return(list(Psi_mrNDE_Astar_condL0 = Psi_mrNDE, Psi_mrNIE_Astar_condL0 = Psi_mrNIE,
              Psi_mrNDE_Astar = Psi_mrNDE_bis, Psi_mrNIE_Astar = Psi_mrNIE_bis))
}

# tableau pour enregistrer les résultats
results.gcomp <- matrix(nrow = 1000, ncol = 4)
colnames(results.gcomp) <- c("Psi_mrNDE_Astar_condL0", "Psi_mrNIE_Astar_condL0",
                             "Psi_mrNDE_Astar", "Psi_mrNIE_Astar")

# simu 
set.seed(1234)
for (j in 1:1000){
  data.sim <- gen.data.time.varying.L(N=10000, A.M.inter=1)
  
  results <- g.comp.direct.indirect.estimation(data.sim)
  
  results.gcomp[j,"Psi_mrNDE_Astar_condL0"] <- results$Psi_mrNDE_Astar_condL0 - true.marg.random2.with.inter$mrNDE.death
  results.gcomp[j,"Psi_mrNIE_Astar_condL0"] <- results$Psi_mrNIE_Astar_condL0 - true.marg.random2.with.inter$mrNIE.death
  results.gcomp[j,"Psi_mrNDE_Astar"] <- results$Psi_mrNDE_Astar - true.marg.random2.with.inter$mrNDE.death
  results.gcomp[j,"Psi_mrNIE_Astar"] <- results$Psi_mrNIE_Astar - true.marg.random2.with.inter$mrNIE.death  
}

boxplot(results.gcomp) 
abline(h = 0)
plot(results.gcomp[,1] ~ results.gcomp[,3])
plot(results.gcomp[,2] ~ results.gcomp[,4])
# les deux méthodes donne des résultats quasiment identiques +++




#############################################################################################################################################
#############################################################################################################################################
### 5) Estimation en utilisant le package stremr
#############################################################################################################################################
#############################################################################################################################################
library(stremr)
library(data.table)

# base à analyser : 
names(data.sim)
# [1] "L0_male"               "L0_parent_low_educ_lv" "A0_ace"                "L1"                    "M_smoking"             "Y_death"              
# [7] "Y_qol"

# stremr est plutôt développer pour des analyses de survies (je ne suis pas sur que l'on puisse l'utiliser pour l'outcome quantitatif)
DATA_stremr <-  subset(data.sim, 
                       select = - c(Y_death) ) # on exclue la variable d'outcome continu des données, car elle est inutile
names(DATA_stremr)

### structure causale hypothétique : 
#
#                                 L1               Y_death
#                               /   \            / 
#                         A0_ace      M_smoking    
#                       /        
# L0_male, L0_parent_low_educ_lv

### on a besoin d'une base de données qui a ce format :

# ID | L0  | t   | L_t   |   A.t           | Y.tplus1 | A.tminus1
# -----------------------------------------------------------------------
#    | L0  | 0   |       |   A0_ace        | Y_death  | 
#    | L0  | 1   | L1    |   M_smoking     | Y_death  | A0_ace


DATA_stremr$ID <- 1:nrow(DATA_stremr)
DATA_stremr$t0 <- rep(0, nrow(DATA_stremr))
DATA_stremr$t1 <- rep(1, nrow(DATA_stremr))
DATA_stremr$alive <- rep(0, nrow(DATA_stremr))

dataDT0 <- data.frame(ID = DATA_stremr$ID, L0_male = DATA_stremr$L0_male, L0_parent_low_educ_lv = DATA_stremr$L0_parent_low_educ_lv,
                      t = DATA_stremr$t0, 
                      L1  = rep(NA, nrow(DATA_stremr)), 
                      A.t = DATA_stremr$A0_ace,
                      Y.tplus1 = DATA_stremr$alive, # il sont tous vivants avant l'exposition à M_smoking, pas de censure dans cet exemple
                      A.tminus1 = rep(NA, nrow(DATA_stremr)) )
                      
dataDT1 <- data.frame(ID = DATA_stremr$ID, L0_male = DATA_stremr$L0_male, L0_parent_low_educ_lv = DATA_stremr$L0_parent_low_educ_lv,
                      t = DATA_stremr$t1, 
                      L1  = DATA_stremr$L1, 
                      A.t = DATA_stremr$M_smoking,
                      Y.tplus1 = (DATA_stremr$Y_qol - min(DATA_stremr$Y_qol))/(max(DATA_stremr$Y_qol) - min(DATA_stremr$Y_qol)) ,
                      A.tminus1 = DATA_stremr$A0_ace )

dataDT <- rbind(dataDT0, dataDT1)

# on fait une data.table
dataDT <- data.table(dataDT)
sort.dataDT <- dataDT[order(ID, t)]
head(sort.dataDT, 20)

# les modèles \bar{Q} à appliquer à chaque temps des régressions itératives (on peut préciser une formule pour chaque temps t)
Qforms <- c("Qkplus1 ~ L0_male + L0_parent_low_educ_lv + A.t ",
            "Qkplus1 ~ L0_male + L0_parent_low_educ_lv + A.t + L1 + A.tminus1")

# colonnes indiquant les expositions contrefactuelles à A0 et M_smoking
# on a 3 expositions contrefactuelles à A0 et M pour refaire l'analyse précédente : 
#   pred_A1_Gamma_A0
#   pred_A0_Gamma_A0
#   pred_A1_Gamma_A1
sort.dataDT[, ("p_A1_Gamma_A0") := NA]
sort.dataDT[, ("p_A0_Gamma_A0") := NA]
sort.dataDT[, ("p_A1_Gamma_A1") := NA]


# pour l'exposition à A0_ace, il s'agit d'une exposition contrefactuelle "statique" : do(A0 = 0) ou do(A0 = 1)
# pour l'exposition à M comme "random draw in the distribution Gamma_A0 or Gamma_A1", 
#              on peut refaire les étapes a.1) et a.2) vue dans la méthode manuelle précédente

#      a.i) Il s'agit de distributions marginales, conditionnellement à L0 (ce sont des scores de propension, on les notes g() par convention)
g_M_model <- glm(M_smoking ~ L0_male + L0_parent_low_educ_lv + A0_ace, family = "binomial", data = data.sim, )
# note : ici, on n'ajuste pas sur L1 car on veut l'effet marginal de A -> M, conditionnellement à L0

#      a.ii) Ensuite prédire la distribution attendue si toute la population est exposée à A=1
data.sim.A1 <- data.sim
data.sim.A1$A0_ace <- 1
g_M_1 <- predict(g_M_model, newdata = data.sim.A1, type="response")

#                      la distribution attendue si toute la population est exposée à A=0
data.sim.A0 <- data.sim
data.sim.A0$A0_ace <- 0
g_M_0 <- predict(g_M_model, newdata = data.sim.A0, type="response")

sort.dataDT$p_A1_Gamma_A0[sort.dataDT$t==0] <- rep(1, nrow(data.sim))
sort.dataDT$p_A1_Gamma_A0[sort.dataDT$t==1] <- g_M_0

sort.dataDT$p_A0_Gamma_A0[sort.dataDT$t==0] <- rep(0, nrow(data.sim))
sort.dataDT$p_A0_Gamma_A0[sort.dataDT$t==1] <- g_M_0

sort.dataDT$p_A1_Gamma_A1[sort.dataDT$t==0] <- rep(1, nrow(data.sim))
sort.dataDT$p_A1_Gamma_A1[sort.dataDT$t==1] <- g_M_1

head(sort.dataDT, 20)
# OK

# ici, il n'y a pas de censure avant l'exposition à M1 (pas de décès intermédiaire), donc le tableau reste tel quel

#### import des données pour stremr
OData <- stremr::importData(sort.dataDT, ID = "ID", t_name = "t",  
                            covars = c("L0_male", "L0_parent_low_educ_lv", "L1"),
                            TRT = "A.t", OUTCOME = "Y.tplus1",
                            CENS = NULL, MONITOR = NULL)
get_data(OData)

# g-compuation pour estimer E(Y_{a=1, \Gamma_a=0})
gcomp_Y_a1_G0 <- fit_GCOMP(OData, tvals = c(0:1), intervened_TRT = "p_A1_Gamma_A0", Qforms = Qforms) 
gcomp_Y_a1_G0

# g-compuation pour estimer E(Y_{a=0, \Gamma_a=0})
gcomp_Y_a0_G0 <- fit_GCOMP(OData, tvals = c(0:1), intervened_TRT = "p_A0_Gamma_A0", Qforms = Qforms) 

# g-compuation pour estimer E(Y_{a=1, \Gamma_a=1})
gcomp_Y_a1_G1 <- fit_GCOMP(OData, tvals = c(0:1), intervened_TRT = "p_A1_Gamma_A1", Qforms = Qforms) 


###########################################
### Calcul des effets directs et indirects
Psi_mrNDE_GCOMP <- gcomp_Y_a1_G0$estimates$cum.inc[2] - gcomp_Y_a0_G0$estimates$cum.inc[2]
Psi_mrNDE_GCOMP
# [1] 0.0625841

### Psi_mrNIE = E(Y_{A=1,m=Gamma_1}) - E(Y_{A=1,m=Gamma_0})
Psi_mrNIE_GCOMP <- gcomp_Y_a1_G1$estimates$cum.inc[2] - gcomp_Y_a1_G0$estimates$cum.inc[2] 
Psi_mrNIE_GCOMP
# [1] 0.009845864

# On retrouve bien la même chose que le calcul manuel précédent 




#############################################################################################################################################
#############################################################################################################################################
### 6) Estimation en utilisant un MSM estimé par IPTW
#############################################################################################################################################
#############################################################################################################################################
set.seed(1234)
data.sim <- gen.data.time.varying.L(N=10000, A.M.inter=1)


iptw.direct.indirect <- function(data) {
  # 1. g functions
  g.A <- glm(A0_ace ~ 1, family = "binomial", data = data)
  g.A.L0 <- glm(A0_ace ~ L0_male + L0_parent_low_educ_lv, family = "binomial", data = data)
  g.M.A <- glm(M_smoking ~ A0_ace, family = "binomial", data = data)
  g.M.L <- glm(M_smoking ~ L0_male + L0_parent_low_educ_lv + A0_ace + L1, family = "binomial", data = data)
  
  pred.g1.A <- predict(g.A, type="response")
  pred.g0.A <- 1 - pred.g1.A
  
  pred.g1.A.L0 <- predict(g.A.L0, type="response")
  pred.g0.A.L0 <- 1 - pred.g1.A.L0
  
  pred.g1.M.A <- predict(g.M.A, type="response")
  pred.g0.M.A <- 1 - pred.g1.M.A
  
  pred.g1.M.L <- predict(g.M.L, type="response")
  pred.g0.M.L <- 1 - pred.g1.M.L
  
  gA <- gM.A <- gA.L <- gM.AL <- rep(NA, nrow(data))
  gA[data$A0_ace==1] <- pred.g1.A[data$A0_ace==1]
  gA[data$A0_ace==0] <- pred.g0.A[data$A0_ace==0]
  gA.L[data$A0_ace==1] <- pred.g1.A.L0[data$A0_ace==1]
  gA.L[data$A0_ace==0] <- pred.g0.A.L0[data$A0_ace==0]
  gM.A[data$M_smoking==1] <- pred.g1.M.A[data$M_smoking==1]
  gM.A[data$M_smoking==0] <- pred.g0.M.A[data$M_smoking==0]
  gM.AL[data$M_smoking==1] <- pred.g1.M.L[data$M_smoking==1]
  gM.AL[data$M_smoking==0] <- pred.g0.M.L[data$M_smoking==0]
  
  sw <- (gA * gM.A) / (gA.L * gM.AL)
  
  msm_y <- glm(Y_death ~ A0_ace + M_smoking + A0_ace:M_smoking, family = "gaussian", data = data, weights = sw)
  msm_m <- glm(M_smoking ~ A0_ace + L0_male + L0_parent_low_educ_lv, family = "binomial", data = data)
  
  iptw.EDN.L0 <- msm_y$coefficients["A0_ace"] + 
    (msm_y$coefficients["A0_ace:M_smoking"] * plogis(rep(msm_m$coefficients["(Intercept)"], nrow(data)) + 
                                                      msm_m$coefficients["L0_male"] * data$L0_male + 
                                                      msm_m$coefficients["L0_parent_low_educ_lv"] * data$L0_parent_low_educ_lv))
  
  iptw.EIN.L0 <- (msm_y$coefficients["M_smoking"] + msm_y$coefficients["A0_ace:M_smoking"]) *
    (plogis(rep(msm_m$coefficients["(Intercept)"], nrow(data)) + 
              msm_m$coefficients["A0_ace"] + 
              msm_m$coefficients["L0_male"] * data$L0_male + 
              msm_m$coefficients["L0_parent_low_educ_lv"] * data$L0_parent_low_educ_lv) - 
       plogis(rep(msm_m$coefficients["(Intercept)"], nrow(data)) + 
                msm_m$coefficients["L0_male"] * data$L0_male + 
                msm_m$coefficients["L0_parent_low_educ_lv"] * data$L0_parent_low_educ_lv))

  iptw.EDN <- mean(iptw.EDN.L0)
  iptw.EIN <- mean(iptw.EIN.L0)
  
  return(list(iptw.EDN = iptw.EDN, iptw.EIN = iptw.EIN))
}

# tableau pour enregistrer les résultats
results.iptw <- matrix(nrow = 1000, ncol = 2)
colnames(results.iptw) <- c("iptw.EDN", "iptw.EIN")

# simu 
set.seed(1234)
for (j in 1:1000){
  data.sim <- gen.data.time.varying.L(N=10000, A.M.inter=1)
  
  results <- iptw.direct.indirect(data.sim)
  
  results.iptw[j,"iptw.EDN"] <- results$iptw.EDN - true.marg.random2.with.inter$mrNDE.death
  results.iptw[j,"iptw.EIN"] <- results$iptw.EIN - true.marg.random2.with.inter$mrNIE.death

}

boxplot(results.iptw) 
abline(h = 0)
# intéressant, ici, sans problème de positivité, l'iptw fonctionne aussi bien que la g-computation

file_path <- "../Data/"
data <- data.frame(read.csv(paste(file_path, "data_sim.csv", sep = "")))
data <- subset(data, select = -c(Y_qol)) # remove Y_qol
head(data)

manip_res <- iptw.direct.indirect(data)
manip_res

# data[data$L0_male == 0 & data$L0_parent_low_educ_lv == 1 & data$A0_ace == 0 & 
# data$L1 == 1 & data$M_smoking == 1 & data$Y_death == 1, ] = 0

# En changeant certaines cases, les estimations varient
# mais restent assez proches de la valeur trouvée

##### 7.1) on modifie les paramètres de sorte à ce que les groupes d'exposition à A et à M aient des probabilités faibles dans certaines strates de L0
sim.param.time.varying.L <- function(A.M.interaction = NULL) {
  # L0
  p_L0_male <- 0.5
  p_L0_parent_low_educ_lv <- 0.65
  
  # A: A0_ace <- rbinom( 0.05 + 0.04 * L0_male + 0.06 * L0_parent_low_educ_lv ) 
  b_A <- 1/10000   # MODIFICATION ICI : baseline 1/10000
  b_male_A <- 0.04  # + 0.04 for the effect of L0_male -> A0_ace
  b_parent_educ_A <- 0.06  # +0.06 for the effect of L0_parent_low_educ_lv -> A0_ace
  
  # L1: L1 <- rbinom( 0.30 - 0.05 * L0_male + 0.08 * L0_parent_low_educ_lv + 0.2 * A0_ace ) 
  b_L1 <- 0.30   # reference prevalence is 30%
  b_male_L1 <- -0.05  # - 0.05 for the effect of L0_male -> L1
  b_parent_L1 <- +0.08 # + 0.08 for the effect of L0_parent_low_educ_lv -> L1
  b_A_L1 <- +0.2 # +0.2 for the effect of A0_ace -> L1
  
  # M: M_smoking <- rbinom( 0.2 + 0.05 * L0_male + 0.06 * L0_parent_low_educ_lv + 0.2 * L1 + 0.1 * A0_ace ) 
  b_M <- 1/10000 # MODIFICATION ICI : baseline 1/10000
  b_male_M <- 0.05 # +0.05 for the effect of L0_male -> M_smoking
  b_parent_educ_M <- 0.06 # +0.06 for the effect of L0_parent_low_educ_lv -> M_smoking
  b_A_M <- 0.1 # +0.10 for the effect of A0_ace -> M_smoking
  b_L1_M <- 0.2 # +0.2 for the effect of L1 -> M_smoking
  
  # Y binary: rbinom( 0.10 + 0.06 * L0_male + 0.04 * L0_parent_low_educ_lv + 0.05 * A0_ace + 0.07 * L1 + 0.08 * M_smoking +
  #                   0.03 * A0_ace * M_smoking * A.M.inter ) 
  b_Y <- 0.1 # reference prevalence is 10%
  b_male_Y <- 0.06 # +0.06 for the effect of L0_male -> Y
  b_parent_educ_Y <- 0.04 # +0.04 for the effect of L0_parent_low_educ_lv -> Y
  b_A_Y <- 0.05 # 0.05 for the effect of A0_ace -> Y
  b_L1_Y <- 0.07 # +0.07 for the effect of L1 -> Y
  b_M_Y <- 0.08 # 0.08 for the effect of M_smoking -> Y
  b_AM_Y <- 0.03 # 0.03 for the interaction effect A0_ace * M_smoking -> Y
  
  # Y continuous: (75 - 1 * L0_male - 3 * L0_parent_low_educ_lv - 4 * A0_ace -3.5 * L1 - 9 * M_smoking + 
  #             -5 * A0_ace * M_smoking * A.M.inter ) + rnorm(N, mean = 0, sd = 10)
  mu_Y <- 75 # reference mean for QoL
  c_male_Y <- -1 # -1 for the effect of L0_male -> Y
  c_parent_educ_Y <- -3 # -3 for the effect of L0_parent_low_educ_lv -> Y
  c_A_Y <- -4 # -4 for the effect of A0_ace -> Y
  c_L1_Y <- -5 # -5 for the effect of L1 -> Y
  c_M_Y <- -9 # -9 for the effect of M_smoking -> Y
  c_AM_Y <- -5  # - 5 for the interaction effect A0_ace * M_smoking  -> Y
  sd_Y <- 10 # standard deviation of the residuals
  
  # A*M interaction ?
  A.M.inter <- A.M.interaction
  
  coef <- c( p_L0_male = p_L0_male, p_L0_parent_low_educ_lv = p_L0_parent_low_educ_lv, 
             b_A = b_A, b_male_A = b_male_A, b_parent_educ_A = b_parent_educ_A, 
             b_L1 = b_L1, b_male_L1 = b_male_L1, b_parent_L1 = b_parent_L1, b_A_L1 = b_A_L1,
             b_M = b_M, b_male_M = b_male_M, b_parent_educ_M = b_parent_educ_M, b_L1_M = b_L1_M, b_A_M = b_A_M,
             b_Y = b_Y, b_male_Y = b_male_Y, b_parent_educ_Y = b_parent_educ_Y, b_A_Y = b_A_Y, b_L1_Y = b_L1_Y, b_M_Y = b_M_Y, b_AM_Y = b_AM_Y,
             mu_Y = mu_Y, c_male_Y = c_male_Y, c_parent_educ_Y = c_parent_educ_Y, c_A_Y = c_A_Y, c_L1_Y = c_L1_Y, c_M_Y = c_M_Y, c_AM_Y = c_AM_Y, 
             sd_Y = sd_Y, A.M.inter = A.M.inter)
  
  return(coef)
}

##### 7.2) les vraies valeurs attendues sont alors : 
# sans terme d'interaction
true.ATE2.no.inter <- true.ATE.time.var.conf(interaction = 0)
true.ATE2.no.inter
# $ATE.death
# [1] 0.0752
# 
# $ATE.qol
# [1] -6.26

# avec un terme d'interaction
true.ATE2.with.inter <- true.ATE.time.var.conf(interaction = 1)
true.ATE2.with.inter
# $ATE.death
# [1] 0.089282
# 
# $ATE.qol
# [1] -8.607
# ça ne change rien pour les effets attendus, dans ces modèles linéaires


true.marg.random2.no.inter <- true.marg.random.time.var(interaction = 0)
true.marg.random2.no.inter
# $mrNDE.death
# [1] 0.064
# 
# $mrNIE.death
# [1] 0.0112
# 
# $mrNDE.qol
# [1] -5
# 
# $mrNIE.qol
# [1] -1.26


# avec terme d'interaction
true.marg.random2.with.inter <- true.marg.random.time.var(interaction = 1)
true.marg.random2.with.inter
# $mrNDE.death
# [1] 0.067912
# 
# $mrNIE.death
# [1] 0.0154
# 
# $mrNDE.qol
# [1] -5.652
# 
# $mrNIE.qol
# [1] -1.96


##### 7.3) estimation dans une base et vérification de la distribution des poids
set.seed(1234)
data.sim <- gen.data.time.varying.L(N=10000, A.M.inter=1)

results <- iptw.direct.indirect(data.sim)

boxplot(results$sw)
# la distribution des poids est un peu plus dispersée


##### 7.4) simulation dans 1000 bases

# tableau pour enregistrer les résultats
results.iptw <- matrix(nrow = 1000, ncol = 2)
colnames(results.iptw) <- c("iptw.EDN", "iptw.EIN")


# simu 
set.seed(1234)
for (j in 1:1000) {
  data.sim <- gen.data.time.varying.L(N=10000, A.M.inter=1)
  
  results <- iptw.direct.indirect(data.sim)
  
  results.iptw[j,"iptw.EDN"] <- results$iptw.EDN - true.marg.random2.with.inter$mrNDE.death
  results.iptw[j,"iptw.EIN"] <- results$iptw.EIN - true.marg.random2.with.inter$mrNIE.death
  
}

boxplot(results.iptw) 
abline(h = 0)
# on voit que les résultats sont un peu plus biaisés et plus dispersés