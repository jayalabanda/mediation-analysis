library(tmle3)
library(sl3)

file_path <- "../../Data/"
data <- data.frame(read.csv(paste(file_path, "data_sim.csv", sep = "")))
data <- subset(data, select = -c(Y_qol)) # remove Y_qol
head(data)

npsem <- list(
  define_node("W", c("L0_male", "L0_parent_low_educ_lv")),
  define_node("A", c("A0_ace"), c("W")),
  define_node("Z", c("L1"), c("W", "A")),
  define_node("M", c("M_smoking"), c("A", "W")),
  define_node("Y", c("Y_death"), c("A", "W"))
)

tmle_task <- tmle3_Task$new(data, npsem = npsem)
head(tmle_task$get_tmle_node("Y"))
tmle_task$get_regression_task("Y")

lrnr_glm_fast <- make_learner(Lrnr_glm_fast)
lrnr_mean <- make_learner(Lrnr_mean)

factor_list <- list(
  define_lf(LF_emp, "W"),
  define_lf(LF_fit, "A", lrnr_glm_fast),
  define_lf(LF_fit, "Z", lrnr_glm_fast),
  define_lf(LF_fit, "M", lrnr_glm_fast),
  define_lf(LF_fit, "Y", lrnr_glm_fast, type = "mean")
)

likelihood_def <- Likelihood$new(factor_list)
likelihood <- likelihood_def$train(tmle_task)
print(likelihood)

likelihood_values <- likelihood$get_likelihoods(tmle_task, "Y")
head(likelihood_values)

intervention <- define_lf(LF_static, "A", value = 1)

cf_likelihood <- make_CF_Likelihood(likelihood, intervention)

cf_likelihood_values <- cf_likelihood$get_likelihoods(tmle_task, "A")
head(cf_likelihood_values)

cf_likelihood_tasks <- cf_likelihood$enumerate_cf_tasks(tmle_task)
head(cf_likelihood_tasks[[1]]$data)

tsm <- define_param(Param_TSM, likelihood, intervention)

updater <- tmle3_Update$new()
targeted_likelihood <- Targeted_Likelihood$new(likelihood, updater)

tsm <- define_param(Param_TSM, likelihood, intervention)
updater$tmle_params <- tsm

tmle_fit <- fit_tmle3(tmle_task, targeted_likelihood, tsm, updater)

print(tmle_fit)
#    type      param init_est tmle_est         se     lower     upper
# 1:  TSM E[Y_{A=1}] 0.277099 0.277099 0.01364148 0.2503622 0.3038358
#    psi_transformed lower_transformed upper_transformed
# 1:        0.277099         0.2503622         0.3038358

#######################################
library(data.table)
library(sl3)
library(tmle3)
library(tmle3mediate)

new_data <- data[sample(nrow(data), 2000), ]

node_list <- list(
  W = c("L0_male", "L0_parent_low_educ_lv"),
  A = c("A0_ace"),
  Z = c("L1"),
  M = c("M_smoking"),
  Y = c("Y_death")
)

processed <- process_missing(new_data, node_list)
new_data <- processed$data
node_list <- processed$node_list

# SL learners used for continuous data (the nuisance parameter Z)
enet_contin_learner <- Lrnr_glmnet$new(
  alpha = 0.5, family = "gaussian", nfolds = 3
)
lasso_contin_learner <- Lrnr_glmnet$new(
  alpha = 1, family = "gaussian", nfolds = 3
)
fglm_contin_learner <- Lrnr_glm_fast$new(family = gaussian())
mean_learner <- Lrnr_mean$new()
contin_learner_lib <- Stack$new(
  enet_contin_learner, lasso_contin_learner, fglm_contin_learner, mean_learner
)
sl_contin_learner <- Lrnr_sl$new(learners = contin_learner_lib)

# SL learners used for binary data (nuisance parameters G and E in this case)
enet_binary_learner <- Lrnr_glmnet$new(
  alpha = 0.5, family = "binomial", nfolds = 3
)
lasso_binary_learner <- Lrnr_glmnet$new(
  alpha = 1, family = "binomial", nfolds = 3
)
fglm_binary_learner <- Lrnr_glm_fast$new(family = binomial())
binary_learner_lib <- Stack$new(
  enet_binary_learner, lasso_binary_learner, fglm_binary_learner, mean_learner
)
sl_binary_learner <- Lrnr_sl$new(learners = binary_learner_lib)

# create list for treatment and outcome mechanism regressions
learner_list <- list(
  Y = sl_binary_learner,
  A = sl_binary_learner
)

tmle_spec_nie <- tmle_NIE(
  e_learners = Lrnr_cv$new(lasso_binary_learner, full_fit = TRUE),
  psi_Z_learners = Lrnr_cv$new(lasso_contin_learner, full_fit = TRUE),
  max_iter = 1
)

data_nie <- tmle3(
  tmle_spec_nie, new_data, node_list, learner_list
)
data_nie
#    type                  param   init_est   tmle_est       se     lower
# 1:  NIE NIE[Y_{A=1} - Y_{A=0}] 0.01371743 0.01373467 2.259678 -4.415152
#       upper psi_transformed lower_transformed upper_transformed
# 1: 4.442621      0.01373467         -4.415152          4.442621