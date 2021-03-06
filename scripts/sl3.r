## ----setup--------------------------------------------------------------------
library(data.table)
library(dplyr)
library(readr)
library(ggplot2)
library(SuperLearner)
library(origami)
library(sl3)
library(knitr)
library(kableExtra)

# load data set and take a peek
washb_data <- fread(
  paste0(
    "https://raw.githubusercontent.com/tlverse/tlverse-data/master/",
    "wash-benefits/washb_data.csv"
  ),
  stringsAsFactors = TRUE
)


## ----sl3_washb_example_table1, echo=FALSE-------------------------------------
if (knitr::is_latex_output()) {
  head(washb_data) %>%
    kable(format = "latex")
} else if (knitr::is_html_output()) {
  head(washb_data) %>%
    kable() %>%
    kable_styling(fixed_thead = TRUE) %>%
    scroll_box(width = "100%", height = "300px")
}


## ----task, warning=TRUE-------------------------------------------------------
# specify the outcome and covariates
outcome <- "whz"
covars <- colnames(washb_data)[-which(names(washb_data) == outcome)]

# create the sl3 task
washb_task <- make_sl3_Task(
  data = washb_data,
  covariates = covars,
  outcome = outcome
)


## ----task-examine-------------------------------------------------------------
washb_task


## ----task-folds-examine-------------------------------------------------------
length(washb_task$folds) # how many folds?

head(washb_task$folds[[1]]$training_set) # row indexes for fold 1 training
head(washb_task$folds[[1]]$validation_set) # row indexes for fold 1 validation

any(
  washb_task$folds[[1]]$training_set %in%
    washb_task$folds[[1]]$validation_set
)


## ----list-properties----------------------------------------------------------
sl3_list_properties()


## ----list-learners------------------------------------------------------------
sl3_list_learners("continuous")


## ----baselearners-------------------------------------------------------------
# choose base learners
lrn_glm <- make_learner(Lrnr_glm)
lrn_mean <- Lrnr_mean$new()


## ----extra-lrnr-awesome-------------------------------------------------------
lrn_lasso <- make_learner(Lrnr_glmnet) # alpha default is 1
lrn_ridge <- Lrnr_glmnet$new(alpha = 0)
lrn_enet.5 <- make_learner(Lrnr_glmnet, alpha = 0.5)

lrn_polspline <- Lrnr_polspline$new()

lrn_ranger100 <- make_learner(Lrnr_ranger, num.trees = 100)

lrn_hal_faster <- Lrnr_hal9001$new(max_degree = 2, reduce_basis = 0.05)

xgb_fast <- Lrnr_xgboost$new() # default with nrounds = 20 is pretty fast
xgb_50 <- Lrnr_xgboost$new(nrounds = 50)


## ----interaction-learner------------------------------------------------------
interactions <- list(c("elec", "tr"), c("tr", "hfiacat"))
# main terms as well as the interactions above will be included
lrn_interaction <- make_learner(Lrnr_define_interactions, interactions)


## ----interaction-pipe---------------------------------------------------------
# we already instantiated a linear model learner above, no need to do it again
lrn_glm_interaction <- make_learner(Pipeline, lrn_interaction, lrn_glm)
lrn_glm_interaction


## ----extra-lrnr-woah----------------------------------------------------------
lrn_bayesglm <- Lrnr_pkg_SuperLearner$new("SL.bayesglm")


## ----extra-lrnr-mindblown-svm, eval = FALSE-----------------------------------
## # I like to crock pot my SLs
## grid_params <- list(
##   cost = c(0.01, 0.1, 1, 10, 100, 1000),
##   gamma = c(0.001, 0.01, 0.1, 1),
##   kernel = c("polynomial", "radial", "sigmoid"),
##   degree = c(1, 2, 3)
## )
## grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
## svm_learners <- apply(grid, MARGIN = 1, function(tuning_params) {
##   do.call(Lrnr_svm$new, as.list(tuning_params))
## })


## ----extra-lrnr-mindblown-xgboost---------------------------------------------
grid_params <- list(
  max_depth = c(2, 4, 6),
  eta = c(0.001, 0.1, 0.3),
  nrounds = 100
)
grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
grid

xgb_learners <- apply(grid, MARGIN = 1, function(tuning_params) {
  do.call(Lrnr_xgboost$new, as.list(tuning_params))
})
xgb_learners


## ----carotene, eval=FALSE-----------------------------------------------------
## # Unlike xgboost, I have no idea how to tune a neural net or BART machine, so
## # I let caret take the reins
## lrnr_caret_nnet <- make_learner(Lrnr_caret, algorithm = "nnet")
## lrnr_caret_bartMachine <- make_learner(Lrnr_caret,
##   algorithm = "bartMachine",
##   method = "boot", metric = "Accuracy",
##   tuneLength = 10
## )


## ----stack--------------------------------------------------------------------
stack <- make_learner(
  Stack, lrn_glm, lrn_polspline, lrn_enet.5, lrn_ridge, lrn_lasso, xgb_50
)
stack


## ----alt-stack----------------------------------------------------------------
# named vector of learners first
learners <- c(lrn_glm, lrn_polspline, lrn_enet.5, lrn_ridge, lrn_lasso, xgb_50)
names(learners) <- c(
  "glm", "polspline", "enet.5", "ridge", "lasso", "xgboost50"
)
# next make the stack
stack <- make_learner(Stack, learners)
# now the names are pretty
stack


## ----alt-stack-cv-------------------------------------------------------------
cv_stack <- Lrnr_cv$new(stack)
cv_stack


## ----screener-----------------------------------------------------------------
miniforest <- Lrnr_ranger$new(
  num.trees = 20, write.forest = FALSE,
  importance = "impurity_corrected"
)

# learner must already be instantiated, we did this when we created miniforest
screen_rf <- Lrnr_screener_importance$new(learner = miniforest, num_screen = 5)
screen_rf

# which covariates are selected on the full data?
screen_rf$train(washb_task)


## ----screener-augment---------------------------------------------------------
keepme <- c("aged", "momage")
# screener must already be instantiated, we did this when we created screen_rf
screen_augment_rf <- Lrnr_screener_augment$new(
  screener = screen_rf, default_covariates = keepme
)
screen_augment_rf


## ----screener-coefs-----------------------------------------------------------
# we already instantiated a lasso learner above, no need to do it again
screen_lasso <- Lrnr_screener_coefs$new(learner = lrn_lasso, threshold = 0)
screen_lasso


## ----screener-pipe------------------------------------------------------------
screen_rf_pipe <- make_learner(Pipeline, screen_rf, stack)
screen_lasso_pipe <- make_learner(Pipeline, screen_lasso, stack)


## ----screeners-stack----------------------------------------------------------
# pretty names again
learners2 <- c(learners, screen_rf_pipe, screen_lasso_pipe)
names(learners2) <- c(names(learners), "randomforest_screen", "lasso_screen")

fancy_stack <- make_learner(Stack, learners2)
fancy_stack


## ----make-sl------------------------------------------------------------------
sl <- make_learner(Lrnr_sl, learners = fancy_stack)


## ----make-sl-discrete---------------------------------------------------------
discrete_sl_metalrn <- Lrnr_cv_selector$new()
discrete_sl <- Lrnr_sl$new(
  learners = fancy_stack,
  metalearner = discrete_sl_metalrn
)


## ----make-sl-plot, eval=FALSE-------------------------------------------------
## dt_sl <- delayed_learner_train(sl, washb_task)
## plot(dt_sl, color = FALSE, height = "400px", width = "90%")


## ----sl-----------------------------------------------------------------------
set.seed(4197)
sl_fit <- sl$train(washb_task)


## ----sl-predictions-----------------------------------------------------------
# we did it! now we have SL predictions
sl_preds <- sl_fit$predict()
head(sl_preds)
