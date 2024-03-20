library(doParallel)
library(dplyr) 
library(gridExtra)
library(lubridate)
library(progress)
library(purrr)
library(ranger)
library(stacks)
library(stringr)
library(tidyverse)
library(tidymodels)

## Load data
#df_ar <- read.csv("~/Downloads/CogDeclineRawData/df_gonogoten.csv")

do.one <- function(seed, df_ar, nhidden, nepoch, nnpenalty, ndeg, nknots){
  
  mAge <- mean(df_ar$age)
  sAge <- sd(df_ar$age) 
  df_ar <- df_ar %>% 
    mutate(zAge = (age - mAge)/sAge, 
           zAgeF = zAge * Female, 
           zAgeAdv=zAge*AdvancedDeg,
           zAgeCol=zAge*CollegeOrLess,
           FCollegeOrLess=Female*CollegeOrLess,
           FAdvancedDeg=Female*AdvancedDeg
    ) 
  
  
  set.seed(seed)
  df_split <- initial_split(df_ar, prop = 0.8)
  df_train <- training(df_split)
  df_test <- testing(df_split)
  
  metric <- metric_set(rmse) #metric to evaluate
  ctrl_res <- control_stack_resamples()
  
  ## Must be done AFTER you standardize
  folds <- vfold_cv(df_train, v = 5) 
  
  #Define models (pre-tuned elsewhere) 
  #general recipe. High school is base case
  gen_rec <- 
    recipe(scaled_score ~ zAge + Female + zAgeF + CollegeOrLess + 
             AdvancedDeg + zAgeAdv + zAgeCol + FCollegeOrLess + FAdvancedDeg, data = df_train) 
  
  
  #neural net
  #model definition 
  nn_spec <- mlp(
    hidden_units = nhidden,
    penalty = nnpenalty,
    epochs = nepoch) %>%  
    set_engine("nnet") %>% 
    set_mode("regression") 
  
  
  #recipe extension (no extension, interaction captured)
  nn_rec <-
    gen_rec 
  
  #workflow
  nn_wflow <- 
    workflow() %>%
    add_model(nn_spec) %>%
    add_recipe(nn_rec)
  
  # fit to the 5-fold cv
  nn_res <- 
    fit_resamples(
      nn_wflow,
      resamples = folds,
      metrics = metric,
      control = ctrl_res
    )
  
  nn.res <- fit(nn_wflow, data = df_train)
  
  #polynomial 
  #create model definition 
  lin_reg_spec <-
    linear_reg(penalty = 0.0) %>%
    set_engine("glmnet")
  
  
  #extend recipe
  lin_reg_rec <-
    gen_rec %>%
    step_poly(zAge, degree = ndeg)  %>%
    step_poly(zAgeF, degree = ndeg) %>%
    step_poly(zAgeAdv, degree = ndeg) %>%
    step_poly(zAgeCol, degree = ndeg)
  
  
  # add both to a workflow
  lin_reg_wflow <- 
    workflow() %>%
    add_model(lin_reg_spec) %>%
    add_recipe(lin_reg_rec)
  
  # fit to the 5-fold cv
  lin_reg_res <- 
    fit_resamples(
      lin_reg_wflow,
      resamples = folds,
      metrics = metric,
      control = ctrl_res
    )
  
  lin.res <- fit(lin_reg_wflow, data = df_train)
  
  ## Spline 
  spline_spec <-
    linear_reg() %>%
    set_engine(engine = 'lm') %>%
    set_mode('regression')
  
  
  spline_rec <- gen_rec %>%
    step_ns(zAge, deg_free = nknots) %>%
    step_ns(zAgeF, deg_free = nknots) %>%
    step_ns(zAgeCol, deg_free = nknots) %>%
    step_ns(zAgeAdv, deg_free = nknots) 
  
  spline_wflow <- workflow() %>%
    add_model(spline_spec) %>%
    add_recipe(spline_rec)
  
  spline_res <- fit_resamples(
    spline_wflow,
    resamples = folds,
    metrics = metric,
    control = ctrl_res
  )
  
  spline.res <- fit(spline_wflow, data = df_train)
  
  ## Stacked Model
  arith_st <- 
    stacks() %>%
    add_candidates(lin_reg_res) %>%
    add_candidates(nn_res) %>%
    add_candidates(spline_res) 
  
  st.res <-
    arith_st %>%
    blend_predictions() %>%
    fit_members()
  
  #keeping track of weights
  lin <- collect_parameters(st.res, "lin_reg_res")
  nn <- collect_parameters(st.res, "nn_res")
  spl <- collect_parameters(st.res, "spline_res")
  weights <- rbind(lin, nn, spl)
  
  preds <- data.frame(
    NN = predict(nn.res, new_data = df_test)$`.pred`, 
    PL = predict(lin.res, new_data = df_test)$`.pred`,
    SP = predict(spline.res, new_data = df_test)$`.pred`,
    ST = predict(st.res, new_data = df_test)$`.pred`,  
    ObsScore = df_test$scaled_score)
  
  preds$Avg <- apply(preds[,c( "NN","PL", "SP")], 1, mean)
  
  ## MSE on scaled test data 
  pred.rmse <- apply(preds, 2, FUN = function(y) {
    sqrt(sum((preds$ObsScore - y)^2/nrow(preds)))
  })
  
  ## Set up plotting data
  plot.dat <- data.frame(age = rep(25:80, 6), 
                         Female = rep(rep(c(0,1), each = 168)),
                         CollegeOrLess = c(rep(0, 56), rep(1,56), rep(0,56)),
                         AdvancedDeg = c(rep(0,56),rep(0,56), rep(1,56))) %>% 
    mutate(zAge = (age - mAge) / sAge, 
           zAgeF = zAge * Female,
           zAgeAdv=zAge*AdvancedDeg,
           zAgeCol=zAge*CollegeOrLess,
           FCollegeOrLess=Female*CollegeOrLess,
           FAdvancedDeg=Female*AdvancedDeg)
  
  ## Get standard predictions 
  plot.preds <- data.frame(
    NN = predict(nn.res, new_data = plot.dat)$`.pred`, 
    PL = predict(lin.res, new_data = plot.dat)$`.pred`, 
    SP = predict(lin.res, new_data = plot.dat)$`.pred`, 
    ST = predict(st.res, new_data = plot.dat)$`.pred`)
  
  plot.out <- cbind(plot.dat, plot.preds)
  
  return(list(RMSE = pred.rmse, 
              PlotData = plot.out, 
              Weights = weights))
}

do.one.boot <- function(seed, df, nhidden, nepoch, nnpenalty, ndeg, nknots) {
  set.seed(seed)
  this.test <- df[sample(1:nrow(df), replace = TRUE), ]
  res <- do.one(seed, this.test, 
                nhidden = nhidden, 
                nepoch = nepoch, 
                nnpenalty = nnpenalty,
                ndeg = ndeg, 
                nknots = nknots)
  return(res)
}

startseed = 1
bigB = 1000

bigres.dividedres <- mclapply(1:bigB, FUN = function(i) {
  do.one.boot(startseed + i - 1, df = df_ar,  
              nhidden = 10, nepoch = 20, nnpenalty = 0.1 , ndeg = 2, nknots = 2)
}) 

