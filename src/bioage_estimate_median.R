
################################################################
# This script contains all prediction functions: GrimAge, AFT-mrl Weibull, AFT-mrl Gompertz, PH-mrl and AFT-mrl semiparametric
#
#
#
################################################################

library(eha)
library(survival)
source("src/gompertz_draw.R")
source("src/weibull_draw.R")

################################################################
# Method A: AFT-mrl Gompertz

GetB_AFTmrl_Gompertz <- function(df_train, df_test, params, lt){
  
  M <- params$M
  
  # fit the AFT model
  fit_aft <- eha::aftreg(formula = Surv(c, age_end, status) ~ V1 + V2, 
                         data = df_train, dist = "gompertz") 
  
  # get estimated parameters
  sigma_est <- exp(fit_aft$coefficients["log(scale)"])
  tau_est <- exp(fit_aft$coefficients["log(shape)"])
  a_est <- tau_est / sigma_est
  b_est <- 1 / sigma_est
  
  # obtain linear predictors
  betas_est_aft <- fit_aft$coefficients[1:M]
  linpred_aft <- rowSums(sweep(df_test[,1:M], 2, betas_est_aft, "*"))
  
  # obtain predicted MRL
  pred_mrl<- vector(length = nrow(df_test))
  
  for (i in 1:nrow(df_test)){
    surv_prob <- gomp_baseline_surv(df_test$c[i] * exp(linpred_aft[i]), a = a_est, b = b_est)
    med_prob <- surv_prob / 2
    t_unadj <- inverse_gomp_baseline_surv(med_prob, a=a_est, b=b_est)
    t_adj <- t_unadj / exp(linpred_aft[i])
    pred_mrl[i] <-  t_adj - df_test$c[i]
  }
  
  # obtain biological age, assuming true lifetable is known
  aft_b <- vector(length = length(pred_mrl))
  
  for (i in 1:length(pred_mrl)){
    aft_b[i] <- lt$t[which.min(abs(lt$medrl - pred_mrl[i]))]
  }
  # obtain RMSE
  rmse <- sqrt(mean((aft_b - df_test$b)^2))
  
  return(list(b_pred = aft_b, mrl_pred = pred_mrl, rmse = rmse, linpred_pred = linpred_aft, betas_est = betas_est_aft,
              sigma_est = sigma_est, tau_est = tau_est, a_est = a_est, b_est = b_est))
  
}


################################################################
# Method B: GrimAge

GetB_GrimAge <- function(df_train, df_test){
  
  # fit Cox PH model
  fit_PH <- coxph(Surv(follow_up_time, status) ~ c + V1 + V2, 
                  data = df_train)
  
  # Get linear predictors 
  linpred_train <- predict(fit_PH, newdata = df_train, type = "lp", reference = "zero")
  linpred_test <- predict(fit_PH, newdata = df_test, type = "lp", reference = "zero")
  
  # scale linear predictors to same mean and sd as chronological age 
  slope <- mean(df_train$c) - (sd(df_train$c)/ sd(linpred_train))* mean(linpred_train)
  intc <- (sd(df_train$c) / sd(linpred_train))
  
  ga_b <- slope + intc * linpred_test
  
  # obtain RMSE
  rmse <- sqrt(mean((ga_b - df_test$b)^2))
  
  return(list(b_pred = ga_b, rmse = rmse, linpred_pred = linpred_test, betas_est = fit_PH$coefficients,
              slope = slope, intc = intc))
  
}

################################################################
# Method C: Cox PH-mrl 

GetB_CoxPHmrl <- function(df_train, df_test, params, lt){
  
  M <- params$M
  
  # fit Cox PH model
  fit_ph <- coxph(Surv(c, age_end, status) ~ V1 + V2, data = df_train)
  
  # obtain linear predictors
  betas_est_ph <- fit_ph$coefficients[1:M]
  linpred_ph <- rowSums(sweep(df_test[,1:M], 2, betas_est_ph, "*"))
  
  # obtain estimate of baseline survival
  base_surv <- survfit(fit_ph, newdata = data.frame(V1 = 0, V2 = 0))
  approx_base_surv <- approxfun(base_surv$time, base_surv$surv, yleft = 1, yright = 0, n = 100)
  inverse_approx_base_surv <- approxfun(base_surv$surv, base_surv$time, yleft = max(base_surv$time), yright = min(base_surv$time), n = 100)
  
  # obtain predicted MRL
  pred_mrl <- vector(length = nrow(df_test))
  
  for (i in 1:nrow(df_test)){                                          
    surv_prob <- approx_base_surv(df_test$c[i])^(exp(linpred_ph[i]))
    med_prob_adj <- (surv_prob / 2)^(1/exp(linpred_ph[i]))
    t_adj <- inverse_approx_base_surv(med_prob_adj)
    pred_mrl[i] <-  t_adj - df_test$c[i]
  }
  
  # obtain biological age, assuming true lifetable is known
  ph_b <- vector(length = length(pred_mrl))
  
  for (i in 1:length(pred_mrl)){
    ph_b[i] <- lt$t[which.min(abs(lt$medrl - pred_mrl[i]))]
  }
  
  # obtain RMSE
  rmse <- sqrt(mean((ph_b - df_test$b)^2))
  
  return(list(b_pred = ph_b, mrl_pred = pred_mrl, rmse = rmse, linpred_pred = linpred_ph, betas_est = betas_est_ph))
  
}

################################################################
# Method D: AFT-mrl Weibull

GetB_AFTmrl_Weibull <- function(df_train, df_test, params, lt){
  
  M <- params$M
  
  # fit the AFT model
  fit_aft <- eha::aftreg(formula = Surv(c, age_end, status) ~ V1 + V2, 
                         data = df_train, dist = "weibull") 
  
  scale <- exp(fit_aft$coefficients["log(scale)"])
  shape <- exp(fit_aft$coefficients["log(shape)"])
  lambda_est <- scale^(-shape)
  nu_est <- shape
  
  # obtain linear predictors 
  betas_est <- fit_aft$coefficients[1:M]
  linpred_aft <- rowSums(sweep(df_test[,1:M], 2, betas_est, "*"))
  
  # obtain predicted MRL
  pred_mrl<- vector(length = nrow(df_test))
  
  for (i in 1:nrow(df_test)){
    surv_prob <- weib_baseline_surv(df_test$c[i] * exp(linpred_aft[i]), lambda = lambda_est, nu = nu_est)
    med_prob <- surv_prob / 2
    t_unadj <- inverse_weib_baseline_surv(med_prob, lambda = lambda_est, nu = nu_est)
    t_adj <- t_unadj / exp(linpred_aft[i])
    pred_mrl[i] <-  t_adj - df_test$c[i]
  }
  
  # obtain biological age, assuming true lifetable is known
  aft_b <- vector(length = length(pred_mrl))
  
  for (i in 1:length(pred_mrl)){
    aft_b[i] <- lt$t[which.min(abs(lt$medrl - pred_mrl[i]))]
  }
  
  # obtain RMSE
  rmse <- sqrt(mean((aft_b - df_test$b)^2))
  
  return(list(b_pred = aft_b, mrl_pred = pred_mrl, rmse = rmse, linpred_pred = linpred_aft, betas_est = betas_est,
              scale_est = scale, shape_est = shape, lambda_est = lambda_est, nu_est = nu_est))
  
}


################################################################
# Method E: Semiparametric AFT-mrl

GetB_AFTmrl_semipar <- function(df_train, df_test, params, lt){
  
  # df_train <- list_df_iter[[6]]
  M <- params$M
  
  # get the KM plot for the weights
  km <- survfit(Surv(c, age_end, status) ~ 1, data = df_train, timefix = F) # taking left truncation into account
  
  # get the weights
  adj_weight <- diff(c(1, km$surv)) * -1  # add 1 because survival starts at 1, * -1 to make weights positive
  
  # sort by T
  df_train <- df_train[order(df_train$age_end),]
  df_train$adj_weight <- adj_weight
  
  # fit weighted regression model
  lm.1 <- lm(log(age_end) ~ V1 + V2, data = df_train, weights = adj_weight)
  
  # obtain baseline survival 
  eps <- log(df_train$age_end) - fitted(lm.1)
  base_surv <- as.data.frame(cbind(time = exp(eps)[order(eps)] * exp(lm.1$coefficients[1]), 
                                   surv =  1 - cumsum(adj_weight[order(eps)]) /sum(adj_weight)))
  approx_base_surv <- approxfun(base_surv$time, base_surv$surv, yleft = 1, yright = 0, n = 100)
  inverse_approx_base_surv <- approxfun(base_surv$surv, base_surv$time, yleft = max(base_surv$time), yright = min(base_surv$time), n = 100)
  
  # obtain linear predictors for test data
  linpred_lm <- predict(lm.1, newdata = df_test)
  linpred_aft <- (linpred_lm - lm.1$coefficients[1]) * -1 # subtract intercept and multiply by -1 to get to same linear predictors as aftreg (take exponent for acceleration factor)
  
  # obtain predicted MRL
  pred_mrl <- vector(length = nrow(df_test))
  
  for (i in 1:nrow(df_test)){
    surv_prob <- approx_base_surv(df_test$c[i] * exp(linpred_aft[i]))
    med_prob <- surv_prob / 2
    t_unadj <- inverse_approx_base_surv(med_prob)
    t_adj <- t_unadj / exp(linpred_aft[i])
    pred_mrl[i] <-  t_adj - df_test$c[i]
  }
  
  # obtain biological age, assuming true lifetable is known
  aft_b <- vector(length = length(pred_mrl))
  
  for (i in 1:length(pred_mrl)){
    aft_b[i] <- lt$t[which.min(abs(lt$medrl - pred_mrl[i]))]
  }
  
  # obtain RMSE
  rmse <- sqrt(mean((aft_b - df_test$b)^2))
  
  return(list(b_pred = aft_b, mrl_pred = pred_mrl, rmse = rmse, linpred_pred = linpred_aft, betas_est = lm.1$coefficients))
  
}
