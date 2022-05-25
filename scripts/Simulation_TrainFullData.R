rm(list = ls())
library(MASS)
library(lme4)
library(tidyverse)
library(randomForestSRC)
library(purrr)
library(caret)
library(xgboost)
source('scripts/functions.R')

# Import data -------------------------------------------------------------
# import simulated training and test data objects
train_obj <- readRDS("data/simulated_data/case2_train.RData")
test_obj <- readRDS("data/simulated_data/case2_test.RData")

# extract training data, training extended data, test data, test extended data and true mixed effect
# train
train <- lapply(train_obj, '[[', "train")
train_ex <- lapply(train_obj, '[[', "train_extend")
# test data
test.ext <- lapply(test_obj, '[[', "test")
true_theta <- lapply(test_obj, '[[', "new_theta")  # true test data theta
test_extend1 <- lapply(test_obj, '[[', "test_extend_left") # test data from extended model
test_extend2 <- lapply(test_obj, '[[', "test_extend_right")

# T test ------------------------------------------------------------------
tscore<-c()
for (i in 1:100){
  ttest <- t.test(train[[i]]$y, test.ext[[i]]$y)
  tscore[i] <- ttest$statistic
}
mean(tscore)
sd(tscore)


###########################################################################
# Training with full data -------------------------------------------------
###########################################################################
# CMMP, lmer, CMMP-extend, rp, rf, boosting, 
# CMMP regroup, lmer regroup, rp regroup

# CMMP ------------------------------------------------
# CMMP train with the full data, test for extreme data
# training data model
objs <- 
  train %>%
  purrr::map( ~ lmer(y ~ X.1 + X.2 + X.3 + X.4 + X.5 + X.6-1 + (Z-1 | id), data = .x, REML = TRUE)) 

# test data group mean
test_mean <-
  test.ext %>%
  purrr::map( ~ group_mean(data = .x, nn = 5))

# obtain CMMP estimator 
cmmp_estimator <-
  map2(test_mean, objs, 
     ~ cmmpmixed(fixed_ef = fixef(.y), random_ef =  ranef(.y)$id, group_mean = .x))

# CMMP mixed effect MSE
mse.cmmp <-
  map2(cmmp_estimator, true_theta,
     ~ mean((.x$cmmpmixed -.y)^2)) %>%
  unlist()
# CMMP ------------------------------------------------

# lmer ------------------------------------------------
# lmer prediction 
nn <- 5
lmer_predicted <-
  map2(objs, test.ext, 
       ~ predict(.x, newdata = .y) %>%
         as_tibble() %>%
         mutate(index = rep(1: (dim(test.ext[[1]])[1]/nn), each = nn))) %>%
  purrr::map(~.x %>%
        group_by(index) %>%
        summarise_all(mean) %>%
        select("value"))

mse.lmer <-
  map2(lmer_predicted, true_theta,
     ~ mean((.x$value -.y)^2)) %>%
  unlist()
# lmer ------------------------------------------------

# cmmp estimator extended----------------------------------
# Extended model
# training extended model
objs_ex <- 
  train_ex %>%
  purrr::map( ~ lmer(y ~ X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10+X.11-1 + (Z-1 | id), data = .x, REML = TRUE)) 

# test data group mean
test_mean_exl <-
  test_extend1 %>%
  purrr::map( ~ group_mean(data = .x, nn = 5))

test_mean_exr <-
  test_extend2 %>%
  purrr::map( ~ group_mean(data = .x, nn = 5))

# obtain CMMP estimator for two test data 
cmmp_estimator_exl <-
  map2(test_mean_exl, objs_ex, 
       ~ cmmpmixed(fixed_ef = fixef(.y), random_ef =  ranef(.y)$id, group_mean = .x))

cmmp_estimator_exr <-
  map2(test_mean_exr, objs_ex, 
       ~ cmmpmixed(fixed_ef = fixef(.y), random_ef =  ranef(.y)$id, group_mean = .x))

# CMMP mixed effect MSE
mse.cmmp.exl <-
  map2(cmmp_estimator_exl, true_theta,
       ~ mean((.x$cmmpmixed -.y)^2)) %>%
  unlist()

mse.cmmp.exr <-
  map2(cmmp_estimator_exr, true_theta,
       ~ mean((.x$cmmpmixed -.y)^2)) %>%
  unlist()

# lmer prediction 
# first extend data 
lmer_predicted_exl <-
  map2(objs_ex, test_mean_exl, 
       ~ predict(.x, newdata = .y) %>%
         as_tibble() %>%
         mutate(index = rep(1: (dim(test_mean_exl[[1]])[1]/nn), each = nn))) %>%
  purrr::map(~.x %>%
        group_by(index) %>%
        summarise_all(mean) %>%
        select("value"))

mse.lmer.exl <-
  map2(lmer_predicted_exl, true_theta,
       ~ mean((.x$value -.y)^2)) %>%
  unlist()

# second extend data
lmer_predicted_exr <-
  map2(objs_ex, test_mean_exr, 
       ~ predict(.x, newdata = .y) %>%
         as_tibble() %>%
         mutate(index = rep(1: (dim(test_mean_exr[[1]])[1]/nn), each = nn))) %>%
  purrr::map(~.x %>%
        group_by(index) %>%
        summarise_all(mean) %>%
        select("value"))

mse.lmer.exr <-
  map2(lmer_predicted_exr, true_theta,
       ~ mean((.x$value -.y)^2)) %>%
  unlist()
# cmmp estimator extended----------------------------------

# linear model -------------------------------------
# least square estimator 
lm_objs <- 
  train %>%
  purrr::map( ~ lm(y ~ X.1 + X.2 + X.3 + X.4 + X.5 + X.6 - 1 , data = .x)) 

lm_predicted <-
  map2(lm_objs, test.ext, 
       ~ predict(.x, newdata = .y) %>%
         as_tibble() %>%
         mutate(index = rep(1: (dim(test.ext[[1]])[1]/nn), each = nn))) %>%
  purrr::map(~.x %>%
        group_by(index) %>%
        summarise_all(mean) %>%
        select("value"))

mse.rp  <-
  map2(lm_predicted, true_theta,
       ~ mean((.x$value -.y)^2)) %>%
  unlist()


# linear model -------------------------------------

# random forest -------------------------------------
set.seed(100)
rf_objs <- 
  train %>%
  purrr::map( ~ rfsrc(y ~ X.1 + X.2 + X.3 + X.4 + X.5 + X.6 - 1,data = .x, importance = FALSE, proximity = TRUE))

rf_predicted <-
  map2(rf_objs, test.ext, 
       ~ randomForestSRC::predict.rfsrc(.x, .y)$predicted %>%
         as_tibble() %>%
         mutate(index = rep(1: (dim(test.ext[[1]])[1]/nn), each = nn))) %>%
  purrr::map(~.x %>%
        group_by(index) %>%
        summarise_all(mean) %>%
        select("value"))

mse.rf  <-
  map2(rf_predicted, true_theta,
       ~ mean((.x$value -.y)^2)) %>%
  unlist()

# random forest -------------------------------------

# boosting-------------------------------------
set.seed(100)
xgb_objs <- 
  train %>%
  purrr::map( ~ train(
    y ~ X.1 + X.2 + X.3 + X.4 + X.5 + X.6 -1, data = .x, method = "xgbTree",
    trControl = trainControl("cv", number = 10)))

xgb_predicted <-
  map2(xgb_objs, test.ext, 
       ~ predict(.x, .y) %>%
         as_tibble() %>%
         mutate(index = rep(1: (dim(test.ext[[1]])[1]/nn), each = nn))) %>%
  purrr::map(~.x %>%
        group_by(index) %>%
        summarise_all(mean) %>%
        select("value"))

mse.boosting   <-
  map2(xgb_predicted, true_theta,
       ~ mean((.x$value -.y)^2)) %>%
  unlist()

round(mean(mse.boosting), 2)
# 321.13
round(sd(mse.boosting), 2)
# 119.13

# boosting -------------------------------------


# CMMP regroup------------------------------------------------
train_regroup <- train

# regroup new id is based on the distribution of y 
kk <- 30
newids <- 
train_regroup %>%
  purrr::map(~hclust(dist(.x$y),"ave") %>%
        cutree(., kk))

# assign the new id to training data
train_regroup <- 
train_regroup %>%
  map2(.y = newids, ~.x %>%
        mutate(id =  as.factor(paste0("Group",.y) )))

# training data model
objs_frg <- 
  train_regroup %>%
  purrr::map( ~ lmer( y ~ X.1 + X.2 + X.3 + X.4 + X.5 + X.6 - 1 + (Z-1 | id), data = .x, REML = TRUE))

# test data group mean
test_mean <-
  test.ext %>%
  purrr::map( ~ group_mean(data = .x, nn = 5))

# obtain CMMP estimator 
cmmp_estimator <-
  map2(test_mean, objs_frg , 
       ~ cmmpmixed(fixed_ef = fixef(.y), random_ef =  ranef(.y)$id, group_mean = .x))

# CMMP mixed effect MSE
mse.cmmp.frg <-
  map2(cmmp_estimator, true_theta,
       ~ mean((.x$cmmpmixed -.y)^2)) %>%
  unlist()
# CMMP regroup------------------------------------------------

# lmer regroup------------------------------------------------
# lmer prediction --This is group level prediction fixed + random
lmpred.frg <-
  objs_frg %>%
  map2(.y = test_mean, ~predict(.x,  newdata =  .y, allow.new.levels = TRUE))

mse.lmer.frg <-
  map2(lmpred.frg, true_theta,
       ~ mean((.x -.y)^2)) %>%
  unlist()
# lmer regroup------------------------------------------------

# rp regroup------------------------------------------------
ls_objs_frg <- 
  train_regroup %>%
  purrr::map( ~ lm( y ~ X.1 + X.2 + X.3 + X.4 + X.5 + X.6 - 1 , data = .x))

ls.frg <-
  ls_objs_frg %>%
  map2(.y = test_mean, ~ predict(.x,  newdata =  .y))

mse.rp.frg <-
  map2(ls.frg, true_theta,
       ~ mean((.x -.y)^2)) %>%
  unlist()
# rp regroup------------------------------------------------

# Table 2 Estimated G and R for toy example row 1------------------------------------------------
# obtain estimator of the random intercept variance VarD, and estimator of the residual variance VarR
VarDR <- 
  objs %>%
  purrr::map( ~ tibble(VarD =  as.data.frame(VarCorr(.x))[1,4],
                VarR =  as.data.frame(VarCorr(.x))[2,4])) %>%
  bind_rows()

# obtain estimator of the random intercept variance VarD, and estimator of the residual variance VarR
VarDR.frg <- 
  objs_frg %>%
  purrr::map( ~ tibble(VarD =  as.data.frame(VarCorr(.x))[1,4],
                VarR =  as.data.frame(VarCorr(.x))[2,4])) %>%
  bind_rows()


T2r1c1 <-
  VarDR %>%
  summarise(VarD_bar = round(mean(VarD),2),
            VarD_sd = round(sd(VarD)/10,2),
            VarR_bar = round(mean(VarR),2),
            VarR_sd = round(sd(VarR)/10,2),
            Ratio_GR = round(mean(VarD/VarR), 2),
            Ratio_GR_sd = round(sd(VarD/VarR)/10,2)) %>%
  mutate("Estimated G" = paste0(VarD_bar,"(", VarD_sd, ")"),
         "Estimated R" = paste0(VarR_bar,"(", VarR_sd, ")"),
         "Ration G/R" =  paste0(Ratio_GR,"(", Ratio_GR_sd, ")")) %>%
  select( `Estimated G`, `Estimated R` ,`Ration G/R`)%>%
  pivot_longer(cols = everything(), names_to = "Train with full data", values_to = "CMMP & Lmer*")


T2r1c2 <-
  VarDR.frg %>%
  summarise(VarD_bar = round(mean(VarD),2),
            VarD_sd = round(sd(VarD)/10,2),
            VarR_bar = round(mean(VarR),2),
            VarR_sd = round(sd(VarR)/10,2),
            Ratio_GR = round(mean(VarD/VarR), 2),
            Ratio_GR_sd = round(sd(VarD/VarR)/10,2)) %>%
  mutate("Estimated G" = paste0(VarD_bar,"(", VarD_sd, ")"),
         "Estimated R" = paste0(VarR_bar,"(", VarR_sd, ")"),
         "Ration G/R" =  paste0(Ratio_GR,"(", Ratio_GR_sd, ")")) %>%
  select( `Estimated G`, `Estimated R` ,`Ration G/R`) %>%
  pivot_longer(cols = everything(), names_to = "Train with full data", values_to = "Regrouped CMMP & Regrouped Lmer*")

d_cov <- 
T2r1c1 %>%
  left_join(T2r1c2, c("Train with full data" = "Train with full data")) 

d_cov %>%
  knitr::kable(caption = "Estimated G and R and (SEs in parentheses) for toy example over 100 runs")

writexl::write_xlsx(d_cov, "data/simulated_data/case2_cov_full.xlsx")

# Table 2 Estimated G and R for toy example row 1------------------------------------------------

# Training with full data end -------------------------------------------------


# second row of Table - rp, cmmp and lmer- minority data
average <- (mse.cmmp.exl+ mse.cmmp.exr)/2
dats <- 
  list(mse.rp, mse.rf, mse.boosting, mse.cmmp, mse.lmer, mse.lmer.frg,mse.cmmp.frg,
       mse.cmmp.exl,
       mse.cmmp.exr,
       average) %>%
  setNames(c("RP", "Random Forest", "Boosted Reg Tree", "CMMP", "Lmer*", "Regrouped Lmer*", "Regrouped CMMP",
            "ULMM CMMP left", "ULMM CMMP right", "Average"))

d_mse <- 
dats %>%
  purrr::imap(~ .x %>%
                as_tibble() %>%
                summarise(mean = round(mean(value), 2) , sd = mean(sd(value)/10, 2)) %>%
                mutate(method = .y)) %>%
  bind_rows() %>%
  mutate("Train with full data (n = 200)" = mean,
         "SE" = paste0("(", round(sd, 2), ")")) %>%
  mutate(RP_baseline = mean(mse.rp))%>%
  mutate("% Decrease relative to RP" = round((RP_baseline - mean)/RP_baseline*100, 2))%>%
  select(method, "Train with full data (n = 200)",SE, "% Decrease relative to RP") %>%
  column_to_rownames(var = "method") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column()

d_mse %>%
  knitr::kable(caption = "Summary of MSPE simulation results for prediction of mixed effect for case2")

writexl::write_xlsx(d_mse, "data/simulated_data/case2_mse_full.xlsx")

