rm(list = ls())
library(mgcv)
library(mclust)
library(purrr)
library(MASS)
library(lme4)
library(tidyverse)
library(randomForestSRC)
library(caret)
library(xgboost)
source('scripts/functions.R')
# Import data -------------------------------------------------------------
# import simulated training and test data objects
train_obj <- readRDS("Data/simulated_data/case1_train.RData")
test_obj <- readRDS("Data/simulated_data/case1_test.RData")

# extract training data, training extended data, test data, test extended data and true mixed effect
# train
train <- lapply(train_obj, '[[', "train")
train_ex <- lapply(train_obj, '[[', "train_extend")
# test data
test.ext <- lapply(test_obj, '[[', "test")
true_theta <- lapply(test_obj, '[[', "new_theta")  # true test data theta
test_extend1 <- lapply(test_obj, '[[', "test_extend_left") # test data from extended model
test_extend2 <- lapply(test_obj, '[[', "test_extend_right")

###################################################################
############# Random effects extrapolation algorithm  #############
################## Step 1 #########################################
################## GIVH ###########################################

# gIVH 
# training data # eliminate x.1 z.1, usually keep z.2, but we do not have intercepts 
##
max.var.obs <- c()
varprednum <- list()
out <- list()
set.seed(100);
for (i in 1:100){
  rm1 <- gam(y ~ s(id, bs="re") + s(X.2) + s(X.3) + s(X.4) + s(X.5) + s(X.6), data = train[[i]],method="REML") 
  
  # delta method
  bread <-  predict(rm1, type="lpmatrix")
  varmu <- bread %*% (vcov(rm1) %*% t(bread))
  max.var.obs[i] <- max(diag(varmu))
  
  ### Variance for predictions
  # L_p matrix for predictions
  breadpred <- predict(rm1, newdata=test.ext[[i]], exclude = "s(id)", type="lpmatrix")
  # delta method
  varpred <- breadpred %*% (vcov(rm1) %*% t(breadpred))
  # indicator for whether prediction outside gIVH 
  varprednum[[i]]<-diag(varpred) 
  out[[i]]<- varprednum[[i]]> max.var.obs[i]
}

lapply(out,table) 
sum(unlist(map(lapply(out,table), length))==1)
# the majority are outside of the GIVH

################## Step 2 #########################################
# Model Based Clustering for majority and minority
set.seed(101);
for (i in 1:100){
  trainXZ <- train[[i]][ ,c("X.2", "X.3","X.4", "X.5", "X.6")]
  mcfit <- Mclust(trainXZ)
  train[[i]]$cluster <- mcfit$classification
}

sapply(train, function(x) length(unique(x$cluster))) # every training data were classified into 2 groups

################## Step 3 determine which group is the minority data #########################################
## for each clusterd training data # check train_clus_1 and train_clus_2 
max.var.obs1 <- c()
max.var.obs2 <- c()
varpred1 <- list()
varpred2 <- list()
out1 <- list()
out2 <- list()
train2<-train

set.seed(101);
for (i in 1:100){
  train_clus_1 <- train2[[i]][train2[[i]]$cluster == 1,]
  train_clus_2 <- train2[[i]][train2[[i]]$cluster == 2,]
  rmc1 <- gam(y ~ s(id,bs="re") + s(X.2) + s(X.3)+s(X.4)+s(X.5)+s(X.6),data=train_clus_1,method="REML") 
  rmc2 <- gam(y ~ s(id,bs="re") + s(X.2)+s(X.3)+s(X.4)+s(X.5)+s(X.6),data=train_clus_2,method="REML") 
  
  # the following deals with the case predict.gam cannot handle new levels of groups in the test data
  # note we are not using the random effects for the prediction, we skip the random part by exclude="s(id)"
  # however, we still need a group variable to feed in the gam() in order for the program to proceed
  # any group works as long as it is in the training group, the results won't change since we are not using it.
  if (all(test.ext[[i]]$id %in%unique(train_clus_1$id)==TRUE)==FALSE){
    newgroupid1<-which(test.ext[[i]]$id %in%unique(train_clus_1$id)==FALSE)
    test.ext[[i]]$id[newgroupid1] <- sample(test.ext[[i]]$id[-newgroupid1],1)
  }
  if (all(test.ext[[i]]$id %in%unique(train_clus_2$id)==TRUE)==FALSE){
    newgroupid<-which(test.ext[[i]]$id %in%unique(train_clus_2$id)==FALSE)
    test.ext[[i]]$id[newgroupid] <- sample(test.ext[[i]]$id[-newgroupid],1)
  }
  
  # delta method for train_clus_1
  bread1  <- predict(rmc1, type="lpmatrix")
  varmu1 <- bread1 %*% (vcov(rmc1) %*% t(bread1))
  max.var.obs1[i] <- max(diag(varmu1))
  breadpred.ext1 <- predict(rmc1, newdata=test.ext[[i]],exclude="s(id)", type="lpmatrix") #L_p matrix for predictions
  varpred.ext1 <- breadpred.ext1%*% (vcov(rmc1) %*% t(breadpred.ext1))
  varpred1[[i]] <- diag(varpred.ext1)
  
  # delta method for train_clus_2
  bread2<- predict(rmc2, type="lpmatrix")
  varmu2 <- bread2 %*% (vcov(rmc2) %*% t(bread2))
  max.var.obs2[i] <- max(diag(varmu2))
  breadpred.ext2  <- predict(rmc2, newdata = test.ext[[i]],exclude="s(id)", type="lpmatrix")# L_p matrix for predictions
  varpred.ext2 <- breadpred.ext2 %*% (vcov(rmc2) %*% t(breadpred.ext2))
  varpred2[[i]] <- diag(varpred.ext2)
  
  out1[[i]]<- varpred1[[i]]> max.var.obs1[i]
  out2[[i]]<- varpred2[[i]]> max.var.obs2[i]
}

lapply(out1,table) # majority out of hull again
lapply(out2,table) # majority out of hull again
sum(unlist(map(lapply(out1,table), length))==1)
sum(unlist(map(lapply(out2,table), length))==1)
iddd<-unlist(map(lapply(out2,table), length))==1&unlist(map(lapply(out1,table), length))==1
# interquantile range
# IQR1 <- sapply(varpred1,function(x) IQR(x,na.rm = TRUE))
# IQR2 <- sapply(varpred2,function(x) IQR(x,na.rm = TRUE))
# first quantile
# first.Qu1 <- sapply(varpred1, function(x) as.numeric(quantile(x, c(0.25),na.rm=TRUE)))
# first.Qu2 <- sapply(varpred2, function(x) as.numeric(quantile(x, c(0.25),na.rm=TRUE)))

# median
med1 <- sapply(varpred1, median)
med2 <- sapply(varpred2, median)

dis1 <- (med1 - max.var.obs1)#/IQR1 
dis2 <- (med2 - max.var.obs2)#/IQR2 

# rule1 <- sign(dis1 -dis2)
# rule2 <- sign(dis1/max.var.obs1- dis2/max.var.obs2) 
# sum(rule1 == rule2)

# data_diff<- data.frame(cbind(cluster1 = abs(dis1), cluster2 = abs(dis2))) # for observe only
# this is the rule
diffe = abs(dis1 - dis2); b = apply(cbind(dis1,dis2), 1, min)
ab <- data.frame(b, diffe)
ab %>% mutate(c= b*0.05, d = diffe > c) ->ab2
table(ab2$d)

# cpick is the pick of the cluster
cpick <- ifelse(abs(dis1) < abs(dis2),1,2)

min(abs(abs(dis1) - abs(dis2)))  

# decide the minority data
train_minority <- list()
train_majority <- list()
for (i in 1:100){
  train_minority[[i]] <- train[[i]][train[[i]]$cluster == cpick[i],]
  train_majority[[i]] <- train[[i]][train[[i]]$cluster != cpick[i],]
}

sapply(train_minority, function(x){dim(x)[1]})
sapply(train_majority, function(x){dim(x)[1]})
################################################################################################################
################## poor performance of training with minority data(before regrouping) ##########################
objs <- 
  train_minority %>%
  purrr::map( ~ lmer(y ~ X.1 + X.2 + X.3 + X.4 + X.5 + X.6-1 + (Z-1 | id), data = .x, REML = TRUE)) 

# test data group mean
test_mean <-
  test.ext %>%
  purrr::map( ~ group_mean(data = .x, nn = 5))

# obtain CMMP estimator 
cmmp_estimator <-
  map2(test_mean, objs, 
       ~ cmmpmixed(fixed_ef = fixef(.y), random_ef =  ranef(.y)$id, group_mean = .x))

# CMMP MSE-------------------------------------------
mse.cmmpmin <-
  map2(cmmp_estimator, true_theta,
       ~ mean((.x$cmmpmixed -.y)^2)) %>%
  unlist()

# lmer MSE-------------------------------------------
nn<- 5
lmer_predicted <-
  map2(objs, test.ext, 
       ~ predict(.x, newdata = .y) %>%
         as_tibble() %>%
         mutate(index = rep(1: (dim(test.ext[[1]])[1]/nn), each = nn))) %>%
  purrr::map(~.x %>%
               group_by(index) %>%
               summarise_all(mean) %>%
               select("value"))

mse.lmermin <-
  map2(lmer_predicted, true_theta,
       ~ mean((.x$value -.y)^2)) %>%
  unlist()

VarDR <- 
  objs %>%
  purrr::map( ~ tibble(VarD =  as.data.frame(VarCorr(.x))[1,4],
                       VarR =  as.data.frame(VarCorr(.x))[2,4])) %>%
  bind_rows()

T2r2c1 <-
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
  pivot_longer(cols = everything(), names_to = "Train with minority data", values_to = "CMMP & Lmer*")

T2r2c1 
# linear model -------------------------------------
lm_objs <- 
  train_minority %>%
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

mse.rpmin <-
  map2(lm_predicted, true_theta,
       ~ mean((.x$value -.y)^2)) %>%
  unlist()
# linear model -------------------------------------

# random forest -------------------------------------
set.seed(100)
rf_objs <-
  train_minority %>%
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

mse.rfmin  <-
  map2(rf_predicted, true_theta,
       ~ mean((.x$value -.y)^2)) %>%
  unlist()

# random forest -------------------------------------

# boosting-------------------------------------
set.seed(100)
xgb_objs <-
  train_minority %>%
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

mse.boostmin   <-
  map2(xgb_predicted, true_theta,
       ~ mean((.x$value -.y)^2)) %>%
  unlist()


# boosting -------------------------------------

############################ Step 4   Regrouping of the minority training ######################################
################################################################################################################
# We don't need Step 2.1 since we are forcing m to be the same
# Step 2.1 : Estimate k by CV_a
# mm <- round(length(train_minority$y)/3,0) 
# kk <- estimateK(y = train_minority$y, m = mm, K = mm, C = 1000,  method="ave")  
train_minority_regroup <- train_minority

# regroup new id is based on the distribution of y 
kk <- 30 # n # force m not change
# newids <- 
#   train_minority_regroup %>%
#   purrr::map( ~ hclust(dist(.x$y),"ave")  %>%
#                 cutree(., kk))


newids <- list()
for (i in 1:100){
  trainXZ <- train_minority_regroup[[i]][ ,c("X.2", "X.3","X.4", "X.5", "X.6")]
  mcfit <- Mclust(trainXZ, G = 30)
  newids[[i]]<- mcfit$classification
}



#####


# assign the new id to training data
train_minority_regroup  <- 
  train_minority_regroup  %>%
  map2(.y = newids, ~.x %>%
         mutate(id =  as.factor(paste0("Group",.y) )))

# training data model
objs_frg <- 
  train_minority_regroup %>%
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
mse.cmmp.rg <-
  map2(cmmp_estimator, true_theta,
       ~ mean((.x$cmmpmixed -.y)^2)) %>%
  unlist()

# lmer regroup------------------------------------------------
# lmer prediction --This is group level prediction fixed + random
lmpred.frg <-
  objs_frg %>%
  map2(.y = test_mean, ~predict(.x,  newdata =  .y, allow.new.levels = TRUE))

mse.lmer.rg <-
  map2(lmpred.frg, true_theta,
       ~ mean((.x -.y)^2)) %>%
  unlist()

# obtain estimator of the random intercept variance VarD, and estimator of the residual variance VarR
VarDR.frg <- 
  objs_frg %>%
  purrr::map( ~ tibble(VarD =  as.data.frame(VarCorr(.x))[1,4],
                       VarR =  as.data.frame(VarCorr(.x))[2,4])) %>%
  bind_rows()


T2r2c2 <-
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
  pivot_longer(cols = everything(), names_to = "Train with minority data", values_to = "Regrouped CMMP & Regrouped Lmer*")
T2r2c2

# Table 2 Estimated G and R for toy example row 2------------------------------------------------
d_cov <- 
  T2r2c1 %>%
  left_join(T2r2c2, c("Train with minority data" = "Train with minority data")) 

d_cov %>%
  knitr::kable(caption = "Estimated G and R and (SEs in parentheses) for toy example over 100 runs")
writexl::write_xlsx(d_cov, "data/simulated_data/m_based/case1_cov_min_mbased.xlsx")

# Table 2 Estimated G and R for toy example row 2------------------------------------------------


###################### Table Results Entries ########################################
# second row of Table - rp, cmmp and lmer- minority data
dats <- 
  list(mse.rpmin, mse.rfmin, mse.boostmin, mse.cmmpmin, mse.lmermin, mse.lmer.rg,mse.cmmp.rg) %>%
  setNames(c("RP", "Random Forest", "Boosted Reg Tree", "CMMP", "Lmer*", "Regrouped Lmer*","PURE"))

d_mse <-
  dats %>%
  purrr::imap(~ .x %>%
                as_tibble() %>%
                summarise(mean = round(mean(value), 2) , sd = mean(sd(value)/10, 2)) %>%
                mutate(method = .y)) %>%
  bind_rows() %>%
  mutate("Train with minority data (n = 120)" = mean,
         "SE" = paste0("(", round(sd, 2), ")")) %>%
  mutate(RP_baseline = mean(mse.rpmin))%>%
  mutate("% Decrease relative to RP" = round((RP_baseline - mean)/RP_baseline*100, 2))%>%
  select(method, "Train with minority data (n = 120)",SE, "% Decrease relative to RP") %>%
  column_to_rownames(var = "method")%>%
  t()%>%
  as.data.frame() %>%
  rownames_to_column()

d_mse %>%
  knitr::kable(caption = "Summary of MSPE simulation results for prediction of mixed effect for case1")

writexl::write_xlsx(d_mse, "data/simulated_data/m_based/case1_mse_min_mbased.xlsx")

