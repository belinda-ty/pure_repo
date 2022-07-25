rm(list = ls())
library(dplyr)
library(tidyverse)
library(class)
library(lme4)
library(mgcv)

# clean data -------------------------------------------------------------
nlms <- read_csv("data/real_data/tobaco_clean.csv")

d <-
nlms %>%
  select(c('agesmk','race','age','sex','ms','wt','hhnum','follow')) %>%
  filter(complete.cases(.)) %>%
  filter(race %in% c(1, 2, 3, 5)) %>%
  mutate(race_orignal = race) %>%
  mutate(race = case_when(race == 1 ~"white",
                          race == 2 ~"black",
                          TRUE ~"nonwhite")) %>%
  mutate(race = factor(race),
         sex = factor(sex),
         ms = factor(ifelse(ms == 1, 1, 0)),
         agesmk = as.numeric(agesmk))

# sampling for 100 times
set.seed(100)
dats <- map(1:5, ~ sample_frac(d, size = 0.02))

# PURE test data --------------------------------------------------------------------
ms <- map(dats, 
          ~.x %>%
            filter(race == "nonwhite") %>%
            nrow()/3)

responses <- dats %>%
  map(~ .x %>%
        filter(race == "nonwhite") %>%
        select(agesmk) %>%
        log(.) %>%
        unlist() %>%
        as.vector())

# estimated k from Wang's cross validation method
hatk <-
responses %>%
  map2(.y = ms , ~ 
      estimateK(y = .x,
                m = round(.y, 0), 
                K = round(.y, 0) , 
                C = 1000,  
                method="ave")  )

# assigning group membership based on k
clustercut.nonwhite <-
responses %>%
  map2(.y = hatk,
      ~ hclust(dist(.x),"ave")%>%
        cutree(., .y)
        )

# nonwhite data with imputed group memberships
dats_nonwhite <- 
  map2(dats, 
       .y = clustercut.nonwhite,
        ~.x %>%
        filter(race == "nonwhite") %>%
        mutate(fakeid = as.factor(paste0("Group",.y)) ))

# obj is the lmer object
# VarR is the estimated variance R
# raneff.test is the random effect coefficients
# theta.test is the predicted mixed effects
# nonwhiteN is the test data
set.seed(100);
lmer_objects <-
  dats_nonwhite %>%
  map(~lmer(log(agesmk) ~ age + sex +  ms + wt + hhnum  +  follow  + (1|fakeid ), data = .x)) %>%
  map(~list(obj = .x,
            VarR = as.data.frame(VarCorr(.x))[2,4],
            raneff.test = unlist(ranef(.x)$fakeid),
            theta.test = predict(.x),
            N = .x@frame %>% nrow()*5
            )) %>%
  map(~list(obj = .x$obj,
            VarR = .x$VarR,
            raneff.test = .x$raneff.test,
            theta.test = .x$theta.test,
            N = .x$N,
            nonwhiteN = as_tibble(.x$obj@frame) %>%
              slice(rep(1:n(), each = 5))%>%
              select(-fakeid) %>%
              mutate(thetaN = rep(as.numeric(.x$theta.test), each = 5),
                     agesmk = ceiling(exp(
                # errorN + thetaN
                rnorm(.x$N, 0, sqrt(.x$VarR)) + thetaN
              )) )
  ))

nonwhiteN  <- lapply(lmer_objects, '[[', "nonwhiteN")
thetaN <- map(nonwhiteN, ~ 
                .x %>%
                select(thetaN)%>%
                unlist())

# PURE training data --------------------------------------------------------------------
# full data training 
fulld_re <- map(dats, 
                ~.x %>%
                  filter(race == "white"|race == "black") )

responses_full <- fulld_re %>%
  map(~ .x %>%
        select(agesmk) %>%
        log(.) %>%
        unlist() %>%
        as.vector())

# assigning group membership based on 20 groups
# full training data with IDs
fulld_re <-
  fulld_re %>%
  map2(.y = responses_full,
       ~ .x %>%
         mutate(id = as.factor(paste0("Group",
                                       hclust(dist(.y),"ave")%>%
                                         cutree(., 20))))
  )

# training on full and CMMP predict 
fullMod <-
  fulld_re  %>%
  map(~lmer(log(agesmk)~ age + sex +  ms + wt + hhnum  +  follow + (1|id ), data = .x))%>%
  map2(.y = nonwhiteN,
       ~list(obj = .x,
             VarR = as.data.frame(VarCorr(.x))[2,4],
             cmmp_randomeffects = unlist(ranef(.x)$id),
             mu.n.f = predict(.x, newdata = .y,re.form=NA),
             yN = .y$agesmk
  )) %>%
  map(~list(obj = .x$obj,
            VarR = .x$VarR,
            cmmp_randomeffects = .x$cmmp_randomeffects,
            mu.n.f = .x$mu.n.f,
            theta.cmmp.f = cmmp.theta(mu.n = .x$mu.n.f, raneff = .x$cmmp_randomeffects,
                                       y = .x$yN)
  ))

theta.cmmp.f <- lapply(fullMod, '[[', "theta.cmmp.f")
spe.cmmp.f <- map2(theta.cmmp.f , thetaN,
                   ~ (.x - .y) ^2)

# training on full and least square predict 
thetals.f <-
  fulld_re  %>%
  map(~lm(log(agesmk)~ age + sex +  ms + wt + hhnum  +  follow , data = .x))%>%
  map2(.y = nonwhiteN,
       ~ predict(.x, newdata = .y,re.form=NA))

spe.ls.f <- map2(thetals.f  , thetaN,
                ~ (.x - .y) ^2)

# training on full and standard cmmp predict 
fullMod_standard <-
  fulld_re %>%
  map(~lmer(log(agesmk)~ age + sex +  ms + wt + hhnum  +  follow + (1| race_orignal), data = .x))%>%
  map2(.y = nonwhiteN,
       ~list(obj = .x,
             standardcmmp_randomeffects= unlist(ranef(.x)$race_orignal),
             mu.n.full_standard = predict(.x, newdata = .y,re.form=NA),
             yN = .y$agesmk
       )) %>%
  map(~list(obj = .x$obj,
            standardcmmp_randomeffects= .x$standardcmmp_randomeffects,
            mu.n.full_standard = .x$mu.n.full_standard,
            theta.cmmp.full_standard = cmmp.theta(mu.n = .x$mu.n.full_standard , 
                                                  raneff = .x$standardcmmp_randomeffects,
                                      y = .x$yN)
  ))

theta.cmmp.full_standard <- lapply(fullMod_standard, '[[', "theta.cmmp.full_standard")
spe.cmmp.full_standard <- map2(theta.cmmp.full_standard , thetaN,
                   ~ (.x - .y) ^2)

# L_p matrix for predictions

funlldXZ <- map(fulld_re,
                ~ .x %>%
                  select(c("age", "sex","ms", "wt", "hhnum", "follow")))
  
full_2cluster <-
  fulld_re %>%
  map2(.y = funlldXZ,
       ~ .x %>%
         mutate(cluster =  hclust(dist(.y),"ave") %>%
                                        cutree(., 2)) %>%
         split(.$cluster)
  )

train_clus_1 <- lapply(full_2cluster , '[[', "1")
train_clus_2 <- lapply(full_2cluster , '[[', "2")

rmc1 <-
train_clus_1 %>%
  map(~  tryCatch(mgcv::gam(log(agesmk) ~ s(id,bs="re") + age  + sex +  ms + wt + hhnum  +  follow  ,
                      data=.x,method="REML")  , 
                  error = function(e) {
                    NA})
)

rmc2 <-
  train_clus_2 %>%
  map(~  tryCatch(mgcv::gam(log(agesmk) ~ s(id,bs="re") + age  + sex +  ms + wt + hhnum  +  follow  ,
                      data=.x,method="REML")  , 
                  error = function(e) {
                    NA})
  )

# we only create this id column, which is not used, since the following line can work only when we have the value of the column, but the value is not take into account
intersect_grp <- 
map2(train_clus_1, train_clus_2, ~
       intersect(.x$id, .y$id)[1])

nonwhiteN <-
nonwhiteN %>%
  map2(.y = intersect_grp , 
       ~ .x %>%
         mutate(id = rep(.y, each= dim(.x)[1])))
# delta method
# object new data

delta_method <- function(object, data ){
  bread  <- tryCatch(predict(object, type="lpmatrix"),
                     error = function(e){
                       NA
                     } )
  
  if(is.na(bread)) return(NA) 
  
    varmu <- bread %*% (vcov(object) %*% t(bread))
    max.var.obs <- max(diag(varmu))
    breadpred.ext <- predict(object, newdata = data, exclude="s(id)", type="lpmatrix") #L_p matrix for predictions
    varpred.ext <- breadpred.ext%*% (vcov(object) %*% t(breadpred.ext))
    varpred <- diag(varpred.ext)
    return(varpred)

}

varpred1 <- list()
varpred2 <- list()
for (i in seq_along(nonwhiteN)){
  varpred1[[i]] <- delta_method(rmc1[[i]], nonwhiteN[[i]])
  varpred2[[i]] <- delta_method(rmc2[[i]], nonwhiteN[[i]])
}

varpred1 <-
map2(rmc1, nonwhiteN, ~
       delta_method(.x, .y))

varpred2<-
  map2(rmc2, nonwhiteN, ~
         delta_method(.x, .y))
