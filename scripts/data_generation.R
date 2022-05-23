rm(list = ls())
library(purrr)
source('Scripts/PURE_2022/functions.R')

# training data parameter setup ---------------------------------------------------------
# Specify the true covariance structure
true.D = matrix(c(20), nrow=1, ncol=1) 
R = 1
q <- ncol(true.D)

# 30 groups 10 per group
n = 30; m =10 

# percetage of majority population is 60%
majorp <- 0.60

# parameters
# for normal people (majority people)
beta1 <- c(0.5,0.4,0.3,0.2,0.1,0)
p <- length(beta1)
Sigma1 <- diag(c(40,40,40,1,1))
mu1 <- c(30, 150,105, 0, 0) 

# for minority people'
# beta2 <- c(0.52,0.42,0.32,0.2,0.1,0) # for case 2
# beta2<-beta1 #for case 3
# beta2 <- beta1 # for case 6
beta2 <- c(0.55,0.45,0.35,0.2,0.1,0) # for case 1
Sigma2 <- diag(c(40,40,40,1,1))
mu2 <- c(90, 250,140, 0, 0)
mu3 <- c(100, 280,150, 0, 0)

# obtain group label and random effect 
id = gl(n = n,k = m, labels = paste0("Group",1:(n)))
set.seed(100); true.b <- get.trueb(true.D, n = (n)) 


# simulate training data --------------------------------------------------
# Simulate 100 sets of training data
set.seed(101)
out <-
  purrr::map(1:100, ~ sim_X(majorp, mu1, mu2, Sigma1, Sigma2, n = 30, m = 10, true.b, R = 1, beta1, beta2)) 

saveRDS(out, file="Data/simulated_data/Case1_train.RData")

# test data parameter setup ---------------------------------------------------------

# construct a extreme test population; 100 of them 
beta3 <- c(0.65,0.20,0.35,0.2,0.1,0) # more aggressive beta parameters
# beta3 <- c(0.60,0.50,0.4,0.2,0.1,0) # This is case 4 beta_3
# beta3 <- beta1  # This is case 6 beta_3
# beta3 <- beta2 # case 5
# md = 100; nn = 5
Sigma3 <- Sigma2
nn = 5

# simulate test data --------------------------------------------------
set.seed(101)
out_test <-
  purrr::map(1:100, ~ sim_test(pop1_prop = majorp, mu3, Sigma3,beta3, n = 30, md =100, nn = 5, true.b, R = 1))

saveRDS(out_test, file="Data/simulated_data/Case1_test.RData")


