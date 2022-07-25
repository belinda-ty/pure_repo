library(mvtnorm)

####################################################
########### some necessary functions ###############
# function to get random effects
get.trueb <- function(true.D, n){
  get.nonzero.D <- which(diag(true.D) > 0)
  true.b <- matrix(0, nrow = n, ncol = dim(true.D)[1])
  true.b[, get.nonzero.D] <- rmvnorm(n, rep(0, length(get.nonzero.D)),as.matrix(true.D[get.nonzero.D, get.nonzero.D]) )
  true.b
}

# function to simulate the training data 
sim_X <- function(pop1_prop, pop1_mu, pop2_mu, pop1_sigma, pop2_sigma, n, m, true.b, R, beta1, beta2){
  # simulate fixed effect X for majority and minority
  X1 <- cbind(1, as.matrix(mvrnorm(n = n*m*pop1_prop, pop1_mu ,pop1_sigma)), 1)
  X2 <- cbind(1, as.matrix(mvrnorm(n = n*m*(1-pop1_prop), pop2_mu, pop2_sigma)),2)
  fixed_X1 <- X1[,1:(ncol(X1)-1)]%*% beta1
  fixed_X2 <- X2[,1:(ncol(X2)-1)]%*% beta2
  fixed_temp <- rbind(fixed_X1, fixed_X2)
  order <- sample(nrow(fixed_temp))
  fixed <- fixed_temp[order, ]
  X <- rbind(X1, X2)[order, ]
  # obtain Z
  Z <- rep(1, n*m)
  id <- gl(n=n,k=m,labels=paste0("Group",1:(n)))
  theta <- fixed + (Z * true.b[id, ])
  y <- rnorm(n*m, mean = theta, sd = sqrt(R))
  # obtain training data
  train <- data.frame(id= id, X = X, Z = Z, y =  y)
  # add extend x
  Xextend = t(apply(X, 1, function(x){if(x[7]==1){cbind(x[1:6],rep(0,6))}else{cbind(rep(0,6),x[1:6])}}))
  Xextend[,1]<-1
  train_extend <- data.frame(id= id, X = Xextend[, -7], Z = Z, y = y)
  return(list(train = train, train_extend = train_extend))
}

# simulate test data
sim_test <- function(pop1_prop, pop3_mu, pop3_sigma, beta3, n, md, nn, true.b, R){
  # test data has the same random effects as in training data
  true.I3 <- sample(1:n, md*(1-pop1_prop), replace = TRUE) 
  new.b3 <- true.b[true.I3, ]
  teid3 <- as.factor(rep(paste0("Group",true.I3), each = nn))
  # test data fixed effect has its own parameters
  Xt3  <- cbind(1, mvrnorm(n = md*(1 - pop1_prop), pop3_mu, pop3_sigma))
  Zt3 <- as.matrix(rep(1, md*(1-majorp)))
  new_X3 <- Xt3[rep(1:nrow(Xt3), each = nn), ]
  new_Z3 <- Zt3[rep(1:nrow(Zt3), each = nn), ]
  new_theta3 <- Xt3 %*% beta3  +  (Zt3* new.b3)
  new_error3 <- rnorm(nn*md*(1-pop1_prop), mean = 0, sd = sqrt(R))
  new_y3 <- rep(new_theta3, each = nn) + new_error3
  test_data <- data.frame(id = teid3, X = new_X3, Z = new_Z3 , y = new_y3)
  
  # get the extended data
  Xexleft = cbind(new_X3, matrix(0, nrow =dim(new_X3)[1] , ncol = 5))
  Xexright = cbind(1,matrix(0, nrow =dim(new_X3)[1] , ncol = 5), new_X3[,2:ncol(new_X3)])
  test_extend_left <- data.frame(id = teid3, X = Xexleft, Z = new_Z3 , y = new_y3)
  test_extend_right <- data.frame(id = teid3, X = Xexright, Z = new_Z3 , y = new_y3)
  
  return(list(test = test_data, new_theta = new_theta3, test_extend_left = test_extend_left , 
              test_extend_right = test_extend_right ))
}

# fixed = fixef(obj.full); nn= nn; newX = test3[ ,2:(p+1)]; 
# newZ = as.matrix(test3[ ,(p+2)]); ynew = test3$y; 
# md =  dim(test3)[1]/nn ; random = ranef(obj.full)$id 

# obtain the group mean
group_mean <- function(data, nn, md){
  md <- dim(data)[1]/nn
  index <- seq(from = 1, to= (md), by = nn)  
  
  group_mean <-
    data %>%
    mutate(index = rep(1:(nrow(data)/nn), each = nn) )%>%
    group_by(index) %>%
    summarise_if(is.numeric, mean )%>%
    left_join(
      data %>%
        mutate(index = rep(1:(nrow(data)/nn), each = nn) )%>%
        group_by(index) %>%
        summarise( id = unique(id)),
      c("index" = "index")
    )%>%
    select(-"index")
  
  return(group_mean)

}




# obtain mixed effect
cmmpmixed <- function(fixed_ef, random_ef, group_mean){
  md <- dim(group_mean)[1]
  xbar <-
    group_mean %>%
    select_at(vars(starts_with("X"))) %>%
    as.matrix() 
  
  zbar <- 
    group_mean %>%
    select_at(vars(starts_with("Z"))) %>%
    as.matrix() 
  
  ybar <- 
    group_mean %>%
    select_at(vars(starts_with("y")))  %>%
    as.matrix() 
  
  fixed.temp <-  xbar %*% as.matrix(fixed_ef)
  randomall <-  zbar %*% t(random_ef)
  
  # CMMP
  cmmpmixed <- matrix(0, md) 
  I <- matrix(0, md)
  for (k in 1 : md){
    theta.temp <- fixed.temp[k]+ randomall[k,]   # assume match = true
    sp.temp <- theta.temp^2-2 * theta.temp * ybar[k]
    I.temp <- match(min(sp.temp),sp.temp)
    I[k,]<- I.temp
    cmmpmixed[k,]<- theta.temp[I.temp]
  }

  return(list(cmmpmixed = cmmpmixed, I= I))
  
}



# The Mode function return the mode, if ties occurs, the largest one is chosen as suggested by WANG(2010)
Mode <- function(y) {
  x <- sort(y)
  ux <- unique(x)
  index <- which(tabulate(match(x, ux))==max(tabulate(match(x, ux))))
  ux[index[length(index)]]
}

# function that estimate instability measure s for k = 2, ..., K (Step 1-3 of the algorithm)
estimateS <- function(y , m , K){
  n<-length(y)
  # step 1. data permutation
  yc <- sample(y, size = length(y), replace = FALSE)
  # step 2. Split the permuted data into three parts
  yc1<- yc[1:m]
  yc2<- yc[(m+1):(2*m)]
  yc3<- yc[(2*m+1):n]
  # step 3. Construct hierarchal clustering based on yc1 and yc2 for k from 2 to K
  clusters1 <- hclust(dist(yc1),"ave")
  clusters2 <- hclust(dist(yc2),"ave")
  s<-c()
  for (k in 2:K){
    clusterCut1 <- cutree(clusters1, k)
    clusterCut2 <- cutree(clusters2, k)
    knn1 <- knn(train = as.matrix(yc1), test = as.matrix(yc3), cl = clusterCut1, k = 10, prob=TRUE)
    testcl1<-as.numeric(knn1) # prediction test label through nearest neighbour algorithm
    knn2 <- knn(train = as.matrix(yc2), test = as.matrix(yc3), cl = clusterCut2, k = 10, prob=TRUE)
    testcl2<-as.numeric(knn2)
    v<-NULL
    for (i in 1: (length(yc3)-1)){
      for (j in (i+1): length(yc3)){
        v<-c(v,(testcl1[i]==testcl1[j])+(testcl2[i]==testcl2[j]))
      }
    }
    s[k]<-sum(v==1)
  }
  # step 4. Compute khat as the argmin of s and return s
  smin <- min(s, na.rm=TRUE)
  khat <- which(s==smin)
  return(list(kvote=khat[length(khat)], s=s))
}

# function that estimate group k for given s with voting or averaging 
estimateK <- function(y , m , K , C = 50,  method=c("vote","ave")){
  obj<-replicate(C,estimateS(y = y, m = m, K = K))
  if (method=="vote"){
    k <- Mode(as.numeric(unlist(obj[1,])))
  }else if(method=="ave"){
    smatrix <- matrix(unlist(obj[2,]), nrow=C, byrow = TRUE)
    save <- colMeans(smatrix)
    sminave <- min(save, na.rm=TRUE)
    k <- match(sminave,save)
  }
  return(k)
}


cmmp.theta <- function(mu.n, raneff,y){
  theta.cmmp <- theta.ls <- stage.cmmp <- c()	
  N = length(mu.n)
  for (i in 1: N){
    mu.est <- mu.n[i]
    # CMMP 	
    theta.temp <- mu.est + raneff
    sp.temp <- theta.temp^2-2*theta.temp*y[i]
    I.temp <- match(min(sp.temp),sp.temp)
    theta.cmmp[i] <- theta.temp[I.temp]
    stage.cmmp[i] <- I.temp
  }
  return(theta.cmmp)
}

########## some necessary functions end ############
####################################################