#### Homework 2
#-------------------------------------------------------------------------------
install.packages("mnormt")
install.packages("rmarkdown")
install.packages("moments")
install.packages("tidyverse")

## Load library
library(mnormt)
library(mvtnorm)
library(moments)
library(tidyverse)

#-------------------------------------------------------------------------------
#### Question 1
#-------------------------------------------------------------------------------
## Function that returns the p-value for a two-sample t-test
twosample.t <- function(X1, X2, n.perm)
{
  X <- c(X1, X2)
  group1 <- rep(1, length(X1))
  group2 <- rep(2, length(X2))
  data <- data.frame("X"=X, "group"=c(group1,group2))
  data$group <- as.factor(data$group)
  t.obs <- t.test(X~group, data=data)$statistic
  perm.data <- data
  t.star <- replicate(n.perm,{
    perm.data$group <- sample(data$group)
    t.perm <- t.test(X~group, data=perm.data)$statistic
  })
  mean(abs(t.obs) <= abs(t.star))
}

## Function that returns the p-value for a paired t-test
paired.t <- function(X1, X2, n.perm)
{
  data <- data.frame("X1"=X1,"X2"=X2)
  t.obs <- t.test(X1, X2, paired=TRUE)$statistic
  perm.data <- data
  x <- as.factor(c(1, 2))
  t.star <- replicate(n.perm,{
    x <- replicate(length(X1),{
      x <- sample(2)
    })
    for(i in 1:length(X1)){
      perm.data[i,] <- perm.data[i,x[,i]]
    }
    t.perm <- t.test(perm.data[,1], perm.data[,2], paired=TRUE)$statistic
  })
  mean(abs(t.obs) <= abs(t.star)) 
}

#-------------------------------------------------------------------------------
#### Question 2: Write a function to simulate the data
#-------------------------------------------------------------------------------
## Function to transform normal into skewed distribution
dist.skewed <- function(rho_skewed, delta, n){
  Y <- rmvnorm(n, mean = c(0,0),
               sigma = matrix(nrow=2, ncol=2, c(1, rho_skewed, rho_skewed, 1)))
  Z <- pnorm(Y)
  X1 <- qexp(Z[,1]) + delta
  X2 <- qexp(Z[,2])
  data<-data.frame("X1" = X1,"X2" = X2)
}

## Function to approximate rho_skewed needed to generate skewed distribution
rho.approximate <- function(rho, delta, n=1000000, eps=0.0001){
  rho_skewed <- rho + 0.01
  X <- dist.skewed(rho_skewed, delta, n)
  rho_true <- cor(X[,1], X[,2])
  while (abs(rho_true - rho) > eps) {
    rho_skewed <- rho_skewed + (rho - rho_true) * 0.1
    X <- dist.skewed(rho_skewed, delta, n)                                                                                                                                                    
    rho_true <- cor(X[,1], X[,2])
    print(c(rho, rho_true, rho_skewed))
  }
  rho_skewed
}

## Approximate rho_skewed for each correlation level
mean.X1 <- 20
correlation <- 1:5
rho <- c(0, 0.1, 0.25, 0.5, 0.7)
rho_skewed <- c()
for (i in correlation){
  delta <- 1 ## Results of approximation are independent of the value of delta
  rho_skewed[i] <- rho.approximate(rho[i], delta)
}

## Function to simulate data
simulate.data <- function(n, delta, correlation){
  rho_skewed <- rho_skewed[correlation]
  dist.skewed(rho_skewed, delta, n)
}

#-------------------------------------------------------------------------------
#### Question 3
#-------------------------------------------------------------------------------
rho <- c(0, 0.1, 0.25, 0.5, 0.7)
N_sim <- 100
n.perm <- 100

# Function that calculates empirical type I for two-sample test
typeI.twosample <- function(N_sim, n, n.perm, delta, rho){
  typeI_twosample <- c(0,0,0,0,0)
  for(i in correlation){
    t.test.count <- replicate(N_sim,{
      X <- simulate.data(n, delta, i)
      X1 <- X[,1]
      X2 <- X[,2]
      pvalue <- twosample.t(X1, X2, n.perm)
      error <- (pvalue < 0.05)
    })
    typeI_twosample[i] <- sum(t.test.count) / length(t.test.count)
  }
  plot(rho, typeI_twosample)
  return(typeI_twosample)
}

# Function that calculates empirical type II for two-sample test
typeII.twosample <- function(N_sim, n, n.perm, delta, rho){
  typeII_twosample <- c(0,0,0,0,0)
  for(i in correlation){
    t.test.count <- replicate(N_sim,{
      X <- simulate.data(n, delta, i)
      X1 <- X[,1]
      X2 <- X[,2]
      pvalue <- twosample.t(X1, X2, n.perm)
      error <- (pvalue > 0.05)                                                 
    })
    typeII_twosample[i] <- sum(t.test.count) / length(t.test.count)
  }
  plot(rho, typeII_twosample)
  return(typeII_twosample)
}

# Function that calculates empirical type I for paired t-test
typeI.paired <- function(N_sim, n, n.perm, delta, rho){
  typeI_paired <- c(0,0,0,0,0)
  for(i in correlation){
    tpair.test.count <- replicate(N_sim,{
      X <- simulate.data(n, delta, i)
      X1 <- X[,1]
      X2 <- X[,2]
      pvalue <- paired.t(X1, X2, n.perm)
      error <- (pvalue < 0.05)
    })
    typeI_paired[i] <- sum(tpair.test.count) / length(tpair.test.count)
  }
  plot(rho, typeI_paired)
  return(typeI_paired)
}

# Function that calculates empirical type II for paired t-test
typeII.paired <- function(N_sim, n, n.perm, delta, rho){
  typeII_paired <- c(0,0,0,0,0)
  for(i in correlation){
    tpair.test.count <- replicate(N_sim,{
      X <- simulate.data(n, delta, i)
      X1 <- X[,1]
      X2 <- X[,2]
      pvalue <- paired.t(X1, X2, 100)
      error <- (pvalue > 0.05)
    })
    typeII_paired[i] <- sum(tpair.test.count) / length(tpair.test.count)
  }
  plot(rho, typeII_paired)
  return(typeII_paired)
}

------------------------------------------------------------
#### Small sample size
n <-5

## Type I: simuleren onder H0
delta <- 0
typeI_twosample_small <- typeI.twosample(N_sim, n, n.perm, delta, rho)
typeI_paired_small <- typeI.paired(N_sim, n, n.perm, delta, rho)

## Type II: simuleren onder Ha
delta <- 2                                                                    
typeII_twosample_small <- typeII.twosample(N_sim, n, n.perm, delta, rho)
typeII_paired_small <- typeII.paired(N_sim, n, n.perm, delta, rho)

## Plot results
errors_n <- data.frame(rho, 
                     typeI_twosample_small, 
                     typeII_twosample_small, 
                     typeI_paired_small, 
                     typeII_paired_small
                    )

shapes<-c("Two sample type I"= 1, "Paired type I" = 2, 
          "Two sample type II" = 16, "Paired type II" = 17)
errors_n %>%                                                                    
  ggplot() +
  geom_point(aes(x=rho, y=typeI_twosample_small,  shape = "Two sample type I"), 
             size=2) +
  geom_point(aes(x=rho, y=typeII_twosample_small, shape = "Two sample type II"), 
             size=2) +
  geom_point(aes(x=rho, y=typeI_paired_small,     shape = "Paired type I"),     
             size=2) +
  geom_point(aes(x=rho, y=typeII_paired_small,   shape = "Paired type II"),     
             size=2) +
  labs(x="correlation", y="proportion", title="Errors for small  n", 
       shape="Legend")+
  scale_shape_manual(values=shapes) +
  theme_bw()

------------------------------------------------------------
#### Large sample size
N <- 100

### Type I: simuleren onder H0
delta <- 0
typeI_twosample_large <- typeI.twosample(N_sim, N, n.perm, delta, rho)
typeI_paired_large <- typeI.paired(N_sim, N, n.perm, delta, rho)

## Type II: simuleren onder Ha
delta <- 0.2                                                                     
typeII_twosample_large <- typeII.twosample(N_sim, N, n.perm, delta, rho)
typeII_paired_large <- typeII.paired(N_sim, N, n.perm, delta, rho)

## Plot results
errors_N <- data.frame(rho, 
                     typeI_twosample_large, 
                     typeII_twosample_large, 
                     typeI_paired_large, 
                     typeII_paired_large
                     )
shapes<-c("Two sample type I"= 1, "Paired type I" = 2, 
          "Two sample type II" = 16, "Paired type II" = 17)
errors_N %>%                                                                    
  ggplot() +
  geom_point(aes(x=rho, y=typeI_twosample_large,  shape = "Two sample type I"), 
             size=2) +
  geom_point(aes(x=rho, y=typeII_twosample_large, shape = "Two sample type II"), 
             size=2) +
  geom_point(aes(x=rho, y=typeI_paired_large,     shape = "Paired type I"),     
             size=2) +
  geom_point(aes(x=rho, y=typeII_paired_large,   shape = "Paired type II"),     
             size=2) +
  labs(x="correlation", y="proportion", title="Errors for large  n",  
       shape="Legend")+
  scale_shape_manual(values=shapes) +
  theme_bw()

#-------------------------------------------------------------------------------
#### Question 6
#-------------------------------------------------------------------------------
alpha <- 0.05
z_alpha <- qnorm(0.975)
exp.variance <- c(1, 1) 
width <- c(0,0,0,0,0)
coverage <- c(0,0,0,0,0)

## Simulate data:
N_sim <- 10000
n.samples <- 100

## Calculate empirical coverage and 95%-confidence intervals
for(cor in correlation){
  for (i in 1:N_sim){
    X <- simulate.data(n.samples, delta, cor)
    diff <- z_alpha * sqrt(sum(exp.variance) - 2*rho[cor]) / sqrt(length(X[,1]))
    conf_interval.left <- mean(X[,1]) - mean(X[,2]) - diff
    conf_interval.right <- mean(X[,1]) - mean(X[,2]) + diff
    coverage[cor] <- coverage[cor] + ((conf_interval.left < delta) & 
                                        (conf_interval.right > delta)) / N_sim
    width[cor] <- width[cor] +  2*diff/ N_sim
  }
}

## Plot the results
shapes<-c("Coverage" = 1, "Width" = 16)
conf.int <- data.frame(rho, coverage, width)
conf.int %>%
  ggplot() +
  geom_point(aes(x=rho, y=coverage, shape = "Coverage") )+
  geom_point(aes(x=rho, y=width, shape = "Width")) +
  labs(x="correlation", y="proportion", shape="Legend")+
  theme_bw()