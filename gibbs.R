#Applied Bayesian Methods
#Homework 04
library(gtools)
library(stats)

#2  Gibbs sampling for 1-dim case: ###############################

set.seed(123)
## Data(X)
x = faithful$eruptions
hist(x)

## Define function
### sample_y <- sampling y
propotion <- function(x){return(x/sum(x))}

sample_y <- function(mu, tau, alpha, x){
  p.y.given.x <- matrix(0, nrow = length(alpha), ncol = length(x))
  for(i in 1:length(alpha)){
    p.y.given.x[i,] = alpha[i] * dnorm(x, mu[i], sqrt(1/tau[i]))
  }
  p.y.given.x = apply(p.y.given.x, 2, propotion) 
  y = rep(0, length(x))
  for(i in 1:length(y)){
    y[i] = sample(1:length(alpha), size = 1, 
                  prob = p.y.given.x[,i], replace=TRUE)
  }
  return(y)
}

### sample_mu <- sampling mu
sample_mu <- function(prior, tau, k, y, x){
  df = data.frame(x = x, y = y)
  mu = rep(0, k)
  for(i in 1:k){
    sample.size = sum(y == i)
    sample.mean = ifelse(sample.size == 0, 0, mean(x[y == i]))
    post.prec = sample.size * tau[i] + prior$prec[i]
    post.mean = (prior$prec[i] * prior$mean[i] + sample.size * tau[i] * sample.mean) / post.prec
    mu[i] = rnorm(1, post.mean, sqrt(1/post.prec))
  }
  return(mu)
}

### sample_tau <- sampling tau
sample_tau <- function(mu, prior, k, y, x){
  df = data.frame(x = x, y = y)
  tau = rep(0, k)
  for(i in 1:k){
    post.a = prior[1] + sum(y == i)/2
    post.b = prior[2] + sum((x[y == i]-mu[i])^2)/2
    tau[i] = rgamma(1, post.a, post.b)
  }
  return(tau)
}

### sample_alpha <- sampling alpha
sample_alpha <- function(prior, y){
  counts = colSums(outer(y, 1:length(prior), FUN = "=="))
  alpha = rdirichlet(1, prior + counts)
  return(alpha)
}

## Combine functions above to stimulate gibbs
gibbs <- function(x, k = 2, niter = 100, mu.prior = list(mean = c(1, 1), prec = c(1, 1)),
                  tau.prior = c(1, 1), alpha.prior = c(1, 1)){
  # set output list
  res = list(y = matrix(nrow = niter, ncol = length(x)), mu = matrix(nrow = niter, ncol = k), 
             tau = matrix(nrow = niter, ncol = k), alpha = matrix(nrow = niter,ncol = k))
  
  # initial value
  mu = rnorm(k, mu.prior$mean, mu.prior$prec)
  tau = rgamma(k, tau.prior[1], tau.prior[2])
  alpha = rdirichlet(1, alpha.prior)
  
  # sampling y
  y = sample_y(mu, tau, alpha, x)
  
  # get the first row of output list
  res$y[1,]=y
  res$mu[1,]=mu
  res$tau[1,]=tau
  res$alpha[1,]=alpha
  
  # get the rest of row of output list
  for(i in 2:niter){
    mu = sample_mu(mu.prior, tau, k, y, x)
    tau = sample_tau(mu, tau.prior, k, y, x)
    alpha = sample_alpha(alpha.prior, y)
    y = sample_y(mu, tau, alpha, x)
    res$y[i,]=y
    res$mu[i,]=mu
    res$tau[i,]=tau
    res$alpha[i,]=alpha
  }
  return(res)
}


res = gibbs(x, k=2, niter = 100)
par(mfrow=c(2, 2))
hist(x)
plot(res$mu[,1], ylim=c(-2,6),type="l", lty = 1, ylab = "mu")
lines(res$mu[,2], lty = 2)
grid()
legend("bottomright", legend = c("Line 1", "Line 2"),
       lty = 1:2, cex = 0.8)
plot(res$tau[,1], ylim=c(0,20),type="l", lty = 1, ylab = "tau")
lines(res$tau[,2], lty = 2)
grid()
legend("bottomright", legend = c("Line 1", "Line 2"),
       lty = 1:2, cex = 0.8)
plot(res$alpha[,1], ylim=c(0,1),type="l", lty = 1, ylab = "alpha")
lines(res$alpha[,2], lty = 2)
grid()
legend("bottomright", legend = c("Line 1", "Line 2"),
       lty = 1:2, cex = 0.8)




