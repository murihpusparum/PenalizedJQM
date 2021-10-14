library(invgamma)

# function for data generation - error~N(0, theta_i), u~N(0,1) 
SimData_NvarN<-function(N, n, a, b, beta, alpha=0.05) {
  # N number of subjects
  # n numbr of observations in each subject, for rand. effects ui
  # a and b: true parameters of inverse gamma
  # beta: intercept in LMM
  # alpha is the confindence level
  
  tau1 = alpha/2
  tau2 = 1-tau1
  
  beta.true = beta
  df<-data.frame(subject=rep(1:N,rep(n,N)),
                 time=rep(1:n,N))
  #simulate f(u)
  u<-rnorm(N,mean=0,sd=sqrt(1))
  
  #simulate epsilon with var ~ IG(a,b)
  theta_sq <- rinvgamma(N, shape = a, rate = b)
  theta_i <- sqrt(theta_sq)
  epsilon <- matrix(NA, nrow = n, ncol = N)
  epsilon <- list()
  for (i in 1:N) {
    epsilon_i <- list(rnorm(n, mean = 0, sd = sqrt(theta_sq[i])))
    epsilon <- append(epsilon, epsilon_i)
  }
  epsilon <- unlist(epsilon)
  
  #simulated data 
  y<- beta.true + rep(u,rep(n,N)) + epsilon
  df$y<-y
  
  return(df)
}

#function for data generation - error~N(0, theta_i), u~chisq(4)
SimData_Nvarchi<-function(N, n, a, b, beta, alpha=0.05) {
  tau1 = alpha/2
  tau2 = 1-tau1
  
  beta.true = beta
  df<-data.frame(subject=rep(1:N,rep(n,N)),
                 time=rep(1:n,N))
  #simulate f(u)
  u<-rchisq(N, df=4)
  #scaled to zero mean and variance=1
  u <- (u-4)/(sqrt(2*4))
  
  #simulate epsilon with var ~ IG(a,b)
  theta_sq <- rinvgamma(N, shape = a, rate = b)
  theta_i <- sqrt(theta_sq)
  epsilon <- matrix(NA, nrow = n, ncol = N)
  epsilon <- list()
  for (i in 1:N) {
    epsilon_i <- list(rnorm(n, mean = 0, sd = sqrt(theta_sq[i])))
    epsilon <- append(epsilon, epsilon_i)
  }
  epsilon <- unlist(epsilon)
  
  #simulated data 
  y<- beta.true + rep(u,rep(n,N)) + epsilon
  df$y<-y
  
  return(df)
}

# function for data generation - error~N(0, 1), u~N(0,1) 
SimData_NN<-function(N, n, a, b, beta, alpha=0.05) {
  tau1 = alpha/2
  tau2 = 1-tau1
  
  beta.true = beta
  df<-data.frame(subject=rep(1:N,rep(n,N)),
                 time=rep(1:n,N))
  #simulate f(u)
  u<-rnorm(N,mean=0,sd=sqrt(1))
  
  #simulate epsilon with N(0, 1)
  epsilon <- rnorm(N*n,mean=0,sd=sqrt(1))
  
  #simulated data 
  y<- beta.true + rep(u,rep(n,N)) + epsilon
  df$y<-y
  
  return(df)
}

# function for data generation - error~N(0, 1), u~chisq(4)
SimData_Nchi<-function(N, n, a, b, beta, alpha=0.05) {
  tau1 = alpha/2
  tau2 = 1-tau1
  
  beta.true = beta
  df<-data.frame(subject=rep(1:N,rep(n,N)),
                 time=rep(1:n,N))
  #simulate f(u)
  u<-rchisq(N, df=4)
  #scaled to zero mean and variance=1
  u<-(u-4)/(sqrt(2*4))
  
  #simulate epsilon with N(0, 1)
  epsilon <- rnorm(N*n,mean=0,sd=sqrt(1))
  
  #simulated data 
  y<- beta.true + rep(u,rep(n,N)) + epsilon
  df$y<-y
  
  return(df)
}

#function for data generation - error~t(3), u~N(0,1)
SimData_tN<-function(N, n, a, b, beta, alpha=0.05) {
  tau1 = alpha/2
  tau2 = 1-tau1
  
  beta.true = beta
  df<-data.frame(subject=rep(1:N,rep(n,N)),
                 time=rep(1:n,N))
  #simulate f(u)
  u<-rnorm(N,mean=0,sd=sqrt(1))
  
  #simulate epsilon with single skewed distribution, t(3)
  epsilon <- rt(N*n, df=3)
  #scaled to zero mean and variance=1
  epsilon <- epsilon/(3/(3-2))
  
  #simulated data 
  y<- beta.true + rep(u,rep(n,N)) + epsilon
  df$y<-y
  
  return(df)
}

#function for data generation - error~t(3), u~chisq(4)
SimData_tchi<-function(N, n, a, b, beta, alpha=0.05) {
  tau1 = alpha/2
  tau2 = 1-tau1
  
  beta.true = beta
  df<-data.frame(subject=rep(1:N,rep(n,N)),
                 time=rep(1:n,N))
  #simulate f(u)
  u<-rchisq(N, df=4)
  #scaled to zero mean and variance=1
  u <- (u-4)/(sqrt(2*4))
  
  #simulate epsilon with single skewed distribution, t(3)
  epsilon <- rt(N*n, df=3)
  #scaled to zero mean and variance=1
  epsilon <- epsilon/(3/(3-2))
  
  #simulated data 
  y<- beta.true + rep(u,rep(n,N)) + epsilon
  df$y<-y
  
  return(df)
}
