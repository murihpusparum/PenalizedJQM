library(lqmm)
library(invgamma)

#data generation for e_ij and u_i follow N
sim_N <- function(N, n, alpha, beta, psiu, a, b, seed){
  set.seed(seed)
  tau1 = alpha/2
  tau2 = 1-tau1
  
  df<-data.frame(subject=rep(1:N,rep(n,N)),
                 time=rep(1:n,N))
  #simulate u
  u<-rnorm(N,mean=0,sd=1)
  #simulate epsilon with N(0,1)
  epsilon = rnorm(N*n, mean = 0, sd = 1)
  #simulated data 
  y <- beta + rep(u,rep(n,N)) + epsilon
  df$y<-y
  
  return(list(d=df, beta1=beta1, beta2=beta2))
}

#data generation for e_ij follow N and u_i follow Chisq
sim_N_chi <- function(N, n, alpha, beta, psiu, a, b, seed){
  set.seed(seed)
  tau1 = alpha/2
  tau2 = 1-tau1
  
  df<-data.frame(subject=rep(1:N,rep(n,N)),
                 time=rep(1:n,N))
  #simulate u
  u<-rchisq(N, df=4)
  #scaled to zero mean and variance=1
  u_i <- (u-4)/(sqrt(2*4))
  #simulate epsilon with N(0,1)
  epsilon = rnorm(N*n, mean = 0, sd = 1)
  #simulated data 
  y <- beta + rep(u_i,rep(n,N)) + epsilon
  df$y<-y
  
  return(list(d=df, beta1=beta1, beta2=beta2))
}

lqmm_call <- function(N, n, alpha, beta, NSim, psiu, a, b, seed, FUN, dist){
  tau1 = alpha/2
  tau2 = 1-tau1
  
  #object initializations
  est <- data.frame(matrix(nrow=NSim, ncol=2))
  warn1a <- NA
  warn_out1 <- NA
  warn2a <- NA
  warn_out2 <- NA
  u1 <- list()
  u2 <- list()
  Results<-matrix(nrow=NSim,ncol=10)
  Results<-as.data.frame(Results)
  names(Results)<-c("N", "n", "alpha", "a", "b", "beta01", "beta02",
                    "EmpCov.subj","EmpCov.time","EmpCov")
  
  #function to exclude dataset with warning
  catch_lqmm <- function(theta, cov_name){
    out <- tryCatch(covHandling(theta = theta, n = 1, cov_name = cov_name, quad_type = "normal"), 
                    error=function(e){
                      message("catch error")
                      print(e)
                    },
                    warning=function(w){
                      message("catch warning")
                      print(w)
                    })
    out <- inherits(out, what = "warning")
    return(out)
  }
  
  for (i in 1:NSim) {
    df.l <- FUN(N=N+1, n=n+1, alpha=alpha, beta=beta, psiu=psiu, a=a, b=b, seed=seed+i)
    df <- df.l$d
    df.fit<-df[(df$time<=n)&(df$subject<=N),]
    fit1 <- lqmm(fixed = y ~ 1, random = ~ 1, group = subject, tau = tau1,
                 nK = 11, type = "normal", data = df.fit, 
                 control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 100, UP_tol = 1e-5,
                                       beta = 0.5, gamma = 1))
    fit2 <- lqmm(fixed = y ~ 1, random = ~ 1, group = subject, tau = tau2,
                 nK = 11, type = "normal", data = df.fit, 
                 control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 100, UP_tol = 1e-5,
                                       beta = 0.5, gamma = 1))
    
    warn_out1[i] <- catch_lqmm(theta = fit1$theta_z, cov_name = fit1$cov_name)
    warn_out2[i] <- catch_lqmm(theta = fit2$theta_z, cov_name = fit2$cov_name)
    
    warn1a[i] <- catch_lqmm(theta = fit1$theta_z, cov_name = fit1$cov_name)
    warn2a[i] <- catch_lqmm(theta = fit2$theta_z, cov_name = fit2$cov_name)
    
    warn1 <- catch_lqmm(theta = fit1$theta_z, cov_name = fit1$cov_name)
    warn2 <- catch_lqmm(theta = fit2$theta_z, cov_name = fit2$cov_name)
    
    j <- 1
    while (warn1[j] == TRUE | warn2[j] == TRUE) {
      df.l <- FUN(N=N+1, n=n+1, alpha=alpha, beta=beta, psiu=psiu, a=a, b=b, seed=seed+(NSim+i+j))
      df <- df.l$d
      df.fit <-df[(df$time<=n)&(df$subject<=N),]
      fit1 <- lqmm(fixed = y ~ 1, random = ~ 1, group = subject, tau = tau1,
                   nK = 11, type = "normal", data = df.fit, 
                   control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 100, UP_tol = 1e-5,
                                         beta = 0.5, gamma = 1))
      fit2 <- lqmm(fixed = y ~ 1, random = ~ 1, group = subject, tau = tau2,
                   nK = 11, type = "normal", data = df.fit, 
                   control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 100, UP_tol = 1e-5,
                                         beta = 0.5, gamma = 1))
      
      j <- j+1
      warn1[j] <- catch_lqmm(theta = fit1$theta_z, cov_name = fit1$cov_name)
      warn2[j] <- catch_lqmm(theta = fit2$theta_z, cov_name = fit2$cov_name)
      
    }
    warn1a <- append(warn1a, warn1)
    warn2a <- append(warn2a, warn2)
    
    est[i,1] <- fit1$theta[1]
    est[i,2] <- fit2$theta[1]
    #random effect
    u_lqmm1 = data.frame(ranef(fit1))
    colnames(u_lqmm1) <- c("u_hat")  
    u1 <- append(u1, u_lqmm1)
    
    #random effect
    u_lqmm2 = data.frame(ranef(fit2))
    colnames(u_lqmm2) <- c("u_hat")  
    u2 <- append(u2, u_lqmm2)
    
    #coverage when there is one new measurement
    coverage<-numeric(N)
    for(s in 1:N) {
      y <- df$y[(df$subject==s) & (df$time==(n+1))]
      coverage[s]<-((y > est[i,1]+u_lqmm1[,1])&
                      (y < est[i,2]+u_lqmm2[,1]))
    }
    coverage<-mean(coverage)
    
    #coverage when there is a new subject
    df.i <- df[df$subject==N+1,]
    fit1.i <- lqmm(fixed = y ~ 1, random = ~ 1, group = subject, tau = tau1,
                   nK = 11, type = "normal", data = df.i, 
                   control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 100, UP_tol = 1e-5,
                                         beta = 0.5, gamma = 1))
    fit2.i <- lqmm(fixed = y ~ 1, random = ~ 1, group = subject, tau = tau2,
                   nK = 11, type = "normal", data = df.i, 
                   control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 100, UP_tol = 1e-5,
                                         beta = 0.5, gamma = 1))
    u1.i <- data.frame(ranef(fit1.i))[1,1]
    u2.i <- data.frame(ranef(fit2.i))[1,1]
    coverage2<-mean((df.i$y > est[i,1]+u1.i)& (df.i$y < est[i,2]+u2.i))
    
    Results[i,1:10] <- c(N, n, alpha, a, b, est[i,1], est[i,2], 
                     coverage2, coverage, mean(c(coverage,coverage2)))
    mean.Results <- colMeans(Results)
  }
  
  save(Results, file=paste("GB_lqmm_",dist, "_", i, "_alpha_",alpha, "_N_",N, "_n_", n, "_beta_", 
                           beta, ".RData",sep = ""))
  save(u1, file=paste("GB_lqmm_u1_",dist, "_", i, "_alpha_",alpha, "_N_",N, "_n_", n, "_beta_", 
                       beta, ".RData",sep = ""))
  save(u2, file=paste("GB_lqmm_u2_",dist, "_", i, "_alpha_",alpha, "_N_",N, "_n_", n, "_beta_", 
                       beta, ".RData",sep = ""))
  return(mean.Results)
}
