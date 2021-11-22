library(lqmm)
library(ald)

lqmm.cov.new.subject<-function(res,df.i,y.i,y.s,alpha,side="both"){
  tau1=alpha/2
  tau2=1-tau1
  
  Z_new<-model.matrix(~1, data.frame(rep(1, nrow(df.i))))
  X_new<-model.matrix(~1, data.frame(rep(1, nrow(df.i))))
  y_new<-y.i
  
  n<-length(unique(df.i$subject))
  id_nomiss<-rep(seq(1:n), each=nrow(df.i))
  
  # for lower bound
  sigma1<-sqrt(varAL(sigma = res$scale1, tau=tau1)) # variance of the model error
  mean.e1<-meanALD(mu=0, sigma=res$scale1, p=tau1)
  betas1<-res$beta01
  Vi_inv<-solve(Z_new %*% tcrossprod(res$D1, Z_new) + 
                  sigma1^2 * diag(sum(id_nomiss))) # Vi_inverse = Z_i D Zi_inv + sigma
  DZtVinv<-tcrossprod(res$D1, Z_new) %*% Vi_inv
  u1<-c(DZtVinv %*% (y_new - X_new %*% betas1 - mean.e1)) #u1 = D Zi_t Vi_inv (y - Xi beta1)
  
  # for upper bound
  sigma2<-sqrt(varAL(sigma = res$scale2, tau=tau2))
  mean.e2<-meanALD(mu=0, sigma=res$scale2, p=tau2)
  betas2<-res$beta02
  Vi_inv<-solve(Z_new %*% tcrossprod(res$D2, Z_new) + 
                  sigma2^2 * diag(sum(id_nomiss)))
  DZtVinv<-tcrossprod(res$D2, Z_new) %*% Vi_inv
  u2<-c(DZtVinv %*% (y_new - X_new %*% betas2 - mean.e2))
  
  if(side=="both"){
    coverage2<-mean((y.s>betas1+u1)& (y.s<betas2+u2))    
  }
  if(side=="up"){
    coverage2<-mean(y.s<betas2+u2)
  }
  if(side=="low"){
    coverage2<-mean(y.s>betas1+u1)
  }
  
  return(coverage2) 
}

#function to exclude dataset with warning
catch.warn.lqmm<-function(theta, cov_name){
  out<-tryCatch(covHandling(theta=theta, n=1, cov_name=cov_name, quad_type="normal"), 
                error=function(e){
                  message("catch error")
                  print(e)
                },
                warning=function(w){
                  message("catch warning")
                  print(w)
                })
  out<-inherits(out, what="warning")
  return(out)
}

lqmm.fit<-function(db,alpha){
  tau1=alpha/2
  tau2=1-tau1
  
  fit1<-lqmm(fixed=y ~ 1, random=~ 1, group=subject, tau=tau1,
             nK=11, type="normal", data=db, 
             control=lqmmControl(LP_max_iter=1000, LP_tol_ll=1e-4, UP_max_iter=1000, UP_tol=1e-4,
                                 beta=0.5, gamma=1))
  fit2<-lqmm(fixed=y ~ 1, random=~ 1, group=subject, tau=tau2,
             nK=11, type="normal", data=db, 
             control=lqmmControl(LP_max_iter=1000, LP_tol_ll=1e-4, UP_max_iter=1000, UP_tol=1e-4,
                                 beta=0.5, gamma=1))
  return(list(theta_z1=fit1$theta_z,
              theta_z2=fit2$theta_z,
              cov_name1=fit1$cov_name, cov_name2=fit2$cov_name,
              beta01=fit1$theta[1], beta02=fit2$theta[1],
              D1=VarCorr(fit1), D2=VarCorr(fit2),
              scale1=fit1$scale,scale2=fit2$scale,
              u.i1=ranef(fit1)[,1], u.i2=ranef(fit2)[,1]))
}

lqmm.fun<-function(alpha, N, n, beta, a, b, seed, si, FUN){
  #generate data that can fit the model without warning
  warn1<-NA
  warn2<-NA
  
  set.seed(seed)
  d<-FUN(N=N+1, n=n+1, alpha=alpha, beta=beta, a=a, b=b)
  subjects<-unique(d$df$subject)
  df.fit<-d$df[(d$df$time<=n),]
  fit<-lqmm.fit(df.fit,alpha = alpha)
  
  warn1<-catch.warn.lqmm(theta=fit$theta_z1, cov_name=fit$cov_name1)
  warn2<-catch.warn.lqmm(theta=fit$theta_z2, cov_name=fit$cov_name2)
  coverage<-numeric(N)
  
  cnt<-1
  for(s in subjects) {
    y<-d$df$y[(d$df$subject==s)&(d$df$time==(n+1))]
    coverage[cnt]<-mean((y>fit$beta01+fit$u.i1[cnt])&
                          (y<fit$beta02+fit$u.i2[cnt]))
    cnt<-cnt+1
  }
  coverage<-mean(coverage)
  
  j<-1
  while (warn1==TRUE | warn2==TRUE) {
    seed<-seed+si+j
    set.seed(seed)
    d<-FUN(N=N+1, n=n+1, alpha=alpha, beta=beta, a=a, b=b)
    df.fit<-d$df[(d$df$time<=n),]
    
    fit<-lqmm.fit(df.fit,alpha = alpha)
    warn1<-catch.warn.lqmm(theta=fit$theta_z1, cov_name=fit$cov_name1)
    warn2<-catch.warn.lqmm(theta=fit$theta_z2, cov_name=fit$cov_name2)
    
    coverage<-numeric(N)
    cnt<-1
    for(s in subjects) {
      y<-d$df$y[(d$df$subject==s)&(d$df$time==(n+1))]
      coverage[cnt]<-mean((y>fit$beta01+fit$u.i1[cnt])&
                            (y<fit$beta02+fit$u.i2[cnt]))
      cnt<-cnt+1
    }
    coverage<-mean(coverage)
    
    j<-j+1
  }
  
  return(list(beta01=fit$beta01,
              beta02=fit$beta02,
              u1=fit$u.i1, u2=fit$u.i2,
              df=d$df, cov=coverage, seed=seed,
              sd=d$sd))
}
