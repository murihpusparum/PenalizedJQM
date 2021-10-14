library(lqmm)
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
              u.i1=ranef(fit1)[,1], u.i2=ranef(fit1)[,1]))
}

lqmm.fun<-function(alpha, N, n, beta, a, b, seed, si, FUN){
  #generate data that can fit the model without warning
  warn1<-NA
  warn2<-NA
  
  set.seed(seed)
  df<-FUN(N=N, n=n, alpha=alpha, beta=beta, a=a, b=b)
  subjects<-unique(df$subject)
  
  fit<-lqmm.fit(df,alpha = alpha)
  
  warn1<-catch.warn.lqmm(theta=fit$theta_z1, cov_name=fit$cov_name1)
  warn2<-catch.warn.lqmm(theta=fit$theta_z2, cov_name=fit$cov_name2)
  
  coverage<-numeric(N)
  cnt<-1
  for(s in subjects) {
    y<-df$y[(df$subject==s)]

    coverage[cnt]<-mean((y>fit$beta01+fit$u.i1[cnt])&
                      (y<fit$beta02+fit$u.i2[cnt]))
    cnt<-cnt+1
  }
  coverage<-mean(coverage)

  j<-1
  while (warn1==TRUE | warn2==TRUE) {
    seed<-seed+j
    set.seed(seed)
    df<-FUN(N=N, n=n, alpha=alpha, beta=beta, a=a, b=b)
    fit<-lqmm.fit(df,alpha = alpha)
    
    warn1<-catch.warn.lqmm(theta=fit$theta_z1, cov_name=fit$cov_name1)
    warn2<-catch.warn.lqmm(theta=fit$theta_z2, cov_name=fit$cov_name2)
    
    coverage<-numeric(N)
    cnt<-1
    for(s in subjects) {
      y<-df$y[(df$subject==s)]
      
      coverage[cnt]<-mean((y>fit$beta01+fit$u.i1[s])&
                            (y<fit$beta02+fit$u.i2[s]))
      cnt<-cnt+1
    }
    coverage<-mean(coverage)
    
    j<-j+1
  }
  
  
  
  return(list(beta01=fit$beta01,
              beta02=fit$beta02,
              u1=fit$u.i1, u2=fit$u.i2,
              df=df, cov=coverage, seed=seed))
}