source('./codes/main/JQM_Function.R')
source('./codes/main/RQPD_Function.R')
source('./codes/main/LQMM_Function.R')

source('./codes/simulation/simdata_sd.R')
source('./codes/simulation/generate_data_sd.R')

SimStudy.length<-function(nr, N, n, beta, a, b, alpha, NSim, seed, FUN) {
  
  res.logit<-matrix(nrow=NSim,ncol=20)
  res.logit<-as.data.frame(res.logit)
  names(res.logit)<-c("iter","N","n","alpha",
                      "EmpCov.jqm","mean.length.jqm","median.length.jqm","min.length.jqm","max.length.jqm",
                      "beta0.logit.jqm","beta1.logit.jqm",
                      "EmpCov.rqpd","length.rqpd","beta0.logit.rqpd","beta1.logit.rqpd",
                      "EmpCov.lqmm","mean.length.lqmm","median.length.lqmm","min.length.lqmm","max.length.lqmm")

  seed.new<-d(alpha, N, n, beta, a, b, seed, NSim, FUN)
  
  len.all<-NA
  for(si in 1:NSim) {
    seed<-seed.new$seed$new[si]
    set.seed(seed)
    
    res<-lqmm.fun(alpha=alpha, N=N, n, beta=beta, a=a, b=b, seed=seed, si=si, FUN=FUN) 
    df<-res$df
    subjects<-unique(df$subject)
    #calculate the empirical coverage for a new subject and a new measurement
    ys.lqmm<-numeric(N)
    coverage.lqmm<-numeric(N)
    length.lqmm<-numeric(N)
    for(s in subjects) {
      y<-df$y[(df$subject==s)&(df$time==(n+1))]
      ys.lqmm[s]<-y
      coverage.lqmm[s]<-((y>res$beta01+res$u1[s])&
                          (y<res$beta02+res$u2[s]))
      length.lqmm[s]<-((res$beta02+res$u2[s])-(res$beta01+res$u1[s]))
    }
    EmpCov.lqmm<-mean(coverage.lqmm)
    
    method<-"LQMM"
    len.lqmm<-cbind.data.frame(si,method,ys.lqmm,coverage.lqmm,length.lqmm)
    colnames(len.lqmm)<-c("iter","method","ys","coverage","length")
    
    set.seed(seed)
    d<-FUN(N=N+1,n=n+1,beta=beta,a=a,b=b,alpha=alpha)
    df<-d$df
    subjects<-unique(df$subject)
    
    #use all subjects, but exclude the last time point
    df.fit<-df[(df$time<=n),]
    
    # for JQM
    res<-jqm(db=df.fit,
             alpha=alpha,
             lambda.u.seq = seq(0.5,4,0.5),
             lambda.z.seq = seq(0.5,5,0.5))
    #calculate the empirical coverage for a new subject and a new measurement
    ys.jqm<-numeric(N)
    coverage.jqm<-numeric(N)
    length.jqm<-numeric(N)
    for(s in subjects) {
      y<-df$y[(df$subject==s)&(df$time==(n+1))]
      ys.jqm[s]<-y
      coverage.jqm[s]<-((y>res$beta0+res$u[s]+res$z[s]*res$beta1)&
                          (y<res$beta0+res$u[s]+res$z[s]*res$beta2))
      length.jqm[s]<-((res$beta0+res$u[s]+res$z[s]*res$beta2)-(res$beta0+res$u[s]+res$z[s]*res$beta1))
    }
    EmpCov.jqm<-mean(coverage.jqm)
    
    method<-"JQM"
    len.jqm<-cbind.data.frame(si,method,ys.jqm,coverage.jqm,length.jqm)
    colnames(len.jqm)<-c("iter","method","ys","coverage","length")
    
    pairs<-cbind.data.frame(coverage.jqm, length.jqm)
    logit.jqm<-glm(coverage.jqm ~ abs(length.jqm), data=pairs, family = binomial)
    beta0.jqm<-logit.jqm$coefficients[1]
    beta1.jqm<-logit.jqm$coefficients[2]
    
    set.seed(seed)
    res<-rqpd.fun(db=df.fit,alpha=alpha,lambda=seq(0.1,5,0.5))
    beta01<-res$beta01[[1]]
    beta02<-res$beta02[[1]]
    alpha.i<-res$alpha.i

    #calculate the empirical coverage for a new measurement (in-sample)
    ys.rqpd<-numeric(N+1)
    coverage.rqpd<-numeric(N+1)
    length.rqpd<-numeric(N+1)
    for(s in subjects) {
      db.y<-df[(df$subject==s),]
      y<-db.y$y[(db.y$time==max(db.y$time))]
      ys.rqpd[s]<-y
      coverage.rqpd[s]<-((y>beta01+alpha.i[[s]])&
                           (y<beta02+alpha.i[[s]]))
      length.rqpd[s]<-((beta02+alpha.i[[s]])-(beta01+alpha.i[[s]]))
    }
    EmpCov.rqpd<-mean(coverage.rqpd)
    
    method<-"RQPD"
    len.rqpd<-cbind.data.frame(si,method,ys.rqpd,coverage.rqpd,length.rqpd)
    colnames(len.rqpd)<-c("iter","method","ys","coverage","length")
    
    len.all<-rbind(len.all,len.lqmm,len.jqm,len.rqpd)
    
    pairs<-cbind(pairs, coverage.rqpd, length.rqpd)
    logit.rqpd<-glm(coverage.rqpd ~ abs(length.jqm), data=pairs, family = binomial)
    beta0.rqpd<-logit.rqpd$coefficients[1]
    beta1.rqpd<-logit.rqpd$coefficients[2]
    
    res.logit[si,1:4]<-c(si,N,n,alpha)
    res.logit[si,5:11]<-c(EmpCov.jqm,mean(length.jqm),median(length.jqm),min(length.jqm),max(length.jqm),beta0.jqm,beta1.jqm)
    res.logit[si,12:15]<-c(EmpCov.rqpd,mean(length.rqpd),beta0.rqpd,beta1.rqpd)
    res.logit[si,16:20]<-c(EmpCov.lqmm,mean(length.lqmm),median(length.lqmm),min(length.lqmm),max(length.lqmm))
  }
  save(res.logit, file=paste("./output/output_final/", nr, "_res.logit.jqm.rqpd_", si, "_alpha_",alpha, "_N_",N, "_n_", n, "_beta_", 
                           beta, "_a_", a, "_b_", b, ".Rdata",sep = ""))
  save(len.all, file=paste("./output/output_final/", nr, "_len.all_", si, "_alpha_",alpha, "_N_",N, "_n_", n, "_beta_", 
                             beta, "_a_", a, "_b_", b, ".Rdata",sep = ""))
}

SimStudy.length(nr=4, N=30, n=5, alpha=0.05, beta=0, NSim=20, a=2, b=0.4, seed=7606663, FUN=SimData_NvarN)