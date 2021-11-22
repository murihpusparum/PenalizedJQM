source("./codes/main/RQPD_Function.R")
source("./codes/simulation/simdata_sd.R")

SimStudy<-function(nr, N, n, beta, a, b, alpha, NSim, seed, FUN) {
  Results<-matrix(nrow=NSim,ncol=16)
  Results<-as.data.frame(Results)
  names(Results)<-c("iter","alpha","beta01","beta02",
                    "lambda",
                    "Cov.subj","Cov.time","Cov",
                    "EmpCov.subj","EmpCov.subj1obs","EmpCov.time","EmpCov", "seed",
                    "sd","avg","median")
  
  #generate seed of data without LQMM warning
  source("./codes/simulation/generate_data_sd.R")
  seed.new<-d(alpha, N, n, beta, a, b, seed, NSim, FUN)
  
  for(si in 1:NSim) {
    seed<-seed.new$seed$new[si]
    set.seed(seed)
    #modify SimData function according to the target scenario
    d<-FUN(N=N+1,n=n+1,beta=beta,a=a,b=b,
           alpha=alpha)
    df<-d$df
    subjects<-unique(df$subject)
    
    #use all subjects, but exclude the last time point
    df.fit<-df[(df$time<=n),]
    
    res<-rqpd.fun(db=df.fit,alpha=alpha,lambda=seq(0.1,5,0.5))
    beta01<-res$beta01[[1]]
    beta02<-res$beta02[[1]]
    alpha.i<-res$alpha.i
    
    Results[si,1:8]<-c(si,alpha,beta01,beta02,res$lambda,
                       res$cov.subj,res$cov.time,res$cov.tot)
    
    #calculate the empirical coverage for a new measurement (in-sample)
    coverage<-numeric(N+1)
    for(s in subjects) {
      db.y<-df[(df$subject==s),]
      y<-db.y$y[(db.y$time==max(db.y$time))]
      coverage[s]<-((y>beta01+alpha.i[[s]])&
                      (y<beta02+alpha.i[[s]]))
    }
    coverage<-mean(coverage)
    
    #estimate the parameters again, without the last subject
    #out-sample, without LOOCV
    df.fit2<-df[df$subject!=(N+1),]
    df.fit2<-df.fit2[df.fit2$time!=max(df.fit2$time),]
    
    res2<-rqpd.fun(db=df.fit2,alpha=alpha,lambda=seq(0.1,5,0.5))
    
    db.y.i<-df[df$subject==(N+1),]
    y.i<-db.y.i$y[db.y.i$time!=max(db.y.i$time)] # Hs of the new subject
    y.s<-db.y.i$y[db.y.i$time==max(db.y.i$time)] # future obs (ns) of the new subject
    coverage2<-rqpd.coverage.new.subject(res=res2,y=y.i,y.s=y.s,alpha=alpha)
    
    #out-sample, without LOOCV, one new observation
    y.i1<-db.y.i$y[db.y.i$time==(max(db.y.i$time)-1)] # Hs of the new subject
    coverage3<-rqpd.coverage.new.subject1obs(res=res2,y=y.i1,y.s=y.s,alpha=alpha)
    
    
    #out-sample, with LOOCV
    # coverage2<-numeric(N+1)
    # for (s in 1:N+1) {
    #   df.fit2<-df[df$subject!=s,]
    #   df.fit2<-df.fit2[df.fit2$time!=max(df.fit2$time),]
    #   
    #   res2<-rqpd.fun(db=df.fit2,alpha=alpha,lambda=seq(0.1,5,0.5))
    #   
    #   db.y.i<-df[df$subject==s,]
    #   y.i<-db.y.i$y[db.y.i$time!=max(db.y.i$time)] # Hs of the new subject
    #   y.s<-db.y.i$y[db.y.i$time==max(db.y.i$time)] # future obs (ns) of the new subject
    #   coverage2[s]<-rqpd.coverage.new.subject(res=res2,y=y.i,y.s=y.s,alpha=alpha)
    # 
    # }
    # coverage2<-mean(coverage2)
    
    Results[si,9:12]<-c(coverage2,coverage3,coverage,mean(c(coverage,coverage2)))
    Results[si,13:16]<-c(seed,sd(d$sd),mean(d$sd),median(d$sd))
    
    #save results for every 10 iterations
    if((si%%10)==0) {
      save(Results, file=paste("./output/", nr, "_RQPDTmp_", si, "_alpha_",alpha, "_N_",N, "_n_", n, "_beta_", 
                               beta, "_a_", a, "_b_", b,  
                               ".Rdata",sep = ""))
    }
    gc(verbose = FALSE)
  }
  save(Results, file=paste("./output/output_final/", nr, "_RQPDTmp_", si, "_alpha_",alpha, "_N_",N, "_n_", n, "_beta_", 
                           beta, "_a_", a, "_b_", b,  
                           ".Rdata",sep = ""))
  
}