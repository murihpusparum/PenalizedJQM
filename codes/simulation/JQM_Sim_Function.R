#function for running the simulation
SimStudy<-function(nr, N, n, beta, a, b, alpha, NSim, seed, FUN) {
  Results<-matrix(nrow=NSim,ncol=18)
  Results<-as.data.frame(Results)
  names(Results)<-c("iter","alpha","beta0","beta1","beta2",
                    "lambda.u","lambda.z",
                    "Cov.subj","Cov.time","Cov",
                    "EmpCov.subj","EmpCov.subj1obs","EmpCov.time","EmpCov", "seed",
                    "sd","avg","median")
  
  #generate seed of data without LQMM warning
  source("generate_data_sd.R")
  seed.new<-d(alpha, N, n, beta, a, b, seed, NSim, FUN)
  
  for(si in 1:NSim) {
    seed<-seed.new$seed$new[si]
    set.seed(seed)
    d<-FUN(N=N+1,n=n+1,beta=beta,a=a,b=b,alpha=alpha)
    df<-d$df
    subjects<-unique(df$subject)
    
    #use all subjects, but exclude the last time point
    df.fit<-df[(df$time<=n),]
    
    res<-jqm(db=df.fit,
             alpha=alpha,
             lambda.u.seq = seq(0.5,4,0.5),
             lambda.z.seq = seq(0.5,5,0.5))
    
    Results[si,1:10]<-c(si,alpha,res$beta0,res$beta1,res$beta2,
                        res$lambda.u,res$lambda.z,
                        res$cov.subj,res$cov.time,res$cov.tot)
    
    #calculate the empirical coverage for a new subject and a new measurement
    coverage<-numeric(N)
    for(s in subjects) {
      y<-df$y[(df$subject==s)&(df$time==(n+1))]
      coverage[s]<-((y>res$beta0+res$u[s]+res$z[s]*res$beta1)&
                      (y<res$beta0+res$u[s]+res$z[s]*res$beta2))
    }
    coverage<-mean(coverage)
    
    # #estimate the parameters again, without the last subject
    df.fit2<-df[df$subject!=(N+1),]
    df.fit2<-df.fit2[df.fit2$time!=max(df.fit2$time),]
    
    res2<-jqm(db=df.fit2,alpha=alpha,
              lambda.u.seq = seq(0.5,4,0.5),
              lambda.z.seq = seq(0.5,5,0.5))
    
    db.y.i<-df[df$subject==(N+1),]
    y.i<-db.y.i$y[db.y.i$time!=max(db.y.i$time)] # Hs of the new subject
    y.s<-db.y.i$y[db.y.i$time==max(db.y.i$time)] # future obs (ns) of the new subject
    coverage2<-jqm.coverage.new.subject(res=res2,y=y.i,y.s=y.s,alpha=alpha)
    
    #out-sample, without LOOCV, one new observation
    y.i1<-db.y.i$y[db.y.i$time==(max(db.y.i$time)-1)] # Hs of the new subject
    coverage3<-jqm.coverage.new.subject1obs(res=res2,y=y.i1,y.s=y.s,alpha=alpha)
    
    
    #using LOOCV
    # coverage2<-numeric(N+1)
    # for (s in 1:(N+1)) {
    #   df.fit2<-df[df$subject!=s,]
    #   df.fit2<-df.fit2[df.fit2$time!=max(df.fit2$time),]
    #   
    #   res2<-jqm(db=df.fit2,alpha=alpha,
    #             lambda.u.seq = seq(0.5,4,0.5),
    #             lambda.z.seq = seq(0.5,5,0.5))
    #   
    #   db.y.i<-df[df$subject==s,]
    #   y.i<-db.y.i$y[db.y.i$time!=max(db.y.i$time)] # Hs of the new subject
    #   y.s<-db.y.i$y[db.y.i$time==max(db.y.i$time)] # future obs (ns) of the new subject
    #   coverage2[s]<-jqm.coverage.new.subject(res=res2,y=y.i,y.s=y.s,alpha=alpha)
    #   
    # }
    # coverage2<-mean(coverage2)
    # 
    Results[si,11:14]<-c(coverage2,coverage3,coverage,mean(c(coverage,coverage2)))
    Results[si,15:18]<-c(seed,sd(d$sd),mean(d$sd),median(d$sd))
    
    #or save results for every 10 iterations
    if((si%%10) ==0) {
      save(Results, file=paste("./output/", nr, "_JQMTmp_", si, "_alpha_",alpha, "_N_",N, "_n_", n, "_beta_",
                               beta, "_a_", a, "_b_", b, ".Rdata",sep = ""))
    }
    gc(verbose = FALSE)
  }
  # #save final results
  save(Results, file=paste("./output/output_final/", nr, "_JQMTmp_", si, "_alpha_",alpha, "_N_",N, "_n_", n, "_beta_", 
                           beta, "_a_", a, "_b_", b, ".Rdata",sep = ""))
}
