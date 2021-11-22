source("./codes/main/LQMM_Function.R")
source("./codes/simulation/simdata_sd.R")

SimStudy<-function(nr, N, n, alpha, beta, NSim, a, b, seed, FUN){
  Results<-matrix(nrow=NSim,ncol=16)
  Results<-as.data.frame(Results)
  names(Results)<-c("iter", "N", "n", "alpha", "a", "b", "beta01", "beta02",
                    "EmpCov.subj","EmpCov.subj1obs","EmpCov.time","EmpCov","seed",
                    "sd","avg","median")
  
  #generate seed of data without LQMM warning
  source("./codes/simulation/generate_data_sd.R")
  seed.new<-d(alpha, N, n, beta, a, b, seed, NSim, FUN)
  
  for (si in 1:NSim) {
    seed<-seed.new$seed$new[si]
    res<-lqmm.fun(alpha=alpha, N=N, n, beta=beta, a=a, b=b, seed=seed, si=si, FUN=FUN) 
    Results[si,1:8]<-c(si,N,n,alpha,a,b,res$beta01,res$beta02)
    
    df<-res$df
    subjects<-unique(df$subject)
    coverage<-numeric(N+1)
    for(s in subjects) {
      y<-df$y[(df$subject==s) & (df$time==(n+1))]
      coverage[s]<-((y>res$beta01+res$u1[s])&
                      (y<res$beta02+res$u2[s]))
    }
    coverage<-mean(coverage)
    
    # #estimate the parameters again, without the last subject
    df.fit2<-df[df$subject!=(N+1),]
    df.fit2<-df.fit2[df.fit2$time!=max(df.fit2$time),]
    
    res2<-lqmm.fit(df.fit2,alpha)
    df.i<-df[df$subject==N+1,]
    y.s<-df.i$y[df.i$time==max(df.i$time)] # future obs (ns) of the new subject
    
    df.i<-df.i[df.i$time!=max(df.i$time),]
    y.i<-df.i$y # Hs of the new subject
    coverage2<-lqmm.cov.new.subject(res2,df.i,y.i,y.s,alpha)
    
    df.i1<-df[df$subject==N+1,]
    df.i1<-df.i[df.i$time==(max(df.i1$time)-1),]
    y.i1<-df.i1$y # Hs of the new subject
    coverage3<-lqmm.cov.new.subject(res2,df.i1,y.i1,y.s,alpha)
    
    Results[si,9:12]<-c(coverage2,coverage3,coverage,mean(c(coverage,coverage2)))
    Results[si,13:16]<-c(seed,sd(res$sd),mean(res$sd),median(res$sd))
    
    #or save results for every 10 iterations
    if((si%%10) ==0) {
      save(Results, file=paste("./output/", nr, "_LQMMTmp_", si, "_alpha_",alpha, "_N_",N, "_n_", n, "_beta_",
                               beta, "_a_", a, "_b_", b, ".Rdata",sep = ""))
    }
  }
  save(Results, file=paste("./output/output_final/", nr, "_LQMMTmp", "_", si, "_alpha_",alpha, "_N_",N, "_n_", n, "_beta_", 
                           beta, "_a_", a, "_b_", b, ".RData",sep=""))
  #  return(Results)
}
