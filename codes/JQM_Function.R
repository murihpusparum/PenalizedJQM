library(quantreg)
library(invgamma)

#function for data generation - error~N(0, theta_i), u~chisq(4)
SimData_chi<-function(N, n, psiu2, a, b, beta, alpha=0.05) {
  # N number of subjects
  # n numbr of observations in each subject, for rand. effects ui
  # psiu2 true variance of u
  # a and b: true parameters of inverse gamma
  # beta: intercept in LMM
  # alpha is the confindence level
  
  tau1 = alpha/2
  tau2 = 1-tau1
  
  beta.true = beta
  df<-data.frame(subject=rep(1:N,rep(n,N)),
                 time=rep(1:n,N))
  #simulate f(u)
  u<-rchisq(N, df=psiu2)
  #scaled to zero mean and variance=1
  u <- (u-psiu2)/(sqrt(2*psiu2))

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

#function for data generation - error~t(3), u~N(0,1)
SimData_eps_t<-function(N, n, psiu2, a, b, beta, alpha=0.05) {
  # N number of subjects
  # n numbr of observations in each subject, for rand. effects ui
  # psiu2 true variance of u
  # a and b: true parameters of inverse gamma
  # beta: intercept in LMM
  # alpha is the confindence level
  
  tau1 = alpha/2
  tau2 = 1-tau1
  
  beta.true = beta
  df<-data.frame(subject=rep(1:N,rep(n,N)),
                 time=rep(1:n,N))
  #simulate f(u)
  u<-rnorm(N,mean=0,sd=sqrt(psiu2))
  
  #simulate epsilon with single skewed distribution, t(3)
  epsilon <- rt(N*n, df=3)
  #scaled to zero mean and variance=1
  epsilon <- epsilon/(3/(3-2))
  
  #simulated data 
  y<- beta.true + rep(u,rep(n,N)) + epsilon
  df$y<-y
  
  return(df)
}

#function for data generation - error~N(0, theta_i), u~N(0,1) 
SimData<-function(N, n, psiu2, a, b, beta, alpha=0.05) {
  # N number of subjects
  # n numbr of observations in each subject, for rand. effects ui
  # psiu2 true variance of u
  # a and b: true parameters of inverse gamma
  # beta: intercept in LMM
  # alpha is the confindence level
  
  tau1 = alpha/2
  tau2 = 1-tau1
  
  beta.true = beta
  df<-data.frame(subject=rep(1:N,rep(n,N)),
                 time=rep(1:n,N))
  #simulate f(u)
  u<-rnorm(N,mean=0,sd=sqrt(psiu2))
  
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

#checker function
check<-function(u,tau) {
  u*(tau-as.numeric(u<=0))
}

jqm.update<-function(db,N,beta0,beta1,beta2,u,z,lambda.u,lambda.z,alpha) {
  tau1 = alpha/2
  tau2 = 1-tau1
  
  subjects<-unique(db$subject)
  # estimation of z
  n<-length(z)
  # estimation of z
  for(i in 1:N) {
    y<-db$y[db$subject==subjects[i]]-beta0-u[i]
    db$y2[db$subject==subjects[i]]<-db$y[db$subject==subjects[i]]-beta0-u[i]
    #z.approx<-mean(c(quantile(y,probs=0.025)/beta1,quantile(y,probs=0.975)/beta2))
    #z.seq<-seq(z.approx*ifelse(z.approx<1,0.5,2),1, length.out=20)
    z.seq<-seq(0.1,10,length.out = 20)
    obj<-sapply(z.seq,function(x) {
      sum(check(y-x*beta1,tau=tau1)+check(y-x*beta2,tau=tau2))+lambda.z*(x-1)^2
    })
    z.i<-z.seq[which.min(obj)[1]]
    z.seq<-seq(z.i-(z.seq[2]-z.seq[1]),z.i+(z.seq[2]-z.seq[1]),length.out = 20)
    obj<-sapply(z.seq,function(x) {
      sum(check(y-x*beta1,tau=tau1)+check(y-x*beta2,tau=tau2))+lambda.z*(x-1)^2
    })
    z.i<-z.seq[which.min(obj)[1]]
    z.seq<-runif(20,min=z.i-(z.seq[2]-z.seq[1]),max=z.i+(z.seq[2]-z.seq[1]))
    obj<-sapply(z.seq,function(x) {
      sum(check(y-x*beta1,tau=tau1)+check(y-x*beta2,tau=tau2))+lambda.z*(x-1)^2
    })
    z[i]<-z.seq[which.min(obj)[1]]
    db$z[db$subject==subjects[i]]<-z[i]
  }
  
  # scaling
  #s.seq<-seq(0.2,5,length.out = 20)
  #obj<-sapply(s.seq,function(x) {
  #  sum(check(db$y2-db$z*beta1,tau=0.025)+check(y-db$z*beta2,tau=0.975))+lambda.z*sum((x*db$z-1)^2)
  #})
  #z<-z*s.seq[which.min(obj)[1]]
  #db$z<-db$z*s.seq[which.min(obj)[1]]
  g<-(sum(z))/(sum(z^2))
  z<-z*g
  db$z<-db$z*g
  
  
  # estimation of beta1 and beta2
  m1<-rq(y2~z-1, data=db, tau=alpha/2)
  m2<-rq(y2~z-1, data=db, tau=1-alpha/2)
  beta1<-coef(m1)
  beta2<-coef(m2)
  
  # estimation of u
  for(i in 1:N) {
    y<-db$y[db$subject==subjects[i]]-beta0
    z.i<-z[i]
    
    u.seq<-seq(1.5*median(y),sign(-u[i])*sd(y),length.out = 100)
    obj<-sapply(u.seq,function(x) {
      sum(check(y-z.i*beta1-x,tau=tau1)+check(y-z.i*beta2-x,tau=tau2))+lambda.u*x^2
    })
    u.i<-u.seq[which.min(obj)[1]]
    u.seq<-runif(20,min=u.i-abs(u.seq[2]-u.seq[1]),max=u.i+abs(u.seq[2]-u.seq[1]))
    obj<-sapply(u.seq,function(x) {
      sum(check(y-z.i*beta1-x,tau=tau1)+check(y-z.i*beta2-x,tau=tau2))+lambda.u*x^2
    })
    u[i]<-u.seq[which.min(obj)[1]]
  }
  
  return(list(beta1=beta1,
              beta2=beta2,
              z=z,
              u=u))
}

jqm.est.fixed<-function(db,lambda.u,lambda.z,alpha) {
  tau1 = alpha/2
  tau2 = 1-tau1
  
  subjects<-unique(db$subject)
  N<-length(subjects)
  beta0<-median(db$y)
  u<-numeric(N)
  
  # estimation of u
  db$y2<-db$y-beta0
  w<-seq(min(db$y2),max(db$y2), length.out = 2*N)
  cnt<-1
  for(s in subjects) {
    y<-db$y2[db$subject==s]
    obj<-sapply(w, function(x) {
      (sum(sign(y-x))-4*lambda.u*x)^2})
    u[cnt]<-w[which.min(obj)[1]]
    #u[cnt]<-quantile(db$y2[db$subject==s],probs=0.5+lambda.u)
    db$y2[db$subject==s]<-db$y2[db$subject==s]-u[cnt]
    cnt<-cnt+1
  }
  
  # initial estimate of beta1 and beta2
  beta1<-quantile(db$y2,probs=alpha/2)
  beta2<-quantile(db$y2,probs=1-alpha/2)
  db$y3<-db$y2
  
  z<-rep(1,N)
  d0<-mean(z*(beta2-beta1))
  d<-10
  cnt<-1
  while((abs(log(d,base=10))>0.01)&(cnt<10)) {
    update<-jqm.update(db=db,N=N,beta0=beta0,
                       beta1=beta1,beta2=beta2,u=u,z=z,
                       lambda.u=lambda.u,lambda.z=lambda.z,alpha=alpha)
    
    
    beta1<-update$beta1
    beta2<-update$beta2
    z<-update$z
    u<-update$u
    d<-abs(mean(z*(beta2-beta1))/d0)
    d0<-mean(z*(beta2-beta1))
    
    for(i in 1:N) {
      db$z[db$subject==subjects[i]]<-z[i]
    }
    
    cnt<-cnt+1
  }
  
  return(list(beta0=beta0,
              beta1=beta1,
              beta2=beta2,
              lambda.u=lambda.u,
              lambda.z=lambda.z,
              z=z,
              u=u,
              cnt=cnt))
}

#compute coverage for IRI with both lower and upper bound
jqm.coverage.new.subject<-function(res,y,l.u,l.z,alpha) {
  l.u<-res$lambda.u
  l.z<-res$lambda.z
  tau1 = alpha/2
  tau2 = 1-tau1
  
  
  y.i<-y-res$beta0
  w<-seq(min(y.i),max(y.i), length.out = 100)
  obj<-sapply(w,function(x) {
    sum(check(y.i-res$beta1-x,tau=tau1)+check(y.i-res$beta2-x,tau=tau2))+l.u*x^2
  })
  u.i<-w[which.min(obj)[1]]
  w<-runif(20,min=u.i-(w[2]-w[1]), max=u.i+(w[2]-w[1]))
  obj<-sapply(w,function(x) {
    sum(check(y.i-res$beta1-x,tau=tau1)+check(y.i-res$beta2-x,tau=tau2))+l.u*x^2
  })
  u.i<-w[which.min(obj)[1]]
  
  y.i<-y.i-u.i
  z.seq<-seq(0.9*min(res$z),1.1*max(res$z), length.out=100)
  obj<-sapply(z.seq,function(x) {
    sum(check(y.i-x*res$beta1,tau=tau1)+check(y.i-x*res$beta2,tau=tau2))+l.z*(x-1)^2
  })
  z.i<-z.seq[which.min(obj)[1]]
  
  y.i<-y-res$beta0
  w<-seq(min(y.i),max(y.i), length.out = 100)
  obj<-sapply(w,function(x) {
    sum(check(y.i-z.i*res$beta1-x,tau=tau1)+check(y.i-z.i*res$beta2-x,tau=tau2))+l.u*x^2
  })
  u.i<-w[which.min(obj)[1]]
  w<-runif(20,min=u.i-(w[2]-w[1]), max=u.i+(w[2]-w[1]))
  obj<-sapply(w,function(x) {
    sum(check(y.i-z.i*res$beta1-x,tau=tau1)+check(y.i-z.i*res$beta2-x,tau=tau2))+l.u*x^2
  })
  u.i<-w[which.min(obj)[1]]
  
  y.i<-y.i-u.i
  z.seq<-seq(0.9*min(res$z),1.1*max(res$z), length.out=100)
  obj<-sapply(z.seq,function(x) {
    sum(check(y.i-x*res$beta1,tau=tau1)+check(y.i-x*res$beta2,tau=tau2))+l.z*(x-1)^2
  })
  z.i<-z.seq[which.min(obj)[1]]
  z.seq<-runif(20,min=z.i-(z.seq[2]-z.seq[1]), max=z.i+(z.seq[2]-z.seq[1]))
  obj<-sapply(z.seq,function(x) {
    sum(check(y.i-x*res$beta1,tau=tau1)+check(y.i-x*res$beta2,tau=tau2))+l.z*(x-1)^2
  })
  z.i<-z.seq[which.min(obj)[1]]
  
  coverage2<-mean((y.i>z.i*res$beta1)&
                    (y.i<z.i*res$beta2))
  
  return(coverage2)
}

#compute coverage for IRI with only upper bound
jqm.coverage.new.subject.up<-function(res,y,l.u,l.z,alpha) {
  l.u<-res$lambda.u
  l.z<-res$lambda.z
  tau1 = alpha/2
  tau2 = 1-tau1
  
  
  y.i<-y-res$beta0
  w<-seq(min(y.i),max(y.i), length.out = 100)
  obj<-sapply(w,function(x) {
    sum(check(y.i-res$beta1-x,tau=tau1)+check(y.i-res$beta2-x,tau=tau2))+l.u*x^2
  })
  u.i<-w[which.min(obj)[1]]
  w<-runif(20,min=u.i-(w[2]-w[1]), max=u.i+(w[2]-w[1]))
  obj<-sapply(w,function(x) {
    sum(check(y.i-res$beta1-x,tau=tau1)+check(y.i-res$beta2-x,tau=tau2))+l.u*x^2
  })
  u.i<-w[which.min(obj)[1]]
  
  y.i<-y.i-u.i
  z.seq<-seq(0.9*min(res$z),1.1*max(res$z), length.out=100)
  obj<-sapply(z.seq,function(x) {
    sum(check(y.i-x*res$beta1,tau=tau1)+check(y.i-x*res$beta2,tau=tau2))+l.z*(x-1)^2
  })
  z.i<-z.seq[which.min(obj)[1]]
  
  y.i<-y-res$beta0
  w<-seq(min(y.i),max(y.i), length.out = 100)
  obj<-sapply(w,function(x) {
    sum(check(y.i-z.i*res$beta1-x,tau=tau1)+check(y.i-z.i*res$beta2-x,tau=tau2))+l.u*x^2
  })
  u.i<-w[which.min(obj)[1]]
  w<-runif(20,min=u.i-(w[2]-w[1]), max=u.i+(w[2]-w[1]))
  obj<-sapply(w,function(x) {
    sum(check(y.i-z.i*res$beta1-x,tau=tau1)+check(y.i-z.i*res$beta2-x,tau=tau2))+l.u*x^2
  })
  u.i<-w[which.min(obj)[1]]
  
  y.i<-y.i-u.i
  z.seq<-seq(0.9*min(res$z),1.1*max(res$z), length.out=100)
  obj<-sapply(z.seq,function(x) {
    sum(check(y.i-x*res$beta1,tau=tau1)+check(y.i-x*res$beta2,tau=tau2))+l.z*(x-1)^2
  })
  z.i<-z.seq[which.min(obj)[1]]
  z.seq<-runif(20,min=z.i-(z.seq[2]-z.seq[1]), max=z.i+(z.seq[2]-z.seq[1]))
  obj<-sapply(z.seq,function(x) {
    sum(check(y.i-x*res$beta1,tau=tau1)+check(y.i-x*res$beta2,tau=tau2))+l.z*(x-1)^2
  })
  z.i<-z.seq[which.min(obj)[1]]
  
  coverage2<-mean(y.i<z.i*res$beta2)
  
  return(coverage2)
}

#compute coverage for IRI with only lower bound
jqm.coverage.new.subject.low<-function(res,y,l.u,l.z,alpha) {
  l.u<-res$lambda.u
  l.z<-res$lambda.z
  tau1 = alpha/2
  tau2 = 1-tau1
  
  
  y.i<-y-res$beta0
  w<-seq(min(y.i),max(y.i), length.out = 100)
  obj<-sapply(w,function(x) {
    sum(check(y.i-res$beta1-x,tau=tau1)+check(y.i-res$beta2-x,tau=tau2))+l.u*x^2
  })
  u.i<-w[which.min(obj)[1]]
  w<-runif(20,min=u.i-(w[2]-w[1]), max=u.i+(w[2]-w[1]))
  obj<-sapply(w,function(x) {
    sum(check(y.i-res$beta1-x,tau=tau1)+check(y.i-res$beta2-x,tau=tau2))+l.u*x^2
  })
  u.i<-w[which.min(obj)[1]]
  
  y.i<-y.i-u.i
  z.seq<-seq(0.9*min(res$z),1.1*max(res$z), length.out=100)
  obj<-sapply(z.seq,function(x) {
    sum(check(y.i-x*res$beta1,tau=tau1)+check(y.i-x*res$beta2,tau=tau2))+l.z*(x-1)^2
  })
  z.i<-z.seq[which.min(obj)[1]]
  
  y.i<-y-res$beta0
  w<-seq(min(y.i),max(y.i), length.out = 100)
  obj<-sapply(w,function(x) {
    sum(check(y.i-z.i*res$beta1-x,tau=tau1)+check(y.i-z.i*res$beta2-x,tau=tau2))+l.u*x^2
  })
  u.i<-w[which.min(obj)[1]]
  w<-runif(20,min=u.i-(w[2]-w[1]), max=u.i+(w[2]-w[1]))
  obj<-sapply(w,function(x) {
    sum(check(y.i-z.i*res$beta1-x,tau=tau1)+check(y.i-z.i*res$beta2-x,tau=tau2))+l.u*x^2
  })
  u.i<-w[which.min(obj)[1]]
  
  y.i<-y.i-u.i
  z.seq<-seq(0.9*min(res$z),1.1*max(res$z), length.out=100)
  obj<-sapply(z.seq,function(x) {
    sum(check(y.i-x*res$beta1,tau=tau1)+check(y.i-x*res$beta2,tau=tau2))+l.z*(x-1)^2
  })
  z.i<-z.seq[which.min(obj)[1]]
  z.seq<-runif(20,min=z.i-(z.seq[2]-z.seq[1]), max=z.i+(z.seq[2]-z.seq[1]))
  obj<-sapply(z.seq,function(x) {
    sum(check(y.i-x*res$beta1,tau=tau1)+check(y.i-x*res$beta2,tau=tau2))+l.z*(x-1)^2
  })
  z.i<-z.seq[which.min(obj)[1]]
  
  coverage2<-mean((y.i>z.i*res$beta1))
  
  return(coverage2)
}

#function for running the joint-IRI estimation, it includes calling the jqm.update, 
#jqm.est.fixed, and jqm.coverage.new.subject functions
jqm<-function(db, alpha=0.05,lambda.u.seq=seq(0.5,4,0.5),
              lambda.z.seq=seq(0.5,5,0.5)) {
    beta0<-median(db$y)
    subjects<-unique(db$subject)
    N<-length(subjects)
    n<-length(unique(db$time)) # assuming n is constant over all subjects
    obj.tot<-0 # objective function to be minized
    obj.check<-0
    obj.pen<-0
    
    u<-numeric(N)
    z<-rep(1,N)
    db$z<-1
    
    w<-seq(min(db$y),max(db$y), length.out = 2*N)-median(db$y)
    w2<-seq(0.1,10,length.out = 100)
    
    db.del<-db[db$time!=max(db$time),]
    all.coverages<-c()
    cv.results<-data.frame(lambda.u=0,lambda.z=0,cov.time=0,cov.subj=0,coverage=0)
    for(l.z in lambda.z.seq) {
      above<-0
      below<-0
      for(l.u in lambda.u.seq) {
        coverage.tot<-c(0,0)
        for(i in 1:N) {
          db.del<-db[db$subject!=subjects[i],]
          db.del<-db.del[db.del$time!=max(db.del$time),]
          res<-jqm.est.fixed(db=db.del,lambda.u=l.u,lambda.z=l.z, alpha=alpha)
          
          coverage<-numeric(N-1)
          cnt<-1
          for(s in subjects[subjects!=i]) {
            y<-db$y[(db$subject==s)&(db$time==max(db$time))]
            coverage[cnt]<-((y>res$beta0+res$u[cnt]+res$z[cnt]*res$beta1)&
                              (y<res$beta0+res$u[cnt]+res$z[cnt]*res$beta2))
            cnt<-cnt+1
          }
          coverage<-mean(coverage)
          
          coverage2<-jqm.coverage.new.subject(res=res,y=db$y[db$subject==i],alpha=alpha)
          
          coverage.tot<-coverage.tot+c(coverage,coverage2)
          coverage.tot.both<-sum(coverage.tot)/(2*N)
        }
        cv.results<-rbind(cv.results,c(l.u,l.z,coverage.tot/N,coverage.tot.both))
        gc(verbose = FALSE)
        if(coverage.tot.both>0.95) {
          above<-above+1
          if((below>0)|(above==4)) {
            below<--1
          }
        }
        if(coverage.tot.both<0.95) {
          below<-below+1
          if((above>0)|(below==4)) {
            above<--1
          }
        }
        if(sign(above*below)==-1) {
          above<-0
          below<-0
          break
        }
      }
    }
    cv.results<-cv.results[-1,]
    optimum<-cv.results[which.min((cv.results$coverage-(1-alpha))^2)[1],]
    res<-jqm.est.fixed(db=db,alpha=alpha,lambda.u=optimum$lambda.u,lambda.z=optimum$lambda.z)
    gc(verbose = FALSE)
  
   plot(cv.results$lambda.u,cv.results$coverage)
   plot(cv.results$lambda.z,cv.results$coverage)
   interaction.plot(response=cv.results$coverage,
                    x.factor=cv.results$lambda.z,trace.factor = cv.results$lambda.u)
   abline(h=0.95,col=2,lty=2)
   interaction.plot(response=cv.results$cov.time,
                    x.factor=cv.results$lambda.z,trace.factor = cv.results$lambda.u)
   abline(h=0.95,col=2,lty=2)
   interaction.plot(response=cv.results$cov.subj,
                    x.factor=cv.results$lambda.z,trace.factor = cv.results$lambda.u)
   abline(h=0.95,col=2,lty=2)


  return(list(beta0=res$beta0,
              beta1=res$beta1,
              beta2=res$beta2,
              u=res$u,
              z=res$z,
              lambda.u=optimum$lambda.u,
              lambda.z=optimum$lambda.z,
              cov.subj=optimum$cov.subj,
              cov.time=optimum$cov.time,
              cov.tot=optimum$coverage))
}

#function for running the simulation
SimStudy<-function(N, n, beta, a, b, psiu2, alpha, NSim, seed) {
  Results<-matrix(nrow=NSim,ncol=14)
  Results<-as.data.frame(Results)
  names(Results)<-c("iter","alpha","beta0","beta1","beta2",
                    "lambda.u","lambda.z",
                    "Cov.subj","Cov.time","Cov",
                    "EmpCov.subj","EmpCov.time","EmpCov", "seed")

  for(si in 1:NSim) {
    #modify SimData function according to the target scenario
    df<-SimData(N=N+1,n=n+1,beta=beta,a=a,b=b,psiu2=psiu2,
              alpha=alpha)
    df.fit<-df[(df$time<=n)&(df$subject<=N),]
    
    res<-jqm(db=df.fit,
           alpha=alpha,
           lambda.u.seq = seq(0.5,4,0.5),
           lambda.z.seq = seq(0.5,5,0.5))
  
    Results[si,1:10]<-c(si,alpha,res$beta0,res$beta1,res$beta2,
                         res$lambda.u,res$lambda.z,
                         res$cov.subj,res$cov.time,res$cov.tot)
  
    #calculate the empirical coverage for a new subject and a new measurement
    coverage<-numeric(N)
    for(s in 1:N) {
      y<-df$y[(df$subject==s)&(df$time==(n+1))]
      coverage[s]<-((y>res$beta0+res$u[s]+res$z[s]*res$beta1)&
                        (y<res$beta0+res$u[s]+res$z[s]*res$beta2))
    }
    coverage<-mean(coverage)
    
    y.i<-df$y[df$subject==N+1]
    coverage2<-jqm.coverage.new.subject(res=res,y=y.i,alpha=alpha)
    
    Results[si,11:13]<-c(coverage2,coverage,mean(c(coverage,coverage2)))
    Results[si,14]<-seed
      
    #save results for every 10 iterations
    if((si%%10) ==0) {
     save(Results, file=paste("./output/PJQM/JQMTmp_", si, "_alpha_",alpha, "_N_",N, "_n_", n, "_beta_", 
                               beta, "_a_", a, "_b_", b, "_Psiu2_", psiu2, 
                               ".Rdata",sep = ""))
    }
   gc(verbose = FALSE)
  }
}