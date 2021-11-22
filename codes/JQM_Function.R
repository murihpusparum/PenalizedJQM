library(quantreg)
library(invgamma)

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

jqm.coverage.new.subject<-function(res,y,y.s,alpha,side="both") {
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
  
  if(side=="both"){
    coverage2<-mean((y.s>res$beta0+u.i+z.i*res$beta1)&
                      (y.s<res$beta0+u.i+z.i*res$beta2))    
  }
  if(side=="up"){
    coverage2<-mean(y.s<res$beta0+u.i+z.i*res$beta2)
  }
  if(side=="low"){
    coverage2<-mean(y.s>res$beta0+u.i+z.i*res$beta1)
  }
  
  return(coverage2)
}

jqm.coverage.new.subject1obs<-function(res,y,y.s,alpha,side="both") {
  l.u<-res$lambda.u
  l.z<-res$lambda.z
  tau1 = alpha/2
  tau2 = 1-tau1
  
  
  y.i<-y-res$beta0
  w<-seq((min(y.i)-l.u-l.z),(max(y.i)+l.u+l.z), length.out = 100)
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
  w<-seq((min(y.i)-l.u-l.z),(max(y.i)+l.u+l.z), length.out = 100)
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
  
  if(side=="both"){
    coverage2<-mean((y.s>res$beta0+u.i+z.i*res$beta1)&
                      (y.s<res$beta0+u.i+z.i*res$beta2))    
  }
  if(side=="up"){
    coverage2<-mean(y.s<res$beta0+u.i+z.i*res$beta2)
  }
  if(side=="low"){
    coverage2<-mean(y.s>res$beta0+u.i+z.i*res$beta1)
  }
  
  return(coverage2)
}

#function for running the joint-IRI estimation, it includes calling the jqm.update, 
#jqm.est.fixed, and jqm.coverage.new.subject functions
jqm<-function(db, alpha=0.05,lambda.u.seq=seq(0.5,4,0.5),
              lambda.z.seq=seq(0.5,5,0.5),side="both") {
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
  
  # if time points are imbalance, max(time) == max(time) - 1
  all.coverages<-c()
  cv.results<-data.frame(lambda.u=0,lambda.z=0,cov.time=0,cov.subj=0,coverage=0)
  for(l.z in lambda.z.seq) {
    above<-0
    below<-0
    for(l.u in lambda.u.seq) {
      # if time points are imbalance, exclude the last time point that is present in all subject
      max.time.del <- ifelse(length(db[db$time==max(db$time),]$subject)==length(unique(db$subject)), 
                             yes=max(db$time), no=(max(db$time)-1))
      db.del<-db[db$time!=max.time.del,]
      res<-jqm.est.fixed(db=db.del,lambda.u=l.u,lambda.z=l.z, alpha=alpha)
      
      coverage<-numeric(N)
      cnt<-1
      for(s in subjects) {
        db.y<-db[(db$subject==s),]
        y<-db.y$y[(db.y$time==max(db.y$time))]
        
        if(side=="both"){
          coverage[cnt]<-((y>res$beta0+res$u[cnt]+res$z[cnt]*res$beta1)&
                            (y<res$beta0+res$u[cnt]+res$z[cnt]*res$beta2))    
        }
        if(side=="up"){
          coverage[cnt]<-(y<res$beta0+res$u[cnt]+res$z[cnt]*res$beta2)
        }
        if(side=="low"){
          coverage[cnt]<-(y>res$beta0+res$u[cnt]+res$z[cnt]*res$beta1)
        }
        
        cnt<-cnt+1
      }
      coverage<-mean(coverage)
      
      # # #estimate the parameters again, without the last subject
      # db.del2<-db[db$subject!=max(db$subject),]
      # N.del<-length(unique(db.del2$subject))
      # max.time<-ifelse(length(db.del2[db.del2$time==max(db.del2$time),]$subject)==N.del, 
      #                  yes=max(db.del2$time), no=(max(db.del2$time)-1))
      # db.del2<-db.del2[db.del2$time!=max.time,]
      # 
      # res2<-jqm.est.fixed(db=db.del2,lambda.u=l.u,lambda.z=l.z, alpha=alpha)
      # # coverage for a new subject minus the last time point
      # db.y<-db[db$subject==s,]
      # y.i<-db.y$y[db.y$time!=max(db.y$time)]
      # y.s<-db.y$y[db.y$time==max(db.y$time)]
      # 
      # coverage2<-jqm.coverage.new.subject(res=res2,y=y.i,y.s=y.s,alpha=alpha)
      # 
      # coverage for the new subject, using LOOCV
      coverage.tot<-c(0,0)
      coverage2<-numeric(N)
      cnt<-1
      for(s in subjects) {
        db.del2<-db[db$subject!=s,]
        N.del<-length(unique(db.del2$subject))
        max.time<-ifelse(length(db.del2[db.del2$time==max(db.del2$time),]$subject)==N.del,
                         yes=max(db.del2$time), no=(max(db.del2$time)-1))
        db.del2<-db.del2[db.del2$time!=max.time,]
        
        res2<-jqm.est.fixed(db=db.del2,lambda.u=l.u,lambda.z=l.z, alpha=alpha)
        # coverage for a new subject minus the last time point
        db.y<-db[db$subject==s,]
        y.i<-db.y$y[db.y$time!=max(db.y$time)]
        y.s<-db.y$y[db.y$time==max(db.y$time)]
        
        coverage2[cnt]<-jqm.coverage.new.subject(res=res2,y=y.i,y.s=y.s,alpha=alpha,side=side)
        cnt<-cnt+1
      }
      coverage2<-mean(coverage2)
      coverage.tot<-c(coverage,coverage2)
      coverage.tot.both<-mean(coverage.tot)
      
      cv.results<-rbind(cv.results,c(l.u,l.z,coverage.tot,coverage.tot.both))
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
  for (i in 1:nrow(cv.results)) {
    cv.results$min[i]<-min(cv.results$cov.time[i],cv.results$cov.subj[i])
  }
  optimum<-cv.results[which.min((cv.results$min-(1-alpha))^2)[1],]
  #    optimum<-cv.results[which.min((cv.results$coverage-(1-alpha))^2)[1],]
  res<-jqm.est.fixed(db=db,alpha=alpha,lambda.u=optimum$lambda.u,lambda.z=optimum$lambda.z)
  gc(verbose = FALSE)
  
  plot(cv.results$lambda.u,cv.results$coverage)
  plot(cv.results$lambda.z,cv.results$coverage)
  interaction.plot(response=cv.results$coverage,
                   x.factor=cv.results$lambda.z,trace.factor = cv.results$lambda.u)
  abline(h=1-alpha,col=2,lty=2)
  interaction.plot(response=cv.results$cov.time,
                   x.factor=cv.results$lambda.z,trace.factor = cv.results$lambda.u)
  abline(h=1-alpha,col=2,lty=2)
  interaction.plot(response=cv.results$cov.subj,
                   x.factor=cv.results$lambda.z,trace.factor = cv.results$lambda.u)
  abline(h=1-alpha,col=2,lty=2)
  
  
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
