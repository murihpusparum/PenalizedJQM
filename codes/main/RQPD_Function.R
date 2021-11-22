#install.packages("rqpd", repos="http://R-Forge.R-project.org")
library(rqpd)

#checker function
check<-function(u,tau) {
  u*(tau-as.numeric(u<=0))
}

#compute coverage for IRI with both lower and upper bound
rqpd.coverage.new.subject<-function(res,y,y.s,alpha,side="both") {
  tau1 = alpha/2
  tau2 = 1-tau1
  
  l<-res$lambda
  beta01<-res$beta01
  beta02<-res$beta02

  #compute alpha for the new subject
  y.i<-y
  w<-seq(min(y.i),max(y.i), length.out = 100)
  obj<-sapply(w,function(x) {
    sum(0.5*check(y.i-beta01-x,tau=tau1)+0.5*check(y.i-beta02-x,tau=tau2))+(l*abs(x))
  })
  alpha.i<-w[which.min(obj)[1]]
  w<-runif(20,min=alpha.i-(w[2]-w[1]), max=alpha.i+(w[2]-w[1]))
  obj<-sapply(w,function(x) {
    sum(0.5*check(y.i-beta01-x,tau=tau1)+0.5*check(y.i-beta02-x,tau=tau2))+l*abs(x)
  })
  alpha.i<-w[which.min(obj)[1]]
  
  if(side=="both"){
    coverage2<-mean((y.s>beta01+alpha.i)&
                      (y.s<beta02+alpha.i))    
  }
  if(side=="up"){
    coverage2<-mean(y.s<beta02+alpha.i)
  }
  if(side=="low"){
    coverage2<-mean(y.s>beta01+alpha.i)
  }

  return(coverage2)
}


#compute coverage for IRI with both lower and upper bound
rqpd.coverage.new.subject1obs<-function(res,y,y.s,alpha,side="both") {
  tau1 = alpha/2
  tau2 = 1-tau1
  
  l<-res$lambda
  beta01<-res$beta01
  beta02<-res$beta02
  
  #compute alpha for the new subject
  y.i<-y
  w<-seq(min(y.i-l),max(y.i+l), length.out = 100)
  obj<-sapply(w,function(x) {
    sum(0.5*check(y.i-beta01-x,tau=tau1)+0.5*check(y.i-beta02-x,tau=tau2))+(l*abs(x))
  })
  alpha.i<-w[which.min(obj)[1]]
  w<-runif(20,min=alpha.i-(w[2]-w[1]), max=alpha.i+(w[2]-w[1]))
  obj<-sapply(w,function(x) {
    sum(0.5*check(y.i-beta01-x,tau=tau1)+0.5*check(y.i-beta02-x,tau=tau2))+l*abs(x)
  })
  alpha.i<-w[which.min(obj)[1]]

  if(side=="both"){
    coverage2<-mean((y.s>beta01+alpha.i)&
                      (y.s<beta02+alpha.i))    
  }
  if(side=="up"){
    coverage2<-mean(y.s<beta02+alpha.i)
  }
  if(side=="low"){
    coverage2<-mean(y.s>beta01+alpha.i)
  }

  return(coverage2)
}

rqpd.fun<-function(db, alpha=0.05,lambda=seq(0.1,5,0.5),side="both") {
  subjects<-unique(db$subject)
  N<-length(subjects)
  n<-length(unique(db$time)) # assuming n is constant over all subjects
  obj.tot<-0 # objective function to be minimized
  obj.check<-0
  obj.pen<-0

  tau1 = alpha/2
  tau2 = 1-tau1
  
  # fit the model in all possible lambda values
  cv.results<-data.frame(lambda=0,cov.time=0,cov.subj=0,coverage=0)
  for(l in lambda) {
    above<-0
    below<-0
    
    # for in-sample, use all subjects except last measurements
    db.del<-db[db$time!=max(db$time),]
    db.del$subject<-as.factor(db.del$subject)
    
    res<-rqpd(y ~ 1 | subject, panel(lambda = l, taus=c(tau1, tau2), tauw=rep(1/2, 2)), data=db.del)
    beta01<-res$coef[1]
    beta02<-res$coef[2]
    alpha.i<-res$coef[-c(1,2)]
    
    coverage<-NA
    cnt<-1
    for(s in subjects) {
      db.y<-db[(db$subject==s),]
      y<-db.y$y[(db.y$time==max(db.y$time))]
      
      if(side=="both"){
        coverage[cnt]<-((y>beta01+alpha.i[cnt])&
                          (y<beta02+alpha.i[cnt]))   
      }
      if(side=="up"){
        coverage[cnt]<-(y<beta02+alpha.i[cnt])
      }
      if(side=="low"){
        coverage[cnt]<-(y>beta01+alpha.i[cnt])
      }

      cnt<-cnt+1
    }
    coverage<-mean(coverage)
    
    # for out-sample, use N-1 subjects and exclude last measurements
    # in leave-one-out cross-validation
    coverage2<-numeric(N)
    cnt<-1
    for(s in subjects) {
      db.del2<-db[db$subject!=s,]
      db.del2<-db.del2[db.del2$time!=max(db.del2$time),]
      db.del2$subject<-as.factor(db.del2$subject)
      
      res2<-rqpd(y ~ 1 | subject, panel(lambda = l, taus=c(tau1, tau2), tauw=rep(1/2, 2)), data=db.del2)
      res2<-list(beta01=res2$coef[[1]],
                 beta02=res2$coef[[2]],
                 lambda=l)
      # coverage for a new subject minus the last time point
      db.y<-db[db$subject==s,]
      y.i<-db.y$y[db.y$time!=max(db.y$time)]
      y.s<-db.y$y[db.y$time==max(db.y$time)]
      
      coverage2[cnt]<-rqpd.coverage.new.subject(res=res2,y=y.i,y.s=y.s,alpha=alpha,side=side)
      cnt<-cnt+1
    }
    coverage2<-mean(coverage2)
    coverage.tot<-c(coverage,coverage2)
    coverage.tot.both<-mean(coverage.tot)
    
    cv.results<-rbind(cv.results,c(l,coverage.tot,coverage.tot.both))
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
  cv.results<-cv.results[-1,]
  for (i in 1:nrow(cv.results)) {
    cv.results$min[i]<-min(cv.results$cov.time[i],cv.results$cov.subj[i])
  }
  optimum<-cv.results[which.min((cv.results$min-(1-alpha))^2)[1],]
#  optimum<-cv.results[which.min((cv.results$coverage-(1-alpha))^2)[1],]
  db$subject<-as.factor(db$subject)
  res<-rqpd(y ~ 1 | subject, panel(lambda = optimum$lambda, taus=c(tau1, tau2), tauw=rep(1/2, 2)), data=db)
  gc(verbose = FALSE)
  
  plot(cv.results$lambda,cv.results$coverage)

  return(list(beta01=res$coef[1],
              beta02=res$coef[2],
              alpha.i=res$coef[-c(1,2)],
              lambda=optimum$lambda,
              cov.subj=optimum$cov.subj,
              cov.time=optimum$cov.time,
              cov.tot=optimum$coverage))
}