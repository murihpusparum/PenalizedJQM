source("./codes/main/RQPD_Function.R")

##coverage for RI with both lower and upper bound
covrg<-function(df,alpha=0.05,side="both"){
  N<-length(unique(df$subject))
  n<-length(unique(df$time))
  subjects<-unique(df$subject)
  
  
  #coverage when there is one new measurement
  df.fit<-df[(df$time<n),]
  #  df.fit$subject<-rep(seq(1:(N)), each=n-1)
  res<-rqpd.fun(db=df.fit,alpha=alpha,lambda=seq(0.1,5,0.5),side=side)
  
  cv<-cbind.data.frame(res$alpha.i)
  cv$low<-res$alpha.i + res$beta01[[1]]
  cv$up<-res$alpha.i + res$beta02[[1]]
  
  y.s<- df$y[(df$time==(n))]
  
  if(side=="both"){
    coverage<-((y.s>cv$low)&(y.s<cv$up))    
  }
  if(side=="up"){coverage<-(y.s<cv$up)}
  if(side=="low"){coverage<-(y.s>cv$low)}
  
  coverage<-mean(coverage)
  
  coverage2<-numeric(N)
  coverage3<-numeric(N)
  cnt<-1
  for(s in subjects) {
    df.del2<-df[df$subject!=s,]
    N.del<-length(unique(df.del2$subject))
    max.time<-ifelse(length(df.del2[df.del2$time==max(df.del2$time),]$subject)==N.del,
                     yes=max(df.del2$time), no=(max(df.del2$time)-1))
    df.del2<-df.del2[df.del2$time!=max.time,]
    
    res2<-rqpd.fun(db=df.del2,alpha=alpha,lambda=seq(0.1,5,0.5),side=side)
    
    # coverage for a new subject minus the last time point
    df.y<-df[df$subject==s,]
    y.i<-df.y$y[df.y$time!=max(df.y$time)]
    y.s<-df.y$y[df.y$time==max(df.y$time)]
    coverage2[cnt]<-rqpd.coverage.new.subject(res=res2,y=y.i,y.s=y.s,alpha=alpha,side=side)
    
    # out-sample, coverage for one new observation in the new subject
    y.i1<-df.y$y[df.y$time==(max(df.y$time)-1)] # Hs of the new subject
    coverage3[cnt]<-rqpd.coverage.new.subject1obs(res=res2,y=y.i1,y.s=y.s,alpha=alpha,side=side)
    cnt<-cnt+1
  }
  coverage2<-mean(coverage2)
  coverage3<-mean(coverage3)
  
  co<-c(coverage2, coverage, mean(c(coverage,coverage2)), coverage3)
  return(co)
}

load(".data/iam/all_data.RData")
all<-all_clin_metabo

# read data with randomized id
load(".data/iam/rand.RData")
ran<-ran

summary_method<-data.frame(matrix(nrow=22, ncol=8))
colnames(summary_method)<-c("par", "beta1", "beta2", "lambda",
                              "cov.subj", "cov.time", "cov.tot", "cov.subj1obs")

#for RI with only upper bound (total cholesterol, LDL cholesterol), lower bound=0
for (j in 1:2) {
  N<-length(unique(all$subject))
  n<-length(unique(all$time))
  
  df.m<-all[,c(1:2)]
  df.m$y<-all[, j+2]
  
  rownames(df.m)<-NULL
  df.m$data<-"metabolomics"
  
  alpha=0.1
  met.iam<-rqpd.fun(db=df.m,alpha=alpha,lambda=seq(0.1,5,0.5),side="up")
  
  uz.m<-cbind.data.frame(met.iam$alpha.i)
  uz.m$low<-0
  uz.m$up <-met.iam$alpha.i + met.iam$beta02[[1]]
  uz.m$data<-"metabolomics"
  uz.m$id<-c(1:N)
  colnames(uz.m)<-c("alpha", "low", "up", "data","id")
  
  #calculate empirical coverage
  coverg<-covrg(df=df.m, alpha=alpha, side="up")

  summary_method[j, 1]<-paste0("met_", colnames(all)[j+2])
  summary_method[j, 2:8]<-c(met.iam$beta01[[1]], met.iam$beta02[[1]], met.iam$lambda,
                               coverg[1], coverg[2], coverg[3], coverg[4])
  
  ##CLINICAL
  df.c<-all[,c(1:2)]
  df.c$y<-all[, j+13]
  
  rownames(df.c)<-NULL
  df.c$data<-"clinical"
  
  clin.iam<-rqpd.fun(db=df.c,alpha=alpha,lambda=seq(0.1,5,0.5),side="up")
  
  uz.c<-cbind.data.frame(clin.iam$alpha.i)
  uz.c$low<-0
  uz.c$up <-clin.iam$alpha.i + clin.iam$beta02[[1]]
  uz.c$data<-"clinical"
  uz.c$id<-c(1:N)
  colnames(uz.c)<-c("alpha", "low", "up", "data", "id")
  
  #calculate empirical coverage
  coverag<-covrg(df=df.c, alpha=alpha, side="up")
  
  summary_method[j+11, 1]<-paste0("clin_", colnames(all)[j+13])
  summary_method[j+11, 2:8]<-c(clin.iam$beta01[[1]], clin.iam$beta02[[1]], clin.iam$lambda,
                              coverag[1], coverag[2], coverag[3], coverag[4])
  
  ####################
  #one combined plot
  subj<-ran[,j+1]
  
  #random subjects
  uz.a<-rbind.data.frame(uz.c, uz.m)
  uz.a$id2<-subj
  
  df.a<-rbind.data.frame(df.c, df.m)
  df.a$id2<-rep(subj, each=n)
  
  save(uz.a, file=paste0("./output/iam/", colnames(all)[j+2], ".uz5.rqpd.RData"))
  save(df.a, file=paste0("./output/iam/", colnames(all)[j+2], ".df5.rqpd.RData"))
}


#for RI with both lower and upper bound
for (j in 3:11) {
  df.m<-all[,c(1:2)]
  df.m$y<-all[, j+2]
  
  rownames(df.m)<-NULL
  df.m$data<-"metabolomics"
  
  alpha=0.05
  met.iam<-rqpd.fun(db=df.m,alpha=alpha,lambda=seq(0.1,5,0.5),side="both")
  
  uz.m<-cbind.data.frame(met.iam$alpha.i)
  uz.m$low<-met.iam$alpha.i + met.iam$beta01[[1]]
  uz.m$up <-met.iam$alpha.i + met.iam$beta02[[1]]
  uz.m$data<-"metabolomics"
  uz.m$id<-c(1:N)
  colnames(uz.m)<-c("alpha", "low", "up", "data","id")
  
  #calculate empirical coverage
  coverg<-covrg(df=df.m, alpha=alpha, side="both")
  
  summary_method[j, 1]<-paste0("met_", colnames(all)[j+2])
  summary_method[j, 2:8]<-c(met.iam$beta01[[1]], met.iam$beta02[[1]], met.iam$lambda,
                              coverg[1], coverg[2], coverg[3], coverg[4])
  
  ##CLINICAL
  df.c<-all[,c(1:2)]
  df.c$y<-all[, j+13]
  
  rownames(df.c)<-NULL
  df.c$data<-"clinical"
  
  clin.iam<-rqpd.fun(db=df.c,alpha=alpha,lambda=seq(0.1,5,0.5),side="both")
  
  uz.c<-cbind.data.frame(clin.iam$alpha.i)
  uz.c$low<-clin.iam$alpha.i + clin.iam$beta01[[1]]
  uz.c$up <-clin.iam$alpha.i + clin.iam$beta02[[1]]
  uz.c$data<-"clinical"
  uz.c$id<-c(1:N)
  colnames(uz.c)<-c("alpha", "low", "up", "data", "id")
  
  #calculate empirical coverage
  coverag<-covrg(df=df.c, alpha=alpha, side="both")
  
  summary_method[j+11, 1]<-paste0("clin_", colnames(all)[j+13])
  summary_method[j+11, 2:8]<-c(clin.iam$beta01[[1]], clin.iam$beta02[[1]], clin.iam$lambda,
                              coverag[1], coverag[2], coverag[3], coverag[4])
  
  ####################
  #one combined plot
  subj<-ran[,j+1]
  
  #random subjects
  uz.a<-rbind.data.frame(uz.c, uz.m)
  uz.a$id2<-subj
  
  df.a<-rbind.data.frame(df.c, df.m)
  df.a$id2<-rep(subj, each=n)
  
  save(uz.a, file=paste0("./output/iam/", colnames(all)[j+2], ".uz5.rqpd.RData"))
  save(df.a, file=paste0("./output/iam/", colnames(all)[j+2], ".df5.rqpd.RData"))
}


save(summary_method, file="./output/iam/sum.rqpd.RData")