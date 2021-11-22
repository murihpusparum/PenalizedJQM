source("./codes/main/LQMM_Function.R")

##coverage for RI with both lower and upper bound
covrg<-function(df,alpha=0.05,side="both"){
  N<-length(unique(df$subject))
  n<-length(unique(df$time))
  subjects<-unique(df$subject)
  
  
  #coverage when there is one new measurement
  df.fit<-df[(df$time<n),]
  #  df.fit$subject<-rep(seq(1:(N)), each=n-1)
  res<-lqmm.fit(df.fit,alpha=alpha)
  
  cv<-cbind.data.frame(res$u.i1,res$u.i2)
  cv$low<-res$beta01 + res$u.i1
  cv$up<-res$beta02 + res$u.i2
  
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
    
    res2<-lqmm.fit(df.del2,alpha=alpha)
    
    # coverage for a new subject minus the last time point
    df.y<-df[df$subject==s,]
    df.i<-df.y[df.y$time!=max(df.y$time),]
    
    y.i<-df.y$y[df.y$time!=max(df.y$time)]
    y.s<-df.y$y[df.y$time==max(df.y$time)]
    
    coverage2[cnt]<-lqmm.cov.new.subject(res2,df.i,y.i,y.s,alpha,side)
    
    # out-sample, coverage for one new observation in the new subject
    df.i1<-df.y[df.y$time==(max(df.y$time)-1),]
    y.i1<-df.y$y[df.y$time==(max(df.y$time)-1)] # Hs of the new subject
    coverage3[cnt]<-lqmm.cov.new.subject(res2,df.i1,y.i1,y.s,alpha,side)
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

summary_method<-data.frame(matrix(nrow=22, ncol=7))
colnames(summary_method)<-c("par", "beta1", "beta2",
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
  met.iam<-lqmm.fit(df.m,alpha=alpha)
  
  uz.m<-cbind.data.frame(met.iam$u.i1,met.iam$u.i2)
  uz.m$low<-0
  uz.m$up <-met.iam$beta02 + met.iam$u.i2
  uz.m$data<-"metabolomics"
  uz.m$id<-c(1:N)
  colnames(uz.m)<-c("u1","u2", "low", "up", "data","id")
  
  #calculate empirical coverage
  coverg<-covrg(df=df.m, alpha=alpha, side="up")
  
  summary_method[j, 1]<-paste0("met_", colnames(all)[j+2])
  summary_method[j, 2:7]<-c(met.iam$beta01, met.iam$beta02, 
                              coverg[1], coverg[2], coverg[3], coverg[4])
  
  ##CLINICAL
  df.c<-all[,c(1:2)]
  df.c$y<-all[, j+13]
  
  rownames(df.c)<-NULL
  df.c$data<-"clinical"
  
  clin.iam<-lqmm.fit(df.c,alpha=alpha)
  
  uz.c<-cbind.data.frame(clin.iam$u.i1,clin.iam$u.i2)
  uz.c$low<-0
  uz.c$up <-clin.iam$beta02 + clin.iam$u.i2
  uz.c$data<-"clinical"
  uz.c$id<-c(1:N)
  colnames(uz.c)<-c("u1","u2",  "low", "up", "data", "id")
  
  #calculate empirical coverage
  coverag<-covrg(df=df.c, alpha=alpha, side="up")
  
  summary_method[j+11, 1]<-paste0("clin_", colnames(all)[j+13])
  summary_method[j+11, 2:7]<-c(clin.iam$beta01, clin.iam$beta02, 
                                 coverag[1], coverag[2], coverag[3], coverag[4])
  
  ####################
  #one combined plot
  subj<-ran[,j+1]
  
  #random subjects
  uz.a<-rbind.data.frame(uz.c, uz.m)
  uz.a$id2<-subj
  
  df.a<-rbind.data.frame(df.c, df.m)
  df.a$id2<-rep(subj, each=n)
  
  save(uz.a, file=paste0("./output/iam/", colnames(all)[j+2], ".uz.lqmm.RData"))
  save(df.a, file=paste0("./output/iam/", colnames(all)[j+2], ".df.lqmm.RData"))
}


#for RI with both lower and upper bound
for (j in 3:11) {
  df.m<-all[,c(1:2)]
  df.m$y<-all[, j+2]
  
  rownames(df.m)<-NULL
  df.m$data<-"metabolomics"
  
  alpha=0.05
  met.iam<-lqmm.fit(df.m,alpha=alpha)
  
  uz.m<-cbind.data.frame(met.iam$u.i1,met.iam$u.i2)
  uz.m$low<-met.iam$beta01 + met.iam$u.i1
  uz.m$up <-met.iam$beta02 + met.iam$u.i2
  uz.m$data<-"metabolomics"
  uz.m$id<-c(1:N)
  colnames(uz.m)<-c("u1","u2", "low", "up", "data","id")
  
  #calculate empirical coverage
  coverg<-covrg(df=df.m, alpha=alpha, side="both")
  
  summary_method[j, 1]<-paste0("met_", colnames(all)[j+2])
  summary_method[j, 2:7]<-c(met.iam$beta01, met.iam$beta02, 
                              coverg[1], coverg[2], coverg[3], coverg[4])
  
  ##CLINICAL
  df.c<-all[,c(1:2)]
  df.c$y<-all[, j+13]
  
  rownames(df.c)<-NULL
  df.c$data<-"clinical"
  
  clin.iam<-lqmm.fit(df.c,alpha=alpha)
  
  uz.c<-cbind.data.frame(clin.iam$u.i1,clin.iam$u.i2)
  uz.c$low<-clin.iam$beta01 + clin.iam$u.i1
  uz.c$up <-clin.iam$beta02 + clin.iam$u.i2
  uz.c$data<-"clinical"
  uz.c$id<-c(1:N)
  colnames(uz.c)<-c("u1","u2",  "low", "up", "data", "id")
  
  #calculate empirical coverage
  coverag<-covrg(df=df.c, alpha=alpha, side="both")
  
  summary_method[j+11, 1]<-paste0("clin_", colnames(all)[j+13])
  summary_method[j+11, 2:7]<-c(clin.iam$beta01, clin.iam$beta02, 
                                 coverag[1], coverag[2], coverag[3], coverag[4])
  
  ####################
  #one combined plot
  subj<-ran[,j+1]
  
  #random subjects
  uz.a<-rbind.data.frame(uz.c, uz.m)
  uz.a$id2<-subj
  
  df.a<-rbind.data.frame(df.c, df.m)
  df.a$id2<-rep(subj, each=n)
  
  save(uz.a, file=paste0("./output/iam/", colnames(all)[j+2], ".uz.lqmm.RData"))
  save(df.a, file=paste0("./output/iam/", colnames(all)[j+2], ".df.lqmm.RData"))
}


save(summary_method, file="./output/iam/sum.lqmm.RData")
