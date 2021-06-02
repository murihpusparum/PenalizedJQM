####produce Table 2 in the main article

source("./codes/JQM_Function.R")

#load all the coverage functions
##coverage for RI with both lower and upper bound
covrg <- function(df,alpha=0.05){
  N <- length(unique(df$subject))
  n <- length(unique(df$time))
  
  #coverage when there is one new measurement
  coverage<-NULL
  coverage2<-NULL
  for(i in 1:N) {
    df.fit2<-df[(df$time!=n)&(df$subject!=i),]
    df.fit2$subject <- rep(seq(1:(N-1)), each=n-1)
    met.cv<-jqm(db=df.fit2,
                alpha=alpha,
                lambda.u.seq = seq(0.02,0.1,0.02),
                lambda.z.seq = seq(0.5,5,0.5))
    
    cv <- cbind.data.frame(met.cv$u, met.cv$z)
    cv$low <- met.cv$beta0 + met.cv$u + met.cv$z*met.cv$beta1
    cv$up  <- met.cv$beta0 + met.cv$u + met.cv$z*met.cv$beta2
    
    y <- df$y[(df$time==(n))&(df$subject!=i)]
    coverage[i]<-mean(((y > cv$low)&(y < cv$up)))
    
    
    y.i<-df$y[df$subject==i]
    coverage2[i]<-jqm.coverage.new.subject(res=met.cv,y=y.i,alpha=alpha)
  }
  coverage<-mean(coverage)
  coverage2<-mean(coverage2)
  
  co <- c(coverage2, coverage, mean(c(coverage,coverage2)))
  return(co)
}

##coverage for RI with upper bound, lower bound = 0
covrg.up <- function(df,alpha=0.1){
  N <- length(unique(df$subject))
  n <- length(unique(df$time))
  
  #coverage when there is one new measurement
  coverage<-NULL
  coverage2<-NULL
  for(i in 1:N) {
    df.fit2<-df[(df$time!=n)&(df$subject!=i),]
    df.fit2$subject <- rep(seq(1:(N-1)), each=n-1)
    met.cv<-jqm(db=df.fit2,
                alpha=alpha,
                lambda.u.seq = seq(0.02,0.1,0.02),
                lambda.z.seq = seq(0.5,5,0.5))
    
    cv <- cbind.data.frame(met.cv$u, met.cv$z)
    cv$low <- met.cv$beta0 + met.cv$u + met.cv$z*met.cv$beta1
    cv$up  <- met.cv$beta0 + met.cv$u + met.cv$z*met.cv$beta2
    
    y <- df$y[(df$time==(n))&(df$subject!=i)]
    coverage[i]<-mean(((y > 0)&(y < cv$up)))
    
    
    y.i<-df$y[df$subject==i]
    coverage2[i]<-jqm.coverage.new.subject.up(res=met.cv,y=y.i,alpha=alpha)
  }
  coverage<-mean(coverage)
  coverage2<-mean(coverage2)
  
  co <- c(coverage2, coverage, mean(c(coverage,coverage2)))
  return(co)
}

##coverage for RI with lower bound
covrg.low <- function(df,alpha=0.1){
  N <- length(unique(df$subject))
  n <- length(unique(df$time))
  
  #coverage when there is one new measurement
  coverage<-NULL
  coverage2<-NULL
  for(i in 1:N) {
    df.fit2<-df[(df$time!=n)&(df$subject!=i),]
    df.fit2$subject <- rep(seq(1:(N-1)), each=n-1)
    met.cv<-jqm(db=df.fit2,
                alpha=alpha,
                lambda.u.seq = seq(0.02,0.1,0.02),
                lambda.z.seq = seq(0.5,5,0.5))
    
    cv <- cbind.data.frame(met.cv$u, met.cv$z)
    cv$low <- met.cv$beta0 + met.cv$u + met.cv$z*met.cv$beta1
    cv$up  <- met.cv$beta0 + met.cv$u + met.cv$z*met.cv$beta2
    
    y <- df$y[(df$time==(n))&(df$subject!=i)]
    coverage[i]<-mean((y > cv$low))
    
    
    y.i<-df$y[df$subject==i]
    coverage2[i]<-jqm.coverage.new.subject.low(res=met.cv,y=y.i,alpha=alpha)
  }
  coverage<-mean(coverage)
  coverage2<-mean(coverage2)
  
  co <- c(coverage2, coverage, mean(c(coverage,coverage2)))
  return(co)
}

################
#JQM
################
load("./data/combined_data.RData")
load("./data/ran.RData")

jqm_summary_method <- data.frame(matrix(nrow = 22, ncol = 9))
colnames(jqm_summary_method) <- c("par","beta0", "beta1", "beta2", "lambda.u", "lambda.z",
                              "cov.subj", "cov.time", "cov.tot")
#start = Sys.time()
#for RI with upper bound, lower bound = 0
for (j in 1:4) {
  df.fit.m <- m.sub3[,c(1:2)]
  df.fit.m$y <- m.sub3[, j+2]
  
  rownames(df.fit.m) <- NULL
  df.fit.m$data <- "metabolomics"
  
  alpha = 0.1
  met.iam<-jqm(db=df.fit.m,
               alpha=alpha,
               lambda.u.seq = seq(0.02,0.1,0.02),
               lambda.z.seq = seq(0.5,5,0.5))
  
  uz.m <- cbind.data.frame(met.iam$u, met.iam$z)
  uz.m$low <- 0
  uz.m$up  <- met.iam$beta0 + met.iam$u + met.iam$z*met.iam$beta2
  uz.m$data <- "metabolomics"
  uz.m$id <- c(1:30)
  colnames(uz.m) <- c("u", "z", "low", "up", "data","id")
  
  #calculate empirical coverage
  coverg <- covrg.up(df = df.fit.m, alpha = alpha)
  
  jqm_summary_method[j, 1] <- paste0("met_", colnames(m.sub3)[j+2])
  jqm_summary_method[j, 2:9] <- c(met.iam$beta0, met.iam$beta1, met.iam$beta2, met.iam$lambda.u,
                             met.iam$lambda.z, coverg[1], coverg[2], coverg[3])

  ##CLINICAL
  df.fit <- m.sub3[,c(1:2)]
  df.fit$y <- m.sub3[, j+13]
  
  rownames(df.fit) <- NULL
  df.fit$data <- "clinical"
  
  res.iam<-jqm(db=df.fit,
               alpha=alpha,
               lambda.u.seq = seq(0.02,0.1,0.02),
               lambda.z.seq = seq(0.5,5,0.5))
  
  uz <- cbind.data.frame(res.iam$u, res.iam$z)
  uz$low <- 0
  uz$up  <- res.iam$beta0 + res.iam$u + res.iam$z*res.iam$beta2
  uz$data <- "clinical"
  uz$id <- c(1:30)
  colnames(uz) <- c("u", "z", "low", "up", "data", "id")
  
  #calculate empirical coverage
  coverag <- covrg.up(df = df.fit, alpha = alpha)
  
  jqm_summary_method[j+11, 1] <- paste0("clin_", colnames(m.sub3)[j+13])
  jqm_summary_method[j+11, 2:9] <- c(met.iam$beta0, met.iam$beta1, met.iam$beta2, met.iam$lambda.u,
                             met.iam$lambda.z, coverag[1], coverag[2], coverag[3])
  
  ####################
  #one combined plot
  subj <- ran[,j+1]
  
  #random subjects
  uz.a <- rbind.data.frame(uz, uz.m)
  uz.a$id2 <- subj
  
  df.a <- rbind.data.frame(df.fit, df.fit.m)
  df.a$id2 <- rep(subj, each=6)
  
  save(uz.a, file=paste0("./output/IAM/", colnames(m.sub3)[j+2], ".uz.jqm.RData"))
}


#for RI with both lower and upper bound
for (j in 5:10) {
  df.fit.m <- m.sub3[,c(1:2)]
  df.fit.m$y <- m.sub3[, j+2]
  
  rownames(df.fit.m) <- NULL
  df.fit.m$data <- "metabolomics"
  
  alpha = 0.05
  met.iam<-jqm(db=df.fit.m,
               alpha=alpha,
               lambda.u.seq = seq(0.02,0.1,0.02),
               lambda.z.seq = seq(0.5,5,0.5))
  
  uz.m <- cbind.data.frame(met.iam$u, met.iam$z)
  uz.m$low <- met.iam$beta0 + met.iam$u + met.iam$z*met.iam$beta1
  uz.m$up  <- met.iam$beta0 + met.iam$u + met.iam$z*met.iam$beta2
  uz.m$data <- "metabolomics"
  uz.m$id <- c(1:30)
  colnames(uz.m) <- c("u", "z", "low", "up", "data","id")
  
  #calculate empirical coverage
  coverg <- covrg(df = df.fit.m, alpha = alpha)
  
  jqm_summary_method[j, 1] <- paste0("met_", colnames(m.sub3)[j+2])
  jqm_summary_method[j, 2:9] <- c(met.iam$beta0, met.iam$beta1, met.iam$beta2, met.iam$lambda.u,
                              met.iam$lambda.z, coverg[1], coverg[2], coverg[3])
  
  ##CLINICAL
  df.fit <- m.sub3[,c(1:2)]
  df.fit$y <- m.sub3[, j+13]
  
  rownames(df.fit) <- NULL
  df.fit$data <- "clinical"
  
  res.iam<-jqm(db=df.fit,
               alpha=alpha,
               lambda.u.seq = seq(0.02,0.1,0.02),
               lambda.z.seq = seq(0.5,5,0.5))
  
  uz <- cbind.data.frame(res.iam$u, res.iam$z)
  uz$low <- res.iam$beta0 + res.iam$u + res.iam$z*res.iam$beta1
  uz$up  <- res.iam$beta0 + res.iam$u + res.iam$z*res.iam$beta2
  uz$data <- "clinical"
  uz$id <- c(1:30)
  colnames(uz) <- c("u", "z", "low", "up", "data", "id")
  
  #calculate empirical coverage
  coverag <- covrg(df = df.fit, alpha = alpha)
  
  jqm_summary_method[j+11, 1] <- paste0("clin_", colnames(m.sub3)[j+13])
  jqm_summary_method[j+11, 2:9] <- c(met.iam$beta0, met.iam$beta1, met.iam$beta2, met.iam$lambda.u,
                                 met.iam$lambda.z, coverag[1], coverag[2], coverag[3])
  
  ####################
  #one combined plot
  subj <- ran[,j+1]
  
  #random subjects
  uz.a <- rbind.data.frame(uz, uz.m)
  uz.a$id2 <- subj
  
  df.a <- rbind.data.frame(df.fit, df.fit.m)
  df.a$id2 <- rep(subj, each=6)
  
  save(uz.a, file=paste0("./output/IAM/", colnames(m.sub3)[j+2], ".uz.jqm.RData"))
}

#for RI with lower bound, upper bound = NA
for (j in 11) {
  df.fit.m <- m.sub3[,c(1:2)]
  df.fit.m$y <- m.sub3[, j+2]
  
  rownames(df.fit.m) <- NULL
  df.fit.m$data <- "metabolomics"
  
  alpha = 0.1
  met.iam<-jqm(db=df.fit.m,
               alpha=alpha,
               lambda.u.seq = seq(0.02,0.1,0.02),
               lambda.z.seq = seq(0.5,5,0.5))
  
  uz.m <- cbind.data.frame(met.iam$u, met.iam$z)
  uz.m$low <- met.iam$beta0 + met.iam$u + met.iam$z*met.iam$beta1
  uz.m$up  <- NA
  uz.m$data <- "metabolomics"
  uz.m$id <- c(1:30)
  colnames(uz.m) <- c("u", "z", "low", "up", "data","id")
  
  #calculate empirical coverage
  coverg <- covrg.low(df = df.fit.m, alpha = alpha)
  
  jqm_summary_method[j, 1] <- paste0("met_", colnames(m.sub3)[j+2])
  jqm_summary_method[j, 2:9] <- c(met.iam$beta0, met.iam$beta1, met.iam$beta2, met.iam$lambda.u,
                              met.iam$lambda.z, coverg[1], coverg[2], coverg[3])
  
  ##CLINICAL
  df.fit <- m.sub3[,c(1:2)]
  df.fit$y <- m.sub3[, j+13]
  
  rownames(df.fit) <- NULL
  df.fit$data <- "clinical"
  
  res.iam<-jqm(db=df.fit,
               alpha=alpha,
               lambda.u.seq = seq(0.02,0.1,0.02),
               lambda.z.seq = seq(0.5,5,0.5))
  
  uz <- cbind.data.frame(res.iam$u, res.iam$z)
  uz$low <- res.iam$beta0 + res.iam$u + res.iam$z*res.iam$beta1
  uz$up  <- NA
  uz$data <- "clinical"
  uz$id <- c(1:30)
  colnames(uz) <- c("u", "z", "low", "up", "data", "id")
  
  #calculate empirical coverage
  coverag <- covrg.low(df = df.fit, alpha = alpha)
  
  jqm_summary_method[j+11, 1] <- paste0("clin_", colnames(m.sub3)[j+13])
  jqm_summary_method[j+11, 2:9] <- c(met.iam$beta0, met.iam$beta1, met.iam$beta2, met.iam$lambda.u,
                                 met.iam$lambda.z, coverag[1], coverag[2], coverag[3])
  
  ####################
  #one combined plot
  subj <- ran[,j+1]
  
  #random subjects
  uz.a <- rbind.data.frame(uz, uz.m)
  uz.a$id2 <- subj
  
  df.a <- rbind.data.frame(df.fit, df.fit.m)
  df.a$id2 <- rep(subj, each=6)
  
  save(uz.a, file=paste0("./output/IAM/", colnames(m.sub3)[j+2], ".uz.jqm.RData"))
}

save(jqm_summary_method, file="./output/jqm_summary_method.RData")