####produce Table 2 in the main article

library(lqmm)
library(tidyverse)
#load all the coverage functions
##coverage for RI with upper bound, lower bound = 0
covrg <- function(df){
  N <- length(unique(df$subject))
  n <- length(unique(df$time))
  
  #coverage when there is one new measurement
  coverage<-NULL
  coverage2<-NULL
  for(i in 1:N) {
    df.fit2<-df[(df$time!=n)&(df$subject!=i),]
    fit1 <- lqmm(fixed = y ~ 1, random = ~ 1, group = subject, tau = c(0.025, 0.975),
                 nK = 11, type = "normal", data = df.fit2, 
                 control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 1000, UP_tol = 1e-5,
                                       beta = 0.5, gamma = 1))
    pred1 <- predict(fit1, level = 1)
    pred1 <- data.frame(pred1)
    a <- cbind.data.frame(df.fit2$subject, pred1)
    pred1 <- data.frame(a %>% distinct())
    ri1 <- pred1$X0.025
    ri2 <- pred1$X0.975
    
    est1 <- fit1$theta_x[1,1]
    est2 <- fit1$theta_x[1,2]
    
    y <- df$y[(df$time==(n))&(df$subject!=i)]
    coverage[i]<-mean(((y > ri1)&(y < ri2)))
    
    #coverage when there is one new subject
    df.i <- df[df$subject==i,]
    fit2 <- lqmm(fixed = y ~ 1, random = ~ 1, group = subject, tau = c(0.025, 0.975),
                 nK = 11, type = "normal", data = df.i, 
                 control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 1000, UP_tol = 1e-5,
                                       beta = 0.5, gamma = 1))
    u1.i <- ranef(fit2)[[1]][1,1]
    u2.i <- ranef(fit2)[[2]][1,1]
    ri.low <- est1+u1.i
    ri.up <- est2+u2.i
    if(ri.up < ri.low){
      ri.low2 <- ri.low
      ri.low <- ri.up
      ri.up <- ri.low2
    }
    
    coverage2[i]<-mean((df.i$y > ri.low)& (df.i$y < ri.up))
    
  }
  coverage<-mean(coverage)
  coverage2<-mean(coverage2)
  
  co <- c(coverage2, coverage, mean(c(coverage,coverage2)))
  return(co)
}

#coverage for only RI with upper bound, lower bound = 0
covrg.up <- function(df){
  N <- length(unique(df$subject))
  n <- length(unique(df$time))
  
  #coverage when there is one new measurement
  coverage<-NULL
  coverage2<-NULL
  for(i in 1:N) {
    df.fit2<-df[(df$time!=n)&(df$subject!=i),]
    fit1 <- lqmm(fixed = y ~ 1, random = ~ 1, group = subject, tau = c(0.05, 0.95),
                 nK = 11, type = "normal", data = df.fit2, 
                 control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 1000, UP_tol = 1e-5,
                                       beta = 0.5, gamma = 1))
    pred1 <- predict(fit1, level = 1)
    pred1 <- data.frame(pred1)
    a <- cbind.data.frame(df.fit2$subject, pred1)
    pred1 <- data.frame(a %>% distinct())
    ri1 <- 0
    ri2 <- pred1$X0.95
    
    est1 <- fit1$theta_x[1,1]
    est2 <- fit1$theta_x[1,2]
    
    y <- df$y[(df$time==(n))&(df$subject!=i)]
    coverage[i]<-mean(((y > ri1)&(y < ri2)))
    
    #coverage when there is one new subject
    df.i <- df[df$subject==i,]
    fit2 <- lqmm(fixed = y ~ 1, random = ~ 1, group = subject, tau = c(0.025, 0.95),
                 nK = 11, type = "normal", data = df.i, 
                 control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 1000, UP_tol = 1e-5,
                                       beta = 0.5, gamma = 1))
    u1.i <- ranef(fit2)[[1]][1,1]
    u2.i <- ranef(fit2)[[2]][1,1]
    ri.low <- est1+u1.i
    ri.up <- est2+u2.i
    if(ri.up < ri.low){
      ri.low2 <- ri.low
      ri.low <- ri.up
      ri.up <- ri.low2
    }
    
    coverage2[i]<-mean((df.i$y > 0)& (df.i$y < ri.up))
    
  }
  coverage<-mean(coverage)
  coverage2<-mean(coverage2)
  
  co <- c(coverage2, coverage, mean(c(coverage,coverage2)))
  return(co)
}

#coverage for only RI with lower bound
covrg.low <- function(df){
  N <- length(unique(df$subject))
  n <- length(unique(df$time))
  
  #coverage when there is one new measurement
  coverage<-NULL
  coverage2<-NULL
  for(i in 1:N) {
    df.fit2<-df[(df$time!=n)&(df$subject!=i),]
    fit1 <- lqmm(fixed = y ~ 1, random = ~ 1, group = subject, tau = c(0.05, 0.95),
                 nK = 11, type = "normal", data = df.fit2, 
                 control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 1000, UP_tol = 1e-5,
                                       beta = 0.5, gamma = 1))
    pred1 <- predict(fit1, level = 1)
    pred1 <- data.frame(pred1)
    a <- cbind.data.frame(df.fit2$subject, pred1)
    pred1 <- data.frame(a %>% distinct())
    ri1 <- pred1$X0.05
    ri2 <- pred1$X0.95
    
    est1 <- fit1$theta_x[1,1]
    est2 <- fit1$theta_x[1,2]
    
    y <- df$y[(df$time==(n))&(df$subject!=i)]
    coverage[i]<-mean((y > ri1))
    
    #coverage when there is one new subject
    df.i <- df[df$subject==i,]
    fit2 <- lqmm(fixed = y ~ 1, random = ~ 1, group = subject, tau = c(0.05, 0.975),
                 nK = 11, type = "normal", data = df.i, 
                 control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 1000, UP_tol = 1e-5,
                                       beta = 0.5, gamma = 1))
    u1.i <- ranef(fit2)[[1]][1,1]
    u2.i <- ranef(fit2)[[2]][1,1]
    ri.low <- est1+u1.i
    ri.up <- est2+u2.i
    if(ri.up < ri.low){
      ri.low2 <- ri.low
      ri.low <- ri.up
      ri.up <- ri.low2
    }
    
    coverage2[i]<-mean((df.i$y > ri.low))
    
  }
  coverage<-mean(coverage)
  coverage2<-mean(coverage2)
  
  co <- c(coverage2, coverage, mean(c(coverage,coverage2)))
  return(co)
}



################
#LQMM
load("./data/combined_data.RData")
load("./data/ran.RData")

lqmm_summary_method <- data.frame(matrix(nrow = 22, ncol = 6))
colnames(lqmm_summary_method) <- c("par","beta01", "beta02",
                              "cov.subj", "cov.time", "cov.tot")
start = Sys.time()
#for RI with upper bound, lower bound = 0
for (j in 1:4) {
  df.fit.m <- m.sub3[,c(1:2)]
  df.fit.m$y <- m.sub3[, j+2]
  
  rownames(df.fit.m) <- NULL
  df.fit.m$data <- "metabolomics"
  met.iam <- lqmm(fixed = y ~ 1, random = ~ 1, group = subject, tau = c(0.05, 0.95),
               nK = 10, type = "normal", data = df.fit.m, 
               control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 10000, UP_tol = 1e-4,
                                     beta = 0.5, gamma = 1))
  pred1 <- predict(met.iam, level = 1)
  pred1 <- data.frame(pred1)
  a <- cbind.data.frame(df.fit.m$subject, pred1)
  pred1 <- data.frame(a %>% distinct())
  pred1[,2] <- 0
  beta1m <- met.iam$theta_x[1,1]
  beta2m <- met.iam$theta_x[1,2]

  uz.m <- pred1[,2:3]
  uz.m$data <- "metabolomics"
  uz.m$id <- c(1:30)
  colnames(uz.m) <- c("low", "up", "data", "id")
  
  cvrg1 <- covrg.up(df.fit.m)
  
  lqmm_summary_method[j, 1] <- paste0("met_", colnames(m.sub3)[j+2])
  lqmm_summary_method[j, 2:6] <- c(beta1m, beta2m, 
                               cvrg1[1], cvrg1[2], cvrg1[3])

  df.fit <- m.sub3[,c(1:2)]
  df.fit$y <- m.sub3[, j+13]
  
  rownames(df.fit) <- NULL
  df.fit$data <- "clinical"
  res.iam <- lqmm(fixed = y ~ 1, random = ~ 1, group = subject, tau = c(0.05, 0.95),
                  nK = 13, type = "normal", data = df.fit, 
                  control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 1000, UP_tol = 1e-4,
                                        beta = 0.5, gamma = 1))
  pred2 <- predict(res.iam, level = 1)
  pred2 <- data.frame(pred2)
  a <- cbind.data.frame(df.fit$subject, pred2)
  pred2 <- data.frame(a %>% distinct())
  pred2[,2] <- 0 
  beta1c <- res.iam$theta_x[1,1]
  beta2c <- res.iam$theta_x[1,2]
  
  uz <- pred2[,2:3]
  uz$data <- "clinical"
  uz$id <- c(1:30)
  colnames(uz) <- c("low", "up", "data", "id")
  
  cvrg2 <- covrg.up(df.fit)
  
  lqmm_summary_method[j+11, 1] <- paste0("clin_", colnames(m.sub3)[j+13])
  lqmm_summary_method[j+11, 2:6] <- c(beta1c, beta2c, 
                                 cvrg2[1], cvrg2[2], cvrg2[3])
  
  ####################
  #one combined plot
  subj <- ran[,j+1]
  
  #random subjects
  uz.a <- rbind.data.frame(uz, uz.m)
  uz.a$id2 <- subj
  
  df.a <- rbind.data.frame(df.fit, df.fit.m)
  df.a$id2 <- rep(rep(subj, each=6),2)
  
  save(uz.a, file=paste0("./output/IAM/", colnames(m.sub3)[j+2], ".uz.lqmm.RData"))
  save(df.a, file=paste0("./output/IAM/", colnames(m.sub3)[j+2], ".df.lqmm.RData"))
}

#for RI with both lower and upper bound
for (j in 5:10) {
  df.fit.m <- m.sub3[,c(1:2)]
  df.fit.m$y <- m.sub3[, j+2]
  
  rownames(df.fit.m) <- NULL
  df.fit.m$data <- "metabolomics"
  met.iam <- lqmm(fixed = y ~ 1, random = ~ 1, group = subject, tau = c(0.025, 0.975),
                  nK = 11, type = "normal", data = df.fit.m, 
                  control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 1000, UP_tol = 1e-4,
                                        beta = 0.5, gamma = 1))
  pred1 <- predict(met.iam, level = 1)
  pred1 <- data.frame(pred1)
  a <- cbind.data.frame(df.fit.m$subject, pred1)
  pred1 <- data.frame(a %>% distinct())
  beta1m <- met.iam$theta_x[1,1]
  beta2m <- met.iam$theta_x[1,2]
  
  uz.m <- pred1[,2:3]
  uz.m$data <- "metabolomics"
  uz.m$id <- c(1:30)
  colnames(uz.m) <- c("low", "up", "data", "id")
  
  cvrg1 <- covrg(df.fit.m)
  
  lqmm_summary_method[j, 1] <- paste0("met_", colnames(m.sub3)[j+2])
  lqmm_summary_method[j, 2:6] <- c(beta1m, beta2m, 
                                   cvrg1[1], cvrg1[2], cvrg1[3])
  
  df.fit <- m.sub3[,c(1:2)]
  df.fit$y <- m.sub3[, j+13]
  
  rownames(df.fit) <- NULL
  df.fit$data <- "clinical"
  res.iam <- lqmm(fixed = y ~ 1, random = ~ 1, group = subject, tau = c(0.025, 0.975),
                  nK = 11, type = "normal", data = df.fit, 
                  control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 1000, UP_tol = 1e-4,
                                        beta = 0.5, gamma = 1))
  pred2 <- predict(res.iam, level = 1)
  pred2 <- data.frame(pred2)
  a <- cbind.data.frame(df.fit$subject, pred2)
  pred2 <- data.frame(a %>% distinct())
  beta1c <- res.iam$theta_x[1,1]
  beta2c <- res.iam$theta_x[1,2]
  
  uz <- pred2[,2:3]
  uz$data <- "clinical"
  uz$id <- c(1:30)
  colnames(uz) <- c("low", "up", "data", "id")
  
  cvrg2 <- covrg(df.fit)
  
  lqmm_summary_method[j+11, 1] <- paste0("clin_", colnames(m.sub3)[j+13])
  lqmm_summary_method[j+11, 2:6] <- c(beta1c, beta2c, 
                                      cvrg2[1], cvrg2[2], cvrg2[3])
  
  ####################
  #one combined plot
  subj <- ran[,j+1]
  
  #random subjects
  uz.a <- rbind.data.frame(uz, uz.m)
  uz.a$id2 <- subj
  
  df.a <- rbind.data.frame(df.fit, df.fit.m)
  df.a$id2 <- rep(rep(subj, each=6),2)
  
  save(uz.a, file=paste0("./output/IAM/", colnames(m.sub3)[j+2], ".uz.lqmm.RData"))
  save(df.a, file=paste0("./output/IAM/", colnames(m.sub3)[j+2], ".df.lqmm.RData"))
}

#for only RI with lower bound
for (j in 11) {
  df.fit.m <- m.sub3[,c(1:2)]
  df.fit.m$y <- m.sub3[, j+2]
  
  rownames(df.fit.m) <- NULL
  df.fit.m$data <- "metabolomics"
  met.iam <- lqmm(fixed = y ~ 1, random = ~ 1, group = subject, tau = c(0.05, 0.95),
                  nK = 11, type = "normal", data = df.fit.m, 
                  control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 1000, UP_tol = 1e-4,
                                        beta = 0.5, gamma = 1))
  pred1 <- predict(met.iam, level = 1)
  pred1 <- data.frame(pred1)
  a <- cbind.data.frame(df.fit.m$subject, pred1)
  pred1 <- data.frame(a %>% distinct())
  pred1[,3] <- NA
  beta1m <- met.iam$theta_x[1,1]
  beta2m <- met.iam$theta_x[1,2]
  
  uz.m <- pred1[,2:3]
  uz.m$data <- "metabolomics"
  uz.m$id <- c(1:30)
  colnames(uz.m) <- c("low", "up", "data", "id")
  
  cvrg1 <- covrg.low(df.fit.m)
  
  lqmm_summary_method[j, 1] <- paste0("met_", colnames(m.sub3)[j+2])
  lqmm_summary_method[j, 2:6] <- c(beta1m, beta2m, 
                                   cvrg1[1], cvrg1[2], cvrg1[3])
  
  df.fit <- m.sub3[,c(1:2)]
  df.fit$y <- m.sub3[, j+13]
  
  rownames(df.fit) <- NULL
  df.fit$data <- "clinical"
  res.iam <- lqmm(fixed = y ~ 1, random = ~ 1, group = subject, tau = c(0.05, 0.95),
                  nK = 11, type = "normal", data = df.fit, 
                  control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 1000, UP_tol = 1e-4,
                                        beta = 0.5, gamma = 1))
  pred2 <- predict(res.iam, level = 1)
  pred2 <- data.frame(pred2)
  a <- cbind.data.frame(df.fit$subject, pred2)
  pred2 <- data.frame(a %>% distinct())
  pred2[,3] <- NA 
  beta1c <- res.iam$theta_x[1,1]
  beta2c <- res.iam$theta_x[1,2]
  
  uz <- pred2[,2:3]
  uz$data <- "clinical"
  uz$id <- c(1:30)
  colnames(uz) <- c("low", "up", "data", "id")
  
  cvrg2 <- covrg.low(df.fit)
  
  lqmm_summary_method[j+11, 1] <- paste0("clin_", colnames(m.sub3)[j+13])
  lqmm_summary_method[j+11, 2:6] <- c(beta1c, beta2c, 
                                      cvrg2[1], cvrg2[2], cvrg2[3])
  
  ####################
  #one combined plot
  subj <- ran[,j+1]
  
  #random subjects
  uz.a <- rbind.data.frame(uz, uz.m)
  uz.a$id2 <- subj
  
  df.a <- rbind.data.frame(df.fit, df.fit.m)
  df.a$id2 <- rep(rep(subj, each=6),2)
  
  save(uz.a, file=paste0("./output/IAM/", colnames(m.sub3)[j+2], ".uz.lqmm.RData"))
  save(df.a, file=paste0("./output/IAM/", colnames(m.sub3)[j+2], ".df.lqmm.RData"))
}
end = Sys.time()
end - start

save(lqmm_summary_method, file="./output/lqmm_summary_method.RData")