#estimate separate-IRI and joint-IRI from VITO IAM Frontier data
source("JQM_Function.R")

combined_data <- read.csv("./data/long_combined_data.csv", sep=";", header = T)

################################
#DATA PREPARATION
################################
#Modify from long to wide format
combined_data$person_id <- as.character(combined_data$person_id)
combined_data$subject <- lapply(combined_data$person_id, parse_number)
combined_data$subject <- as.numeric(combined_data$subject)

combined_data$time <- combined_data$month
combined_data[combined_data$time == 3,]$time <- 2
combined_data[combined_data$time == 5,]$time <- 3
combined_data[combined_data$time == 7,]$time <- 4
combined_data[combined_data$time == 9,]$time <- 5
combined_data[combined_data$time == 11,]$time <- 6

CLINICAL <- combined_data[combined_data$type=="clinical",]
METABOLOMICS <- combined_data[combined_data$type=="metabolomics",]

CLINICAL_wide <- dcast(CLINICAL, subject+time ~ label, value.var = "value", fun.aggregate = mean)
CLINICAL_wide <- CLINICAL_wide[, c(1:3,5:13,4)]
colnames(CLINICAL_wide) <- c("subject", "time", "Albumin", "ApoA1", "ApoB", "Total_C", "Creatinine",
                             "Glucose", "HDL_C", "Clinical_LDL_C", 
                             "non_HDL_C", "Total_TG", "ApoB_by_ApoA1")
METABOLOMICS_wide <- dcast(METABOLOMICS, subject+time ~ label, value.var = "value", fun.aggregate = mean)
METABOLOMICS_wide <- METABOLOMICS_wide[, c(1:5,12,8:10,7,11,13,6)]

all.data_wide <- cbind.data.frame(METABOLOMICS_wide, CLINICAL_wide[,c(3:13)]) 

########################
#compute the separate-IRI
########################

#function for computing the empirical coverage
covrg.lqmm <- function(df){
  N <- length(unique(df$subject))
  n <- length(unique(df$time))
  
  coverage<-NULL
  coverage2<-NULL
  for(i in 1:N) {
    #coverage when there is one new measurement
    df.fit<-df[(df$time!=n)&(df$subject!=i),]
    fit1 <- lqmm(fixed = y ~ 1, random = ~ 1, group = subject, tau = c(0.025, 0.975),
                 nK = 11, type = "normal", data = df.fit, 
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
    coverage2[i]<-mean((df.i$y > est1+u1.i)& (df.i$y < est2+u2.i))
    
  }
  coverage<-mean(coverage)
  coverage2<-mean(coverage2)
  co <- c(coverage2, coverage, mean(c(coverage,coverage2)))
  return(co)
}

summary_lqmm <- data.frame(matrix(nrow = 22, ncol = 6))
colnames(summary_lqmm) <- c("par","beta01", "beta02",
                                   "cov.subj", "cov.time", "cov.tot")
for (j in 1:11) {
  #for metabolomics
  df.m <- all.data_wide[,c(1:2)]
  df.m$y <- all.data_wide[, j+2]
  rownames(df.m) <- NULL
  df.m$data <- "metabolomics"
  df.m$id <- df.m$subject
  
  met.iam <- lqmm(fixed = y ~ 1, random = ~ 1, group = subject, tau = c(0.025, 0.975),
                  nK = 11, type = "normal", data = df.m, 
                  control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 1000, UP_tol = 1e-5,
                                        beta = 0.5, gamma = 1))
  pred1 <- predict(met.iam, level = 1)
  pred1 <- data.frame(pred1)
  a <- cbind.data.frame(df.m$subject, pred1)
  pred1 <- data.frame(a %>% distinct())
  beta1m <- met.iam$theta_x[1,1]
  beta2m <- met.iam$theta_x[1,2]
  
  uz.m <- pred1[,2:3]
  uz.m$data <- "metabolomics"
  uz.m$id <- unique(df.m$subject)
  colnames(uz.m) <- c("low", "up", "data", "id")
  
  #calculate empirical coverage
  cvrg.m <- covrg.lqmm(df = df.m)
  
  summary_lqmm[j, 1] <- paste0("met_", colnames(all.data_wide)[j+2])
  summary_lqmm[j, 2:6] <- c(beta1m, beta2m, cvrg.m[1], cvrg.m[2], cvrg.m[3])
  
  
  #for clinical
  df.c <- all.data_wide[,c(1:2)]
  df.c$y <- all.data_wide[, j+13]
  rownames(df.c) <- NULL
  df.c$data <- "clinical"
  df.c$id <- df.c$subject
  
  clin.iam <- lqmm(fixed = y ~ 1, random = ~ 1, group = subject, tau = c(0.025, 0.975),
                   nK = 11, type = "normal", data = df.c, 
                   control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 100, UP_tol = 1e-5,
                                         beta = 0.5, gamma = 1))
  pred2 <- predict(clin.iam, level = 1)
  pred2 <- data.frame(pred2)
  a <- cbind.data.frame(df.c$subject, pred2)
  pred2 <- data.frame(a %>% distinct())
  beta1c <- clin.iam$theta_x[1,1]
  beta2c <- clin.iam$theta_x[1,2]
  
  uz.c <- pred2[,2:3]
  uz.c$data <- "clinical"
  uz.c$id <- unique(df.c$subject)
  colnames(uz.c) <- c("low", "up", "data", "id")
  
  #calculate empirical coverage
  cvrg.c <- covrg.lqmm(df = df.c)
  
  summary_lqmm[j+11, 1] <- paste0("clin_", colnames(all.data_wide)[j+13])
  summary_lqmm[j+11, 2:6] <- c(beta1c, beta2c, cvrg.c[1], cvrg.c[2], cvrg.c[3])
}

########################
#compute the joint-IRI
########################

#function for computing the empirical coverage
covrg.jqm <- function(df,alpha){
  N <- length(unique(df$subject))
  n <- length(unique(df$time))
  
  #coverage when there is one new measurement
  coverage<-NULL
  coverage2<-NULL
  for (i in 1:N) {
    #coverage when there is one new measurement
    df.fit<-df[(df$time!=n)&(df$subject!=i),]
    
    df.fit2<-df[(df$time!=n)&(df$subject!=N),]
    #df.fit2$subject <- rep(seq(1:(N-1)), each=n-1)
    met.cv<-jqm(db=df.fit,
                alpha=alpha,
                lambda.u.seq = seq(0.02,0.1,0.02),
                lambda.z.seq = seq(0.5,5,0.5))
    
    cv <- cbind.data.frame(met.cv$u, met.cv$z)
    cv$low <- met.cv$beta0 + met.cv$u + met.cv$z*met.cv$beta1
    cv$up  <- met.cv$beta0 + met.cv$u + met.cv$z*met.cv$beta2
    
    y <- df$y[(df$time==(n))&(df$subject!=i)]
    coverage[i]<-mean((y>cv$low[s])&
                        (y<cv$up[s]))
    
    #coverage when there is one new subject
    df.i <- df[df$subject==i,]
    met.cv2<-jqm(db=df.i,
                 alpha=alpha,
                 lambda.u.seq = seq(0.02,0.1,0.02),
                 lambda.z.seq = seq(0.5,5,0.5))
    coverage2[i]<-jqm.coverage.new.subject(res=met.cv2,y=df.i$y)
  }
  coverage<-mean(coverage)
  coverage2<-mean(coverage2)
  co <- c(coverage2, coverage, mean(c(coverage,coverage2)))
  return(co)
}

summary_jqm <- data.frame(matrix(nrow = 22, ncol = 9))
colnames(summary_jqm) <- c("par","beta0", "beta1", "beta2", "lambda.u", "lambda.z",
                              "cov.subj", "cov.time", "cov.tot")
alpha = 0.05
for (j in 1:11) {
  #for metabolomics
  df.m <- all.data_wide[,c(1:2)]
  df.m$y <- all.data_wide[, j+2]
  rownames(df.m) <- NULL
  df.m$data <- "metabolomics"
  df.m$id <- df.m$subject
  
  met.iam<-jqm(db=df.m,
               alpha=alpha,
               lambda.u.seq = seq(0.02,0.1,0.02),
               lambda.z.seq = seq(0.5,5,0.5))
  uz.m <- cbind.data.frame(met.iam$u, met.iam$z)
  uz.m$low <- met.iam$beta0 + met.iam$u + met.iam$z*met.iam$beta1
  uz.m$up  <- met.iam$beta0 + met.iam$u + met.iam$z*met.iam$beta2
  uz.m$data <- "metabolomics"
  uz.m$id <- unique(df.m$subject)
  colnames(uz.m) <- c("u", "z", "low", "up", "data", "id")
  
  #calculate empirical coverage
  cvrg.m <- covrg.jqm(df = df.m, alpha = alpha)
  
  summary_jqm[j, 1] <- paste0("met_", colnames(all.data_wide)[j+2])
  summary_jqm[j, 2:9] <- c(met.iam$beta0, met.iam$beta1, met.iam$beta2, met.iam$lambda.u,
                             met.iam$lambda.z, cvrg.m[1], cvrg.m[2], cvrg.m[3])

  #for clinical
  df.c <- all.data_wide[,c(1:2)]
  df.c$y <- all.data_wide[, j+13]
  rownames(df.c) <- NULL
  df.c$data <- "clinical"
  df.c$id <- df.c$subject
  
  clin.iam<-jqm(db=df.c,
               alpha=alpha,
               lambda.u.seq = seq(0.02,0.1,0.02),
               lambda.z.seq = seq(0.5,5,0.5))
  
  uz.c <- cbind.data.frame(clin.iam$u, clin.iam$z)
  uz.c$low <- clin.iam$beta0 + clin.iam$u + clin.iam$z*clin.iam$beta1
  uz.c$up  <- clin.iam$beta0 + clin.iam$u + clin.iam$z*clin.iam$beta2
  uz.c$data <- "clinical"
  uz.c$id <- unique(df.c$subject)
  colnames(uz.c) <- c("u", "z", "low", "up", "data", "id")
  
  #calculate empirical coverage
  cvrg.c <- covrg.jqm(df = df.c, alpha = alpha)
  
  summary_jqm[j+11, 1] <- paste0("clin_", colnames(all.data_wide)[j+13])
  summary_jqm[j+11, 2:9] <- c(clin.iam$beta0, clin.iam$beta1, clin.iam$beta2, clin.iam$lambda.u,
                             clin.iam$lambda.z, cvrg.c[1], cvrg.c[2], cvrg.c[3])
  
}

