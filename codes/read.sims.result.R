###################
#read Penalized JQM simulation results stored in "./output/PJQM" directory
###################

####load data resulted from JQM_Sims.R in Normal PJQM directory
my_files <- list.files("./output/Normal PJQM", pattern = "\\.Rdata$", full.names=TRUE)
NRow <- length(my_files)

all_res <- lapply(my_files, function(x) mget(load(x)))
mean.cov <- data.frame(matrix(nrow=NRow, ncol=length(colnames(all_res[[1]]$Results))))
colnames(mean.cov) <- colnames(all_res[[1]]$Results)
mean.cov$EmpCov.median <- NA
mean.cov <- mean.cov[,-1]
for (i in 1:NRow) {
  mean.cov[i,1] <- all_res[[i]]$Results$alpha[1]
  mean.cov[i,2:12] <- colMeans(all_res[[i]]$Results[,3:13])
  mean.cov[i,13] <- all_res[[i]]$Results$seed[1]
  mean.cov[i,14] <- median(all_res[[i]]$Results[,13])
}
#get the info about N, n, a, and b
resnames <- strsplit(my_files, split = "_")
mean.cov$N <- 0
mean.cov$n <- NA
mean.cov$a <- NA
mean.cov$b <- NA
for (i in 1:NRow) {
  mean.cov[i,15] <- resnames[[i]][6] #N
  mean.cov[i,16] <- resnames[[i]][8] #n
  mean.cov[i,17] <- resnames[[i]][12] #a
  mean.cov[i,18] <- resnames[[i]][14] #b
}
#only take necessary info
mean.res <- mean.cov[,c(15,16,1,17,18,10:12,14,5,6,13)]
mean.res$N <- as.numeric(mean.res$N)
mean.res$n <- as.numeric(mean.res$n)
mean.res$a <- as.numeric(mean.res$a)
mean.res$b <- as.numeric(mean.res$b)

mean.jqm <- mean.res[,c(1:9)]
mean.jqm$method <- "Normal PJQM"

####load data resulted from JQM_Sims.R in Chi-sq PJQM directory
my_files <- list.files("./output/Chi-sq PJQM", pattern = "\\.Rdata$", full.names=TRUE)
NRow <- length(my_files)

all_res <- lapply(my_files, function(x) mget(load(x)))
mean.cov <- data.frame(matrix(nrow=NRow, ncol=length(colnames(all_res[[1]]$Results))))
colnames(mean.cov) <- colnames(all_res[[1]]$Results)
mean.cov$EmpCov.median <- NA
mean.cov <- mean.cov[,-1]
for (i in 1:NRow) {
  mean.cov[i,1] <- all_res[[i]]$Results$alpha[1]
  mean.cov[i,2:12] <- colMeans(all_res[[i]]$Results[,3:13])
  mean.cov[i,13] <- all_res[[i]]$Results$seed[1]
  mean.cov[i,14] <- median(all_res[[i]]$Results[,13])
}
#get the info about N, n, a, and b
resnames <- strsplit(my_files, split = "_")
mean.cov$N <- 0
mean.cov$n <- NA
mean.cov$a <- NA
mean.cov$b <- NA
for (i in 1:NRow) {
  mean.cov[i,15] <- resnames[[i]][6] #N
  mean.cov[i,16] <- resnames[[i]][8] #n
  mean.cov[i,17] <- resnames[[i]][12] #a
  mean.cov[i,18] <- resnames[[i]][14] #b
}
#only take necessary info
mean.res <- mean.cov[,c(15,16,1,17,18,10:12,14,5,6,13)]
mean.res$N <- as.numeric(mean.res$N)
mean.res$n <- as.numeric(mean.res$n)
mean.res$a <- as.numeric(mean.res$a)
mean.res$b <- as.numeric(mean.res$b)

mean.jqm2 <- mean.res[,c(1:9)]
mean.jqm2$method <- "Chi-sq PJQM"

####load data resulted from JQM_Sims.R in Error_t3 PJQM directory
my_files <- list.files("./output/Error_t3 PJQM", pattern = "\\.Rdata$", full.names=TRUE)
NRow <- length(my_files)

all_res <- lapply(my_files, function(x) mget(load(x)))
mean.cov <- data.frame(matrix(nrow=NRow, ncol=length(colnames(all_res[[1]]$Results))))
colnames(mean.cov) <- colnames(all_res[[1]]$Results)
mean.cov$EmpCov.median <- NA
mean.cov <- mean.cov[,-1]
for (i in 1:NRow) {
  mean.cov[i,1] <- all_res[[i]]$Results$alpha[1]
  mean.cov[i,2:12] <- colMeans(all_res[[i]]$Results[,3:13])
  mean.cov[i,13] <- all_res[[i]]$Results$seed[1]
  mean.cov[i,14] <- median(all_res[[i]]$Results[,13])
}
#get the info about N, n, a, and b
resnames <- strsplit(my_files, split = "_")
mean.cov$N <- 0
mean.cov$n <- NA
mean.cov$a <- NA
mean.cov$b <- NA
for (i in 1:NRow) {
  mean.cov[i,15] <- resnames[[i]][7] #N
  mean.cov[i,16] <- resnames[[i]][9] #n
  mean.cov[i,17] <- resnames[[i]][13] #a
  mean.cov[i,18] <- resnames[[i]][15] #b
}
#only take necessary info
mean.res <- mean.cov[,c(15,16,1,17,18,10:12,14,5,6,13)]
mean.res$N <- as.numeric(mean.res$N)
mean.res$n <- as.numeric(mean.res$n)
mean.res$a <- as.numeric(mean.res$a)
mean.res$b <- as.numeric(mean.res$b)

mean.jqm3 <- mean.res[,c(1:9)]
mean.jqm3$method <- "Error t(3) PJQM"

####merge all PJQM sim results
mean.jqm.all <- rbind.data.frame(mean.jqm, mean.jqm2, mean.jqm3)
mean.jqm.all <- mean.jqm.all[,-9]


##################
#read LQMM simulation results stored in "./output/LQMM" directory
##################

#load data resulted from LQMM_Sims.R
load("./output/LQMM/mean.N.Rdata")
load("./output/LQMM/mean.chi.Rdata")

mean.lqmm.N <- mean.N[, c(1:5, 8:10)]
mean.lqmm.chi <- mean.all.chi[, c(1:5, 8:10)]
mean.lqmm.N$method <- "Normal LQMM"
mean.lqmm.chi$method <- "Chi-sq LQMM"

#combine all penalized JQM and LQMM results
all <- rbind.data.frame(mean.jqm.all, mean.lqmm.N, mean.lqmm.chi)

#create group for fixed a and b parameters
all$fixpar<- NA
for (i in 1:nrow(all)) {
  if(all$a[i]==2 & all$b[i]==0.4 | all$a[i]==0 & all$b[i]==0) {all$fixpar[i] = 1}
  else{all$fixpar[i] = 0}
}
#create group for six different settings
for (i in 1:nrow(all)) {
  if(all$N[i]==30 & all$alpha[i]==0.05) {all$group[i] = 1}
  if(all$N[i]==30 & all$alpha[i]==0.1) {all$group[i] = 2}
  if(all$N[i]==30 & all$alpha[i]==0.15) {all$group[i] = 3}
  if(all$N[i]==50 & all$alpha[i]==0.05) {all$group[i] = 4}
  if(all$N[i]==50 & all$alpha[i]==0.1) {all$group[i] = 5}
  if(all$N[i]==50 & all$alpha[i]==0.15) {all$group[i] = 6}
}
#change to factor
all$group <- as.factor(all$group)
all$N <- as.factor(all$N)
all$n <- as.factor(all$n)
all$alpha <- as.factor(all$alpha)

#create group for different distributions settings
for (i in 1:nrow(all)) {
  if(all$N[i]==30 & all$method[i]=="Normal LQMM") {all$group_alpha[i] = 1}
  if(all$N[i]==30 & all$method[i]=="Chi-sq LQMM") {all$group_alpha[i] = 3}
  if(all$N[i]==30 & all$method[i]=="Normal PJQM") {all$group_alpha[i] = 5}
  if(all$N[i]==30 & all$method[i]=="Chi-sq PJQM") {all$group_alpha[i] = 7}
  if(all$N[i]==30 & all$method[i]=="Error t(3) PJQM") {all$group_alpha[i] = 9}
  if(all$N[i]==50 & all$method[i]=="Normal LQMM") {all$group_alpha[i] = 2}
  if(all$N[i]==50 & all$method[i]=="Chi-sq LQMM") {all$group_alpha[i] = 4}
  if(all$N[i]==50 & all$method[i]=="Normal PJQM") {all$group_alpha[i] = 6}
  if(all$N[i]==50 & all$method[i]=="Chi-sq PJQM") {all$group_alpha[i] = 8}
  if(all$N[i]==50 & all$method[i]=="Error t(3) PJQM") {all$group_alpha[i] = 10}
}
#change to factor
all$group_alpha <- as.factor(all$group_alpha)

#save results in output directory
save(all, file="./output/all_sim_results.RData")

#The values of Table 1 and Table S2 are extracted from this 'all' data frame
