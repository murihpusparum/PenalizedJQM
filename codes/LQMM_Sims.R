source("LQMM_Function.R")

#call input scenarios
#set the working directory
p <- read.csv("./data/input_LQMM_N.csv", header = T, sep = ";")
q <- read.csv("./data/input_LQMM_chi.csv", header = T, sep = ";")

#object initializations
mean.N2 <- matrix(nrow=24, ncol=11)
mean.N2 <- as.data.frame(mean.N2)
colnames(mean.N2) <- c("N", "n", "alpha", "a", "b", "beta01", "beta02",
                      "EmpCov.subj","EmpCov.time","EmpCov", "EmpCov.median")
mean.all.chi2 <- matrix(nrow=24, ncol=11)
mean.all.chi2 <- as.data.frame(mean.all.chi2)
colnames(mean.all.chi2) <- c("N", "n", "alpha", "a", "b", "beta01", "beta02",
                            "EmpCov.subj","EmpCov.time","EmpCov", "EmpCov.median")

#call simulation for separate-IRI with N(0,1) and chi-sq(4) random effects distribution
for (r in 1:24) {
  est1 <- lqmm_call(N=p[r,1], n=p[r,2], alpha=p[r,3], beta=p[r,4], NSim=p[r,5], psiu=p[r,6], 
                    a=p[r,7], b=p[r,8], seed=p[r,9], FUN=sim_N, dist="Normal")
  mean.N2[r,1:11] <- est1
  est2 <- lqmm_call(N=q[r,1], n=q[r,2], alpha=q[r,3], beta=q[r,4], NSim=q[r,5], psiu=q[r,6], 
                    a=q[r,7], b=q[r,8], seed=q[r,9], FUN=sim_N_chi, dist="Chi-sq")
  mean.all.chi2[r,1:11] <- est2
}

save(mean.N2, file = "./output/LQMM/mean.N.Rdata")
save(mean.all.chi2, file = "./output/LQMM/mean.chi.Rdata")
