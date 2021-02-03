#call input scenarios
#set the working directory
p <- read.csv("./data/input_LQMM_N.csv", header = T, sep = ";")
q <- read.csv("./data/input_LQMM_chi.csv", header = T, sep = ";")

#object initializations
mean.N <- matrix(nrow=24, ncol=10)
mean.N <- as.data.frame(mean.all.N)
colnames(mean.N) <- c("N", "n", "alpha", "a", "b", "beta01", "beta02",
                      "EmpCov.subj","EmpCov.time","EmpCov")
mean.all.chi <- matrix(nrow=24, ncol=10)
mean.all.chi <- as.data.frame(mean.all.chi)
colnames(mean.all.chi) <- c("N", "n", "alpha", "a", "b", "beta01", "beta02",
                            "EmpCov.subj","EmpCov.time","EmpCov")

#call simulation for separate-IRI with N(0,1) and chi-sq(4) random effects distribution
for (r in 1:24) {
  est1 <- lqmm_call(N=p[r,1], n=p[r,2], alpha=p[r,3], beta=p[r,4], NSim=p[r,5], psiu=p[r,6], 
                    a=p[r,7], b=p[r,8], seed=p[r,9], FUN=sim_N, dist="Normal")
  mean.N[r,1:10] <- est1
  est2 <- lqmm_call(N=q[r,1], n=q[r,2], alpha=q[r,3], beta=q[r,4], NSim=q[r,5], psiu=q[r,6], 
                    a=q[r,7], b=q[r,8], seed=q[r,9], FUN=sim_N_chi, dist="Chi-sq")
  mean.chi[r,1:10] <- est2
}

save(mean.N, file = "mean.N.Rdata")
save(mean.chi, file = "mean.chi.Rdata")