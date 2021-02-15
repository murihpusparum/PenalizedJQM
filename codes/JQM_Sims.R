setwd("./codes")
#call all functions
source("JQM_Function.R")

setwd("../")
#set the seed and call SimStudy function
seed=123
set.seed(seed)
SimStudy(N=30, n=10, beta=0, a=2.8, b=0.8, psiu2=1, alpha=0.05, NSim=1, seed=seed)