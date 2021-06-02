#call all functions
source("./codes/JQM_Function.R")

#set the seed and call SimStudy function
#seed can be found in ./data/PJQM_scenarios.csv
seed=123
set.seed(seed)
SimStudy(N=30, n=10, beta=0, a=2.8, b=0.8, psiu2=1, alpha=0.05, NSim=1, seed=seed)