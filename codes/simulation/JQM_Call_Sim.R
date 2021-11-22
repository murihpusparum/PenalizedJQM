#call all functions
#also provide simdata.R and generate_data.R in the same directory
source("./codes/main/JQM_Function.R")
source("./codes/simulation/simdata_sd.R")
source("./codes/simulation/generate_data_sd.R")
source("./codes/simulation/JQM_Sim_Function.R")
SimStudy(N=30, n=10, beta=0, a=2.8, b=0.8, alpha=0.05, NSim=2, seed=6892, FUN=SimData_NN)
