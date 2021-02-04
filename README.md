# PenalizedJQM
# Requirements
-To compile all the codes, R 4.0 or later is recommended.
-To get the latest R update:
	1. Install and load from R Gui installr: install.packages("installr") and library(installr)
	2. Call update(R)
	3. From within RStudio, run R.version.string to get detailed information about the version of R running
-Go to "./codes/install.package.R" and run the codes for installing all necessary packages.

# Simulation study
-To run the LQMM simulation, the LQMM_Function.R and the LQMM_Sims.R must be in a same local directory.
-The input files can be found in "./data" directory.
-Run the LQMM_Sims.R for the LQMM simulation, the output of the results will be stored in "./output/LQMM" directory.
-To run the Penalized JQM simulation, the JQM_Function.R and the JQM_Sims.R must be in a same local directory.
-Run the JQM_Sims.R: the list of sample size and parameter scenarios with their corresponding seeds can be found in "./data/PJQM_scenarios.csv".
-The output of the results will be stored in "./output/PJQM" directory.

# Read simulation results
-To read the simulation results and produce Figure 1, Table S2, and Table S3, run the read.sims.result.R

# Case Study: IAM Frontier Data
-To produce the case study results and figures, VITO IAM Frontier data is needed. The data is not provided in this repository because of confidentiality reason, please ask the authors for further information.
-To produce the IRI plots in Figure 2 and 4, IAM.IRI.plots.R and JQM_Function.R must be in a same local directory. Run IAM.IRI.plots.R.
-To produce Table 1 and Table S4, IAM.allpar.estimates.R and JQM_Function.R must be in a same local directory. Run IAM.allpar.estimates.R.
-IAM.circos.plot.R is used to produce the circos plot in Figure 3
-IAM.cosine.similarity is used to produce the cosine similarity matrix and Figure S1.
