library(tidyverse)

# load simulation scenarios
m<-read.csv("./data/simulation/sim_scenarios.csv", sep=";", stringsAsFactors = T)

####load data resulted from JQM
my_files<-list.files("./output/output_final/PJQM", pattern = "\\.Rdata$", full.names=TRUE)
NRow<-length(my_files)

all_res<-lapply(my_files, function(x) mget(load(x)))
mean.cov<-data.frame(matrix(nrow=NRow, ncol=length(colnames(all_res[[1]]$Results))-3))
colnames(mean.cov)<-colnames(all_res[[1]]$Results)[1:15]
colnames(mean.cov)
mean.cov$EmpCov.median<-NA
for (i in 1:NRow) {
  mean.cov[i,1]<-nrow(all_res[[i]]$Results[is.na(all_res[[i]]$Results$beta0)==F,])
  mean.cov[i,2]<-all_res[[i]]$Results$alpha[1]
  mean.cov[i,3:14]<-colMeans(all_res[[i]]$Results[,3:14], na.rm = T)
  mean.cov[i,15]<-mean(all_res[[i]]$Results$seed, na.rm = T)
  mean.cov[i,16]<-median(all_res[[i]]$Results[,14], na.rm = T)
}
#get the info about N, n, a, and b
resnames<-strsplit(my_files, split = "_")
mean.cov$N<-NA
mean.cov$n<-NA
mean.cov$a<-NA
mean.cov$b<-NA
mean.cov$nr<-NA
for (i in 1:NRow) {
  mean.cov[i,17]<-resnames[[i]][8] #N
  mean.cov[i,18]<-resnames[[i]][10] #n
  mean.cov[i,19]<-resnames[[i]][14] #a
  mean.cov[i,20]<-strsplit(resnames[[i]][16], split = ".Rdata")[[1]] #b
  mean.cov[i,21]<-strsplit(resnames[[i]][2], split = "/")[[1]][3]
}

mean.cov$nr<-as.integer(mean.cov$nr)
mean.cov<-mean.cov %>% left_join(., m[,c(1,10)], by="nr")
mean.cov$method<-"JQM"
jqm.res<-mean.cov[,c("method","nr","FUN","alpha","N","n","a","b","EmpCov.subj","EmpCov.subj1obs","EmpCov.time","EmpCov","EmpCov.median")]


####load data resulted from RQPD 
my_files<-list.files("./output/output_final/RQPD", pattern = "\\.Rdata$", full.names=TRUE)
NRow<-length(my_files)

all_res<-lapply(my_files, function(x) mget(load(x)))
mean.cov<-data.frame(matrix(nrow=NRow, ncol=length(colnames(all_res[[1]]$Results))-3))
colnames(mean.cov)<-colnames(all_res[[1]]$Results)[1:13]
colnames(mean.cov)
mean.cov$EmpCov.median<-NA
for (i in 1:NRow) {
  mean.cov[i,1]<-nrow(all_res[[i]]$Results[is.na(all_res[[i]]$Results$beta01)==F,])
  mean.cov[i,2]<-all_res[[i]]$Results$alpha[1]
  mean.cov[i,3:12]<-colMeans(all_res[[i]]$Results[,3:12], na.rm = T)
  mean.cov[i,13]<-mean(all_res[[i]]$Results$seed, na.rm = T)
  mean.cov[i,14]<-median(all_res[[i]]$Results[,12], na.rm = T)
}
#get the info about N, n, a, and b
resnames<-strsplit(my_files, split = "_")
mean.cov$N<-NA
mean.cov$n<-NA
mean.cov$a<-NA
mean.cov$b<-NA
mean.cov$nr<-NA
for (i in 1:NRow) {
  mean.cov[i,15]<-resnames[[i]][8] #N
  mean.cov[i,16]<-resnames[[i]][10] #n
  mean.cov[i,17]<-resnames[[i]][14] #a
  mean.cov[i,18]<-strsplit(resnames[[i]][16], split = ".Rdata")[[1]] #b
  mean.cov[i,19]<-strsplit(resnames[[i]][2], split = "/")[[1]][3]
}

mean.cov$nr<-as.integer(mean.cov$nr)
mean.cov<-mean.cov %>% left_join(., m[,c(1,10)], by="nr")
mean.cov$method<-"RQPD"
rqpd.res<-mean.cov[,c("method","nr","FUN","alpha","N","n","a","b","EmpCov.subj","EmpCov.subj1obs","EmpCov.time","EmpCov","EmpCov.median")]


####load data resulted from LQMM
my_files<-list.files("./output/output_final/LQMM", pattern = "\\.Rdata$", full.names=TRUE)
NRow<-length(my_files)

all_res<-lapply(my_files, function(x) mget(load(x)))
mean.cov<-data.frame(matrix(nrow=NRow, ncol=length(colnames(all_res[[1]]$Results))-3))
colnames(mean.cov)<-colnames(all_res[[1]]$Results)[1:13]
colnames(mean.cov)
mean.cov$EmpCov.median<-NA
colnames(all_res[[i]]$Results)
for (i in 1:NRow) {
  mean.cov[i,1]<-nrow(all_res[[i]]$Results[is.na(all_res[[i]]$Results$beta01)==F,])
  mean.cov[i,4]<-all_res[[i]]$Results$alpha[1]
  mean.cov[i,c(2:3,5:12)]<-colMeans(all_res[[i]]$Results[,c(2:3,5:12)], na.rm = T)
  mean.cov[i,13]<-mean(all_res[[i]]$Results$seed, na.rm = T)
  mean.cov[i,14]<-median(all_res[[i]]$Results[,12], na.rm = T)
}
#get the info about N, n, a, and b
resnames<-strsplit(my_files, split = "_")
mean.cov$N<-NA
mean.cov$n<-NA
mean.cov$a<-NA
mean.cov$b<-NA
mean.cov$nr<-NA
for (i in 1:NRow) {
  mean.cov[i,]$N<-resnames[[i]][8] #N
  mean.cov[i,]$n<-resnames[[i]][10] #n
  mean.cov[i,]$a<-resnames[[i]][14] #a
  mean.cov[i,]$b<-strsplit(resnames[[i]][16], split = ".Rdata")[[1]] #b
  mean.cov[i,]$nr<-strsplit(resnames[[i]][2], split = "/")[[1]][3]
}
mean.cov$nr<-as.integer(mean.cov$nr)
mean.cov<-mean.cov %>% left_join(., m[,c(1,10)], by="nr")
mean.cov$method<-"LQMM"
lqmm.res<-mean.cov[,c("method","nr","FUN","alpha","N","n","a","b","EmpCov.subj","EmpCov.subj1obs","EmpCov.time","EmpCov","EmpCov.median")]


#### combine all data
all<-rbind(jqm.res, rqpd.res, lqmm.res)

# group by fixed a and b
all$fixpar<- NA
for (i in 1:nrow(all)) {
  if(all$a[i]==2 & all$b[i]==0.4 | all$a[i]==0 & all$b[i]==0) {all$fixpar[i] = 1}
  else{all$fixpar[i] = 0}
}

# group for six different settings
for (i in 1:nrow(all)) {
  if(all$N[i]==30 & all$alpha[i]==0.05) {all$group[i] = 1}
  if(all$N[i]==30 & all$alpha[i]==0.1) {all$group[i] = 2}
  if(all$N[i]==30 & all$alpha[i]==0.15) {all$group[i] = 3}
  if(all$N[i]==50 & all$alpha[i]==0.05) {all$group[i] = 4}
  if(all$N[i]==50 & all$alpha[i]==0.1) {all$group[i] = 5}
  if(all$N[i]==50 & all$alpha[i]==0.15) {all$group[i] = 6}
}
#change to factor
all$group<-as.factor(all$group)
all$N<-as.factor(all$N)
all$n<-as.numeric(all$n)
all$alpha<-as.factor(all$alpha)

for (i in 1:nrow(all)) {
  if(all$N[i]==30 & all$method[i]=="JQM" & all$FUN[i]=="SimData_NvarN") {all$group_alpha[i] = 1}
  if(all$N[i]==30 & all$method[i]=="JQM" & all$FUN[i]=="SimData_Nvarchi") {all$group_alpha[i] = 2}
  if(all$N[i]==30 & all$method[i]=="JQM" & all$FUN[i]=="SimData_NN") {all$group_alpha[i] = 3}
  if(all$N[i]==30 & all$method[i]=="JQM" & all$FUN[i]=="SimData_Nchi") {all$group_alpha[i] = 4}
  if(all$N[i]==30 &  all$method[i]=="JQM" & all$FUN[i]=="SimData_tN") {all$group_alpha[i] = 5}
  if(all$N[i]==30 &  all$method[i]=="JQM" & all$FUN[i]=="SimData_tchi") {all$group_alpha[i] = 6}
  if(all$N[i]==30 & all$method[i]=="LQMM" & all$FUN[i]=="SimData_NvarN") {all$group_alpha[i] = 7}
  if(all$N[i]==30 & all$method[i]=="LQMM" & all$FUN[i]=="SimData_Nvarchi") {all$group_alpha[i] = 8}
  if(all$N[i]==30 & all$method[i]=="LQMM" & all$FUN[i]=="SimData_NN") {all$group_alpha[i] = 9}
  if(all$N[i]==30 & all$method[i]=="LQMM" & all$FUN[i]=="SimData_Nchi") {all$group_alpha[i] = 10}
  if(all$N[i]==30 &  all$method[i]=="LQMM" & all$FUN[i]=="SimData_tN") {all$group_alpha[i] = 11}
  if(all$N[i]==30 &  all$method[i]=="LQMM" & all$FUN[i]=="SimData_tchi") {all$group_alpha[i] = 12}
  
  if(all$N[i]==30 & all$method[i]=="RQPD" & all$FUN[i]=="SimData_NvarN") {all$group_alpha[i] = 13}
  if(all$N[i]==30 & all$method[i]=="RQPD" & all$FUN[i]=="SimData_Nvarchi") {all$group_alpha[i] = 14}
  if(all$N[i]==30 & all$method[i]=="RQPD" & all$FUN[i]=="SimData_NN") {all$group_alpha[i] = 15}
  if(all$N[i]==30 & all$method[i]=="RQPD" & all$FUN[i]=="SimData_Nchi") {all$group_alpha[i] = 16}
  if(all$N[i]==30 &  all$method[i]=="RQPD" & all$FUN[i]=="SimData_tN") {all$group_alpha[i] = 17}
  if(all$N[i]==30 &  all$method[i]=="RQPD" & all$FUN[i]=="SimData_tchi") {all$group_alpha[i] = 18}
  
  if(all$N[i]==50 & all$method[i]=="JQM" & all$FUN[i]=="SimData_NvarN") {all$group_alpha[i] = 19}
  if(all$N[i]==50 & all$method[i]=="JQM" & all$FUN[i]=="SimData_Nvarchi") {all$group_alpha[i] = 20}
  if(all$N[i]==50 & all$method[i]=="JQM" & all$FUN[i]=="SimData_NN") {all$group_alpha[i] = 21}
  if(all$N[i]==50 & all$method[i]=="JQM" & all$FUN[i]=="SimData_Nchi") {all$group_alpha[i] = 22}
  if(all$N[i]==50 &  all$method[i]=="JQM" & all$FUN[i]=="SimData_tN") {all$group_alpha[i] = 23}
  if(all$N[i]==50 &  all$method[i]=="JQM" & all$FUN[i]=="SimData_tchi") {all$group_alpha[i] = 24}
  
  if(all$N[i]==50 & all$method[i]=="LQMM" & all$FUN[i]=="SimData_NvarN") {all$group_alpha[i] = 25}
  if(all$N[i]==50 & all$method[i]=="LQMM" & all$FUN[i]=="SimData_Nvarchi") {all$group_alpha[i] = 26}
  if(all$N[i]==50 & all$method[i]=="LQMM" & all$FUN[i]=="SimData_NN") {all$group_alpha[i] = 27}
  if(all$N[i]==50 & all$method[i]=="LQMM" & all$FUN[i]=="SimData_Nchi") {all$group_alpha[i] = 28}
  if(all$N[i]==50 &  all$method[i]=="LQMM" & all$FUN[i]=="SimData_tN") {all$group_alpha[i] = 29}
  if(all$N[i]==50 &  all$method[i]=="LQMM" & all$FUN[i]=="SimData_tchi") {all$group_alpha[i] = 30}
  
  if(all$N[i]==50 & all$method[i]=="RQPD" & all$FUN[i]=="SimData_NvarN") {all$group_alpha[i] = 31}
  if(all$N[i]==50 & all$method[i]=="RQPD" & all$FUN[i]=="SimData_Nvarchi") {all$group_alpha[i] = 32}
  if(all$N[i]==50 & all$method[i]=="RQPD" & all$FUN[i]=="SimData_NN") {all$group_alpha[i] = 33}
  if(all$N[i]==50 & all$method[i]=="RQPD" & all$FUN[i]=="SimData_Nchi") {all$group_alpha[i] = 34}
  if(all$N[i]==50 &  all$method[i]=="RQPD" & all$FUN[i]=="SimData_tN") {all$group_alpha[i] = 35}
  if(all$N[i]==50 &  all$method[i]=="RQPD" & all$FUN[i]=="SimData_tchi") {all$group_alpha[i] = 36}
  
  
}

save(all, file="./output/output_final/all_sim.Rdata")