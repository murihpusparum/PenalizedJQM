source("./codes/main/JQM_Function.R")
db<-read.table("./data/result_edge.txt", header = T)
res<-jqm(db=db,alpha=0.05)
res