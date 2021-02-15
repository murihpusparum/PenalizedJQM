#produce Figure 1, Table S2-S3, of simulation results

library(ggplot2)
###################
#read Penalized JQM simulation results stored in "./output/PJQM" directory
###################
#load data resulted from JQM_Sims.R
my_files <- list.files("./output/PJQM", pattern = "\\.Rdata$", full.names=TRUE)
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

mean.jqm <- mean.res[,c(1:8)]
mean.jqm$method <- "Penalized JQM"

##################
#read LQMM simulation results stored in "./output/LQMM" directory
##################
#load data resulted from LQMM_Sims.R
load("./output/LQMM/mean.N.Rdata")
load("./output/LQMM/mean.chi.Rdata")

mean.lqmm.N <- mean.N[, c(1:5, 8:10)]
mean.lqmm.chi <- mean.chi[, c(1:5, 8:10)]
mean.lqmm.N$method <- "Normal LQMM"
mean.lqmm.chi$method <- "Chi-sq LQMM"

#combine all penalized JQM and LQMM results
all <- rbind.data.frame(mean.jqm, mean.lqmm.N, mean.lqmm.chi)

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
all$EmpCov <- round(all$EmpCov, digits=2)


##################
#Figure 1
##################

#create theme for Figure 1a-1c
My_Theme = theme(
  axis.title.x = element_text(size = 12),
  axis.text.x = element_text(size = 10),
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size = 12),
  panel.border = element_rect(size=1, color="grey"),
  legend.position = "bottom",
  legend.title = element_text(size=12),
  legend.text = element_text(size=11))

#Figure 1a
all1 <- all[all$fixpar==1 & all$method=="Normal LQMM",]
lqmm1 <- ggplot(all1, aes(x=n))+
  geom_line(aes(y=EmpCov, group=group, linetype=group, color=group), size=1)+
  geom_point(aes(x=n, y=EmpCov, color=group), size=2)+
  scale_linetype_manual(name="scenario", values=c("dashed", "solid", "dotted",
                                                  "dashed", "solid", "dotted"),
                        labels=c("95% RI - N30", "90% RI - N30", "85% RI - N30",
                                 "95% RI - N50", "90% RI - N50", "85% RI - N50")) +
  scale_color_manual(name="scenario", values=c("#F5C000", "#F5C000", "#F5C000",
                                               "#34A3DC", "#34A3DC", "#34A3DC"),
                     labels=c("95% RI - N30", "90% RI - N30", "85% RI - N30",
                              "95% RI - N50", "90% RI - N50", "85% RI - N50")) +
  scale_y_continuous(limits = c(0.65, 0.95), breaks = seq(0.65, 0.95, by = 0.05)) + 
  labs(y="OEC") +
  theme_bw() + My_Theme

#Figure 1b
all2 <- all[all$fixpar==1 & all$method=="Chi-sq LQMM",]
lqmm2 <- ggplot(all2, aes(x=n))+
  geom_line(aes(y=EmpCov, group=group, linetype=group, color=group), size=1)+
  geom_point(aes(x=n, y=EmpCov, color=group), size=2)+
  scale_linetype_manual(name="scenario", values=c("dashed", "solid", "dotted",
                                                  "dashed", "solid", "dotted"),
                        labels=c("95% RI - N30", "90% RI - N30", "85% RI - N30",
                                 "95% RI - N50", "90% RI - N50", "85% RI - N50")) +
  scale_color_manual(name="scenario", values=c("#F5C000", "#F5C000", "#F5C000",
                                               "#34A3DC", "#34A3DC", "#34A3DC"),
                     labels=c("95% RI - N30", "90% RI - N30", "85% RI - N30",
                              "95% RI - N50", "90% RI - N50", "85% RI - N50")) +
  scale_y_continuous(limits = c(0.65, 0.95), breaks = seq(0.65, 0.95, by = 0.05)) + 
  labs(y="OEC") +
  theme_bw() + My_Theme

#Figure 1c
all3 <- all[all$fixpar==1 & all$method=="Penalized JQM",]
jqm <- ggplot(all3, aes(x=n))+
  geom_line(aes(y=EmpCov, group=group, linetype=group, color=group), size=1)+
  geom_point(aes(x=n, y=EmpCov, color=group), size=2)+
  scale_linetype_manual(name="scenario", values=c("dashed", "solid", "dotted",
                                                  "dashed", "solid", "dotted"),
                        labels=c("95% RI - N30", "90% RI - N30", "85% RI - N30",
                                 "95% RI - N50", "90% RI - N50", "85% RI - N50")) +
  scale_color_manual(name="scenario", values=c("#F5C000", "#F5C000", "#F5C000",
                                               "#34A3DC", "#34A3DC", "#34A3DC"),
                     labels=c("95% RI - N30", "90% RI - N30", "85% RI - N30",
                              "95% RI - N50", "90% RI - N50", "85% RI - N50")) +
  scale_y_continuous(limits = c(0.65, 0.95), breaks = seq(0.65, 0.95, by = 0.05)) + 
  labs(y="OEC") +
  theme_bw() + My_Theme

#Figure 1d
box <- ggplot(all, aes(x=method, y=EmpCov, color=method))+
  geom_boxplot(size=1)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  scale_color_manual(name="method", values=c("#F5C000", "mediumaquamarine", "#34A3DC"),
                     labels=c("LQMM Chi-sq(4)", "LQMM N(0, 1)", "Penalized JQM")) +
  labs(y="OEC") +
  theme_bw() + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    panel.border = element_rect(size=1, color="grey"),
    legend.position = "none" )