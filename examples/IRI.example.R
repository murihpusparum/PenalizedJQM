source("./codes/simdata.R")
source("./codes/LQMM_Function.R")
source("./codes/JQM_Function.R")
library(ggplot2)

# data generation
# data must not result into warning in LQMM method
seed<-4581621 # random seed
res.lqmm<-lqmm.fun(alpha=0.05, N=30, n=10, beta=0, a=2, b=0.4, seed=seed, FUN=SimData_NvarN) 

# data overview
df<-res.lqmm$df
head(df)


####################################
# to get the IRIs from Penalized JQM
####################################
# run the JQM estimation procedure
res<-jqm(db=df,
         alpha=0.05,
         lambda.u.seq = seq(0.5,4,0.5),
         lambda.z.seq = seq(0.5,5,0.5))

iri.jqm<-cbind.data.frame(res$beta0+res$u+res$z*res$beta1, 
                          res$beta0+res$u+res$z*res$beta2, 
                          seq(1:30))
colnames(iri.jqm)<-c("low","up","id")

jqm.plot <- ggplot(iri.jqm) +
  geom_errorbar(aes(x=as.factor(id), ymin = low, ymax = up), color="darkblue", 
                position=position_dodge(width=0.7)) +
  geom_point(data = df, aes(x = subject, y = y), color="darkred", 
             position=position_dodge(width=0.7), size=1.5) +
  geom_vline(xintercept=seq(1.5, length(unique(df$subject))-0.5, 1), 
             lwd=0.5, colour="grey") +
  labs(x="subjects", y = "measurements", title = "IRI with Penalized JQM") +
  theme_classic()
jqm.plot

# empirical coverage
res$cov

####################################
# to get the IRIs from LQMM
####################################
iri.lqmm<-cbind.data.frame(res.lqmm$beta01+res.lqmm$u1, 
                           res.lqmm$beta02+res.lqmm$u2, 
                          seq(1:30))
colnames(iri.lqmm)<-c("low","up","id")

lqmm.plot <- ggplot(iri.lqmm) +
  geom_errorbar(aes(x=as.factor(id), ymin = low, ymax = up), color="darkblue", 
                position=position_dodge(width=0.7)) +
  geom_point(data = df, aes(x = subject, y = y), color="darkred", 
             position=position_dodge(width=0.7), size=1.5) +
  geom_vline(xintercept=seq(1.5, length(unique(df$subject))-0.5, 1), 
             lwd=0.5, colour="grey") +
  labs(x="subjects", y = "measurements", title = "IRI with LQMM") +
  theme_classic()
lqmm.plot

# empirical coverage
res.lqmm$cov
