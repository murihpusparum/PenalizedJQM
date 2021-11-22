# load from here: load("./output/output_final/")
library(reshape2)
library(ggplot2)
library(ggpubr)

len.all<-len.all[-1,]
len.wide1<-dcast(len.all, iter+ys ~ method, value.var= "length")
len.wide2<-dcast(len.all, iter+ys ~ method, value.var="coverage")
len.wide<-cbind(len.wide1, len.wide2)
len.wide<-len.wide[,-c(6:7)]
colnames(len.wide)<-c(colnames(len.wide)[1:5],"JQM_cov","LQMM_cov","RQPD_cov")

for (i in 1:nrow(len.wide)) {
  if(len.wide$JQM_cov[i]==1 & len.wide$RQPD_cov[i]==0) {len.wide$cov[i]="green"}
  if(len.wide$JQM_cov[i]==0 & len.wide$RQPD_cov[i]==1) {len.wide$cov[i]="red"}
  if(len.wide$JQM_cov[i]==1 & len.wide$RQPD_cov[i]==1) {len.wide$cov[i]="black"}
  if(len.wide$JQM_cov[i]==0 & len.wide$RQPD_cov[i]==0) {len.wide$cov[i]="black"}
}

for (i in 1:nrow(len.wide)) {
  if(len.wide$JQM_cov[i]==1 & len.wide$LQMM_cov[i]==0) {len.wide$cov2[i]="green"}
  if(len.wide$JQM_cov[i]==0 & len.wide$LQMM_cov[i]==1) {len.wide$cov2[i]="red"}
  if(len.wide$JQM_cov[i]==1 & len.wide$LQMM_cov[i]==1) {len.wide$cov2[i]="black"}
  if(len.wide$JQM_cov[i]==0 & len.wide$LQMM_cov[i]==0) {len.wide$cov2[i]="black"}
}

jqm.rqpd<-ggplot(len.wide, aes(x=as.factor(iter), y=JQM, color=cov))+
  geom_boxplot(outlier.shape = NA, color="black")+
  geom_jitter(aes(size=cov),shape=16, position=position_jitter(0.2))+
  geom_point(aes(y=RQPD), color="blue", size=3)+
  labs(y="IRI length",x="Data")+
  theme_bw() + 
  scale_color_manual(values = c("gray30","springgreen3","orange"))+
  scale_size_manual(values = c(1,2,2))+
  theme(legend.position = "none",
        text = element_text(size=16),
        axis.text.y = element_text(face = "bold"))


jqm.lqmm<-ggplot(len.wide, aes(x=LQMM, y=JQM, color=cov2))+
  geom_point(size=2)+
  geom_abline(slope=1, intercept = 0, linetype="dashed")+
  labs(y="Joint-IRI length",x="Separate-IRI length")+
  scale_color_manual(values = c("gray30","springgreen3","orange"))+
  theme_bw() + 
  theme(legend.position = "none",
        text = element_text(size=16),
        axis.text.y = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold"))

figure2 <- ggarrange(jqm.rqpd, jqm.lqmm,
                         common.legend = FALSE, legend = "none", align = c("hv"),
                         font.label = list(size=12), ncol = 2, nrow = 1)
figure2