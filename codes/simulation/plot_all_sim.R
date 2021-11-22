#load("./output/output_final/all_sim.Rdata")
# filter: N=30, alpha=0.05, fixpar=1
all<-all[order(all$n),]
all.use<-all[all$alpha==0.05&all$N==30&all$fixpar==1,]

n<-all.use[all.use$group_alpha==1,]$n

jqm1<-all.use[all.use$group_alpha==1,]$EmpCov.time
jqm2<-all.use[all.use$group_alpha==2,]$EmpCov.time
jqm3<-all.use[all.use$group_alpha==3,]$EmpCov.time
jqm4<-all.use[all.use$group_alpha==4,]$EmpCov.time

lqmm1<-all.use[all.use$group_alpha==7,]$EmpCov.time
lqmm2<-all.use[all.use$group_alpha==8,]$EmpCov.time
lqmm3<-all.use[all.use$group_alpha==9,]$EmpCov.time
lqmm4<-all.use[all.use$group_alpha==10,]$EmpCov.time
rqpd1<-all.use[all.use$group_alpha==13,]$EmpCov.time
rqpd2<-all.use[all.use$group_alpha==14,]$EmpCov.time
rqpd3<-all.use[all.use$group_alpha==15,]$EmpCov.time
rqpd4<-all.use[all.use$group_alpha==16,]$EmpCov.time

par(mfrow = c(2, 4), oma=c(4,0,0,0))

plot(n,jqm1,type="b",col="darkblue", lwd=2, lty=6, pch=19, xlab="n", ylab="TEC", ylim=c(0.85,1),
     cex.lab=1.5,cex.axis=1.5)
lines(n,jqm2,type="b",col="darkblue", lwd=2, lty=1, pch=19, xlab="n", ylab="TEC", ylim=c(0.85,1))
lines(n,jqm3,type="b",col="darkblue", lwd=2, lty=2, pch=19, xlab="n", ylab="TEC", ylim=c(0.85,1))
lines(n,jqm4,type="b",col="darkblue", lwd=2, lty=3, pch=19, xlab="n", ylab="TEC", ylim=c(0.85,1))
title("TEC, N=30", line=1, adj=0,cex.main=1.5)


plot(n,lqmm1,type="b",col="#34A3DC", lwd=2, lty=6, pch=19, xlab="n", ylab="", ylim=c(0.85,1),
     cex.lab=1.5,cex.axis=1.5)
lines(n,lqmm2,type="b",col="#34A3DC", lwd=2, lty=1, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))
lines(n,lqmm3,type="b",col="#34A3DC", lwd=2, lty=2, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))
lines(n,lqmm4,type="b",col="#34A3DC", lwd=2, lty=3, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))
lines(n,rqpd1,type="b",col="#F5C000", lwd=2, lty=6, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))
lines(n,rqpd2,type="b",col="#F5C000", lwd=2, lty=1, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))
lines(n,rqpd3,type="b",col="#F5C000", lwd=2, lty=2, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))
lines(n,rqpd4,type="b",col="#F5C000", lwd=2, lty=3, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))


jqm1<-all.use[all.use$group_alpha==1,]$EmpCov.subj
jqm2<-all.use[all.use$group_alpha==2,]$EmpCov.subj
jqm3<-all.use[all.use$group_alpha==3,]$EmpCov.subj
jqm4<-all.use[all.use$group_alpha==4,]$EmpCov.subj

lqmm1<-all.use[all.use$group_alpha==7,]$EmpCov.subj
lqmm2<-all.use[all.use$group_alpha==8,]$EmpCov.subj
lqmm3<-all.use[all.use$group_alpha==9,]$EmpCov.subj
lqmm4<-all.use[all.use$group_alpha==10,]$EmpCov.subj
rqpd1<-all.use[all.use$group_alpha==13,]$EmpCov.subj
rqpd2<-all.use[all.use$group_alpha==14,]$EmpCov.subj
rqpd3<-all.use[all.use$group_alpha==15,]$EmpCov.subj
rqpd4<-all.use[all.use$group_alpha==16,]$EmpCov.subj

plot(n,jqm1,type="b",col="darkblue", lwd=2, lty=6, pch=19, xlab="n", ylab="SEC", ylim=c(0.85,1),
     cex.lab=1.5,cex.axis=1.5)
lines(n,jqm2,type="b",col="darkblue", lwd=2, lty=1, pch=19, xlab="n", ylab="SEC", ylim=c(0.85,1))
lines(n,jqm3,type="b",col="darkblue", lwd=2, lty=2, pch=19, xlab="n", ylab="SEC", ylim=c(0.85,1))
lines(n,jqm4,type="b",col="darkblue", lwd=2, lty=3, pch=19, xlab="n", ylab="SEC", ylim=c(0.85,1))
title("SEC, N=30", line=1, adj=0,cex.main=1.5)


plot(n,lqmm1,type="b",col="#34A3DC", lwd=2, lty=6, pch=19, xlab="n", ylab="", ylim=c(0.85,1),
     cex.lab=1.5,cex.axis=1.5)
lines(n,lqmm2,type="b",col="#34A3DC", lwd=2, lty=1, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))
lines(n,lqmm3,type="b",col="#34A3DC", lwd=2, lty=2, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))
lines(n,lqmm4,type="b",col="#34A3DC", lwd=2, lty=3, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))
lines(n,rqpd1,type="b",col="#F5C000", lwd=2, lty=6, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))
lines(n,rqpd2,type="b",col="#F5C000", lwd=2, lty=1, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))
lines(n,rqpd3,type="b",col="#F5C000", lwd=2, lty=2, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))
lines(n,rqpd4,type="b",col="#F5C000", lwd=2, lty=3, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))


# filter: N=50, alpha=0.05, fixpar=1
all<-all[order(all$n),]
all.use<-all[all$alpha==0.05&all$N==50&all$fixpar==1,]

n<-all.use[all.use$group_alpha==19,]$n

jqm1<-all.use[all.use$group_alpha==19,]$EmpCov.time
jqm2<-all.use[all.use$group_alpha==20,]$EmpCov.time
jqm3<-all.use[all.use$group_alpha==21,]$EmpCov.time
jqm4<-all.use[all.use$group_alpha==22,]$EmpCov.time

lqmm1<-all.use[all.use$group_alpha==25,]$EmpCov.time
lqmm2<-all.use[all.use$group_alpha==26,]$EmpCov.time
lqmm3<-all.use[all.use$group_alpha==27,]$EmpCov.time
lqmm4<-all.use[all.use$group_alpha==28,]$EmpCov.time
rqpd1<-all.use[all.use$group_alpha==31,]$EmpCov.time
rqpd2<-all.use[all.use$group_alpha==32,]$EmpCov.time
rqpd3<-all.use[all.use$group_alpha==33,]$EmpCov.time
rqpd4<-all.use[all.use$group_alpha==34,]$EmpCov.time


plot(n,jqm1,type="b",col="darkblue", lwd=2, lty=6, pch=19, xlab="n", ylab="TEC", ylim=c(0.85,1),
     cex.lab=1.5,cex.axis=1.5)
lines(n,jqm2,type="b",col="darkblue", lwd=2, lty=1, pch=19, xlab="n", ylab="TEC", ylim=c(0.85,1))
lines(n,jqm3,type="b",col="darkblue", lwd=2, lty=2, pch=19, xlab="n", ylab="TEC", ylim=c(0.85,1))
lines(n,jqm4,type="b",col="darkblue", lwd=2, lty=3, pch=19, xlab="n", ylab="TEC", ylim=c(0.85,1))
title("TEC, N=50", line=1, adj=0,cex.main=1.5)

plot(n,lqmm1,type="b",col="#34A3DC", lwd=2, lty=6, pch=19, xlab="n", ylab="", ylim=c(0.85,1),
     cex.lab=1.5,cex.axis=1.5)
lines(n,lqmm2,type="b",col="#34A3DC", lwd=2, lty=1, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))
lines(n,lqmm3,type="b",col="#34A3DC", lwd=2, lty=2, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))
lines(n,lqmm4,type="b",col="#34A3DC", lwd=2, lty=3, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))
lines(n,rqpd1,type="b",col="#F5C000", lwd=2, lty=6, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))
lines(n,rqpd2,type="b",col="#F5C000", lwd=2, lty=1, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))
lines(n,rqpd3,type="b",col="#F5C000", lwd=2, lty=2, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))
lines(n,rqpd4,type="b",col="#F5C000", lwd=2, lty=3, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))

legend(0.6,0.80,c("joint-IRI","separate-IRI","joint FE-IRI"), col=c("darkblue","#34A3DC","#F5C000"), pch=19,bty="n",pt.cex=2,
       cex=1.5,xpd="NA")

legend(12,0.80,c(expression(paste(epsilon," ~ N(0,", theta[i],"), u ~ N(0,1)")),
                 expression(paste(epsilon," ~ N(0,", theta[i],"), u ~ ", chi[(4)]^{2}))), 
       col=c(rep("gray30",2)), lty=c(6,1),bty="n",lwd=2,
       cex=1.5,xpd="NA")

legend(32,0.80,c(expression(paste(epsilon," ~ N(0,1), u ~ N(0,1)")),
                 expression(paste(epsilon," ~ N(0,1), u ~ ", chi[(4)]^{2}))), 
       col=c(rep("gray30",2)), lty=c(2,3),bty="n",lwd=2,
       cex=1.5,xpd="NA")

jqm1<-all.use[all.use$group_alpha==19,]$EmpCov.subj
jqm2<-all.use[all.use$group_alpha==20,]$EmpCov.subj
jqm3<-all.use[all.use$group_alpha==21,]$EmpCov.subj
jqm4<-all.use[all.use$group_alpha==22,]$EmpCov.subj

lqmm1<-all.use[all.use$group_alpha==25,]$EmpCov.subj
lqmm2<-all.use[all.use$group_alpha==26,]$EmpCov.subj
lqmm3<-all.use[all.use$group_alpha==27,]$EmpCov.subj
lqmm4<-all.use[all.use$group_alpha==28,]$EmpCov.subj
rqpd1<-all.use[all.use$group_alpha==31,]$EmpCov.subj
rqpd2<-all.use[all.use$group_alpha==32,]$EmpCov.subj
rqpd3<-all.use[all.use$group_alpha==33,]$EmpCov.subj
rqpd4<-all.use[all.use$group_alpha==34,]$EmpCov.subj


plot(n,jqm1,type="b",col="darkblue", lwd=2, lty=6, pch=19, xlab="n", ylab="SEC", ylim=c(0.85,1),
     cex.lab=1.5,cex.axis=1.5)
lines(n,jqm2,type="b",col="darkblue", lwd=2, lty=1, pch=19, xlab="n", ylab="SEC", ylim=c(0.85,1))
lines(n,jqm3,type="b",col="darkblue", lwd=2, lty=2, pch=19, xlab="n", ylab="SEC", ylim=c(0.85,1))
lines(n,jqm4,type="b",col="darkblue", lwd=2, lty=3, pch=19, xlab="n", ylab="SEC", ylim=c(0.85,1))
title("SEC, N=50", line=1, adj=0,cex.main=1.5)


plot(n,lqmm1,type="b",col="#34A3DC", lwd=2, lty=6, pch=19, xlab="n", ylab="", ylim=c(0.85,1),
     cex.lab=1.5,cex.axis=1.5)
lines(n,lqmm2,type="b",col="#34A3DC", lwd=2, lty=1, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))
lines(n,lqmm3,type="b",col="#34A3DC", lwd=2, lty=2, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))
lines(n,lqmm4,type="b",col="#34A3DC", lwd=2, lty=3, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))
lines(n,rqpd1,type="b",col="#F5C000", lwd=2, lty=6, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))
lines(n,rqpd2,type="b",col="#F5C000", lwd=2, lty=1, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))
lines(n,rqpd3,type="b",col="#F5C000", lwd=2, lty=2, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))
lines(n,rqpd4,type="b",col="#F5C000", lwd=2, lty=3, pch=19, xlab="n",  ylab="", ylim=c(0.85,1))