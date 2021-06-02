#### produce Figure 2, Figure 3, Figure 4, Figure S3, and Figure 5

library(ggplot2)
library(dplyr)
library(ggpubr)

my_files <- list.files("./output/IAM", pattern = ".RData$", full.names=FALSE)
resnames <- strsplit(my_files, split = "\\.")

NRow <- length(my_files)
all_res <- lapply(my_files, function(x) mget(load(x)))

datalist<-list()
uz.list <- list()
for (i in seq(1,33, by=3)) {
  all_res[[i]]$df.a$par <- resnames[[i]][1]
  datalist[[i]] <- all_res[[i]]$df.a
  
  all_res[[i+1]]$uz.a$par <- resnames[[i+1]][1]
  all_res[[i+1]]$uz.a$method <- "PJQM"
  
  all_res[[i+2]]$uz.a$par <- resnames[[i+2]][1]
  all_res[[i+2]]$uz.a$method <- "LQMM"
  all_res[[i+2]]$uz.a$u <- NA
  all_res[[i+2]]$uz.a$z <- NA
  all_res[[i+2]]$uz.a <- all_res[[i+2]]$uz.a[,c(8:9,1:7)]
  
  uz.list[[i]]<-rbind(all_res[[i+1]]$uz.a, all_res[[i+2]]$uz.a)
}
uz.all <- do.call(rbind, datalist)
df.all <- do.call(rbind,uz.list)

uz.all$len <- uz.all$up - uz.all$low

load("./output/jqm_summary_method.RData")
sum.jqm.met <- jqm_summary_method[1:11,]
sum.jqm.met$method <- "PJQM"
sum.jqm.met$group <- c(3,4,2,9,1,7,8,10,5,6,11) 
sum.jqm.met <- sum.jqm.met[order(sum.jqm.met$group),]

load("./output/lqmm_summary_method.RData")
sum.lqmm.met <- lqmm_summary_method[1:11,]
sum.lqmm.met$method <- "LQMM"
sum.lqmm.met$group <- c(3,4,2,9,1,7,8,10,5,6,11) 
sum.lqmm.met <- sum.lqmm.met[order(sum.lqmm.met$group),]


#### Figure 2
uz.plot <- subset(uz.all, par=="Glucose" & data=="metabolomics")
df.plot <- subset(df.all, par=="Glucose" & data=="metabolomics")
df.plot <- rbind(df.plot, df.plot)
rownames(df.plot) <- NULL
df.plot$method <-c(rep("PJQM", nrow(df.plot)/2), rep("LQMM", nrow(df.plot)/2))

iri.plot.glu <- ggplot(uz.plot, aes(x=as.factor(id2), linetype=method)) +
  geom_errorbar(aes(ymin = low, ymax = up, group=method), position=position_dodge(width=0.7)) +
  scale_color_manual(values = c("grey", "black"))+
  scale_linetype_manual(name="method", values=c("dashed", "solid"), labels=c("separate-IRI", "joint-IRI")) +
  geom_point(data = df.plot, aes(x = id2, y = y), color="#34A3DC", 
             size=1.5) +
  geom_vline(xintercept=seq(1.5, length(unique(df.plot$id2))-0.5, 1), 
             lwd=0.5, colour="grey") +
  labs(x="Participants", y = "Glucose (mmol/l)") +
  annotate("rect", xmin = 0, xmax = 31, ymin = 3.9, ymax = 6.1,
           alpha = .08, fill="#92D050") +
  annotate("segment", x = 0, xend = 31, y = 3.9, yend = 3.9, color = "#746f6e", 
           alpha=1, linetype="dashed") +
  annotate("segment", x = 0, xend = 31, y = 6.1, yend = 6.1, color = "#746f6e", 
           alpha=1, linetype="dashed") +
  annotate("rect", xmin = 0, xmax = 31, ymin = 4.768, ymax = 5.715,
           alpha = .25, fill="#92D050") +
  theme_classic()+
  theme(legend.position = "bottom",
        text = element_text(size=15),
        axis.ticks = element_blank(),
        panel.border = element_rect(NA))  
iri.plot.glu

#### Figure 3
uz.plot <- subset(uz.all, par=="Total_TG" & data=="metabolomics")
df.plot <- subset(df.all, par=="Total_TG" & data=="metabolomics")
df.plot <- rbind(df.plot, df.plot)
rownames(df.plot) <- NULL
df.plot$method <-c(rep("PJQM", nrow(df.plot)/2), rep("LQMM", nrow(df.plot)/2))

iri.plot.trigl <- ggplot(uz.plot, aes(x=as.factor(id2), linetype=method)) +
  geom_errorbar(aes(ymin = low, ymax = up, group=method), position=position_dodge(width=0.7)) +
  scale_color_manual(values = c("grey", "black"))+
  scale_linetype_manual(name="method", values=c("dashed", "solid"), labels=c("separate-IRI", "joint-IRI")) +
  geom_point(data = df.plot, aes(x = id2, y = y), color="#34A3DC", 
             size=1.5) +
  geom_vline(xintercept=seq(1.5, length(unique(df.plot$id2))-0.5, 1), 
             lwd=0.5, colour="grey") +
  labs(x="Participants", y = "Triglyceride (mmol/l)") +
  annotate("rect", xmin = 0, xmax = 31, ymin = 0, ymax = 2.15,
           alpha = .08, fill="#92D050") +
  # annotate("segment", x = 0, xend = 31, y = 0, yend = 0, color = "#746f6e",
  #          alpha=1, linetype="dashed") +
  annotate("segment", x = 0, xend = 31, y = 2.15, yend = 2.15, color = "#746f6e",
           alpha=1, linetype="dashed") +
  # annotate("rect", xmin = 0, xmax = 31, ymin = 0.877, ymax = 1.773,
  #          alpha = .08, fill="darkred") +
  theme_classic()+
  theme(legend.position = "bottom",
        text = element_text(size=15),
        axis.ticks = element_blank(),
        panel.border = element_rect(NA))  
iri.plot.trigl


#### Figure 4
uz.all.metabo.chol <- subset(uz.all.metabo, par=="Clinical_LDL_C" | par=="ApoB_by_ApoA1" | par=="Total_C" | par=="Total_TG" | par=="Glucose" | par=="ApoB")

uz.all.metabo.chol$par_ord <- NA
uz.all.metabo.chol[uz.all.metabo.chol$par=="Glucose",]$par_ord <- 1
uz.all.metabo.chol[uz.all.metabo.chol$par=="Total_TG",]$par_ord <- 2
uz.all.metabo.chol[uz.all.metabo.chol$par=="Total_C",]$par_ord <- 3
uz.all.metabo.chol[uz.all.metabo.chol$par=="Clinical_LDL_C",]$par_ord <- 4
uz.all.metabo.chol[uz.all.metabo.chol$par=="ApoB_by_ApoA1",]$par_ord <- 6
uz.all.metabo.chol[uz.all.metabo.chol$par=="ApoB",]$par_ord <- 5
uz.all.metabo.chol <- uz.all.metabo.chol[order(uz.all.metabo.chol$par_ord),]
uz.all.metabo.chol$par <- factor(uz.all.metabo.chol$par, level=unique(uz.all.metabo.chol$par))

par.labs <- c("Apolipoprotein B", "LDL cholesterol", "Glucose", "Ratio Apolipoprotein B/A1", "Total cholesterol", "Triglyceride")
names(par.labs) <- c("ApoB", "Clinical_LDL_C", "Glucose", "ApoB_by_ApoA1", "Total_C", "Total_TG")

sum.method.met <- rbind(sum.lqmm.met[1:6,c(1,4:8)], sum.jqm.met[1:6,c(1,7:11)])
sum.method.met <- sum.method.met[order(sum.method.met$group),]
sum.method.met$par <- rep(unique(uz.all.metabo.chol$par), each=2)
sum.method.met$yloc[1:2] <- max(uz.all.metabo.chol$len[uz.all.metabo.chol$par=="Glucose"]) + 0.025
sum.method.met$yloc[3:4] <- max(uz.all.metabo.chol$len[uz.all.metabo.chol$par=="Total_TG"]) + 0.4
sum.method.met$yloc[5:6] <- max(uz.all.metabo.chol$len[uz.all.metabo.chol$par=="Total_C"]) + 0.5
sum.method.met$yloc[7:8] <- max(uz.all.metabo.chol$len[uz.all.metabo.chol$par=="Clinical_LDL_C"]) + 0.5
sum.method.met$yloc[9:10] <- max(uz.all.metabo.chol$len[uz.all.metabo.chol$par=="ApoB"]) + 0.01
sum.method.met$yloc[11:12] <- max(uz.all.metabo.chol$len[uz.all.metabo.chol$par=="ApoB_by_ApoA1"]) + 0.005
sum.method.met$cov.tot <- round(sum.method.met$cov.tot, digits=3)


box1 <- ggplot(uz.all.metabo.chol, aes(x=par, y=len, color=method)) +
  geom_boxplot(size=0.75, position = position_dodge(width = 0.9)) +
  labs(x=NULL, y="IRI length") +
  geom_point(position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.1))+
  #scale_linetype_manual(name="method", values=c("dashed", "solid"), labels=c("separate-IRI", "joint-IRI")) +
  scale_color_manual(name = "method", labels = c("separate-IRI", "joint-IRI"), values=c("mediumaquamarine", "#34A3DC"))+
  geom_label(data=sum.method.met, aes(y=yloc, label=cov.tot), 
             position = position_dodge(width = 0.75), fill="white", show.legend = FALSE)+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.text = element_text(size=13),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.ticks = element_blank(),
        strip.text.x = element_text(size = 12)) +
  facet_wrap(~par, scale="free", labeller = labeller(par=par.labs))
box1


#### Figure S3
uz.all.metabo.rest <- subset(uz.all.metabo, par=="ApoA1" | par=="non_HDL_C" | par=="Creatinine" | par=="Albumin")

uz.all.metabo.rest$par_ord <- NA
uz.all.metabo.rest[uz.all.metabo.rest$par=="Albumin",]$par_ord <- 1
uz.all.metabo.rest[uz.all.metabo.rest$par=="Creatinine",]$par_ord <- 2
uz.all.metabo.rest[uz.all.metabo.rest$par=="ApoA1",]$par_ord <- 4
uz.all.metabo.rest[uz.all.metabo.rest$par=="non_HDL_C",]$par_ord <- 3
uz.all.metabo.rest <- uz.all.metabo.rest[order(uz.all.metabo.rest$par_ord),]
uz.all.metabo.rest$par <- factor(uz.all.metabo.rest$par, level=unique(uz.all.metabo.rest$par))

par.labs <- c("Albumin", "Creatinine", "Apolipoprotein A1", "non-HDL cholesterol")
names(par.labs) <- c("Albumin", "Creatinine", "ApoA1", "non_HDL_C")

sum.method.met <- rbind(sum.lqmm.met[7:10,c(1,4:8)], sum.jqm.met[7:10,c(1,7:11)])
sum.method.met <- sum.method.met[order(sum.method.met$group),]
sum.method.met$par <- rep(unique(uz.all.metabo.rest$par), each=2)
sum.method.met$yloc[1:2] <- max(uz.all.metabo.rest$len[uz.all.metabo.rest$par=="Albumin"]) + 2
sum.method.met$yloc[3:4] <- max(uz.all.metabo.rest$len[uz.all.metabo.rest$par=="Creatinine"]) + 10
sum.method.met$yloc[5:6] <- max(uz.all.metabo.rest$len[uz.all.metabo.rest$par=="non_HDL_C"]) + 0.5
sum.method.met$yloc[7:8] <- max(uz.all.metabo.rest$len[uz.all.metabo.rest$par=="ApoA1"]) + 0.01
sum.method.met$cov.tot <- round(sum.method.met$cov.tot, digits=3)

box2 <- ggplot(uz.all.metabo.rest, aes(x=par, y=len, color=method)) +
  geom_boxplot(size=0.75, position = position_dodge(width = 0.9)) +
  labs(x=NULL, y="IRI length") +
  geom_point(position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.1))+
  scale_linetype_manual(name="method", values=c("dashed", "solid"), labels=c("separate-IRI", "joint-IRI")) +
  scale_color_manual(name = "method", labels = c("separate-IRI", "joint-IRI"), values=c("mediumaquamarine", "#34A3DC"))+
  geom_label(data=sum.method.met, aes(y=yloc, label=cov.tot), 
             position = position_dodge(width = 0.75), fill="white", show.legend = FALSE)+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.text = element_text(size=13),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.ticks = element_blank(),
        strip.text.x = element_text(size = 12)) +
  facet_wrap(~par, scale="free", labeller = labeller(par=par.labs))
box2

#### Figure 5
uz.glu.dat <- subset(uz.all, par=="Glucose" & method=="PJQM")
df.glu.dat <- subset(df.all, par=="Glucose")

iri.glu <- ggplot(uz.glu.dat, aes(x=as.factor(id2))) +
  geom_errorbar(aes(ymin = low, ymax = up, group=data), position=position_dodge(width=0.7)) +
  geom_point(data = df.glu.dat, aes(x = id2, y = y, group=data, color=data), 
             position=position_dodge(width=0.7), size=1.5) +
  geom_vline(xintercept=seq(1.5, length(unique(df.glu.dat$subject))-0.5, 1), 
             lwd=0.5, colour="grey") +
  labs(x="Participants", y = "Glucose (mmol/l)") +
  scale_color_manual(values=c("#F5C000", "#34A3DC")) +
  annotate("rect", xmin = 0, xmax = 31, ymin = 3.9, ymax = 6.1,
           alpha = .1, fill="#92D050") +
  annotate("segment", x = 0, xend = 31, y = 3.9, yend = 3.9, color = "#746f6e", 
           alpha=1, linetype="dashed") +
  annotate("segment", x = 0, xend = 31, y = 6.1, yend = 6.1, color = "#746f6e", 
           alpha=1, linetype="dashed") +
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.ticks = element_blank(),
        panel.border = element_rect(NA))  
iri.glu


uz.trigl.dat <- subset(uz.all, par=="Total_TG" & method=="PJQM")
df.trigl.dat <- subset(df.all, par=="Total_TG")

iri.trigl <- ggplot(uz.trigl.dat, aes(x=as.factor(id2))) +
  geom_errorbar(aes(ymin = low, ymax = up, group=data), position=position_dodge(width=0.7)) +
  geom_point(data = df.trigl.dat, aes(x = id2, y = y, group=data, color=data), 
             position=position_dodge(width=0.7), size=1.5) +
  geom_vline(xintercept=seq(1.5, length(unique(df.trigl.dat$subject))-0.5, 1), 
             lwd=0.5, colour="grey") +
  labs(x="Participants", y = "Triglyceride (mmol/l)") +
  scale_color_manual(values=c("#F5C000", "#34A3DC")) +
  annotate("rect", xmin = 0, xmax = 31, ymin = 0, ymax = 2.15,
           alpha = .1, fill="#92D050") +
  # annotate("segment", x = 0, xend = 31, y = 0, yend = 0, color = "#746f6e",
  #          alpha=1, linetype="dashed") +
  annotate("segment", x = 0, xend = 31, y = 2.15, yend = 2.15, color = "#746f6e",
           alpha=1, linetype="dashed") +
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.ticks = element_blank(),
        panel.border = element_rect(NA))  
iri.trigl

figure5 <- ggarrange(iri.glu, iri.trigl,
                     common.legend = TRUE, legend = "bottom",
                     ncol = 2, nrow = 1)
figure5