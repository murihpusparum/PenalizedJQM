#produce the plot of separate-IRI vs joint-IRI
#Figure 2 and Figure 4

library(readr)
library(ggplot2)
library(lqmm)
library(tidyverse)
source("JQM_Function.R")

#set confidence level and colors
alpha = 0.05

#light green, normal values
rgb(red=146, green=208, blue=80, maxColorValue = 255) 
"#92D050"
#orange, clinical
rgb(red=245, green=192, blue=0, maxColorValue = 255) 
"#F5C000"
#blue, metabolomics
rgb(red=52, green=163, blue=220, maxColorValue = 255) 
"#34A3DC"

#data preparation
combined_data <- read.csv("./data/long_combined_data.csv", sep=";", header = T)
combined_data$person_id <- as.character(combined_data$person_id)
combined_data$subject <- lapply(combined_data$person_id, parse_number)
combined_data$subject <- as.numeric(combined_data$subject)

combined_data$time <- combined_data$month
combined_data[combined_data$time == 3,]$time <- 2
combined_data[combined_data$time == 5,]$time <- 3
combined_data[combined_data$time == 7,]$time <- 4
combined_data[combined_data$time == 9,]$time <- 5
combined_data[combined_data$time == 11,]$time <- 6

#prepare metabolomics dataset
#change label to get IRI for other parameters, e.g. label=="Total_TG" for triglyceride
m.sub <- subset(combined_data, type=="metabolomics" & label=="Glucose")
m.sub <- m.sub[,c(10,11,8,9)]
colnames(m.sub) <- c("subject", "time", "y", "data")
m.sub <- m.sub[order(m.sub[,2]), ]
m.sub <- m.sub[order(m.sub[,1]), ]

#fit the Penalized JQM for metabolomics data with 95\% nominal level and sets of lambda.u and lambda.z
met.iam<-jqm(db=m.sub,
             alpha=alpha,
             lambda.u.seq = seq(0.02,0.1,0.02),
             lambda.z.seq = seq(0.5,5,0.5))
uz.m <- cbind.data.frame(met.iam$u, met.iam$z)
uz.m$low <- met.iam$beta0 + met.iam$u + met.iam$z*met.iam$beta1
uz.m$up  <- met.iam$beta0 + met.iam$u + met.iam$z*met.iam$beta2
uz.m$data <- "metabolomics"
uz.m$id <- unique(m.sub$subject)
colnames(uz.m) <- c("u", "z", "low", "up", "data", "id")


#prepare clinical dataset
#change label to get IRI for other parameters, e.g. label=="Triglyceriden" for triglyceride
c.sub <- subset(combined_data, type=="clinical" & label=="Glucose")
c.sub <- c.sub[,c(10,11,8,9)]
colnames(c.sub) <- c("subject", "time", "y", "data")
c.sub <- c.sub[order(c.sub[,2]), ]
c.sub <- c.sub[order(c.sub[,1]), ]

#fit the Penalized JQM for clinical data with 95\% nominal level and sets of lambda.u and lambda.z
clin.iam<-jqm(db=c.sub,
             alpha=alpha,
             lambda.u.seq = seq(0.02,0.1,0.02),
             lambda.z.seq = seq(0.5,5,0.5))

uz.c <- cbind.data.frame(clin.iam$u, clin.iam$z)
uz.c$low <- clin.iam$beta0 + clin.iam$u + clin.iam$z*clin.iam$beta1
uz.c$up  <- clin.iam$beta0 + clin.iam$u + clin.iam$z*clin.iam$beta2
uz.c$data <- "clinical"
uz.c$id <- unique(c.sub$subject)
colnames(uz.c) <- c("u", "z", "low", "up", "data", "id")

#combine data and results
uz <- rbind.data.frame(uz.c, uz.m)
df <- rbind.data.frame(c.sub, m.sub)

#Figure 4
iri.data <- ggplot(uz, aes(x=as.factor(id))) +
  geom_errorbar(aes(ymin = low, ymax = up, group=data), position=position_dodge(width=0.7)) +
  geom_point(data = df, aes(x = subject, y = y, group=data, color=data), 
             position=position_dodge(width=0.7), size=1.5) +
  geom_vline(xintercept=seq(1.5, length(unique(df.a$id2))-0.5, 1), 
             lwd=0.5, colour="grey") +
  labs(x="participants", y = "Glucose (mmol/l)") +
  scale_color_manual(values=c("#F5C000", "#34A3DC")) +
  annotate("rect", xmin = 0, xmax = 31, ymin = 3.885, ymax = 6.105,
           alpha = .1, fill="#92D050") +
  annotate("segment", x = 0, xend = 31, y = 3.885, yend = 3.885, color = "#746f6e", 
           alpha=1, linetype="dashed") +
  annotate("segment", x = 0, xend = 31, y = 6.105, yend = 6.105, color = "#746f6e", 
           alpha=1, linetype="dashed") +
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.ticks = element_blank(),
        panel.border = element_rect(NA))  


#fit IRI using LQMM for metabolomics data
fit1 <- lqmm(fixed = y ~ 1, random = ~ 1, group = subject, tau = c(0.025, 0.975),
                 nK = 11, type = "normal", data = m.sub, 
                 control = lqmmControl(LP_max_iter = 1000, LP_tol_ll = 1e-7, UP_max_iter = 100, UP_tol = 1e-5,
                                       beta = 0.5, gamma = 1))
pred1 <- predict(fit1, level = 1)
pred1 <- data.frame(pred1)
pred1 <- data.frame(pred1 %>% distinct())

uz.l <- data.frame(ranef(fit1),pred1$X0.025, pred1$X0.975)
uz.l$method <- "LQMM"
uz.l$id <- unique(m.sub$subject)
colnames(uz.l) <- c("u", "z", "low", "up", "method", "id")

#combine data and results of Penalized JQM
colnames(uz.m) <- c("u", "z", "low", "up", "method", "id")
uz.m$method <- "Penalized JQM"

uz.method <- rbind.data.frame(uz.m, uz.l)
df.method <- rbind.data.frame(m.sub, m.sub)
df.method$method <-c(rep("LQMM", nrow(m.sub)), rep("Penalized JQM", nrow(m.sub)))


#Figure 2
iri.method <- ggplot(uz.method, aes(x=as.factor(id), linetype=method)) +
  geom_errorbar(aes(ymin = low, ymax = up, group=method), position=position_dodge(width=0.7)) +
  scale_color_manual(values = c("grey", "black"))+
  scale_linetype_manual(name="method", values=c("dashed", "solid"), labels=c("separate-IRI", "joint-IRI")) +
  geom_point(data = df.method, aes(x = subject, y = y, group=method), color="#34A3DC", 
             position=position_dodge(width=0.7), size=2.5) +
  geom_vline(xintercept=seq(1.5, length(unique(df.b$id2))-0.5, 1), 
             lwd=0.5, colour="grey") +
  labs(x="participants", y = "Glucose (mmol/l)") +
  annotate("rect", xmin = 0, xmax = 31, ymin = 3.885, ymax = 6.105,
           alpha = .08, fill="#92D050") +
  annotate("segment", x = 0, xend = 31, y = 3.885, yend = 3.885, color = "#746f6e", 
           alpha=1, linetype="dashed") +
  annotate("segment", x = 0, xend = 31, y = 6.105, yend = 6.105, color = "#746f6e", 
           alpha=1, linetype="dashed") +
  theme_classic()+
  theme(legend.position = "bottom",
        text = element_text(size=15),
        axis.ticks = element_blank(),
        panel.border = element_rect(NA))  

