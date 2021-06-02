#produce Figure 1 in the main article

library(ggplot2)
library(ggpubr)

load("./output/all_sim_results.RData")

##################
#Figure 1
##################
#create theme for Figure 1a-1c
My_Theme = theme(
  axis.title.x = element_text(size = 14),
  axis.text.x = element_text(size = 11),
  axis.text.y = element_text(size = 11),
  axis.title.y = element_text(size = 14),
  panel.border = element_rect(size=1, color="grey"),
  legend.position = "bottom",
  legend.title = element_text(size=14),
  legend.text = element_text(size=12))

labels_name <- c("Normal LQMM - N30", "Normal LQMM - N50", "Chi-sq LQMM - N30",
                 "Chi-sq LQMM - N50", "Normal PJQM - N30", "Normal PJQM - N50",
                 "Chi-sq PJQM - N30", "Chi-sq PJQM - N50",
                 "Error t(3) PJQM - N30", "Error t(3) PJQM - N50")
col_values <- c("mediumaquamarine", "mediumaquamarine", 
                "#F5C000", "#F5C000", "#34A3DC", "#34A3DC",
                "red3", "red3", "#8F9498", "#8F9498")
line_values <- c("solid","dashed",
                 "solid","dashed",
                 "solid","dashed",
                 "solid","dashed",
                 "solid","dashed")

#Figure 1a
all1 <- all[all$fixpar==1 & all$alpha ==0.05,]
lqmm1 <- ggplot(all1, aes(x=n))+
  geom_line(aes(y=EmpCov, group=group_alpha, linetype=group_alpha, color=group_alpha), size=1)+
  geom_point(aes(x=n, y=EmpCov, color=group_alpha), size=2)+
  scale_linetype_manual(name="scenario", values=line_values,
                        labels=labels_name) +
  scale_color_manual(name="scenario", values=col_values,
                     labels=labels_name) +
  scale_y_continuous(limits = c(0.65, 0.95), breaks = seq(0.65, 0.95, by = 0.05)) +
  labs(y="OEC 95% IRI") +
  theme_bw() + My_Theme
lqmm1


#Figure 1b
all2 <- all[all$fixpar==1 & all$alpha ==0.1,]
lqmm2 <- ggplot(all2, aes(x=n))+
  geom_line(aes(y=EmpCov, group=group_alpha, linetype=group_alpha, color=group_alpha), size=1)+
  geom_point(aes(x=n, y=EmpCov, color=group_alpha), size=2)+
  scale_linetype_manual(name="scenario", values=line_values,
                        labels=labels_name) +
  scale_color_manual(name="scenario", values=col_values,
                     labels=labels_name) +
  scale_y_continuous(limits = c(0.65, 0.95), breaks = seq(0.65, 0.95, by = 0.05)) +
  labs(y="OEC 90% IRI") +
  theme_bw() + My_Theme
lqmm2

#Figure 1c
all3 <- all[all$fixpar==1 & all$alpha ==0.15,]
lqmm3 <- ggplot(all3, aes(x=n))+
  geom_line(aes(y=EmpCov, group=group_alpha, linetype=group_alpha, color=group_alpha), size=1)+
  geom_point(aes(x=n, y=EmpCov, color=group_alpha), size=2)+
  scale_linetype_manual(name="scenario", values=line_values,
                        labels=labels_name) +
  scale_color_manual(name="scenario", values=col_values,
                     labels=labels_name) +
  scale_y_continuous(limits = c(0.60, 0.90), breaks = seq(0.60, 0.90, by = 0.05)) +
  labs(y="OEC 85% IRI") +
  theme_bw() + My_Theme
lqmm3

#Figure 1d
all <- all[order(all$group_alpha),]
all$method <- factor(all$method, level=unique(all$method))
box <- ggplot(all, aes(x=method, y=EmpCov, color=method))+
  geom_boxplot(size=1)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  scale_color_manual(name="method", values=c("#F5C000", "mediumaquamarine", "#34A3DC", "red3", "#8F9498"),
                     labels=c("Chi-sq LQMM", "Normal LQMM", "Normal PJQM", "Chi-sq PJQM", "Error t(3) PJQM")) +
  labs(y="OEC") +
  scale_y_continuous(limits = c(0.65, 0.95), breaks = seq(0.65, 0.95, by = 0.05)) +
  theme_bw() + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 14),
    panel.border = element_rect(size=1, color="grey"),
    legend.position = "none" )
box

#Figure 1
figure1.a <- ggarrange(lqmm1, lqmm2,
                       lqmm3, box,
                       common.legend = TRUE, legend = "bottom",
                       labels = c("(a)", "(b)", "(c)", "(d)"), align = c("hv"),
                       font.label = list(size=12), ncol = 2, nrow = 2)
figure1.a
