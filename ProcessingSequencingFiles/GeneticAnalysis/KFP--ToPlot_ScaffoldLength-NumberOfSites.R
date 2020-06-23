### To plot correlation between scaffold lengths and number of SNPs

setwd("/Users/hintze/Desktop/PhD\ Action\ Plan/Side\ Projects/KFP/KFP--Analyses/KFP--Miscellaneous/KFP-SitesInfo/")

library(ggplot2)

a <- read.table("KFP--GoodSamplesReads_NoES1_Kher_SITES--Article--Ultra.ScaffoldInfo_OnlyWithSites.txt")
colnames(a) <- c("Scaffold","ScaffoldLength", "NumberOfSNPs")

Regression <- lm(formula = a$ScaffoldLength ~ a$NumberOfSNPs, data = a)
summary(Regression)

ggplot(a,aes(ScaffoldLength, NumberOfSNPs)) + stat_smooth(method="lm") + geom_point(alpha=0.7, color="#000000", size=1) +
  annotate("text", label = "Multiple R-squared: 0.9311", x = 5000000, y = 14000, size = 5, colour = "#FF0000") +
  annotate("text", label = "P-value: < 2.2e-16", x = 5000000, y = 13200, size = 5, colour = "#FF0000") +
  scale_x_continuous("Scaffold Length",
                     breaks = c(2500000, 5000000, 7500000, 10000000, 12500000),
                     labels = c("2.5e+06", "5e+06", "7.5e+06", "1e+07", "1.25e+07"),
                     limits = c(0, 12800000),
                     expand = c(0,0)) +
  scale_y_continuous("# of Sites",
                     breaks = c(2500, 5000, 7500, 10000, 12500, 15000),
                     labels = c("2.5K", "5K", "7.5K", "10K", "12.5K", "15K"),
                     limits = c(0, 15300),
                     expand = c(0,0)) +
  theme(axis.text.x = element_text(size=10, color="#000000"),
        axis.text.y = element_text(size=10, color="#000000")) +
  theme(axis.title.x = element_text(size = 16, color="#000000", face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, color="#000000", face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(panel.background = element_rect(fill = '#FAFAFA')) +
  theme(axis.ticks.x = element_line(size=0.3, color="#000000"),
        axis.ticks.y = element_line(size=0.3, color="#000000")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor=element_blank(), 
        axis.line=element_line(colour="#000000", size=0.3, color="#000000")) +
  theme(panel.border=element_blank())

ggsave("KFP--GoodSamplesReads_NoES1_Kher_SITES--Article--Ultra.ScaffoldInfo_OnlyWithSites.pdf", height=2, width=4, scale=2.9, dpi = 1000)