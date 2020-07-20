setwd("/Users/hintze/Desktop/PhD\ Action\ Plan/Side\ Projects/KFP/KFP--Analyses/KFP--ngsAdmix/")

library(optparse)
library(grid)
library(ggplot2)
library(RColorBrewer)
library(gtable)

samples <- read.table("KFP--GoodSamplesReads--Article--Ultra_QOPTFiles.txt", stringsAsFactors = FALSE, sep = "\t")

ids <- read.table("KFP--GoodSamplesReads--Article--Ultra.annot", stringsAsFactors = FALSE, sep = "\t", header=TRUE)

ids$Population <- factor(ids$Population, ordered=T, levels=c("ES1_Kher","RJ4_Kher","RJ5_Kher","RJ4_Koce","RJ5_Koce","SC1_Koce","SC2_Koce"))
sampleid="Sample_ID"
target="Population"

data_for_plot <- data.frame()

#c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25),
#c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),
#c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23),
#c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22),
#c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21),
#c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),
#c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19),
#c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18),
#c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),
#c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),
#c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),

#c(1,2,3,4,5,6,7,8,9,10,11,12),
#c(1,2,3,4,5,6,7,8,9,10,11),

x <- list(c(10,4,1,8,6,5,3,2,9,7),
          c(8,4,9,2,3,7,6,1,5),
          c(1,6,7,5,2,3,8,4),
          c(2,6,3,4,5,7,1),
          c(1,2,3,4,6,5),
          c(5,3,4,2,1),
          c(3,4,1,2),
          c(3,1,2),
          c(1,2))

for (j in 1:length(samples[,1])){
  data <- read.table(samples[j,1])[,x[[j]]]
  for (i in 1:dim(data)[2]) { 
    temp <- data.frame(value = data[,i])
    temp$k <- as.factor(rep(i, times = length(temp$value)))
    temp[sampleid] <- as.factor(ids[sampleid][,1])
    temp$k_value <- as.factor(rep(paste("k = ", dim(data)[2], sep =""), times = length(temp$value)))
    temp <- merge(ids, temp)
    data_for_plot <- rbind(data_for_plot, temp)
  }
}

x_lab <- (sampleid)

Plot <- ggplot(data_for_plot, aes(x = get(sampleid), y=value, fill=k)) + labs(x=x_lab) +
  geom_bar(stat="identity", width=0.85) +
  facet_grid(k_value ~ get(target), space="free_x", scales="free_x") +
  #scale_fill_manual(values=c("darkcyan", "darkseagreen", "#ffe119", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#d2f53c",
                             #"#fabebe", "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1",
                             #"#000080",  "#808080")) +
  theme(plot.title = element_blank()) +
  theme(panel.spacing=unit(0.25, "lines")) + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0), breaks = NULL) +
  theme(axis.text=element_blank()) +
  theme(axis.text.x=element_text(colour="#000000", size=10, angle=90, vjust=0.5, hjust=1),
        axis.text.y=element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  theme(strip.background=element_rect(colour="#000000", fill='#E6E6E6', size = 0.05)) +
  #theme(strip.text.x=element_text(colour="#000000", face="bold", size=7, margin=margin(0.1, 0, 0.1, 0, "cm"))) +
        #strip.text.y=element_text(colour="#000000", face="bold", size=7, margin=margin(0, 0.1, 0, 0.1, "cm"))) +
  theme(strip.text.x=element_text(colour="#000000", face="bold", size=11, margin=margin(0.25, 0, 0.25, 0, "cm")),
        strip.text.y=element_text(colour="#000000", face="bold", size=9, angle=270, margin=margin(0, 0.175, 0, 0.175, "cm"))) +
  theme(panel.grid.minor=element_blank()) + theme(panel.background=element_rect(fill="#000000")) +
  theme(axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor.x = element_blank(), legend.position="none")

ggsave(file = "KFP--GoodSamplesReads--Article--Ultra.pdf", width=15, height=9, dpi = 1000)

Plot_G <- ggplotGrob(Plot)

Plot_G <- gtable_add_rows(Plot_G, unit(0.85, "cm"), pos = 5)

Plot_G <- gtable_add_grob(Plot_G, list(rectGrob(gp = gpar(col = "#000000", fill = "#E0E0E0", size = .5, lwd = 0.25)), textGrob("Kher", gp = gpar(cex = .85, fontface = 'bold', col = "black"))), t=6, l=4, b=6, r=9, name = c("a", "b"))
Plot_G <- gtable_add_grob(Plot_G, list(rectGrob(gp = gpar(col = "#000000", fill = "#636363", size = .5, lwd = 0.25)), textGrob("Koce", gp = gpar(cex = .85, fontface = 'bold', col = "black"))), t=6, l=11, b=6, r=17, name = c("a", "b"))

Plot_G <- gtable_add_rows(Plot_G, unit(2/10, "line"), 6)

grid.newpage()
grid.draw(Plot_G)

ggsave(Plot_G, file = "PGP--GoodSamples--Article--Ultra_SNPCalling--Article--Ultra_PopLabels_RColours_Colours_Final.pdf", width=15, height=10, dpi = 1000)