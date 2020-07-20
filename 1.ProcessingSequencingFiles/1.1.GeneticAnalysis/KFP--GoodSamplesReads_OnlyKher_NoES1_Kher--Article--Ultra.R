setwd("/Users/hintze/Desktop/PhD\ Action\ Plan/Side\ Projects/KFP/KFP--Analyses/KFP--MDS/")

library(optparse)
library(ggplot2)
library(plyr)
library(RColorBrewer)
library(ggrepel)

option_list <- list(make_option(c('-i','--in_file'), action='store', type='character', default="stdin", help='Input file'),
                    make_option(c('--no_header'), action='store_true', type='logical', default=FALSE, help='Input file has no header'),
                    make_option(c('--var_excl'), action='store', type='character', default=NULL, help='Variables to exclude from analysis'),
                    make_option(c('-a','--annot'), action='store', type='character', default=NULL, help='File with indiv annotations'),
                    make_option(c('--id_column'), action='store', type='numeric', default=1, help='Column to use as ID'),
                    make_option(c('-L','--in_maj_labels'), action='store', type='character', default=NULL, help='Column from annotation file to use as MAJOR label'),
                    make_option(c('-l','--in_min_labels'), action='store', type='character', default=NULL, help='Column from annotation file to use as MINOR label'),
                    make_option(c('-c','--in_colors'), action='store', type='character', default=NULL, help='Column from input file to use as individual colors'),
                    make_option(c('-s','--plot_size'), action='store', type='numeric', default=1, help='Plot size'),
                    make_option(c('-t','--plot_title'), action='store', type='character', default=NULL, help='Plot title'),
                    make_option(c('-x', '--plot_x_limits'), action='store', type='character', default=NULL, help='Comma-sepparated values for plot X-axis limits (eg: "-1,1")'),
                    make_option(c('-y', '--plot_y_limits'), action='store', type='character', default=NULL, help='Comma-sepparated values for plot Y-axis limits (eg: "-1,1")'),
                    make_option(c('-o','--out_file'), action='store', type='character', default=NULL, help='Output file'),
                    make_option(c('--debug'), action='store_true', type='logical', default=FALSE, help='Debug mode')
)
opt <- parse_args(OptionParser(option_list = option_list))
opt$in_file="KFP--GoodSamplesReads_OnlyKher_NoES1_Kher--Article--Ultra.mds"
opt$annot="KFP--GoodSamplesReads_OnlyKher_NoES1_Kher--Article--Ultra.annot"
opt$id_column=1
opt$in_maj_labels="Population"
opt$out_file="KFP--GoodSamplesReads_OnlyKher_NoES1_Kher--Article--Ultra_Auto.pdf"

# Print parameters
cat('# Input file:', opt$in_file, fill=TRUE)
cat('# Input file has header:', !opt$no_header, fill=TRUE)
cat('# Excluded variables:', opt$var_excl, fill=TRUE)
cat('# Annotations file:', opt$annot, fill=TRUE)
cat('# ID column:', opt$id_column, fill=TRUE)
cat('# Major label variable:', opt$in_maj_labels, fill=TRUE)
cat('# Minor label variable:', opt$in_min_labels, fill=TRUE)
cat('# Individual colors:', opt$in_colors, fill=TRUE)
cat('# Plot size:', opt$plot_size, fill=TRUE)
cat('# Plot title:', opt$plot_title, fill=TRUE)
cat('# Plot X-axis limits:', opt$plot_x_limits, fill=TRUE)
cat('# Plot Y-axis limits:', opt$plot_y_limits, fill=TRUE)
cat('# Out file:', opt$out_file, fill=TRUE)
cat('# Debug mode?:', opt$debug, fill=TRUE)

### Check plot X-axis limits
if(!is.null(opt$plot_x_limits))
  opt$plot_x_limits <- as.numeric( gsub("\\(|\\)|\"", "", strsplit(opt$plot_x_limits, ",", fixed=TRUE)[[1]]) )

### Check plot Y-axis limits
if(!is.null(opt$plot_y_limits))
  opt$plot_y_limits <- as.numeric( gsub("\\(|\\)|\"", "", strsplit(opt$plot_y_limits, ",", fixed=TRUE)[[1]]) )

### Read data
cat("# \tReading input file...", fill=TRUE)
data <- read.table(opt$in_file, row.names=1, sep="\t", header=!opt$no_header, stringsAsFactors=FALSE, check.names=FALSE)
n <- ncol(data)
if(opt$debug)
  print(data)

# Read annotation file
if(!is.null(opt$annot)){
  cat("# \tReading annotations file...", fill=TRUE)
  annot <- read.table(opt$annot, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  if(opt$debug)
    print(annot)
  data <- merge(data, annot, by.x=0, by.y=opt$id_column)
  # Get rownames back into place
  rownames(data) <- data[,1]; data <- data[,-1]
  data[colnames(annot)[opt$id_column]] <- rownames(data)
}

# Exclude variables
if( !is.null(opt$var_excl) ) {
  cat("# \tExcluding variables...", fill=TRUE)
  opt$var_excl <- unlist(strsplit(opt$var_excl, ","))
  data <- data[!(rownames(data) %in% opt$var_excl),]
}

# Set plot title
if(is.null(opt$plot_title))
  opt$plot_title <- basename(opt$in_file)

# Get Major labels mean location
colors <- NULL
if(!is.null(opt$in_maj_labels)){
  cat("# Calculating Major labels...", fill=TRUE)
  # Merge Major labels
  in_maj_labels <- unlist(strsplit(opt$in_maj_labels, ",", fixed=TRUE))
  tmp_data <- data[,in_maj_labels[1]]
  data[in_maj_labels[1]] <- NULL
  if(length(in_maj_labels) > 1){
    for (cnt in 2:length(in_maj_labels)){
      tmp_data <- paste(tmp_data, data[,in_maj_labels[cnt]], sep="/")
      data[in_maj_labels[cnt]] <- NULL
    }
    opt$in_maj_labels <- "MERGE"
  }
  
  # Make sure Major label column is after data
  data <- data.frame(data, tmp_data)
  colnames(data)[ncol(data)] <- opt$in_maj_labels
  # Convert to factor, in case there is a Major label with just numbers
  data[,opt$in_maj_labels] <- factor(data[,opt$in_maj_labels])
  # If label was in input file, decrease number of data columns
  if(is.null(opt$annot) || !opt$in_maj_labels %in% colnames(annot))
    n = n - 1
  # Get mean value for Major label
  data_mean <- ddply(data, opt$in_maj_labels, function(x){colMeans(x[, 1:n], na.rm=TRUE)})
  colors <- as.character(opt$in_maj_labels)
}
# If color variable provided, override previous definitions
if (!is.null(opt$in_colors))
  colors <- as.character(opt$in_colors)

### Plot

pdf(opt$out_file, width=opt$plot_size*8, height=opt$plot_size*6)
for(i in 1:(n-1)){
  for(j in (i+1):n){
    plot <- ggplot(data, aes_string(x=colnames(data)[i], y=colnames(data)[j], colour=colors))
    plot <- plot + 
      labs(x = "Dimension (.%)", y = "Dimension (.%)", color = "Population") +
      theme_bw() +
      theme(axis.title.x = element_text(size = 12, color="#000000", face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
            axis.title.y = element_text(size = 12, color="#000000", face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0))) +
      theme(legend.title=element_text(size=10, face="bold")) +
      theme(legend.text=element_text(size=8)) +
      theme(panel.background = element_rect(fill = '#FAFAFA')) +
      theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
            plot.title=element_text(size=10)) +
      theme(axis.line = element_line(colour = "#000000", size = 0.3)) +
      theme(panel.border = element_blank()) +
      guides(colour=guide_legend(override.aes=list(alpha=1, size=3), ncol=1)) +
      coord_cartesian(xlim=opt$plot_x_limits, ylim=opt$plot_y_limits) +
      theme(legend.background = element_rect(fill="#FAFAFA")) +
      theme(axis.text.x = element_text(color="#000000", size=7),
            axis.text.y = element_text(color="#000000", size=7)) +
      theme(axis.ticks.x = element_line(color="#000000", size=0.3), axis.ticks.y = element_line(color="#000000", size=0.3))
    
    if(nrow(data) < 8){
      plot <- plot + scale_colour_brewer(palette="Set1")
    }else{
      plot <- plot + scale_colour_discrete()
    }
    
    # Minor labels
    if(is.null(opt$in_min_labels)){
      plot <- plot + geom_point(alpha=0.7, size=1) +
        theme(legend.position="right", 
              legend.key = element_rect(fill=NA),
              legend.title=element_blank())
           }else{
      plot <- plot + geom_text(aes_string(label=opt$in_min_labels), alpha=0.1, size=1.5)
    }
    
    # Major labels
    if(!is.null(opt$in_maj_labels))
      if(!is.null(opt$in_colors)){
      }else{
      }
    print(plot)
  }
}

x <- dev.off()

#######################################

data$Population <- factor(data$Population, ordered=T, levels=c("RJ4_Kher","RJ5_Kher"))

ggplot(data, aes_string(x="D1_15.0511416716791", y="D2_9.10542166111844", colour="Population")) + geom_point(alpha = 0.6, size = 2.2) + geom_text_repel(aes(label=data$Sample_ID), size = 2.2) +
  
   scale_fill_manual(values=c("#ca0020", "#f4a582"), drop=FALSE) +
  
   scale_colour_manual(values=c("#ca0020", "#f4a582"), drop=FALSE) +
  
  scale_x_continuous("Dimension 1 (35.3%)") +
                     #breaks = c(-0.75, 0, 0.75),
                     #labels = c("-0.75", "0", "0.75"),
                     #expand = c(0,0),
                     #limits = c(-1.25, 1.25)) +
  scale_y_continuous("Dimension 2 (1.9%)") +
                     #breaks = c(-0.5, 0),
                     #labels = c("-0.75", "0")) +
                     #expand = c(0,0)) +
                     #limits = c(-0.085, 0.077)) +
  
  theme(legend.key = element_blank()) +
  theme(legend.title=element_blank()) +
  theme(axis.title.x = element_text(size = 16, color="#000000", face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, color="#000000", face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(legend.text=element_text(size=11)) +
  theme(panel.background = element_rect(fill = '#FAFAFA')) +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  theme(axis.line = element_line(colour = "#000000", size = 0.3)) +
  theme(panel.border = element_blank()) +
  #guides(colour=guide_legend(override.aes=list(alpha=0.6, size=2.2, shape = c(NA, 15, 15, NA, NA, 19, 19, 19, 19)), ncol=1)) +
  coord_cartesian(xlim=opt$plot_x_limits, ylim=opt$plot_y_limits) +
  theme(legend.background = element_rect(fill="#FAFAFA", colour = "#000000", size = 0.3)) +
  theme(axis.text.x = element_text(color="#000000", size=11),
        # angle=270, vjust=0.5, hjust=1
        axis.text.y = element_text(color="#000000", size=12)) +
  theme(axis.ticks.x = element_line(color="#000000", size=0.3), axis.ticks.y = element_line(color="#000000", size=0.3))
  #theme(legend.position=c(.94, 0.20))

ggsave("KFP--GoodSamplesReads_OnlyKher_NoES1_Kher--Article--Ultra_D1-D2_Labels.pdf", height=25, width=28, scale=0.4, dpi=1000)
