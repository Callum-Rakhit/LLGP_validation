# TODO(Callum)
#  - Plot showing ROI % covered to 30x versus total reads for each sample
#  - Separate out GIAB from clinical samples
#  - Raw reads vs processed reads
#  - Identify low coverage regions

#### Load/Install relevant packages ####

# Also need libssl-dev and libxml2-dev on Ubuntu 18.04 (if starting from scratch)
GetPackages <- function(required.packages) {
  packages.not.installed <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if(length(packages.not.installed)){install.packages(packages.not.installed, dependencies = T)}
  suppressMessages(lapply(required.packages, require, character.only = T))
}

GetPackages(c("ggplot2", "reshape2", "wesanderson", "tidyverse", "scales", "doParallel", "reshape2", 
              "devtools", "dplyr", "gtable", "grid", "gridExtra", "data.table", "rlist", "ggrepel"))

# Developmental packages
install_github("kassambara/easyGgplot2")  # Need devtools to use this function
library(easyGgplot2)

#### Load the data ####
exoncoverage <- read.table(file = "~/LLGP_validation/exoncoverage.summaryplus", header = F)
colnames(exoncoverage) <- read.table(file = "~/LLGP_validation/exoncoverage.header", header = F)
readdepth <- read.table(file = "~/LLGP_validation/READS.summaryplus", header = F)
colnames(readdepth) <- read.table(file = "~/LLGP_validation/READS.header", header = F)
duplication <- read.table(file = "~/LLGP_validation/duplication.summaryplus", header = F)
passedreads <- data.frame(do.call('rbind', strsplit(as.character(readdepth$PASSED_READS), ',', fixed = T)))
rawreads <- data.frame(do.call('rbind', strsplit(as.character(readdepth$TOTAL_READS), ',', fixed = T)))

exoncoverage$READ_DEPTH_F <- as.numeric(passedreads$X1)
exoncoverage$READ_DEPTH_R <- as.numeric(passedreads$X2)
exoncoverage$RAW_DEPTH_F <- as.numeric(rawreads$X1)
exoncoverage$RAW_DEPTH_R <- as.numeric(rawreads$X2)
exoncoverage$duplication <- duplication$V1

filenames <- list.files(path = "~/LLGP_validation/LLGP_metrics_files_plus/", pattern = "*.json", )
filenames <- tools::file_path_sans_ext(filenames)
exoncoverage$filenames <- filenames
filenames <- data.frame(do.call('rbind', strsplit(as.character(exoncoverage$filenames), '-', fixed = T)))
exoncoverage$filenames <- filenames$X1

#### Create the melt datasets ####
rm(exoncoveragepassedmelt)
exoncoveragepassedmelt <- exoncoverage[,F]
exoncoveragepassedmelt$READ_DEPTH <- exoncoverage$READ_DEPTH_F + exoncoverage$READ_DEPTH_R
exoncoveragepassedmelt$COVERED_20X_PCT <- exoncoverage$COVERED_20X_PCT
exoncoveragepassedmelt$COVERED_30X_PCT <- exoncoverage$COVERED_30X_PCT
exoncoveragepassedmelt$COVERED_40X_PCT <- exoncoverage$COVERED_40X_PCT
exoncoveragepassedmelt <- melt(exoncoveragepassedmelt, id.vars = "READ_DEPTH")

rm(exoncoveragerawmelt)
exoncoveragerawmelt <- exoncoverage[,F]
exoncoveragerawmelt$READ_DEPTH <- as.numeric(exoncoverage$RAW_DEPTH_F) + as.numeric(exoncoverage$RAW_DEPTH_R)
exoncoveragerawmelt$COVERED_20X_PCT <- as.numeric(exoncoverage$COVERED_20X_PCT)
exoncoveragerawmelt$COVERED_30X_PCT <- as.numeric(exoncoverage$COVERED_30X_PCT)
exoncoveragerawmelt$COVERED_40X_PCT <- as.numeric(exoncoverage$COVERED_40X_PCT)
exoncoveragerawmelt <- melt(exoncoveragerawmelt, id.vars = "READ_DEPTH")

rm(exoncoveragerawnamesmelt)
exoncoveragerawnamesmelt <- exoncoverage[,F]
exoncoveragerawnamesmelt$READ_DEPTH <- as.numeric(exoncoverage$RAW_DEPTH_F) + as.numeric(exoncoverage$RAW_DEPTH_R)
exoncoveragerawnamesmelt$COVERED_30X_PCT <- as.numeric(exoncoverage$COVERED_30X_PCT)
exoncoveragerawnamesmelt$filenames <- exoncoverage$filenames
# exoncoveragerawnamesmelt <- melt(exoncoveragerawnamesmelt, id.vars = "filenames")

# Input amounts if alpha-numeric (had to get these manually from the filenames)
exoncoverage$Input <- c(
  200, 200, 200, 200, 125, 225, 150, 200, 200, 200, 200, 200, 150, 200, 200, 200, 
  200, 200, 200, 200, 200, 125, 200, 225, 200, 200, 200, 200, 200, 200, 200, 75,
  200, 200, 200, 200, 10, 200, 200, 200, 200, 10, 75, 200, 200, 200, 200, 200
)

readdepth$Input <- c(
  200, 200, 200, 200, 125, 225, 150, 200, 200, 200, 200, 200, 150, 200, 200, 200, 
  200, 200, 200, 200, 200, 125, 200, 225, 200, 200, 200, 200, 200, 200, 200, 75,
  200, 200, 200, 200, 10, 200, 200, 200, 200, 10, 75, 200, 200, 200, 200, 200
)

readdepth_200 <- readdepth1 %>% filter(Input > 199)
rawreads_200 <- data.frame(do.call('rbind', strsplit(as.character(readdepth_200$TOTAL_READS), ',', fixed = T)))
passedreads_200 <- data.frame(do.call('rbind', strsplit(as.character(readdepth_200$PASSED_READS), ',', fixed = T)))
rawreads_200$Forward <- as.numeric(rawreads_200$X1)
rawreads_200$Reverse <- as.numeric(rawreads_200$X2)
summary(rawreads_200$Forward)
summary(rawreads_200$Reverse)

passedreads_200$Forward <- as.numeric(passedreads_200$X1)
passedreads_200$Reverse <- as.numeric(passedreads_200$X2)
summary(passedreads_200$Forward)
summary(passedreads_200$Reverse)

#### Plot the data ####

# Reads vs % ROI covered to 200x
ggplot(exoncoverage, aes(READ_DEPTH_F, COVERED_200X_PCT)) +
  geom_point() +
  geom_smooth(method = 'lm', color = "#666666") +
  xlab("Number of Reads (Passed)") +
  ylab("% of ROI at 200x depth") +
  scale_x_continuous() +
  scale_y_continuous() +
  ggtitle("Total Number of Reads Versus Percentage of Amplicons Achieving 200x Coverage") +
  theme(# Lengends to the top
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    # Remove the y-axis
    # axis.title.y = element_blank(),
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major.x = element_blank(),
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line(size = .25, color = "black"),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Rotate the x-axis labels 0 degrees
    axis.text.x = element_text(angle = 45, hjust = 1))

# Reads vs % ROI covered to 200x including sub-samples
ggplot(exoncoverageplus, aes(READ_DEPTH_F, COVERED_200X_PCT)) +
  geom_point() +
  xlab("Number of Reads (Passed)") +
  ylab("% of ROI at 200x Depth") +
  scale_x_continuous() +
  ggtitle("Total Number of Reads (inc. subsamples) Versus Percentage of Amplicons Achieving 200x Coverage") +
  geom_hline(yintercept = 0.95, linetype="dashed", color = "red", size = 0.5) +
  theme(# Lengends to the top
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    # Remove the y-axis
    # axis.title.y = element_blank(),
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major.x = element_blank(),
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line(size = .25, color = "black"),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Rotate the x-axis labels 0 degrees
    axis.text.x = element_text(angle = 45, hjust = 1))

# Plot reads vs % ROI covered for 20, 30, and 40x - line at 95 %
ggplot(exoncoverageplusmelt, aes(READ_DEPTH/1000000, value, col = variable)) + 
  geom_point() +
  geom_label_repel(aes(label = ifelse(value > 0.94 & value < 0.96, as.character(READ_DEPTH), '')),
                   box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
  xlab("Number of Passed Reads (Millions)") +
  ylab("ROI Coverage (%)") +
  scale_x_log10() +
  scale_y_continuous() +
  ggtitle("Number of Reads Versus ROI Coverage (%)") +
  geom_hline(yintercept = 0.95, linetype="dashed", color = "red", size = 0.5) +
  theme(
    # Lengends to the top
    plot.title = element_text(hjust = 0.5),
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major.x = element_blank(),
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line(size = .25, color = "black"),
    #panel.grid.minor = element_line(size = .25, color = "black"),
    # Remove panel background
    panel.background = element_blank(),
    # Rotate the x-axis labels 0 degrees
    axis.text.x = element_text(angle = 45, hjust = 1))

pdf("~/LLGP_validation/95_ROI_vs_Total_Reads_20_30_40.pdf", width = 13.33, height = 7.5)

print(
  ggplot(exoncoverageplusmelt, aes(READ_DEPTH/1000000, value, col = variable)) + 
        geom_point() +
        geom_label_repel(aes(label = ifelse(value > 0.94 & value < 0.96, as.character(READ_DEPTH), '')),
                         box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
        xlab("Number of Passed Reads (Millions)") +
        ylab("ROI Coverage (%)") +
        scale_x_log10() +
        scale_y_continuous() +
        ggtitle("Number of Reads Versus ROI Coverage (%)") +
        geom_hline(yintercept = 0.95, linetype="dashed", color = "red", size = 0.5) +
        theme(
          # Lengends to the top
          plot.title = element_text(hjust = 0.5),
          # Remove panel border
          panel.border = element_blank(),
          # Remove panel grid lines
          panel.grid.major.x = element_blank(),
          # explicitly set the horizontal lines (or they will disappear too)
          panel.grid.major.y = element_line(size = .25, color = "black"),
          #panel.grid.minor = element_line(size = .25, color = "black"),
          # Remove panel background
          panel.background = element_blank(),
          # Rotate the x-axis labels 0 degrees
          axis.text.x = element_text(angle = 45, hjust = 1))
      )

graphics.off()

# 98

pdf("~/LLGP_validation/98_ROI_vs_Total_Reads_20_30_40.pdf", width = 13.33, height = 7.5)

print(
  ggplot(exoncoverageplusmelt, aes(READ_DEPTH/1000000, value, col = variable)) + 
    geom_point() +
    geom_label_repel(aes(label = ifelse(value > 0.975 & value < 0.985, as.character(READ_DEPTH), '')),
                     box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
    xlab("Number of Passed Reads (Millions)") +
    ylab("ROI Coverage (%)") +
    scale_x_log10() +
    scale_y_continuous() +
    ggtitle("Number of Reads Versus ROI Coverage (%)") +
    geom_hline(yintercept = 0.98, linetype="dashed", color = "red", size = 0.5) +
    theme(
      # Lengends to the top
      plot.title = element_text(hjust = 0.5),
      # Remove panel border
      panel.border = element_blank(),
      # Remove panel grid lines
      panel.grid.major.x = element_blank(),
      # explicitly set the horizontal lines (or they will disappear too)
      panel.grid.major.y = element_line(size = .25, color = "black"),
      #panel.grid.minor = element_line(size = .25, color = "black"),
      # Remove panel background
      panel.background = element_blank(),
      # Rotate the x-axis labels 0 degrees
      axis.text.x = element_text(angle = 45, hjust = 1))
)

graphics.off()

# Raw reads 98
options(scipen=999)

pdf("~/LLGP_validation/98_ROI_vs_Raw_Reads_20_30_40.pdf", width = 13.33, height = 7.5)

print(
  ggplot(exoncoverageplusrawmelt, aes(READ_DEPTH/1000000, value, col = variable)) + 
    geom_point() +
    geom_label_repel(aes(label = ifelse(value > 0.975 & value < 0.985, as.character(READ_DEPTH), '')),
                     box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
    xlab("Number of Raw Reads (Millions)") +
    ylab("ROI Coverage (%)") +
    scale_x_log10() +
    scale_y_continuous() +
    ggtitle("Number of Reads Versus ROI Coverage (%)") +
    geom_hline(yintercept = 0.98, linetype="dashed", color = "red", size = 0.5) +
    theme(
      # Lengends to the top
      plot.title = element_text(hjust = 0.5),
      # Remove panel border
      panel.border = element_blank(),
      # Remove panel grid lines
      panel.grid.major.x = element_blank(),
      # explicitly set the horizontal lines (or they will disappear too)
      panel.grid.major.y = element_line(size = .25, color = "black"),
      #panel.grid.minor = element_line(size = .25, color = "black"),
      # Remove panel background
      panel.background = element_blank(),
      # Rotate the x-axis labels 0 degrees
      axis.text.x = element_text(angle = 45, hjust = 1))
)

graphics.off()

# Raw reads, 30x ROI % coverage and including filenames
pdf("~/LLGP_validation/98_ROI_vs_Total_Reads_30_with_filenames.pdf", width = 13.33, height = 7.5)

print(
  ggplot(exoncoveragerawnamesmelt, aes(READ_DEPTH/1000000, COVERED_30X_PCT, col = filenames)) + 
    geom_point() +
    geom_label_repel(aes(label = ifelse(COVERED_30X_PCT > 0.975 & COVERED_30X_PCT < 0.9856, 
                                        as.character(READ_DEPTH), '')),
                     box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
    xlab("Number of Passed Reads (Millions)") +
    ylab("ROI Coverage (%)") +
    scale_x_log10() +
    scale_y_continuous() +
    ggtitle("Number of Reads Versus ROI Coverage (%)") +
    geom_hline(yintercept = 0.98, linetype="dashed", color = "red", size = 0.5) +
    theme(
      # Lengends to the top
      plot.title = element_text(hjust = 0.5),
      # Remove panel border
      panel.border = element_blank(),
      # Remove panel grid lines
      panel.grid.major.x = element_blank(),
      # explicitly set the horizontal lines (or they will disappear too)
      panel.grid.major.y = element_line(size = .25, color = "black"),
      #panel.grid.minor = element_line(size = .25, color = "black"),
      # Remove panel background
      panel.background = element_blank(),
      # Rotate the x-axis labels 0 degrees
      axis.text.x = element_text(angle = 45, hjust = 1))
)

graphics.off()

options(scipen = 999)
exoncoveragerawnamesmelt$giab_status <- with(exoncoveragerawnamesmelt, ifelse(filenames == "gb1", 1, 0))

ggplot(exoncoveragerawnamesmelt, aes(READ_DEPTH/1000000, COVERED_30X_PCT, fill = as.factor(giab_status))) +
  geom_point(data = subset(exoncoveragerawnamesmelt, giab_status == 0), color = "grey") +
  geom_point(data = subset(exoncoveragerawnamesmelt, giab_status == 1), color = "red") +
  geom_label_repel(aes(label = ifelse(
    COVERED_30X_PCT > 0.975 & COVERED_30X_PCT < 0.9856, as.character(READ_DEPTH), '')),
    box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
  scale_fill_manual(values = c("grey", "red"),
                    name = "GIAB Status",
                    breaks = c("0", "1"),
                    labels = c("Clinical", "GIAB")) +
  ggtitle("Number of Reads Versus ROI Coverage at 30x (%)") +
  xlab("Number of Passed Reads (Millions)") +
  ylab("ROI Coverage at 30x (%)") +
  scale_x_log10() +
  scale_y_continuous() +
  geom_hline(yintercept = 0.98, linetype = "dashed", color = "red", size = 0.5) +
  theme(
    # Lengends to the top
    plot.title = element_text(hjust = 0.5),
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major.x = element_blank(),
    # Explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line(size = .25, color = "black"),
    # Remove panel background
    panel.background = element_blank(),
    # Rotate the x-axis labels 0 degrees
    axis.text.x = element_text(angle = 45, hjust = 1))

# Plot reads vs % ROI covered for 20, 30, and 40x - line at 99 %
ggplot(exoncoverageplusmelt, aes(READ_DEPTH/1000000, value, col = variable)) + 
  geom_point() +
  geom_label_repel(aes(label = ifelse(value > 0.985 & value < 0.991, as.character(READ_DEPTH), '')),
                   box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
  xlab("Number of Passed Reads (Millions)") +
  ylab("ROI Coverage (%)") +
  scale_x_log10() +
  scale_y_continuous() +
  ggtitle("Number of Reads Versus ROI Coverage (%)") +
  geom_hline(yintercept = 0.98, linetype="dashed", color = "red", size = 0.5) +
  theme(
    # Lengends to the top
    plot.title = element_text(hjust = 0.5),
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major.x = element_blank(),
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line(size = .25, color = "black"),
    #panel.grid.minor = element_line(size = .25, color = "black"),
    # Remove panel background
    panel.background = element_blank(),
    # Rotate the x-axis labels 0 degrees
    axis.text.x = element_text(angle = 45, hjust = 1))

# Reads vs % duplication rate
ggplot(exoncoverage, aes(READ_DEPTH_F, duplication*100)) +
  geom_point() +
  geom_smooth(method = 'lm', color = "#666666") +
  xlab("Number of Reads (Passed)") +
  ylab("Duplication (%)") +
  scale_x_continuous() +
  scale_y_continuous() +
  ggtitle("Total Number of Reads and Duplication Percentage") +
  theme(
    # Lengends to the top
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    # Remove the y-axis
    # axis.title.y = element_blank(),
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major.x = element_blank(),
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line(size = .25, color = "black"),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Rotate the x-axis labels 0 degrees
    axis.text.x = element_text(angle = 45, hjust = 1))

# DNA Input vs % duplication rate
ggplot(exoncoverage, aes(Input, duplication*100)) +
  geom_point() +
  geom_smooth(method = 'lm', color = "#666666") +
  xlab("DNA Input (ng)") +
  ylab("Duplication (%)") +
  scale_x_continuous() +
  scale_y_continuous() +
  ggtitle("DNA Input (ng) and Duplication Percentage") +
  theme(
    # Lengends to the top
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    # Remove the y-axis
    # axis.title.y = element_blank(),
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major.x = element_blank(),
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line(size = .25, color = "black"),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Rotate the x-axis labels 0 degrees
    axis.text.x = element_text(angle = 45, hjust = 1))

# DNA Input vs % ROI at 200x
ggplot(exoncoverage, aes(Input, COVERED_200X_PCT)) +
  geom_point() +
  geom_smooth(method = 'lm', color = "#666666") +
  xlab("DNA Input (ng)") +
  ylab("% of ROI at 200x depth") +
  scale_x_continuous() +
  scale_y_continuous() +
  ggtitle("DNA Input (ng) and % ROI Covered to 200x Depth") +
  theme(
    # Lengends to the top
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    # Remove the y-axis
    # axis.title.y = element_blank(),
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major.x = element_blank(),
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line(size = .25, color = "black"),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Rotate the x-axis labels 0 degrees
    axis.text.x = element_text(angle = 45, hjust = 1))
