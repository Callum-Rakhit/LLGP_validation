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
  suppressMessages(lapply(required.packages, require, character.only = T))}

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

# Add everything to exoncoverage
exoncoverage$READ_DEPTH <- (as.numeric(passedreads$X1) + as.numeric(passedreads$X2))/2
exoncoverage$RAW_DEPTH <- (as.numeric(rawreads$X1) + as.numeric(rawreads$X2))/2
exoncoverage$duplication <- duplication$V1
fragmentsize <- data.frame(do.call('rbind', strsplit(as.character(readdepth$LENGTH_MEDIAN), ',', fixed = T)))
exoncoverage$fragmentsize <- (as.numeric(fragmentsize$X1) + as.numeric(fragmentsize$X2))/2

# filenames <- Sys.glob(paths = "~/LLGP_validation/LLGP_metrics_files_plus/*.json")
# list.files(path = "~/LLGP_validation/LLGP_metrics_files_plus/", pattern = "*.json", full.names = F, no.)

#### Load the subsampled data ####
# exoncoverageplus <- read.table(file = "~/LLGP_validation/exoncoverage.summaryplus", header = F)
# colnames(exoncoverageplus) <- read.table(file = "~/LLGP_validation/exoncoverage.header", header = F)
# readdepthplus <- read.table(file = "~/LLGP_validation/READS.summaryplus", header = F)
# colnames(readdepthplus) <- read.table(file = "~/LLGP_validation/READS.header", header = F)
# passedreadsplus <- data.frame(do.call('rbind', strsplit(as.character(readdepthplus$PASSED_READS), ',', fixed = T)))

rm(exoncoverageplusmelt)
exoncoverageplusmelt <- exoncoverage[,F]
exoncoverageplusmelt$READ_DEPTH <- exoncoverage$READ_DEPTH
exoncoverageplusmelt$COVERED_20X_PCT <- exoncoverage$COVERED_20X_PCT
exoncoverageplusmelt$COVERED_30X_PCT <- exoncoverage$COVERED_30X_PCT
exoncoverageplusmelt$COVERED_40X_PCT <- exoncoverage$COVERED_40X_PCT
exoncoverageplusmelt <- melt(exoncoverageplusmelt, id.vars = "READ_DEPTH")

rm(exoncoverageplusrawmelt)
exoncoverageplusrawmelt <- exoncoverage[,F]
exoncoverageplusrawmelt$READ_DEPTH <- as.numeric(exoncoverage$RAW_DEPTH)
exoncoverageplusrawmelt$COVERED_20X_PCT <- as.numeric(exoncoverage$COVERED_20X_PCT)
exoncoverageplusrawmelt$COVERED_30X_PCT <- as.numeric(exoncoverage$COVERED_30X_PCT)
exoncoverageplusrawmelt$COVERED_40X_PCT <- as.numeric(exoncoverage$COVERED_40X_PCT)
exoncoverageplusrawmelt <- melt(exoncoverageplusrawmelt, id.vars = "READ_DEPTH")

# Input amounts if alpha-numeric (had to get these manually from the filenames)
exoncoverageNoSubsamples <- exoncoverage[exoncoverage$RAW_DEPTH > 750000,]
exoncoverageNoSubsamples$Input <- c(
  200, 200, 200, 200, 125, 225, 150, 200, 200, 200, 200, 200, 150, 200, 200, 200, 
  200, 200, 200, 200, 200, 125, 200, 225, 200, 200, 200, 200, 200, 200, 200, 75,
  200, 200, 200, 200, 10, 200, 200, 200, 200, 10, 75, 200, 200, 200, 200, 200
)
filenames <- read.csv("~/LLGP_validation/filenames_hashes.csv", header = T, sep = "\t")
exoncoverageNoSubsamples$filenames <- filenames$Name
exoncoverageNoSubsamples200ng <- exoncoverageNoSubsamples[exoncoverageNoSubsamples$Input > 199,]

#### Plot the data ####
options(scipen = 999)

# 98% ROI and 20, 30, 40x for raw paired reads
pdf("~/LLGP_validation/PDF_Plots/98_ROI_vs_Raw_Reads_20_30_40.pdf", width = 13.33, height = 7.5)
print(
  ggplot(exoncoverageplusrawmelt, aes(READ_DEPTH/1000000, value, col = variable)) + 
    geom_point() +
    geom_label_repel(aes(label = ifelse(
      value >= 0.972 & value <= 0.98 & variable == "COVERED_30X_PCT", as.character(READ_DEPTH), '')),
      box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
    xlab("Number of Raw Paired Reads (Millions)") +
    ylab("ROI Coverage (%)") +
    scale_x_log10() +
    scale_y_continuous() +
    ggtitle("Number of Raw Paired Reads Versus ROI Coverage (%)") +
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

# 98% ROI and 20, 30, 40x for processed paired reads
pdf("~/LLGP_validation/PDF_Plots/98_ROI_vs_Processed_Reads_20_30_40.pdf", width = 13.33, height = 7.5)
print(
  ggplot(exoncoverageplusmelt, aes(READ_DEPTH/1000000, value, col = variable)) + 
    geom_point() +
    geom_label_repel(aes(label = ifelse(
      value >= 0.972 & value <= 0.98 & variable == "COVERED_30X_PCT", as.character(READ_DEPTH), '')),
      box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
    xlab("Number of Processed Paired Reads (Millions)") +
    ylab("ROI Coverage (%)") +
    scale_x_log10() +
    scale_y_continuous() +
    ggtitle("Number of Processed Paired Reads Versus ROI Coverage (%)") +
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

# GIAB samples separated out
filenames <- list.files(path = "~/LLGP_validation/LLGP_metrics_files_plus/", pattern = "*.json", )
filenames <- tools::file_path_sans_ext(filenames)
exoncoverage$hashnames <- filenames
exoncoverage$giab_status <- with(
  exoncoverage, ifelse(grepl("gb1", exoncoverage$hashnames, fixed = T), 1, 0))
exoncoverage$RAW_DEPTH
exoncoverage$READ_DEPTH
# GIAB raw
pdf("~/LLGP_validation/PDF_Plots/GIAB_Paired_Raw.pdf", width = 13.33, height = 7.5)
print(
ggplot(exoncoverage, aes(RAW_DEPTH/1000000, COVERED_30X_PCT, fill = as.factor(giab_status))) +
  geom_point(data = subset(exoncoverage, giab_status == 0), color = "grey") +
  geom_point(data = subset(exoncoverage, giab_status == 1), color = "red") +
  geom_label_repel(aes(label = ifelse(
    COVERED_30X_PCT > 0.975 & COVERED_30X_PCT < 0.9856, as.character(RAW_DEPTH), '')),
    box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
  scale_fill_manual(values = c("grey", "red"),
                    name = "GIAB Status",
                    breaks = c("0", "1"),
                    labels = c("Clinical", "GIAB")) +
  ggtitle("Number of Raw Paired Reads Versus ROI 30x (%)") +
  xlab("Number of Raw Paired Reads (Millions)") +
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
)
graphics.off()

# GIAB processed
pdf("~/LLGP_validation/PDF_Plots/GIAB_Paired_Processed.pdf", width = 13.33, height = 7.5)
print(
  ggplot(exoncoverage, aes(READ_DEPTH/1000000, COVERED_30X_PCT, fill = as.factor(giab_status))) +
    geom_point(data = subset(exoncoverage, giab_status == 0), color = "grey") +
    geom_point(data = subset(exoncoverage, giab_status == 1), color = "red") +
    geom_label_repel(aes(label = ifelse(
      COVERED_30X_PCT > 0.975 & COVERED_30X_PCT < 0.9856, as.character(READ_DEPTH), '')),
      box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
    scale_fill_manual(values = c("grey", "red"),
                      name = "GIAB Status",
                      breaks = c("0", "1"),
                      labels = c("Clinical", "GIAB")) +
    ggtitle("Number of Paired Processed Reads Versus ROI at 30x (%)") +
    xlab("Number of Passed Paired Reads (Millions)") +
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
)
graphics.off()

# Reads vs % ROI covered to 200x including sub-samples
pdf("~/LLGP_validation/PDF_Plots/ROI_vs_Total_Reads_200.pdf", width = 13.33, height = 7.5)
print(
  ggplot(exoncoverageNoSubsamples200ng, aes(READ_DEPTH, COVERED_200X_PCT)) +
    geom_point() +
    xlab("Number of Paired Reads (Passed)") +
    ylab("% of ROI at 200x Depth") +
    scale_x_continuous() +
    ggtitle("Total Paired Processed Reads Versus % of ROI at  200x Coverage") +
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
)
graphics.off()

# Raw Reads vs % duplication rate
pdf("~/LLGP_validation/PDF_Plots/Duplication_vs_Raw_Reads.pdf", width = 13.33, height = 7.5)
print(
  ggplot(exoncoverageNoSubsamples, aes(RAW_DEPTH, duplication*100)) +
    geom_point() +
    geom_label_repel(aes(label = ifelse(
    duplication > 0.5, filenames, '')),
    box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
    xlab("Number of Paired Reads (Raw)") +
    ylab("Duplication (%)") +
    labs(caption = "Figure 5. Total paired raw reads need versus duplication rate") +
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
        # Centre the labels
        plot.caption = element_text(hjust = 0.5),
        # Rotate the x-axis labels 0 degrees
        axis.text.x = element_text(angle = 45, hjust = 1))
)
graphics.off()

# Processed Reads vs % duplication rate
pdf("~/LLGP_validation/PDF_Plots/Duplication_vs_Processed_Reads.pdf", width = 13.33, height = 7.5)
print(
  ggplot(exoncoverageNoSubsamples, aes(READ_DEPTH, duplication*100)) +
geom_point() +
geom_label_repel(aes(label = ifelse(
duplication > 0.5, filenames, '')),
box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
ggtitle("Total Processed Paired Reads and Duplication Percentage") +
xlab("Number of Paired Reads (Passed)") +
ylab("Duplication (%)") +
labs(caption = "Figure 6. Total processed reads need versus duplication rate") +
scale_x_continuous() +
scale_y_continuous() +
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
# Centre the labels
plot.caption = element_text(hjust = 0.5),
# Rotate the x-axis labels 0 degrees
axis.text.x = element_text(angle = 45, hjust = 1))
)
graphics.off()

# DNA Input vs % duplication rate
pdf("~/LLGP_validation/PDF_Plots/Input_vs_Duplication.pdf", width = 13.33, height = 7.5)
print(ggplot(exoncoverageNoSubsamples, aes(Input, duplication*100)) +
  geom_point() +
  ggtitle("DNA Input (ng) and Duplication Percentage") +
  xlab("DNA Input (ng)") +
  ylab("Duplication (%)") +
  labs(caption = "Figure 7. Input DNA versus duplication rate") +
  geom_label_repel(aes(label = ifelse(
  duplication > 0.5, filenames, '')),
  box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
  scale_x_continuous() +
  scale_y_continuous() +
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
    # Centre the labels
    plot.caption = element_text(hjust = 0.5),
    # Rotate the x-axis labels 0 degrees
    axis.text.x = element_text(angle = 45, hjust = 1))
)
graphics.off()

# Summary Table
sampleQCtable <- exoncoverageNoSubsamples200ng[,F]
sampleQCtable$Filenames <- exoncoverageNoSubsamples200ng$filenames
sampleQCtable$Input <- exoncoverageNoSubsamples200ng$Input
sampleQCtable$Raw_Reads <- exoncoverageNoSubsamples200ng$RAW_DEPTH
sampleQCtable$Processed_Reads <- exoncoverageNoSubsamples200ng$READ_DEPTH
sampleQCtable$Duplication_Percentage <- exoncoverageNoSubsamples200ng$duplication
sampleQCtable$Median_Fragment_Size <- exoncoverageNoSubsamples200ng$fragmentsize
sampleQCtable$ROI_to_30x <- exoncoverageNoSubsamples200ng$COVERED_30X_PCT
sampleQCtable$ROI_to_40x <- exoncoverageNoSubsamples200ng$COVERED_40X_PCT
sampleQCtable$ROI_to_200x <- exoncoverageNoSubsamples200ng$COVERED_200X_PCT

write.csv(sampleQCtable, file = "~/LLGP_validation/summary_table.csv")
pdf("~/LLGP_validation/sampleQCtable.pdf", height = 12, width = 14)
grid.table(sampleQCtable, rows = NULL)
dev.off()

