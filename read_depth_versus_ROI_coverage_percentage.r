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

# Exon coverage
exoncoverage <- read.table(file = "~/LLGP_Validation/exoncoverage.summary", header = F)
exon_coverage_header <- read.table(file = "~/LLGP_Validation/exoncoverage.header", header = F)
colnames(exoncoverage) <- as.character(unlist(exon_coverage_header[1,]))
rm(exon_coverage_header)

# Read depth
readdepth <- read.table(file = "~/LLGP_Validation/READS.summary", header = F)
read_depth_header <- read.table(file = "~/LLGP_Validation/READS.header", header = F)
colnames(readdepth) <- as.character(unlist(read_depth_header[1,]))
rm(read_depth_header)

# Duplication rates and raw/passed read counts
duplication <- read.table(file = "~/LLGP_Validation/duplication.summary", header = F)
passedreads <- data.frame(do.call('rbind', strsplit(as.character(readdepth$PASSED_READS), ',', fixed = T)))
rawreads <- data.frame(do.call('rbind', strsplit(as.character(readdepth$TOTAL_READS), ',', fixed = T)))

filenames <- Sys.glob(paths = "/home/callum/LLGP_Validation/LLGP_metrics_files/LLGP_MXTHS*/*/watchdog/metrics/*.json")
split_path <- function(path) {
  if (dirname(path) %in% c(".", path)) return(basename(path))
  return(c(basename(path), split_path(dirname(path))))
}

filename_df <- exoncoverage[NULL]

for(i in 1:length(filenames)){
  filename_df$Sample_Name[i] <- split_path(filenames[i])[4]
  filename_df$Hash[i] <- tools::file_path_sans_ext(split_path(filenames[i])[1])
}

# Add everything to exoncoverage
exoncoverage$READ_DEPTH <- (as.numeric(levels(passedreads$X1))[passedreads$X1] + 
                              as.numeric(levels(passedreads$X2))[passedreads$X2])/2
exoncoverage$RAW_DEPTH <-  (as.numeric(levels(rawreads$X1))[rawreads$X1] + 
                              as.numeric(levels(rawreads$X2))[rawreads$X2])/2
exoncoverage$duplication <- duplication$V1*100
fragmentsize <- data.frame(do.call('rbind', strsplit(as.character(readdepth$LENGTH_MEDIAN), ',', fixed = T)))
exoncoverage$fragmentsize <- (as.numeric(levels(fragmentsize$X1))[fragmentsize$X1] + 
                                as.numeric(levels(fragmentsize$X2))[fragmentsize$X2])/2

# Add filenames
exoncoverage$Sample_Name <- filename_df$Sample_Name
exoncoverage$Hash <- filename_df$Hash

# Also produce a melt (raw)
rm(exoncoveragerawmelt)
exoncoveragerawmelt <- exoncoverage[,F]
exoncoveragerawmelt$RAW_DEPTH <- exoncoverage$RAW_DEPTH
exoncoveragerawmelt$COVERED_20X_PCT <- exoncoverage$COVERED_20X_PCT
exoncoveragerawmelt$COVERED_30X_PCT <- exoncoverage$COVERED_30X_PCT
exoncoveragerawmelt$COVERED_40X_PCT <- exoncoverage$COVERED_40X_PCT
exoncoveragerawmelt <- melt(exoncoveragerawmelt, id.vars = "RAW_DEPTH")

# Also produce a melt (filtered)
rm(exoncoveragefilteredmelt)
exoncoveragefilteredmelt <- exoncoverage[,F]
exoncoveragefilteredmelt$READ_DEPTH <- exoncoverage$READ_DEPTH
exoncoveragefilteredmelt$COVERED_20X_PCT <- exoncoverage$COVERED_20X_PCT
exoncoveragefilteredmelt$COVERED_30X_PCT <- exoncoverage$COVERED_30X_PCT
exoncoveragefilteredmelt$COVERED_40X_PCT <- exoncoverage$COVERED_40X_PCT
exoncoveragefilteredmelt <- melt(exoncoveragefilteredmelt, id.vars = "READ_DEPTH")

#### Plot the data ####
options(scipen = 999)

# 98% ROI and 20, 30, 40x for raw paired reads
pdf("~/LLGP_Validation/PDF_Plots/ROI_Coverage_vs_Raw_Reads_20_30_40.pdf", width = 13.33, height = 7.5)
print(
  ggplot(exoncoveragerawmelt, aes(RAW_DEPTH/1000000, value, col = variable)) + 
    geom_point() +
    geom_smooth(se = F) +
    geom_label_repel(aes(label = ifelse(
      value <= 0.990 & variable == "COVERED_40X_PCT", as.character(RAW_DEPTH), '')),
      box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
    labs(caption = "Figure 1. Number of Raw Paired Reads Versus ROI Coverage (%). 99.5% of the ROI was covered to 40x in all samples"
         ) +
    xlab("Number of Raw Paired Reads (Millions)") +
    ylab("ROI Coverage (%)") +
    scale_x_log10() +
    scale_y_continuous() +
    ggtitle("Number of Raw Paired Reads Versus ROI Coverage (%)") +
    # geom_hline(yintercept = 0.98, linetype="dashed", color = "red", size = 0.5) +
    theme(
      # Legends to the top
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      plot.caption = element_text(hjust = 0.5),
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
pdf("~/LLGP_Validation/PDF_Plots/98_ROI_vs_Processed_Reads_20_30_40.pdf", width = 13.33, height = 7.5)
print(
  ggplot(exoncoveragefilteredmelt, aes(READ_DEPTH/1000000, value, col = variable)) + 
    geom_point() +
    geom_smooth(se = F) +
    geom_label_repel(aes(label = ifelse(
      value >= 0.972 & value <= 0.98 & variable == "COVERED_30X_PCT", as.character(READ_DEPTH), '')),
      box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
    xlab("Number of Processed Paired Reads (Millions)") +
    ylab("ROI Coverage (%)") +
    scale_x_log10() +
    scale_y_continuous() +
    ggtitle("Number of Processed Paired Reads Versus ROI Coverage (%)") +
    labs(caption = "Figure 2. Number of Processed Paired Reads Versus ROI Coverage (%). 99.5% of the ROI was covered to 40x in all samples") +
    # geom_hline(yintercept = 0.98, linetype="dashed", color = "red", size = 0.5) +
    theme(
      # Legends to the top
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      plot.caption = element_text(hjust = 0.5),
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

# Reads vs % ROI covered to 200x
pdf("~/LLGP_Validation/PDF_Plots/ROI_vs_Total_Reads_200.pdf", width = 13.33, height = 7.5)
print(
  ggplot(exoncoverage, aes(READ_DEPTH, COVERED_200X_PCT)) +
    geom_point() +
    geom_smooth() +
    geom_label_repel(aes(label = ifelse(
      COVERED_200X_PCT <= 0.965, as.character(Sample_Name), '')),
      box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
    xlab("Number of Paired Reads (Passed)") +
    ylab("% of ROI at 200x Depth") +
    scale_x_continuous() +
    ggtitle("Total Paired Processed Reads Versus % of ROI at  200x Coverage") +
    labs(caption = "Figure 3. Number of Raw Paired Reads Versus ROI Coverage (%). 99.5% of the ROI was covered to 40x in all samples") +
    theme(
      # Legends to the top
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      plot.caption = element_text(hjust = 0.5),
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
pdf("~/LLGP_Validation/PDF_Plots/Duplication_vs_Raw_Reads.pdf", width = 13.33, height = 7.5)
print(
  ggplot(exoncoverage, aes(RAW_DEPTH, duplication*100)) +
    geom_point() +
    geom_smooth() +
    geom_label_repel(aes(label = ifelse(
    duplication > 0.08, Sample_Name, '')),
    box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
    xlab("Number of Paired Reads (Raw)") +
    ylab("Duplication (%)") +
    labs(caption = "Figure 4. Total paired raw reads need versus duplication rate") +
    scale_x_continuous() +
    scale_y_continuous() +
    ggtitle("Total Number of Reads and Duplication Percentage") +
      theme(
        # Legends to the top
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
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
pdf("~/LLGP_Validation/PDF_Plots/Duplication_vs_Processed_Reads.pdf", width = 13.33, height = 7.5)
print(
  ggplot(exoncoverage, aes(READ_DEPTH, duplication*100)) +
    geom_point() +
    geom_smooth() +
    geom_label_repel(aes(label = ifelse(
      duplication > 0.08, Sample_Name, '')),
    box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
    ggtitle("Total Processed Paired Reads and Duplication Percentage") +
    xlab("Number of Paired Reads (Passed)") +
    ylab("Duplication (%)") +
    labs(caption = "Figure 5. Total processed reads need versus duplication rate") +
    scale_x_continuous() +
    scale_y_continuous() +
    theme(
      # Legends to the top
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
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

library(data.table)
library(ggpubr)

Magnis_Summer <- exoncoverage[(exoncoverage$Sample_Name) %like% "MS", ]
Magnis_Winter <- exoncoverage[(exoncoverage$Sample_Name) %like% "MW", ]
Magnis_Summer$Robot_Name <- as.factor("Summer")
Magnis_Winter$Robot_Name <- as.factor("Winter")
Magnis <- rbind(Magnis_Summer, Magnis_Winter)

# Processed Reads vs % duplication rate
pdf("~/LLGP_Validation/PDF_Plots/Robot_Comparison_Total_Reads.pdf", width = 13.33, height = 7.5)
print(
  ggboxplot(Magnis, x = "Robot_Name", y = "READ_DEPTH",
            color = "Robot_Name", palette = "jco",
            add = "jitter") +
    stat_compare_means(method = "t.test", label.x = 1.5) +
    ggtitle("Comparison of Total Processed Paired Reads Between Robots") +
    xlab("Robot Name") +
    ylab("Number of Paired Reads (Passed)") +
    labs(caption = "Figure 6. Comparison of the total processed reads for each robot. There was no significant difference.") +
    theme(
      # Legends to the top
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      plot.caption = element_text(hjust = 0.5)
    )
)
graphics.off()

pdf("~/LLGP_Validation/PDF_Plots/Robot_Comparison_Fragment_Size.pdf", width = 13.33, height = 7.5)
print(
  ggboxplot(Magnis, x = "Robot_Name", y = "fragmentsize",
            color = "Robot_Name", palette = "jco",
            add = "jitter") +
    stat_compare_means(method = "t.test", label.x = 1.5) +
    ggtitle("Comparison of Average Fragment Sizes Between Robots") +
    xlab("Robot Name") +
    ylab("Fragment Size (bp)") +
    labs(caption = "Figure 7. Comparison of the average fragment sizes for samples produced by each robot. There was no significant difference.") +
    theme(
      # Legends to the top
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      plot.caption = element_text(hjust = 0.5)
    )
)
graphics.off()

pdf("~/LLGP_Validation/PDF_Plots/Robot_Comparison_Duplication.pdf", width = 13.33, height = 7.5)
print(
  ggboxplot(Magnis, x = "Robot_Name", y = "duplication",
            color = "Robot_Name", palette = "jco",
            add = "jitter") +
    stat_compare_means(method = "t.test", label.x = 1.5) +
    ggtitle("Comparison of Duplication Rates Between Robots") +
    xlab("Robot Name") +
    ylab("Duplication Rate (%)") +
    labs(caption = "Figure 8. Comparison of duplication rates for samples produced by each robot. There was no significant difference.") +
    theme(
      # Legends to the top
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      plot.caption = element_text(hjust = 0.5)
    )
)
graphics.off()



# Summary Table
sampleQCtable <- exoncoverage[,F]
sampleQCtable$Filenames <- exoncoverage$filenames
sampleQCtable$Input <- exoncoverage$Input
sampleQCtable$Raw_Reads <- exoncoverage$RAW_DEPTH
sampleQCtable$Processed_Reads <- exoncoverage$READ_DEPTH
sampleQCtable$Duplication_Percentage <- exoncoverage$duplication
sampleQCtable$Median_Fragment_Size <- exoncoverage$fragmentsize
sampleQCtable$ROI_to_30x <- exoncoverage$COVERED_30X_PCT
sampleQCtable$ROI_to_40x <- exoncoverage$COVERED_40X_PCT
sampleQCtable$ROI_to_200x <- exoncoverage$COVERED_200X_PCT

write.csv(sampleQCtable, file = "~/LLGP_Validation/summary_table.csv")
pdf("~/LLGP_Validation/PDF_Plots/sampleQCtable.pdf", height = 12, width = 14)
grid.table(sampleQCtable, rows = NULL)
dev.off()



################################################################################


# GIAB samples separated out
filenames <- list.files(path = "~/LLGP_Validation/LLGP_metrics_files_plus/", pattern = "*.json", )
filenames <- tools::file_path_sans_ext(filenames)
exoncoverage$hashnames <- filenames
exoncoverage$giab_status <- with(
  exoncoverage, ifelse(grepl("gb1", exoncoverage$hashnames, fixed = T), 1, 0))
exoncoverage$RAW_DEPTH
exoncoverage$READ_DEPTH

# GIAB raw
pdf("~/LLGP_Validation/PDF_Plots/GIAB_Paired_Raw.pdf", width = 13.33, height = 7.5)
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
pdf("~/LLGP_Validation/PDF_Plots/GIAB_Paired_Processed.pdf", width = 13.33, height = 7.5)
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


# DNA Input vs % duplication rate
pdf("~/LLGP_Validation/PDF_Plots/Input_vs_Duplication.pdf", width = 13.33, height = 7.5)
print(ggplot(exoncoverage, aes(Input, duplication*100)) +
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
