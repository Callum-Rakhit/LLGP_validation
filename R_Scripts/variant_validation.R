#!/usr/bin/env Rscript

# Function to install packages
GetPackages <- function(required.packages) {
  packages.not.installed <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if(length(packages.not.installed)){install.packages(packages.not.installed, dependencies = T)}
  suppressMessages(lapply(required.packages, require, character.only = T))
}

# Install required packages
GetPackages(c("tidyverse", "dplyr", "reshape", "reshape2", "ggrepel", "ggpubr",
              "devtools", "ggplot2", "ggraph", "data.table", "stringr"))

args <- commandArgs(trailingOnly = T)
directory <- "/home/callum//LLGP_Validation/llgp_validation_reference/output/"
args <- list.files(path = directory, pattern = "*.vcf")

vcf_data_merged_allsamples <- vcf_data_merged[NULL]

for(i in 1:length(args)){
  
  # Get the file names
  basename <- basename(args[i])
  basename <- tools::file_path_sans_ext(basename)
  
  # read two times the vcf file, first for the columns names, second for the data
  tmp_vcf <- readLines(paste(directory, args[i], sep = ""))
  tmp_vcf_data <- read.table(paste(directory, args[i], sep = ""), stringsAsFactors = F)
  
  # filter for the columns names
  tmp_vcf <- tmp_vcf[-(grep("#CHROM", tmp_vcf)+1):-(length(tmp_vcf))]
  vcf_names <- unlist(strsplit(tmp_vcf[length(tmp_vcf)], "\t"))
  names(tmp_vcf_data) <- vcf_names
  vcf_data <- na_if(tmp_vcf_data, "./.:.:.:.:.")
  vcf_data$variant_info <- apply(vcf_data[c(10:127)], 1, function(x) paste(x[!is.na(x) & x != "No"], collapse = ", "))
  vcf_data <- vcf_data[-c(10:127)]
  rm(tmp_vcf, tmp_vcf_data)
  
  # Extract information on variant
  vcf_data$actual_AF <- str_split_fixed(vcf_data$variant_info, ":", 5)[,2]
  alt <- as.numeric(str_split_fixed(vcf_data$actual_AF, ",", 2)[,2])
  ref <- as.numeric(str_split_fixed(vcf_data$actual_AF, ",", 2)[,1])
  vcf_data$actual_AF <- alt/(ref+alt)
  
  # Get the truth set
  truth_set <- fread(input = "~/LLGP_Validation/llgp_validation_reference/giab_truthset")
  colnames(truth_set) <- c("expected", "CHROM", "START", "END", "REF", "ALT", "AF", "READS", "SAMPLE", "GT", "INFO")
  
  # Add the truth set
  vcf_data_merged <- left_join(vcf_data, truth_set, by = c("POS" = "START"))
  vcf_data_merged$samplename <- basename
  vcf_data_merged_allsamples <- rbind(vcf_data_merged_allsamples, vcf_data_merged)
}

pdf("~/LLGP_Validation/PDF_Plots/Robot_runs_variants.pdf", width = 13.33, height = 7.5)
print(
  ggplot(vcf_data_merged_allsamples, aes(actual_AF, AF, col = samplename)) + 
    geom_point() +
    geom_label_repel(aes(label = ifelse(
      actual_AF <= 0.75 & AF >= 0.75, as.character(INFO.y), '')),
      box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
    geom_label_repel(aes(label = ifelse(
      actual_AF >= 0.75 & AF <= 0.75, as.character(INFO.y), '')),
      box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50',
      max.overlaps = Inf) +
    # labs(caption = "") +
    xlab("Actual allelic frequency") +
    ylab("GIAB expected frequency") +
    scale_x_log10() +
    scale_y_continuous() +
    ggtitle("Comparing actual to expected GIAB frequencies") +
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

# Plot for average divergence
predicted_AF <- vcf_data_merged_allsamples$AF*100
observed_AF <- vcf_data_merged_allsamples$actual_AF*100

vcf_data_merged_allsamples$deviation <- sqrt(predicted_AF*predicted_AF)-sqrt(observed_AF*observed_AF)

pdf("~/LLGP_Validation/PDF_Plots/Robot_runs_variants_significance.pdf", width = 13.33, height = 7.5)
print(
  ggboxplot(vcf_data_merged_allsamples, x = "samplename", y = "deviation",
            color = "samplename", palette = "jco",
            add = "jitter") +
    stat_compare_means(method = "anova", label.x = 1.5) +
    ggtitle("Comparing actual to expected GIAB allelic frequencies") +
    xlab("Sample") +
    ylab("Deviation of observed to expected") +
    labs(caption = "No signiicant difference between runs (p <= 1)") +
    theme(
      # Legends to the top
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      plot.caption = element_text(hjust = 0.5)
    )
)
graphics.off()

