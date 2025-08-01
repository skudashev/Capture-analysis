###############################################
# IsoformSwitchAnalyzeR Brain Region Analysis
# Author: <Your Name>
# Date: <Date>
###############################################

# ================================
#         CONFIGURATION
# ================================
# Clear environment
rm(list = ls())

# --- Set Paths ---
base_dir <- "/Volumes/Projects/8/8289c66d-2d56-4706-a307-5a9a3eb3747e"

files_path      <- file.path(base_dir, "data/capture_quantification/paths_to_quant_files_final.txt")
samples_path    <- file.path(base_dir, "data/Capture_experiment_data/sample_description.txt")
probes_names    <- file.path(base_dir, "data/Capture_experiment_data/all_probes_names.txt")
lincRNA_probes  <- file.path(base_dir, "data/Capture_experiment_data/Capture_probes_list/lincRNA_probes.txt")

# Annotation
gtf_path        <- file.path(base_dir, "data/Nicola_aorta_out/isoquant_all/Brain/sqanti3_report/filter_out/IsoQuant_Brain.filtered.gtf")
fasta_path      <- file.path(base_dir, "data/Nicola_aorta_out/isoquant_all/Brain/sqanti3_report/filter_out/IsoQuant_Brain.filtered.fasta")
cds_gtf_path    <- file.path(base_dir, "data/Nicola_aorta_out/isoquant_all/Brain/sqanti3_report/Brain.extended_annotation.transcripts_corrected.gtf.cds.gtf")

# Output directory
output_dir      <- file.path(base_dir, "results/isoform_switch_analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ================================
#         LIBRARIES
# ================================
library(readr)
library(dplyr)
library(IsoformSwitchAnalyzeR)
library(DTUrtle)
library(tidyverse)
library(data.table)
library(tibble)
library(ggpmisc)

# ================================
#         DATA IMPORT
# ================================
files <- readLines(files_path)
conditions <- read_table(samples_path)
conditions$Donor <- paste0(conditions$Sex, "_", conditions$AgeDeath, "_", conditions$Race)
ids <- paste0(conditions$Death_group, "_", conditions$Region, "_", conditions$Condition, "_", conditions$Sample_ID)
names(files) <- ids

salmonQuant <- importIsoformExpression(sampleVector = files, ignoreAfterPeriod = FALSE)

# Remove sequin counts
for (slot in c("abundance", "counts")) {
  salmonQuant[[slot]] <- salmonQuant[[slot]][!grepl("^R", salmonQuant[[slot]]$isoform_id), ]
}

# ================================
#         EXPRESSIONS FILTER
# ================================
expressed_isoforms <- salmonQuant$abundance %>%
  mutate(across(-isoform_id, as.numeric)) %>%
  filter(rowSums(across(-isoform_id, ~ . != 0)) >= 3)

write.table(expressed_isoforms$isoform_id, 
            file.path(output_dir, "expressed_isoform_ids.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

ProjectDesign <- data.frame(
  sampleID = colnames(salmonQuant$abundance)[-1],
  condition = factor(paste0(conditions[[1]], "_", conditions[[2]], "_", conditions[[3]])),
  barcode = conditions$Barcode,
  stringsAsFactors = FALSE
)
write.table(ProjectDesign, file.path(output_dir, "project_design.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# ================================
#         REGION ANALYSIS
# ================================
regions <- unique(conditions$Region)

run_region_analysis <- function(region) {
  message("Processing region: ", region)
  
  # Filter comparisons
  comps <- expand.grid(condition_1 = levels(ProjectDesign$condition), condition_2 = levels(ProjectDesign$condition))
  comps <- comps[comps$condition_1 != comps$condition_2, ]
  comps <- comps[grepl(region, comps$condition_1) & grepl(region, comps$condition_2), ]
  comps <- comps[grepl("Control", comps$condition_1) | grepl("Control", comps$condition_2), ]
  
  # Abundance and counts
  abundance <- salmonQuant$abundance[, grepl(paste(c(region, "isoform"), collapse='|'), colnames(salmonQuant$abundance))]
  counts    <- salmonQuant$counts[,    grepl(paste(c(region, "isoform"), collapse='|'), colnames(salmonQuant$counts))]
  
  rownames(abundance) <- abundance$isoform_id; abundance <- abundance %>% select(-isoform_id)
  rownames(counts)    <- counts$isoform_id;    counts    <- counts %>% select(-isoform_id)
  
  design <- ProjectDesign[grepl(region, ProjectDesign$sampleID), ]
  
  # Run IsoformSwitchAnalyzeR
  sw <- importRdata(isoformCountMatrix   = counts,
                    isoformRepExpression = abundance,
                    designMatrix         = design,
                    isoformExonAnnoation = gtf_path,
                    isoformNtFasta       = fasta_path,
                    comparisonsToMake    = comps,
                    showProgress         = TRUE,
                    detectUnwantedEffects= TRUE)
  
  # Gene name assignment if missing
  if(all(is.na(sw$isoformFeatures$gene_name))){
    all_probes <- read.delim(probes_names, header=FALSE)
    colnames(all_probes) <- c("gene_name","gene_id")
    sw$isoformFeatures <- left_join(sw$isoformFeatures, all_probes, by='gene_id')
  }
  
  sw <- preFilter(sw, geneExpressionCutoff = 0.5, isoformExpressionCutoff = 0, removeSingleIsoformGenes = TRUE)
  sw <- isoformSwitchTestDEXSeq(sw, reduceToSwitchingGenes=TRUE)
  sw <- isoformSwitchTestSatuRn(sw, reduceToSwitchingGenes=TRUE, alpha = 0.05, dIFcutoff = 0.1)
  
  top_switch <- extractTopSwitches(sw, filterForConsequences=FALSE, extractGenes=FALSE, alpha=0.05, dIFcutoff = 0.1, n=Inf, sortByQvals=TRUE)
  write.table(top_switch, file.path(output_dir, paste0(region, "_TopSwitch.txt")), sep = "\t", quote = FALSE, row.names = FALSE)
  
  return(list(SwitchList=sw, TopSwitch=top_switch))
}

results <- lapply(regions, run_region_analysis)
names(results) <- regions

message("Analysis complete. Results saved in: ", output_dir)
