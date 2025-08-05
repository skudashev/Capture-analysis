rm(list=ls())
#Set the working directory to my file
getwd()
# files.path <- args[1]
# samples.path <- args[2]
# gtf.path <- args[3]
# fasta.path <- args[4]

files.path <- "/Volumes/Projects/8/8289c66d-2d56-4706-a307-5a9a3eb3747e/data/Capture/capture_quantification/paths_to_quant_files_final.txt"
samples.path <- "/Volumes/Projects/8/8289c66d-2d56-4706-a307-5a9a3eb3747e/data/Capture_experiment_data/sample_description.txt"
# gtf.path <- "/Volumes/Projects/8/8289c66d-2d56-4706-a307-5a9a3eb3747e/data/Nicola_aorta_out/brain_ipsc_heart_report/filter_out/report/brain_ipsc_heart.renamed.fsm_corrected.gtf.cds.gtf"
# fasta.path <- "/Volumes/Projects/8/8289c66d-2d56-4706-a307-5a9a3eb3747e/data/Nicola_aorta_out/brain_ipsc_heart_report/filter_out/report/brain_ipsc_heart.renamed.fsm_corrected.fasta"
CDSgtf.path <- "/Volumes/Projects/8/8289c66d-2d56-4706-a307-5a9a3eb3747e/data/Nicola_aorta_out/isoquant_all/Brain/sqanti3_report/Brain.extended_annotation.transcripts_corrected.gtf.cds.gtf"
gtf.path <- "/Volumes/Projects/8/8289c66d-2d56-4706-a307-5a9a3eb3747e/data/Nicola_aorta_out/isoquant_all/Brain/sqanti3_report/filter_out/IsoQuant_Brain.filtered.gtf"
fasta.path <- "/Volumes/Projects/8/8289c66d-2d56-4706-a307-5a9a3eb3747e/data/Nicola_aorta_out/isoquant_all/Brain/sqanti3_report/filter_out/IsoQuant_Brain.filtered.fasta"
probes_names <- "/Volumes/Projects/8/8289c66d-2d56-4706-a307-5a9a3eb3747e/data/Capture_experiment_data/all_probes_names.txt"
lincRNA_probes <- "/Volumes/Projects/8/8289c66d-2d56-4706-a307-5a9a3eb3747e/data/Capture_experiment_data/Capture_probes_list/lincRNA_probes.txt"

lincRNA_probes <- read.delim(lincRNA_probes, header=FALSE)
colnames(lincRNA_probes) <- c("gene_id")

#Import libraries 
library(readr)
library(dplyr)
library(IsoformSwitchAnalyzeR)
library(DTUrtle)
library(tidyverse)
library(readr)
library(tximport) 
library(data.table)
library(tibble)
library(ggpmisc)

# import the data
files <- readLines(files.path)
conditions <- read_table(samples.path)
# add column "Donor" if Sex, AgeDeath, Race are the same for a sample
conditions$Donor <- paste0(conditions$Sex, "_", conditions$AgeDeath, "_", conditions$Race)
#how many unique donors are there?
length(unique(conditions$Donor))

id_cols <- colnames(conditions)
# ids <- apply(conditions[,id_cols], 1, paste, collapse = "_")
ids <- paste0(conditions$Death_group, "_", conditions$Region, "_", conditions$Condition, "_", conditions$Sample_ID)
names(files) <- ids
salmonQuant <- importIsoformExpression(sampleVector = files, ignoreAfterPeriod = FALSE)
# remove sequin counts from analysis
salmonQuant$abundance <- salmonQuant$abundance[!grepl("^R", salmonQuant$abundance$isoform_id),]
salmonQuant$counts <- salmonQuant$counts[!grepl("^R", salmonQuant$counts$isoform_id),]

# saveRDS(salmonQuant, "salmonQuant.rds")

# expressed_isoforms are those where abundance is not 0 in at least 3 samples (columns )
expressed_isoforms <- salmonQuant$abundance %>% 
  mutate_at(vars(-isoform_id), as.numeric) %>% 
  filter(rowSums(.[, -1] != 0) >= 3)

  
write.table(expressed_isoforms$isoform_id, "expressed_isoform_ids.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

groups <- factor(paste0(conditions[[1]], "_", conditions[[2]], "_", conditions[[3]]))
group_levels <- levels(groups)
ProjectDesign <- data.frame(sampleID = colnames(salmonQuant$abundance)[-1], condition = groups, barcode=conditions$Barcode, stringsAsFactors = FALSE) 
write.table(ProjectDesign, "project_design.txt", sep = "\t", quote = FALSE, row.names = FALSE)

comparisons <- data.frame(
  condition_1 = rep(group_levels, each = length(group_levels)),
  condition_2 = rep(group_levels, length(group_levels))
)

comparisons <- comparisons[comparisons$condition_1 != comparisons$condition_2, ]

# Filter rows where at least 2 of AgeDeath_group, Region, Condition are the same
matching_rows <- apply(comparisons, 1, function(row) {
  conditions_count <- sum(
    unlist(strsplit(row["condition_1"], "_")) == unlist(strsplit(row["condition_2"], "_"))
  )
  conditions_count >= 2
})

comparisons <- comparisons[matching_rows, ]
# Remove rows if they are mirror duplicates (e.g. remove A vs B if there already is a B vs A)

# https://community.rstudio.com/t/remove-duplicates-based-on-pairs/36608/2
comparisons <- comparisons %>%
  dplyr::mutate(normalized = purrr::map2_chr(condition_1, condition_2, ~paste(sort(c(.x, .y)), collapse = ""))) %>%
  dplyr::group_by(normalized) %>%
  dplyr::summarise(condition_1 = dplyr::first(condition_1),
                   condition_2 = dplyr::first(condition_2)) %>%
  dplyr::select(-normalized)


# IDIOT - should break up comparisons by brain region 
# Create dataframes named after each region "Caudate", "DLPFC", "Hippocampus"
regions <- unique(conditions$Region)

# Initialize empty lists to store dataframes for each region
region_comparisons_list <- list()
region_abundance_list <- list()
region_counts_list <- list()
region_design_list <- list()

# Loop through each region to create input for IsoformSwitchAnalyzeR
for (region in regions) {
  # Extract comparisons for the current region
  region_comparisons <- comparisons[grepl(paste0("Adult_",region), comparisons$condition_1) & grepl(paste0("Adult_",region), comparisons$condition_2),]
  region_comparisons <- region_comparisons[grepl("Control", region_comparisons$condition_1) | grepl("Control", region_comparisons$condition_2),]
  
  # Extract abundance and counts for the current region
  region_abundance <- salmonQuant$abundance[, grepl(paste(c(paste0("Adult_",region),"isoform"),collapse='|'), colnames(salmonQuant$abundance))]
  region_counts <- salmonQuant$counts[, grepl(paste(c(paste0("Adult_",region),"isoform"),collapse='|'), colnames(salmonQuant$counts))]
  
  # Set row names for abundance and counts dataframes
  rownames(region_abundance) <- region_abundance$isoform_id
  rownames(region_counts) <- region_counts$isoform_id
  region_abundance <- region_abundance %>% dplyr::select(-c('isoform_id'))
  region_counts <- region_counts %>% dplyr::select(-c('isoform_id'))
  
  # design 
  region_Design <- ProjectDesign[grepl(paste(c("Adult",region),collapse ='_'),  ProjectDesign$sampleID),]
  
  # Store dataframes in lists
  region_comparisons_list[[region]] <- region_comparisons
  region_abundance_list[[region]] <- region_abundance
  region_counts_list[[region]] <- region_counts
  region_design_list[[region]] <- region_Design
}

# 1. DLPFC
# comparisons have to include Adult_DLPFC_Control in at least one condition column but both must contain DLPFC 
DLPFC_comparisons <- region_comparisons_list[["DLPFC"]]
DLPFC_abundance <- region_abundance_list[["DLPFC"]]
DLPFC_counts <- region_counts_list[["DLPFC"]]
DLPFC_Design <- region_design_list[["DLPFC"]]

DLPFC_SwitchList <- importRdata(
  isoformCountMatrix   = DLPFC_counts,
  isoformRepExpression = DLPFC_abundance,
  designMatrix         = DLPFC_Design,
  isoformExonAnnoation = gtf.path,
  isoformNtFasta = fasta.path,
  comparisonsToMake = DLPFC_comparisons,
  showProgress = TRUE,
  detectUnwantedEffects = FALSE)

if(all(is.na(DLPFC_SwitchList$isoformFeatures$gene_name))){
    all_probes_names <- read.delim(probes_names, header=FALSE)
    colnames(all_probes_names) <- c("gene_name","gene_id")
    gene_id <- as.data.frame(DLPFC_SwitchList$isoformFeatures$gene_id) 
    colnames(gene_id) <- "gene_id"
    gene_names <- left_join(gene_id,all_probes_names, by='gene_id')
    DLPFC_SwitchList$isoformFeatures$gene_name <- gene_names$gene_name
} 

DLPFC_SwitchList <- preFilter(
  DLPFC_SwitchList,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE)

DLPFC_SwitchList <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = DLPFC_SwitchList,
  reduceToSwitchingGenes=TRUE
)

DLPFC_SwitchList <- isoformSwitchTestSatuRn(
  switchAnalyzeRlist = DLPFC_SwitchList,
  reduceToSwitchingGenes=TRUE,
  alpha = 0.05,
  dIFcutoff = 0.1
)

DLPFC_SwitchSummary <- extractSwitchSummary(DLPFC_SwitchList)
n=as.numeric(max(DLPFC_SwitchSummary$nrIsoforms))
DLPFC_TopSwitch <- extractTopSwitches(
  DLPFC_SwitchList, 
  filterForConsequences=FALSE,
  extractGenes=FALSE,
  alpha=0.05,
  dIFcutoff = 0.1,
  n=n,
  sortByQvals=TRUE,
  inEachComparison = TRUE
)
DLPFC_TopSwitch <- DLPFC_TopSwitch %>% drop_na(isoform_id) %>% dplyr::select(-c('gene_ref','iso_ref','Rank'))

# Now run this for Caudate_Adult samples 

Caudate_comparisons <- region_comparisons_list[["Caudate"]]
Caudate_abundance <- region_abundance_list[["Caudate"]]
Caudate_counts <- region_counts_list[["Caudate"]]
Caudate_Design <- region_design_list[["Caudate"]]

Caudate_Design <- ProjectDesign[grepl("Caudate", ProjectDesign$sampleID),]
Caudate_SwitchList <- importRdata(
  isoformCountMatrix   = Caudate_counts,
  isoformRepExpression = Caudate_abundance,
  designMatrix         = Caudate_Design,
  isoformExonAnnoation = gtf.path,
  isoformNtFasta = fasta.path,
  comparisonsToMake = Caudate_comparisons,
  showProgress = TRUE,
  detectUnwantedEffects = TRUE)

if(all(is.na(Caudate_SwitchList$isoformFeatures$gene_name))){
    all_probes_names <- read.delim(probes_names, header=FALSE)
    colnames(all_probes_names) <- c("gene_name","gene_id")
    gene_id <- as.data.frame(Caudate_SwitchList$isoformFeatures$gene_id) 
    colnames(gene_id) <- "gene_id"
    gene_names <- left_join(gene_id,all_probes_names, by='gene_id')
    Caudate_SwitchList$isoformFeatures$gene_name <- gene_names$gene_name
}

Caudate_SwitchList <- preFilter(
  Caudate_SwitchList,
  geneExpressionCutoff = 0.5,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE)

Caudate_SwitchList <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = Caudate_SwitchList,
  reduceToSwitchingGenes=TRUE
)

Caudate_SwitchSummary <- extractSwitchSummary(Caudate_SwitchList)
n=as.numeric(max(Caudate_SwitchSummary$nrIsoforms))
Caudate_TopSwitch <- extractTopSwitches(
  Caudate_SwitchList, 
  filterForConsequences=FALSE,
  extractGenes=FALSE,
  alpha=0.05,
  dIFcutoff = 0.1,
  n=n,
  sortByQvals=TRUE,
  inEachComparison = TRUE
)

Caudate_TopSwitch <- Caudate_TopSwitch %>% drop_na(isoform_id) %>% dplyr::select(-c('gene_ref','iso_ref','Rank'))

# now run this for Hippocampus_Adult samples
Hippocampus_comparisons <- region_comparisons_list[["Hippocampus"]]
Hippocampus_abundance <- region_abundance_list[["Hippocampus"]]
Hippocampus_counts <- region_counts_list[["Hippocampus"]]
Hippocampus_Design <- region_design_list[["Hippocampus"]]

Hippocampus_SwitchList <- importRdata(
  isoformCountMatrix   = Hippocampus_counts,
  isoformRepExpression = Hippocampus_abundance,
  designMatrix         = Hippocampus_Design,
  isoformExonAnnoation = gtf.path,
  isoformNtFasta = fasta.path,
  comparisonsToMake = Hippocampus_comparisons,
  showProgress = TRUE,
  detectUnwantedEffects = TRUE)

Hippocampus_SwitchList <- preFilter(
  Hippocampus_SwitchList,
  geneExpressionCutoff = 0.5,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE)

if(all(is.na(Hippocampus_SwitchList$isoformFeatures$gene_name))){
  all_probes_names <- read.delim(probes_names, header=FALSE)
  colnames(all_probes_names) <- c("gene_name","gene_id")
  gene_id <- as.data.frame(Hippocampus_SwitchList$isoformFeatures$gene_id) 
  colnames(gene_id) <- "gene_id"
  gene_names <- left_join(gene_id,all_probes_names, by='gene_id')
  Hippocampus_SwitchList$isoformFeatures$gene_name <- gene_names$gene_name
}

Hippocampus_SwitchList <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = Hippocampus_SwitchList,
  reduceToSwitchingGenes=TRUE
)

Hippocampus_SwitchSummary <- extractSwitchSummary(Hippocampus_SwitchList)
n=as.numeric(max(Hippocampus_SwitchSummary$nrIsoforms))
Hippocampus_TopSwitch <- extractTopSwitches(
  Hippocampus_SwitchList, 
  filterForConsequences=FALSE,
  extractGenes=FALSE,
  alpha=0.05,
  dIFcutoff = 0.1,
  n=n,
  sortByQvals=TRUE,
  inEachComparison = TRUE
)

Hippocampus_TopSwitch <- Hippocampus_TopSwitch %>% drop_na(isoform_id) %>% dplyr::select(-c('gene_ref','iso_ref','Rank'))

# save top switches to file
for (i in ls(pattern = "TopSwitch")) {
  currentList <- get(i)
  write.table(currentList, file = paste0(i, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}

# calculate percentage of isoform_ids in TopSwitch that do not start with ENST
for (i in ls(pattern = "TopSwitch")) 
{
  novel <- get(i) %>% filter(!grepl("ENST",isoform_id))
  print(paste0(round(nrow(novel)/nrow(get(i))*100, digits=2), "% of isoforms are novel in ", i))
}


for (i in ls(pattern = "SwitchList")) {
  region <- strsplit(i, "_")[[1]][1]
  assign(i, addORFfromGTF(get(i), CDSgtf.path, overwriteExistingORF=TRUE), envir = .GlobalEnv)
  assign(i, extractSequence(get(i), onlySwitchingGenes = TRUE, pathToOutput = paste0(region) ,outputPrefix = paste0(i, "_ORF"), writeToFile = TRUE, alsoSplitFastaFile = FALSE), envir = .GlobalEnv)
}

# http://cpc2.gao-lab.org/batch.php
for (i in ls(pattern = "SwitchList")) {
  region <- strsplit(i, "_")[[1]][1]
  assign(i, analyzeCPC2(
    get(i),
    pathToCPC2resultFile = paste0(region,"/",region,"_result_cpc2.txt"),
    codingCutoff = 0.25,
    removeNoncodinORFs = TRUE,
    quiet=FALSE
  ))
}

for (i in ls(pattern = "SwitchList")) {
  currentList <- get(i)
  currentList$isoformFeatures$codingPotential <- NA
  currentList$isoformFeatures$codingPotential[which(currentList$isoformFeatures$codingPotentialValue <= 0.25)] <- FALSE
  currentList$isoformFeatures$codingPotential[which(currentList$isoformFeatures$codingPotentialValue >= 0.50)] <- TRUE
  assign(i, currentList, envir = .GlobalEnv)
}

# ### run PFAM
# run_pfam.sh

pfam_list <- list()

for (i in ls(pattern = "SwitchList")) {
    currentList <- get(i)
    region <- strsplit(i, "_")[[1]][1]
    pfam <- dir(pattern = paste0(region,"_SwitchList_ORF_AA.pfam.txt"), full.names = TRUE, recursive = TRUE, all.files = TRUE)
    pfam_list[[i]] <- pfam
  }

# restart R if h(simpleError(msg, call)) 

for (i in ls(pattern = "SwitchList")) {
  print(i)
  assign(i, analyzePFAM(get(i), 
                        pathToPFAMresultFile = pfam_list[[i]], 
                        showProgress=TRUE))
}

# can use seqkit split -s 500 if server busy  
# biolib run --local DTU/DeepTMHMM --fasta isoformSwitchAnalyzeR_ORF_AA_complete.fasta
# or https://dtu.biolib.com/DeepTMHMM 

for (i in ls(pattern = "SwitchList")) {
  print(i)
  region <- strsplit(i, "_")[[1]][1]
  path_to_gff <- list.files(pattern = paste0(region,"_DeepTMHMM.gff3"), full.names = TRUE, recursive = TRUE, all.files = TRUE)
  assign(i, analyzeDeepTMHMM(
    get(i),
    pathToDeepTMHMMresultFile = path_to_gff,
    showProgress=TRUE))
}

# filter TopSwitches where gene_name is NA
for (i in ls(pattern = "TopSwitch")) {
  currentTable <- get(i)
  currentTable <- currentTable %>% drop_na(gene_name)
  assign(i, currentTable)
}

for (x in Hippocampus_TopSwitch$gene_name) {
  dir.create(file.path("Hippocampus_SwitchPlots"), showWarnings = FALSE)
  pdf(paste0("Hippocampus_SwitchPlots/",x, ".pdf"), width = 11, height = 6,onefile=FALSE)
  switchPlot(Hippocampus_SwitchList, gene = x, condition1 = "Adult_Hippocampus_Control", condition2 = "Adult_Hippocampus_Schizo",
             IFcutoff = 0.05, reverseMinus = FALSE, localTheme = theme_bw(base_size = 10))
  dev.off()
}

for (i in ls(pattern = "SwitchList")) {
  currentList <- get(i)
  switchPlotTopSwitches(currentList, 
                        IFcutoff = 0.05, 
                        n = Inf,
                        splitFunctionalConsequences = TRUE,
                        splitComparison = TRUE)
}

# from fetal top switches get number of those where isoform_id starts with "transcript" and dIF is positive
# and then those where isoform_id starts with "transcript" and dIF is negative and calculate ratio

# get number of isoforms where dIF is positive
num_positive <- nrow(Fetal_TopSwitch[Fetal_TopSwitch$isoform_id %like% "transcript" & Fetal_TopSwitch$dIF > 0,])
num_negative <- nrow(Fetal_TopSwitch[Fetal_TopSwitch$isoform_id %like% "transcript" & Fetal_TopSwitch$dIF < 0,])
ratio <- num_positive/num_negative
print(paste0(ratio, " more novel isoforms are upregulated than downregulated in fetal samples"))



# change the volcano plot so that point shape is  by novel vs known (ie starts with "transcript" vs starts with "ENST")
# add gene_ids to the plot for points where -log10(isoform_switch_q_value) > 20
ggplot(data=Fetal_SwitchList$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
  geom_point(
    aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05, shape=isoform_id %like% "transcript" ), # default cutoff
    size=2
  ) +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
  geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
  facet_wrap( ~ condition_2) +
  #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
  scale_color_manual('Signficant\nIsoform Switch', values = c('lightgrey','darkred')) +
  scale_shape_manual('Novel', values = c(16, 1)) +
  labs(x='dIF (Fetal/Adult)', y='-log10(FDR)') +
  theme_bw(base_size = 18) +
  geom_text(data=Fetal_SwitchList$isoformFeatures[abs(Fetal_SwitchList$isoformFeatures$dIF) > 0.5 & -log10(Fetal_SwitchList$isoformFeatures$isoform_switch_q_value) > 10  ,], aes(label=gene_name), size=3.5, nudge_x=0.1, nudge_y=0.15)

# do the same for Hippocampus switchlist
ggplot(data=DLPFC_SwitchList$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
  geom_point(
    aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05, shape=isoform_id %like% "transcript" ), # default cutoff
    size=1.5
  ) +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
  geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
  facet_wrap( ~ condition_2) +
  #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
  scale_color_manual('Signficant\nIsoform Switch', values = c('lightgrey','darkred')) +
  scale_shape_manual('Novel', values = c(16, 1)) +
  labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
  theme_bw()

ggplot(data=Control_SwitchList$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
  geom_point(
    aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05, shape=isoform_id %like% "transcript" ), # default cutoff
    size=2
  ) +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
  geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
  facet_wrap( ~ condition_1) +
  # facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
  scale_color_manual('Signficant\nIsoform Switch', values = c('lightgrey','darkred')) +
  scale_shape_manual('Novel', values = c(16, 1)) +
  labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
  theme_bw(base_size = 18) +
  geom_text(data=Control_SwitchList$isoformFeatures[abs(Control_SwitchList$isoformFeatures$dIF) > 0.5 & -log10(Control_SwitchList$isoformFeatures$isoform_switch_q_value) > 10  ,], aes(label=gene_name), size=3.5, nudge_x=0.1, nudge_y=0.15)


# Now run for Adult_DLPFC_Control vs Fetal_DLPFC_Control
condition_1 <- c("Adult_DLPFC_Control")
condition_2 <- c("Fetal_DLPFC_Control")
Fetal_comparisons <- data.frame(condition_1,condition_2)
# create abundance and count matrices, grep isoform, Fetal, Adult_DLPFC_Control 
Fetal_abundance <- salmonQuant$abundance[, grepl(paste(c("Fetal","isoform", "Adult_DLPFC_Control"),collapse='|'), colnames(salmonQuant$abundance))]
Fetal_counts <- salmonQuant$counts[, grepl(paste(c("Fetal","isoform", "Adult_DLPFC_Control"),collapse='|'), colnames(salmonQuant$counts))]
rownames(Fetal_abundance) <- Fetal_abundance$isoform_id
rownames(Fetal_counts) <- Fetal_counts$isoform_id
Fetal_abundance <- Fetal_abundance %>% dplyr::select(-c('isoform_id'))
Fetal_counts <- Fetal_counts %>% dplyr::select(-c('isoform_id'))
Fetal_Design <- ProjectDesign[grepl(paste(c("Fetal", "Adult_DLPFC_Control"),collapse='|'),  ProjectDesign$sampleID),]
Fetal_Design <- Fetal_Design %>% dplyr::select(-c('race'))
  
Fetal_SwitchList <- importRdata(
  isoformCountMatrix   = Fetal_counts,
  isoformRepExpression = Fetal_abundance,
  designMatrix         = Fetal_Design,
  isoformExonAnnoation = gtf.path,
  isoformNtFasta = fasta.path,
  comparisonsToMake = Fetal_comparisons,
  showProgress = TRUE,
  detectUnwantedEffects = TRUE)

Fetal_SwitchList <- preFilter(
  Fetal_SwitchList,
  geneExpressionCutoff = 0.5,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE)

if(all(is.na(Fetal_SwitchList$isoformFeatures$gene_name))){
  all_probes_names <- read.delim(probes_names, header=FALSE)
  colnames(all_probes_names) <- c("gene_name","gene_id")
  gene_id <- as.data.frame(Fetal_SwitchList$isoformFeatures$gene_id) 
  colnames(gene_id) <- "gene_id"
  gene_names <- left_join(gene_id,all_probes_names, by='gene_id')
  Fetal_SwitchList$isoformFeatures$gene_name <- gene_names$gene_name
}

Fetal_SwitchList <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = Fetal_SwitchList,
  reduceToSwitchingGenes=TRUE
)

Fetal_SwitchSummary <- extractSwitchSummary(Fetal_SwitchList)
n=as.numeric(max(Fetal_SwitchSummary$nrIsoforms))
Fetal_TopSwitch <- extractTopSwitches(
  Fetal_SwitchList, 
  filterForConsequences=FALSE,
  extractGenes=FALSE,
  alpha=0.05,
  dIFcutoff = 0.1,
  n=n,
  sortByQvals=TRUE,
  inEachComparison = TRUE
)

Fetal_TopSwitch <- Fetal_TopSwitch %>% drop_na(isoform_id) %>% dplyr::select(-c('gene_ref','iso_ref','Rank'))

# Surrogate variable analysis
# analysis failing with Error in estimateDispersionsGeneEst(x, maxit = maxit, quiet = quiet, modelMatrix = modelMatrix, : the number of samples and the number of model coefficients are equal,
  # i.e., there are no replicates to estimate the dispersion. use an alternate design formula
  # run sva to remove batch effects
# 1. make isoform_id column rownames for Control_counts and remove isoform_id column

rownames(Control_counts) <- Control_counts$isoform_id 
Control_counts <- Control_counts %>% dplyr::select(-c('isoform_id'))

# create table pheno from conditions and make rownames ids and remove Sample_ID column

Control_counts <- as.matrix(Control_counts)
# remove rows in Control_counts where values in over 50% of samples are 0
Control_counts <- Control_counts[rowSums(Control_counts != 0) >= ncol(Control_counts)/2,]
rownames(Control_counts) <- Control_counts$isoform_id
Control_counts <- Control_counts %>% dplyr::select(-c('isoform_id'))
Control_pheno <- as.data.frame(conditions)
rownames(Control_pheno) <- ids
Control_pheno$logCounts <- log10(colSums(Control_counts))
# filter Control_pheno to only include samples where Condition = Control but keep rownames
Control_pheno <- Control_pheno[Control_pheno$Condition == "Control",]
# remove Condition column
Control_pheno <- Control_pheno %>% dplyr::select(-c('Condition'))
mod = model.matrix(~as.factor(Age)+as.factor(Region), data=Control_pheno)
mod0 = model.matrix(~1, data=Control_pheno)

Control_abundance <- salmonQuant$abundance[,grepl("Control|isoform", colnames(salmonQuant$abundance))] 
rownames(Control_abundance) <- Control_abundance$isoform_id
Control_abundance <- Control_abundance %>% dplyr::select(-c('isoform_id'))
Control_abundance <- as.matrix(Control_abundance)
Control_abundance <- Control_abundance[rowSums(Control_abundance != 0) >= ncol(Control_abundance)/5,]

fit <- svaseq(Control_abundance,mod,mod0,n.sv=9)
# plot the svaseq results
sv_df <- as.data.frame(fit$sv)
p <- ggplot(sv_df, aes(x = V1, y = V2, color = Control_pheno$Region, pch=Control_pheno$Age)) +  geom_point(size = 2)
p + theme_bw()

mod = model.matrix(~AgeDeath+as.factor(Region), data=Control_pheno)
mod0 = model.matrix(~1, data=Control_pheno)
n.sv = num.sv(Control_abundance,mod,method="be", seed=1234)
fit <- svaseq(Control_abundance,mod,mod0,n.sv=2)
sv_df <- as.data.frame(fit$sv)
p <- ggplot(sv_df, aes(x = V1, y = V2, color = as.factor(Control_pheno$Barcode), pch=Control_pheno$Sex)) +  geom_point(size = 2)
p + theme_bw()
p <- ggplot(sv_df, aes(x = V1, y = V2, color = Control_pheno$AgeDeath, pch=Control_pheno$Sex)) +  geom_point(size = 2)
p + theme_bw()
race.colors <- c("AA" = "#7EA1C4", "CAUC" = "#C51A17")

