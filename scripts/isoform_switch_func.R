###############################################
# Helper Functions for PFAM, CPC2, DeepTMHMM & Plots
###############################################

# ================================
# PFAM Annotation
# ================================
annotate_pfam <- function(results, pfam_dir, overwriteExisting = TRUE) {
  message("[INFO] Running PFAM annotation...")
  
  for (region in names(results)) {
    sw <- results[[region]]$SwitchList
    
    pfam_file <- list.files(pfam_dir, 
                            pattern = paste0(region, "_SwitchList_ORF_AA.pfam.txt"), 
                            full.names = TRUE)
    
    if (length(pfam_file) == 0) {
      warning("PFAM file not found for region: ", region)
      next
    }
    
    sw <- analyzePFAM(sw,
                      pathToPFAMresultFile = pfam_file,
                      showProgress = TRUE)
    
    results[[region]]$SwitchList <- sw
  }
  
  return(results)
}

# ================================
# CPC2 Annotation
# ================================
annotate_cpc2 <- function(results, cpc2_results_dir, codingCutoff = 0.25) {
  message("[INFO] Running CPC2 annotation...")
  
  for (region in names(results)) {
    sw <- results[[region]]$SwitchList
    
    cpc2_file <- file.path(cpc2_results_dir, paste0(region, "_result_cpc2.txt"))
    
    if (!file.exists(cpc2_file)) {
      warning("CPC2 results not found for region: ", region)
      next
    }
    
    sw <- analyzeCPC2(sw,
                      pathToCPC2resultFile = cpc2_file,
                      codingCutoff = codingCutoff,
                      removeNoncodinORFs = TRUE,
                      quiet = FALSE)
    
    results[[region]]$SwitchList <- sw
  }
  
  return(results)
}

# ================================
# DeepTMHMM Annotation
# ================================
annotate_deeptmhmm <- function(results, deeptmhmm_dir) {
  message("[INFO] Running DeepTMHMM annotation...")
  
  for (region in names(results)) {
    sw <- results[[region]]$SwitchList
    
    deeptmhmm_file <- list.files(deeptmhmm_dir, 
                                 pattern = paste0(region, "_DeepTMHMM.gff3"), 
                                 full.names = TRUE)
    
    if (length(deeptmhmm_file) == 0) {
      warning("DeepTMHMM results not found for region: ", region)
      next
    }
    
    sw <- analyzeDeepTMHMM(sw,
                           pathToDeepTMHMMresultFile = deeptmhmm_file,
                           showProgress = TRUE)
    
    results[[region]]$SwitchList <- sw
  }
  
  return(results)
}

# ================================
# Plotting Switches
# ================================
plot_switches <- function(results, output_plot_dir, IFcutoff = 0.05) {
  message("[INFO] Generating switch plots...")
  dir.create(output_plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (region in names(results)) {
    sw <- results[[region]]$SwitchList
    top_switch <- results[[region]]$TopSwitch
    
    # Plot each gene in TopSwitch
    gene_names <- unique(na.omit(top_switch$gene_name))
    
    region_plot_dir <- file.path(output_plot_dir, region)
    dir.create(region_plot_dir, recursive = TRUE, showWarnings = FALSE)
    
    for (gene in gene_names) {
      pdf(file.path(region_plot_dir, paste0(gene, ".pdf")), width = 11, height = 6, onefile = FALSE)
      switchPlot(sw,
                 gene = gene,
                 condition1 = unique(top_switch$condition_1)[1],
                 condition2 = unique(top_switch$condition_2)[1],
                 IFcutoff = IFcutoff,
                 reverseMinus = FALSE,
                 localTheme = theme_bw(base_size = 10))
      dev.off()
    }
    
    # Batch plot top switches
    pdf(file.path(region_plot_dir, paste0(region, "_TopSwitches.pdf")), width = 11, height = 8)
    switchPlotTopSwitches(sw,
                          IFcutoff = IFcutoff,
                          n = Inf,
                          splitFunctionalConsequences = TRUE,
                          splitComparison = TRUE)
    dev.off()
  }
}

###############################################
# Example Usage
###############################################
# results <- lapply(regions, run_region_analysis)
# results <- annotate_pfam(results, pfam_dir = "/path/to/pfam/results")
# results <- annotate_cpc2(results, cpc2_results_dir = "/path/to/cpc2/results")
# results <- annotate_deeptmhmm(results, deeptmhmm_dir = "/path/to/deeptmhmm/results")
# plot_switches(results, output_plot_dir = file.path(output_dir, "SwitchPlots"))
