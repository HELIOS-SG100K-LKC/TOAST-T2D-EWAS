library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)


# Load CpG and binding site data
Sentinel_cpg_data <- fread("~/path_to_sentinel_CpG/Sentinel_CpG.txt")
Sentinel_cpg_data$BP <- Sentinel_cpg_data$Position
names(Sentinel_cpg_data) <- c("CpG","Gene","Chr","BP","Position")
Sentinel_cpg_data[, Chr := paste0("chr",as.character(Chr))]


All_CpGs <- fread("~/Path_to_full_summary_statistics/T2D_EWAS_Meta-analysis.txt",header = T)
All_CpGs <- All_CpGs[,c(1,2,4:7)]
names(All_CpGs) <- c("CpG","Chr","Position","Beta","SE","P")
All_CpGs[, Chr := paste0("chr",as.character(Chr))]


# List all TF BED files in a directory
bed_files <- list.files(path = "~/Path_to_TFBS_annotation_file/TF_BED_files/", pattern = "*.bed", full.names = TRUE)

### TF Bed files downloaded from reMAP 2022 Homo sapiens database (CRM peaks)

# Initialize a summary results table
summary_results <- data.table(TF = character(), Sentinel_CpGs = integer(), All_CpGs = integer())

# Loop through each BED file and perform the overlap analysis
for (bed_file in bed_files) {
  # Extract TF name from the file name
  tf_name <- gsub(".*/|\\.bed$", "", bed_file)
  
  # Load the TF binding data
  tf_binding_data <- fread(bed_file, col.names = c("Chr", "Start", "End"))
  tf_binding_data[, Chr := as.character(Chr)]
  
  # Perform overlap analysis
  results_sentinel <- Sentinel_cpg_data[tf_binding_data, on = .(Chr, Position >= Start, Position <= End), nomatch = 0]
  results_all <- All_CpGs[tf_binding_data, on = .(Chr, Position >= Start, Position <= End), nomatch = 0]
  
  # Add the results to the summary table
  summary_results <- rbind(summary_results, data.table(
    TF = tf_name,
    Sentinel_CpGs = length(unique(results_sentinel$CpG)),
    All_CpGs = length(unique(results_all$CpG))
  ))
}



# Calculate Fold Enrichment and Hypergeometric p-values
summary_results <- summary_results %>%
  mutate(
    # Fold enrichment
    Fold_Enrichment = (Sentinel_CpGs / total_sentinel_cpgs) /
      (All_CpGs / total_all_cpgs),
    
    # Hypergeometric p-value
    Hypergeom_P = phyper(
      Sentinel_CpGs - 1,          # Observed overlaps (successes in subset)
      All_CpGs,                   # Total overlaps in background (successes in population)
      total_all_cpgs - All_CpGs,  # Background without overlaps (failures)
      total_sentinel_cpgs,        # Total subset size (trials)
      lower.tail = FALSE          # Probability of observing this or more overlaps
    )
  )

summary_results_filt <- summary_results %>% filter(Sentinel_CpGs>0)
summary_results_filt$P_Adjusted <- p.adjust(summary_results_filt$Hypergeom_P,method = "fdr")
