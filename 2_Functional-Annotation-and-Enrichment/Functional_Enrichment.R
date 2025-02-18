library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggnewscale)


# Load CpG and binding site data
Sentinel_cpg_data <- fread("Path_to_sentinel_CpG/Sentinel_CpG.txt")
Sentinel_cpg_data$BP <- Sentinel_cpg_data$Position
names(Sentinel_cpg_data) <- c("CpG","Gene","Chr","BP","Position")
Sentinel_cpg_data[, Chr := paste0("chr",as.character(Chr))]


All_CpGs <- fread("Path_to_full_sumstats/T2D_EWAS_Meta-analysis.txt",header = T)
All_CpGs <- All_CpGs[,c(1,2,4:7)]
names(All_CpGs) <- c("CpG","Chr","Position","Beta","SE","P")
All_CpGs[, Chr := paste0("chr",as.character(Chr))]
All_CpGs <- na.omit(All_CpGs)


########
cell_types <- read.table("Path_to_func_annot_files/Func_Annot_bed_Files/cell_types.txt",header = F)
cell_types <- cell_types$V1
Marks <- c("1_TssA","2_TssAFlnk","3_TxFlnk","4_Tx","5_TxWk","6_EnhG","7_Enh","8_ZNF_Rpts","9_Het","10_TssBiv",
           "11_BivFlnk","12_EnhBiv","13_ReprPC","14_ReprPCWk","15_Quies")

summary_results <- data.table(Cell_Acc = character(),Mark = character(),Sentinel_CpGs = integer(), All_CpGs = integer())

for (i in cell_types){
  Func_annot_bed <- read.delim(paste0("Path_to_func_annot_files/Func_Annot_bed_Files/",i,"_15_coreMarks_hg38lift_mnemonics.bed"),
                               col.names = c("Chr","Start","End","Mark"))
  for (j in Marks){
    M1 <- Func_annot_bed %>% filter(Mark %in% paste0(j)) %>% select(Chr, Start, End)
    M1 <- as.data.table(M1)  # Convert to data.table
    M1[, Chr := as.character(Chr)]
    
    results_sentinel <- Sentinel_cpg_data[M1, on = .(Chr, Position >= Start, Position <= End), nomatch = 0]
    results_all <- All_CpGs[M1, on = .(Chr, Position >= Start, Position <= End), nomatch = 0]
    
    summary_results <- rbind(summary_results, data.table(
      Cell_Acc = paste0(i),
      Mark = paste0(j),
      Sentinel_CpGs = length(unique(results_sentinel$CpG)),
      All_CpGs = length(unique(results_all$CpG))
    ))
  }
}


total_all_cpgs <- nrow(na.omit(All_CpGs))  # Total background CpGs
total_sentinel_cpgs <- nrow(Sentinel_cpg_data)  # Total observed CpGs (subset)


summary_results <- summary_results %>%
  mutate(
    # Fold enrichment
    Fold_Enrichment = (Sentinel_CpGs / total_sentinel_cpgs) /
      (All_CpGs / total_all_cpgs),
    
    # Enrichment hypergeometric p-value
    Hypergeom_P = phyper(
      Sentinel_CpGs - 1,          # Observed overlaps (successes in subset)
      All_CpGs,                   # Total overlaps in background (successes in population)
      total_all_cpgs - All_CpGs,  # Background without overlaps (failures)
      total_sentinel_cpgs,        # Total subset size (trials)
      lower.tail = FALSE          # Probability of observing this or more overlaps
    ),
    
    # Depletion hypergeometric p-value
    Depletion_P = phyper(
      Sentinel_CpGs,              # Observed overlaps (successes in subset)
      All_CpGs,                   # Total overlaps in background (successes in population)
      total_all_cpgs - All_CpGs,  # Background without overlaps (failures)
      total_sentinel_cpgs,        # Total subset size (trials)
      lower.tail = TRUE           # Probability of observing this or fewer overlaps
    )
  ) %>%
  mutate(
    # Multiple testing correction for enrichment p-values
    P_FDR = p.adjust(Hypergeom_P, method = "fdr"),
    
    # Multiple testing correction for depletion p-values
    Depletion_FDR = p.adjust(Depletion_P, method = "fdr")
  )
  
  
summary_results$P_Bonf <- summary_results$Depletion_Bonf <- NULL



##########################################


### DHS and H3 in GrCh37 ####


All_CpGs_hg19 <- fread("Path_to_full_summary-statistics/TOAST/T2D_EWAS_Meta-analysis.txt",header = T)
All_CpGs_hg19 <- All_CpGs_hg19[,c(1,2,3,5:7)]
names(All_CpGs_hg19) <- c("CpG","Chr","Position","Beta","SE","P")
All_CpGs_hg19[, Chr := paste0("chr",as.character(Chr))]
All_CpGs_hg19 <- na.omit(All_CpGs_hg19)


Sentinel_cpg_hg19 <- All_CpGs_hg19 %>% filter(CpG %in% Sentinel_cpg_data$CpG)



H3_marks <- c("H3K4me1","H3K4me3","H3K9me3","H3K27me3","H3K36me3")

j <- "H3K4me1"

summary_results_reg <- data.table(Cell_Acc = character(),Mark = character(),Sentinel_CpGs = integer(), All_CpGs = integer())

for (i in cell_types){
  for (j in H3_marks) {
  Func_annot_bed <- read.delim(paste0("Path_to_func_annot_files/Func_Annot_bed_Files/",i,"-",j,".broadPeak"),
                               col.names = c("Chr","Start","End","Rank","V1","V2","V3","V4","P"))
    M1 <- Func_annot_bed %>% filter(P>=2) %>% select(Chr, Start, End)
    M1 <- as.data.table(M1)  # Convert to data.table
    M1[, Chr := as.character(Chr)]
    
    results_sentinel <- Sentinel_cpg_hg19[M1, on = .(Chr, Position >= Start, Position <= End), nomatch = 0]
    results_all <- All_CpGs_hg19[M1, on = .(Chr, Position >= Start, Position <= End), nomatch = 0]
    summary_results_reg <- rbind(summary_results_reg, data.table(
      Cell_Acc = paste0(i),
      Mark = paste0(j),
      Sentinel_CpGs = length(unique(results_sentinel$CpG)),
      All_CpGs = length(unique(results_all$CpG))
    ))
  }
}


for (i in cell_types){
  Func_annot_bed <- read.delim(paste0("Path_to_function_annotation_files/Func_Annot_bed_Files/",i,"-DNase.hotspot.fdr0.01.broad.bed"),
                                 col.names = c("Chr","Start","End","V1","Rank"))
    M1 <- Func_annot_bed %>% select(Chr, Start, End)
    M1 <- as.data.table(M1)  # Convert to data.table
    M1[, Chr := as.character(Chr)]
    
    results_sentinel <- Sentinel_cpg_hg19[M1, on = .(Chr, Position >= Start, Position <= End), nomatch = 0]
    results_all <- All_CpGs_hg19[M1, on = .(Chr, Position >= Start, Position <= End), nomatch = 0]
    summary_results_reg <- rbind(summary_results_reg, data.table(
      Cell_Acc = paste0(i),
      Mark = "DHS",
      Sentinel_CpGs = length(unique(results_sentinel$CpG)),
      All_CpGs = length(unique(results_all$CpG))
    ))
}

summary_results_reg <- summary_results_reg %>%
  mutate(
    # Fold enrichment
    Fold_Enrichment = (Sentinel_CpGs / total_sentinel_cpgs) /
      (All_CpGs / total_all_cpgs),
    
    # Enrichment hypergeometric p-value
    Hypergeom_P = phyper(
      Sentinel_CpGs - 1,          # Observed overlaps (successes in subset)
      All_CpGs,                   # Total overlaps in background (successes in population)
      total_all_cpgs - All_CpGs,  # Background without overlaps (failures)
      total_sentinel_cpgs,        # Total subset size (trials)
      lower.tail = FALSE          # Probability of observing this or more overlaps
    ),
    
    # Depletion hypergeometric p-value
    Depletion_P = phyper(
      Sentinel_CpGs,              # Observed overlaps (successes in subset)
      All_CpGs,                   # Total overlaps in background (successes in population)
      total_all_cpgs - All_CpGs,  # Background without overlaps (failures)
      total_sentinel_cpgs,        # Total subset size (trials)
      lower.tail = TRUE           # Probability of observing this or fewer overlaps
    )
  ) %>%
  mutate(
    # Multiple testing correction for enrichment p-values
    P_FDR = p.adjust(Hypergeom_P, method = "fdr"),
    
    # Multiple testing correction for depletion p-values
    Depletion_FDR = p.adjust(Depletion_P, method = "fdr")
  )


summary_results <- summary_results[,c(1:6,8,7,9)]

Enrichment_reg_all <- rbind.data.frame(summary_results_reg,summary_results)


Cell_type_info <- read.delim("Cell_Information.txt",header = T) ## Info about cell type, and accession codes from roadmap epigenomic consortium

Enrichment_reg_all <- merge.data.frame(Cell_type_info,Enrichment_reg_all,by.x = "Accession",by.y = "Cell_Acc")

