## SMR code 


top_snps_all <- read.delim("toast_sentinel_trans_mqtl_enrichment.txt",header = T)
Signif_snps <- top_snps_all %>% dplyr::filter(P_adj<=0.05)

Coloc_res <- read.delim("SNP_trans_CpG_cis_eQTl_coloc.txt",header = T)
Coloc_res <- Coloc_res %>% dplyr::filter(RSID %in% Signif_snps$SNP)
Coloc_res$SNP <- ifelse(Coloc_res$RSID=="chr7:21927361:T:C","rs1581694",Coloc_res$SNP) ### fixing missing rsid


Coloc_res <- Coloc_res[,c(1,2,3,4,7,8)]


## 
Trans_meQTL_per_numCpG_lowest_pVal_gt_500Kb_assoc <- readxl::read_excel("Trans-meQTL per numCpG lowest pVal gt 500Kb assoc.xlsx")

Cpg_SNP_Assoc <- Trans_meQTL_per_numCpG_lowest_pVal_gt_500Kb_assoc[,c(1,2,9,10,13:16)]
names(Cpg_SNP_Assoc) <- c("RSID","CpG","CpG_A1","CpG_A2","CpG_freq","CpG_Beta","CpG_SE","CpG_P")

Analysis_table <- merge.data.frame(Coloc_res,Cpg_SNP_Assoc,by=c("CpG","RSID"))

Gene_SNP_Assoc <- read.table("snp_gene_assoc.txt",header = T) ## From eQTLgen
Gene_SNP_Assoc$Gene_SE = 1

Analysis_table <- merge.data.frame(Analysis_table,Gene_SNP_Assoc,by=c("SNP","Gene"))
Analysis_table <- Analysis_table %>% mutate(
  Gene_Beta = ifelse(CpG_A1==Gene_A2 & CpG_A2==Gene_A1,-1*Gene_Beta,Gene_Beta)
)

Analysis_table$CpG_Z <- Analysis_table$CpG_Beta/Analysis_table$CpG_SE
Analysis_table$CpG_SE2 <- 1

# Load the dplyr package
library(dplyr)

### SMR using single SNP model ###
Analysis_table <- Analysis_table %>%
  # Filter out rows where Gene_Beta is 0 to avoid division errors
  filter(Gene_Beta != 0) %>%
  # Calculate the Wald's ratio (causal estimate), SE, Z, and P-value
  mutate(
    Beta_causal = CpG_Z / Gene_Beta,  # Causal estimate
    SE_causal = sqrt((CpG_SE2^2 / Gene_Beta^2) + 
                       (CpG_Z^2 * Gene_SE^2 / Gene_Beta^4)),  # Standard error
    Z_causal = Beta_causal / SE_causal,  # Z-statistic
    P_causal = 2*pnorm(abs(Z_causal),lower.tail = F)  # P-value
  )

# Display the updated data frame
head(Analysis_table)


Analysis_table_results <- Analysis_table[,c(1,2,3,4,7,8,18,12,15,20:23)]

########


#Coloc Sample code

library(coloc)

dat1 <- read.table("CpG.summary.txt",header=T)
dat2 <- read.table("Gene.summary.txt",header=T,fill=T)

### Can be any 2 pair of summary statistics (CpG, Gene, Cardio-metabolic traits, etc)
### required columns : SNP, Pos, A1,A2,Freq, Beta, SE, P, N


## harmonize sumstats
Common_table <- merge.data.frame(dat1,dat2,by="SNP")

numeric_cols <- c("freq","b","se","p","freq_gene","b_gene","se_gene","p_gene")
int_cols <- c("Pos","N","N_gene")

Common_table[numeric_cols] <- lapply(Common_table[numeric_cols], as.numeric)
Common_table[int_cols] <- lapply(Common_table[int_cols], as.integer)

Common_table$b_gene <- ifelse(Common_table$A1==Common_table$A1_gene & Common_table$A2==Common_table$A2_gene,Common_table$b_gene,
				ifelse(Common_table$A1==Common_table$A2_gene & Common_table$A2==Common_table$A1_gene,-1*Common_table$b_gene,NA))
Common_table <- na.omit(Common_table)
Common_table <- Common_table[!duplicated(Common_table$SNP), ]
Common_table$varbeta <- (Common_table$se)^2
Common_table$varbeta_gene <- (Common_table$se_gene)^2

## Prepare coloc files
D1 <- as.list(Common_table[,c(1,2,5,6,17,9)])
D1$type <- "quant" ## for binary - change to "cc"
names(D1) <- c("snp","position","MAF","beta","varbeta","N","type")
D2 <- as.list(Common_table[,c(1,2,12,13,18,16)])
D2$type <- "quant"
names(D2) <- c("snp","position","MAF","beta","varbeta","N","type")

## Perform coloc
my_res <- coloc.abf(dataset1=D1,dataset2=D2)
