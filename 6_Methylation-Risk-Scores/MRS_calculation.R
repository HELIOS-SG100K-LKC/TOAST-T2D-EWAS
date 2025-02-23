

#Read phenotype
pheno <- read.csv("GSA_clinical_RS_final.csv")

#Read methylation data
GSA_methylation <- read.csv("GSA_methylation_data.csv")

#Merge methylation data with phenotype 
pheno2 <- merge(pheno,GSA_methylation,by="Sample")

#read beta from lasso regression 43 cpg
lasso1se_beta <- read.csv("coef_1selambda_cleaned.csv")

#scaled methylation data fopr all samples
GSA_methylation_z = GSA_methylation
for (i in 3:ncol(GSA_methylation_z)) {
  GSA_methylation_z[,i] <- scale(GSA_methylation_z[,i])
}
GSA_methylation_zv2 <- GSA_methylation_z[,-c(1:2)]
rownames(GSA_methylation_zv2) <- GSA_methylation_z$Sample

scaled_GSA_methv2_t <- as.data.frame(t(GSA_methylation_zv2))
colnames(scaled_GSA_methv2_t) <- rownames(GSA_methylation_zv2)
scaled_GSA_methv2_tv2 <- scaled_GSA_methv2_t
scaled_GSA_methv2_tv2$CpG <- rownames(scaled_GSA_methv2_tv2)

#calculate MRS for cor glm lasso1se
merge_lasso1se_MRS <- merge(lasso1se_beta,scaled_GSA_methv2_tv2, by = "CpG")

#Calculate score for lasso1se
score_lasso1se_sum_all = data.frame()

for (i in 1:nrow(merge_lasso1se_MRS)){
  
  score_lasso1se = NULL
  score_lasso1se_sum_each = NULL
  
  for (z in seq(from = 3, to = (ncol(merge_lasso1se_MRS)))){
    #for (z in seq(from = 13, to = 16)){  
    
    score_lasso1se = merge_lasso1se_MRS[i,2] * (merge_lasso1se_MRS[i,(z)])
    score_lasso1se_sum_each = cbind(score_lasso1se_sum_each,score_lasso1se)
    
  }
  
  score_lasso1se_sum_all = rbind(score_lasso1se_sum_all,score_lasso1se_sum_each)
  
} 

colnames(score_lasso1se_sum_all) <- GSA_methylation$Sample

MRS_lasso1se_each = NULL
MRS_lasso1se_all=data.frame()

for (x in 1:ncol(score_lasso1se_sum_all)){
  
  MRS_lasso1se_each = sum(score_lasso1se_sum_all[x])
  MRS_lasso1se_all = rbind(MRS_lasso1se_all,MRS_lasso1se_each)
}

rownames(MRS_lasso1se_all) <- colnames(score_lasso1se_sum_all)
colnames(MRS_lasso1se_all) <- c("MRS_lasso1se")
MRS_lasso1se_all$Sample <- rownames(MRS_lasso1se_all)

#write methylation lasso RS into csv
write.csv(MRS_lasso1se_all, file="GSA_MRS_Lasso1se_all.csv", row.names = FALSE)

#_________________


