library(dplyr)
library(ggplot2)
library(MASS)

# regression
genetic <- read.csv("Sensitive_analysis/LOLIPOP_Mich_323snp_final.csv")

ctrlprobPC <- read.table("ctrlprobes_PC1to30.txt")
phe <- read.table("EPIC_LOLIPOP_clinpath.txt", header = TRUE)
phe_genetic <- merge(phe,genetic, by="SampleID")
houseman <- read.table("Houseman_unconstrained.txt")

houseman_ctrl <- merge(houseman, ctrlprobPC, by="row.names")
colnames(houseman_ctrl)[1] <- "SentrixID"

cpg_312 <- read.csv("Sensitive_analysis/finallist_312_cg_st.csv")
colnames(cpg_312)[1] <- "CpG"
cpg_312_list <- as.character(cpg_312$CpG)

load("beta_QN.RData")
beta_312cg <- beta[rownames(beta) %in% c(cpg_312_list), ]
beta_312cg_t <- data.frame(t(beta_312cg))
beta_312cg_t$SentrixID <- row.names(beta_312cg_t)

#filter for 1819 samples in the list of 2000 samples phe file
phe2 <- merge(phe_genetic, houseman_ctrl, by = "SentrixID")

#merge beta with pheno
phe5 <- merge(phe2, beta_312cg_t, by="SentrixID")

#define cpg and snp dosage pairs
analysis <- read.csv("Sensitive_analysis/finallist_312_cg_st.csv")

#define the dimension of matrix
nvar=nrow(analysis)
res=matrix(ncol=4, nrow=nvar)
colnames(res) =c('Estimate','Std.Error', 'z_value','P')  

for (i in 1:nrow(analysis)){

timestamp()
print(i)

cpg <- as.character(analysis[i, 1])  # Convert to character if needed
snp <- as.character(analysis[i, 2])

# Dynamically construct the formula
lfla <- as.formula(sprintf(
  "as.factor(phe5$Incident_T2D) ~ (phe5$%s) + as.factor(phe5$%s) + phe5$Age + phe5$Gender + phe5$CD8T + phe5$CD4T + phe5$NK + phe5$Bcell + phe5$Mono + phe5$Gran + phe5$PC1_cp + phe5$PC2_cp + phe5$PC3_cp + phe5$PC4_cp + phe5$PC5_cp + phe5$PC6_cp + phe5$PC7_cp + phe5$PC8_cp + phe5$PC9_cp + phe5$PC10_cp + phe5$PC11_cp + phe5$PC12_cp + phe5$PC13_cp + phe5$PC14_cp + phe5$PC15_cp + phe5$PC16_cp + phe5$PC17_cp + phe5$PC18_cp + phe5$PC19_cp + phe5$PC20_cp + phe5$PC21_cp + phe5$PC22_cp + phe5$PC23_cp + phe5$PC24_cp + phe5$PC25_cp + phe5$PC26_cp + phe5$PC27_cp + phe5$PC28_cp + phe5$PC29_cp + phe5$PC30_cp",
  cpg, snp
))

fit = summary(glm(lfla, family = binomial))
res[i,] = fit$coefficients[2,]
rm(fit)
    
}

timestamp()
rownames(res) <- analysis$Sentinel_CpG
write.table(res, file="Sensitive_analysis/results_LOLIPOP850K_324cg_adjusted_genetic_unscaled.txt")

#________

#Sensitivity analysis, adjust for bmi, glu

#define formula
lfla_bmi_glu=as.formula('as.factor(phe5$Incident_T2D) ~ scale(beta_324cg[i, ]) + phe5$Age + phe5$Gender + phe5$BMI + phe5$Glucose + phe5$CD8T + phe5$CD4T + phe5$NK + phe5$Bcell + phe5$Mono + phe5$Gran + phe5$PC1_cp + phe5$PC2_cp + phe5$PC3_cp + phe5$PC4_cp + phe5$PC5_cp + phe5$PC6_cp + phe5$PC7_cp + phe5$PC8_cp + phe5$PC9_cp + phe5$PC10_cp + phe5$PC11_cp + phe5$PC12_cp + phe5$PC13_cp + phe5$PC14_cp + phe5$PC15_cp + phe5$PC16_cp + phe5$PC17_cp + phe5$PC18_cp + phe5$PC19_cp + phe5$PC20_cp + phe5$PC21_cp + phe5$PC22_cp + phe5$PC23_cp + phe5$PC24_cp + phe5$PC25_cp + phe5$PC26_cp + phe5$PC27_cp + phe5$PC28_cp + phe5$PC29_cp + phe5$PC30_cp')

nvar = nrow(beta_324cg)

#define the dimension of matrix
res_bmi_glu=matrix(ncol=4, nrow=nvar)

timestamp()

colnames(res_bmi_glu) =c('Estimate','Std.Error', 'z_value','P')
if(!is.factor(phe5$Incident_T2D)){
  phe5$Incident_T2D = as.factor(as.character(phe5$Incident_T2D))
}
for(i in 1:nvar) {
  tryCatch({fit = summary(glm(lfla_bmi_glu, family = binomial))}, error = function(error) {return(NA)})
  if(!exists("fit")){
    res_bmi_glu[i,] = rep(NA,4)
  }else{
    res_bmi_glu[i,] = tryCatch({fit$coefficients[2,]},error = function(error) {return(rep(NA,4))})
    rm(fit)
  }
}

timestamp()
rownames(res_bmi_glu) <- rownames(beta_324cg)
write.table(res_bmi_glu, file="results_LOLIPOP850K_324cg_adjusted_bmi_glu.txt")

#______________

#Sensitivity analysis, adjust for bmi, hba1c

#define formula
lfla_bmi_hba1c=as.formula('as.factor(phe5$Incident_T2D) ~ scale(beta_324cg[i, ]) + phe5$Age + phe5$Gender + phe5$BMI + phe5$HbA1c + phe5$CD8T + phe5$CD4T + phe5$NK + phe5$Bcell + phe5$Mono + phe5$Gran + phe5$PC1_cp + phe5$PC2_cp + phe5$PC3_cp + phe5$PC4_cp + phe5$PC5_cp + phe5$PC6_cp + phe5$PC7_cp + phe5$PC8_cp + phe5$PC9_cp + phe5$PC10_cp + phe5$PC11_cp + phe5$PC12_cp + phe5$PC13_cp + phe5$PC14_cp + phe5$PC15_cp + phe5$PC16_cp + phe5$PC17_cp + phe5$PC18_cp + phe5$PC19_cp + phe5$PC20_cp + phe5$PC21_cp + phe5$PC22_cp + phe5$PC23_cp + phe5$PC24_cp + phe5$PC25_cp + phe5$PC26_cp + phe5$PC27_cp + phe5$PC28_cp + phe5$PC29_cp + phe5$PC30_cp')

nvar = nrow(beta_324cg)

#define the dimension of matrix
res_bmi_hba1c=matrix(ncol=4, nrow=nvar)

timestamp()

colnames(res_bmi_hba1c) =c('Estimate','Std.Error', 'z_value','P')
if(!is.factor(phe5$Incident_T2D)){
  phe5$Incident_T2D = as.factor(as.character(phe5$Incident_T2D))
}
for(i in 1:nvar) {
  tryCatch({fit = summary(glm(lfla_bmi_hba1c, family = binomial))}, error = function(error) {return(NA)})
  if(!exists("fit")){
    res_bmi_hba1c[i,] = rep(NA,4)
  }else{
    res_bmi_hba1c[i,] = tryCatch({fit$coefficients[2,]},error = function(error) {return(rep(NA,4))})
    rm(fit)
  }
}

timestamp()
rownames(res_bmi_hba1c) <- rownames(beta_324cg)
write.table(res_bmi_hba1c, file="results_LOLIPOP850K_324cg_adjusted_bmi_hba1c.txt")

#______________

#Sensitivity analysis, adjust for bmi

#define formula
lfla_bmi=as.formula('as.factor(phe5$Incident_T2D) ~ scale(beta_324cg[i, ]) + phe5$Age + phe5$Gender + phe5$BMI + phe5$CD8T + phe5$CD4T + phe5$NK + phe5$Bcell + phe5$Mono + phe5$Gran + phe5$PC1_cp + phe5$PC2_cp + phe5$PC3_cp + phe5$PC4_cp + phe5$PC5_cp + phe5$PC6_cp + phe5$PC7_cp + phe5$PC8_cp + phe5$PC9_cp + phe5$PC10_cp + phe5$PC11_cp + phe5$PC12_cp + phe5$PC13_cp + phe5$PC14_cp + phe5$PC15_cp + phe5$PC16_cp + phe5$PC17_cp + phe5$PC18_cp + phe5$PC19_cp + phe5$PC20_cp + phe5$PC21_cp + phe5$PC22_cp + phe5$PC23_cp + phe5$PC24_cp + phe5$PC25_cp + phe5$PC26_cp + phe5$PC27_cp + phe5$PC28_cp + phe5$PC29_cp + phe5$PC30_cp')

nvar = nrow(beta_324cg)

#define the dimension of matrix
res_bmi=matrix(ncol=4, nrow=nvar)

timestamp()

colnames(res_bmi) =c('Estimate','Std.Error', 'z_value','P')  
if(!is.factor(phe5$Incident_T2D)){
  phe5$Incident_T2D = as.factor(as.character(phe5$Incident_T2D))
}
for(i in 1:nvar) {
  tryCatch({fit = summary(glm(lfla_bmi, family = binomial))}, error = function(error) {return(NA)})
  if(!exists("fit")){
    res_bmi[i,] = rep(NA,4)
  }else{
    res_bmi[i,] = tryCatch({fit$coefficients[2,]},error = function(error) {return(rep(NA,4))})
    rm(fit)
  }
}

timestamp()
rownames(res_bmi) <- rownames(beta_324cg)
write.table(res_bmi, file="results_LOLIPOP850K_324cg_adjusted_bmi.txt")

#__________________________

#define formula
lfla_glu=as.formula('as.factor(phe5$Incident_T2D) ~ scale(beta_324cg[i, ]) + phe5$Age + phe5$Gender + phe5$Glucose + phe5$CD8T + phe5$CD4T + phe5$NK + phe5$Bcell + phe5$Mono + phe5$Gran + phe5$PC1_cp + phe5$PC2_cp + phe5$PC3_cp + phe5$PC4_cp + phe5$PC5_cp + phe5$PC6_cp + phe5$PC7_cp + phe5$PC8_cp + phe5$PC9_cp + phe5$PC10_cp + phe5$PC11_cp + phe5$PC12_cp + phe5$PC13_cp + phe5$PC14_cp + phe5$PC15_cp + phe5$PC16_cp + phe5$PC17_cp + phe5$PC18_cp + phe5$PC19_cp + phe5$PC20_cp + phe5$PC21_cp + phe5$PC22_cp + phe5$PC23_cp + phe5$PC24_cp + phe5$PC25_cp + phe5$PC26_cp + phe5$PC27_cp + phe5$PC28_cp + phe5$PC29_cp + phe5$PC30_cp')

nvar = nrow(beta_324cg)

#define the dimension of matrix
res_glu=matrix(ncol=4, nrow=nvar)

timestamp()

colnames(res_glu) =c('Estimate','Std.Error', 'z_value','P')
if(!is.factor(phe5$Incident_T2D)){
  phe5$Incident_T2D = as.factor(as.character(phe5$Incident_T2D))
}
for(i in 1:nvar) {
  tryCatch({fit = summary(glm(lfla_glu, family = binomial))}, error = function(error) {return(NA)})
  if(!exists("fit")){
    res_glu[i,] = rep(NA,4)
  }else{
    res_glu[i,] = tryCatch({fit$coefficients[2,]},error = function(error) {return(rep(NA,4))})
    rm(fit)
  }
}

timestamp()
rownames(res_glu) <- rownames(beta_324cg)
write.table(res_glu, file="results_LOLIPOP850K_324cg_adjusted_glu.txt")

#______________

#Sensitivity analysis, adjust for bmi, hba1c

#define formula
lfla_hba1c=as.formula('as.factor(phe5$Incident_T2D) ~ scale(beta_324cg[i, ]) + phe5$Age + phe5$Gender + phe5$HbA1c + phe5$CD8T + phe5$CD4T + phe5$NK + phe5$Bcell + phe5$Mono + phe5$Gran + phe5$PC1_cp + phe5$PC2_cp + phe5$PC3_cp + phe5$PC4_cp + phe5$PC5_cp + phe5$PC6_cp + phe5$PC7_cp + phe5$PC8_cp + phe5$PC9_cp + phe5$PC10_cp + phe5$PC11_cp + phe5$PC12_cp + phe5$PC13_cp + phe5$PC14_cp + phe5$PC15_cp + phe5$PC16_cp + phe5$PC17_cp + phe5$PC18_cp + phe5$PC19_cp + phe5$PC20_cp + phe5$PC21_cp + phe5$PC22_cp + phe5$PC23_cp + phe5$PC24_cp + phe5$PC25_cp + phe5$PC26_cp + phe5$PC27_cp + phe5$PC28_cp + phe5$PC29_cp + phe5$PC30_cp')

nvar = nrow(beta_324cg)

#define the dimension of matrix
res_hba1c=matrix(ncol=4, nrow=nvar)

timestamp()

colnames(res_hba1c) =c('Estimate','Std.Error', 'z_value','P')
if(!is.factor(phe5$Incident_T2D)){
  phe5$Incident_T2D = as.factor(as.character(phe5$Incident_T2D))
}
for(i in 1:nvar) {
  tryCatch({fit = summary(glm(lfla_hba1c, family = binomial))}, error = function(error) {return(NA)})
  if(!exists("fit")){
    res_hba1c[i,] = rep(NA,4)
  }else{
    res_hba1c[i,] = tryCatch({fit$coefficients[2,]},error = function(error) {return(rep(NA,4))})
    rm(fit)
  }
}

timestamp()
rownames(res_hba1c) <- rownames(beta_324cg)
write.table(res_hba1c, file="results_LOLIPOP850K_324cg_adjusted_hba1c.txt")

#______________

#Sensitivity analysis, adjust for bmi, hba1c

#define formula
lfla_glu_hba1c=as.formula('as.factor(phe5$Incident_T2D) ~ scale(beta_324cg[i, ]) + phe5$Age + phe5$Gender + phe5$Glucose + phe5$HbA1c + phe5$CD8T + phe5$CD4T + phe5$NK + phe5$Bcell + phe5$Mono + phe5$Gran + phe5$PC1_cp + phe5$PC2_cp + phe5$PC3_cp + phe5$PC4_cp + phe5$PC5_cp + phe5$PC6_cp + phe5$PC7_cp + phe5$PC8_cp + phe5$PC9_cp + phe5$PC10_cp + phe5$PC11_cp + phe5$PC12_cp + phe5$PC13_cp + phe5$PC14_cp + phe5$PC15_cp + phe5$PC16_cp + phe5$PC17_cp + phe5$PC18_cp + phe5$PC19_cp + phe5$PC20_cp + phe5$PC21_cp + phe5$PC22_cp + phe5$PC23_cp + phe5$PC24_cp + phe5$PC25_cp + phe5$PC26_cp + phe5$PC27_cp + phe5$PC28_cp + phe5$PC29_cp + phe5$PC30_cp')

nvar = nrow(beta_324cg)

#define the dimension of matrix
res_glu_hba1c=matrix(ncol=4, nrow=nvar)

timestamp()

colnames(res_glu_hba1c) =c('Estimate','Std.Error', 'z_value','P')
if(!is.factor(phe5$Incident_T2D)){
  phe5$Incident_T2D = as.factor(as.character(phe5$Incident_T2D))
}
for(i in 1:nvar) {
  tryCatch({fit = summary(glm(lfla_glu_hba1c, family = binomial))}, error = function(error) {return(NA)})
  if(!exists("fit")){
    res_glu_hba1c[i,] = rep(NA,4)
  }else{
    res_glu_hba1c[i,] = tryCatch({fit$coefficients[2,]},error = function(error) {return(rep(NA,4))})
    rm(fit)
  }
}

timestamp()
rownames(res_glu_hba1c) <- rownames(beta_324cg)
write.table(res_glu_hba1c, file="results_LOLIPOP850K_324cg_adjusted_glu_hba1c.txt")

#______________

#define formula
lfla_noadjustment_age_sex_6wbc_30cp=as.formula('as.factor(phe5$Incident_T2D) ~ scale(beta_324cg[i, ]) + phe5$Age + phe5$Gender + phe5$CD8T + phe5$CD4T + phe5$NK + phe5$Bcell + phe5$Mono + phe5$Gran + phe5$PC1_cp + phe5$PC2_cp + phe5$PC3_cp + phe5$PC4_cp + phe5$PC5_cp + phe5$PC6_cp + phe5$PC7_cp + phe5$PC8_cp + phe5$PC9_cp + phe5$PC10_cp + phe5$PC11_cp + phe5$PC12_cp + phe5$PC13_cp + phe5$PC14_cp + phe5$PC15_cp + phe5$PC16_cp + phe5$PC17_cp + phe5$PC18_cp + phe5$PC19_cp + phe5$PC20_cp + phe5$PC21_cp + phe5$PC22_cp + phe5$PC23_cp + phe5$PC24_cp + phe5$PC25_cp + phe5$PC26_cp + phe5$PC27_cp + phe5$PC28_cp + phe5$PC29_cp + phe5$PC30_cp')

nvar = nrow(beta_324cg)

#define the dimension of matrix
res_noadjustment_age_sex_6wbc_30cp=matrix(ncol=4, nrow=nvar)

timestamp()

colnames(res_noadjustment_age_sex_6wbc_30cp) =c('Estimate','Std.Error', 'z_value','P')  
if(!is.factor(phe5$Incident_T2D)){
  phe5$Incident_T2D = as.factor(as.character(phe5$Incident_T2D))
}
for(i in 1:nvar) {
  tryCatch({fit = summary(glm(lfla_noadjustment_age_sex_6wbc_30cp, family = binomial))}, error = function(error) {return(NA)})
  if(!exists("fit")){
    res_noadjustment_age_sex_6wbc_30cp[i,] = rep(NA,4)
  }else{
    res_noadjustment_age_sex_6wbc_30cp[i,] = tryCatch({fit$coefficients[2,]},error = function(error) {return(rep(NA,4))})
    rm(fit)
  }
}

timestamp()
rownames(res_noadjustment_age_sex_6wbc_30cp) <- rownames(beta_324cg)
write.table(res_noadjustment_age_sex_6wbc_30cp, file="results_LOLIPOP850K_324cg_adjusted_noadjustment_age_sex_6wbc_30cp.txt")

#__________________________


#sensitivity analysis in nonprediabetic participant

#grab the sample name of participants that are nonprediabetic
samples_nonprediabetic <- phe7$Row.names

#filter for methylation data of nonprediabetic participants
beta_nonprediabetic_324cg <- beta_324cg[, colnames(beta_324cg) %in% samples_nonprediabetic]

#define formula
lfla_nonprediabetic=as.formula("phe7$Incident_T2D ~ beta_nonprediabetic_324cg[i, ] + phe7$Age + phe7$Gender + phe7$CD8T + phe7$CD4T + phe7$NK + phe7$Bcell + phe7$Mono + phe7$Gran + phe7$PC1_cp + phe7$PC2_cp + phe7$PC3_cp + phe7$PC4_cp + phe7$PC5_cp + phe7$PC6_cp + phe7$PC7_cp + phe7$PC8_cp + phe7$PC9_cp + phe7$PC10_cp + phe7$PC11_cp + phe7$PC12_cp + phe7$PC13_cp + phe7$PC14_cp + phe7$PC15_cp + phe7$PC16_cp + phe7$PC17_cp + phe7$PC18_cp + phe7$PC19_cp + phe7$PC20_cp + phe7$PC21_cp + phe7$PC22_cp + phe7$PC23_cp + phe7$PC24_cp + phe7$PC25_cp + phe7$PC26_cp + phe7$PC27_cp + phe7$PC28_cp + phe7$PC29_cp + phe7$PC30_cp")

nvar = nrow(beta_nonprediabetic_324cg)
#define dimension of martix
res2=matrix(ncol=4, nrow=nvar)

timestamp()

  colnames(res2) =c('Estimate','Std.Error', 'z_value','P')
  if(!is.factor(phe7$Incident_T2D)){
    phe7$Incident_T2D = as.factor(as.character(phe7$Incident_T2D))
  }
  for(i in 1:nvar) {
    tryCatch({fit = summary(glm(lfla_nonprediabetic, family = binomial))}, error = function(error) {return(NA)})
    if(!exists("fit")){
      res2[i,] = rep(NA,4)
    }else{
      res2[i,] = tryCatch({fit$coefficients[2,]},error = function(error) {return(rep(NA,4))})
      rm(fit)
    }
  }

timestamp()
rownames(res2) <- rownames(beta_nonprediabetic_324cg)

write.table(res2, file="results_LOLIPOP850K_324cg_nonprediabetic.txt")

#_______sensitivity analysis done__________
