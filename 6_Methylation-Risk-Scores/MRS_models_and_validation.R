library(corrplot)
library(pscl)# for pR2 function to calculate pseudo R2
library(patchwork)

#load housemen
housemen <- read.table("~/OneDrive - Nanyang Technological University/Marie_TOAST_analysis/Methylation_analysis/LOLIPOP_850K/output/Houseman_unconstrained.txt", header = T)

controlpc <- read.table("~/OneDrive - Nanyang Technological University/Marie_TOAST_analysis/Methylation_analysis/LOLIPOP_850K/output/ctrlprobes_PC1to30.txt", header = T)

housemen_ctrl <- merge(housemen, controlpc, by = "row.names")

sentrixID <- read.table("~/OneDrive - Nanyang Technological University/Marie_TOAST_analysis/Methylation_analysis/LOLIPOP_850K/output/EPIC_LOLIPOP_clinpath.txt", header = T)

sentrixID_housemen <- merge(housemen_ctrl,sentrixID[,c(1:2)], by.x = "Row.names", by.y ="SentrixID")

#PRS
PRS <- read.csv("LOLIPOP_GSA_PRS_MRS_MetRS.csv")
#PRS2 <- PRS[,c(1:8,17,39,40,9,10,13,25,27,28,30:38)]
PRStemp <- PRS %>%
  mutate(
    Sex2 = case_when(
      Sex == 2 ~ 1, 
      Sex == 1 ~ 0
    ),
    Smoker1 = case_when(
      Smoker == 2 ~ 0, 
      Smoker == 1 ~ 1, 
      Smoker == 0 ~ 0
    ),
    Smoker2 = case_when(
      Smoker == 2 ~ 1, 
      Smoker == 1 ~ 0,
      Smoker == 0 ~ 0
    )
    )

PRS3 <- merge(PRStemp,sentrixID_housemen, by.x="Sample", by.y="SampleID")

#MRS
MRS_lasso1se <- read.csv("GSA_MRS_Lasso1se_all.csv")

#Merge all riskscore
clinical_v2 <- merge(PRS3, MRS_lasso1se, by = "Sample")

#scale the riskscore
clinical_v3 = clinical_v2
clinical_v3$PRS_Ge_et_al <- scale(clinical_v3$PRS_Ge_et_al)
clinical_v3$MRS_lasso1se <- scale(clinical_v3$MRS_lasso1se)
clinical_v3$T2D_status <- as.factor(clinical_v3$T2D_status)
clinical_v3$Sex2 <- as.factor(clinical_v3$Sex2)


#_______________

#quartlie the individual RS base on control
#merge sample name and RS value
#or merge sample name and scale MRS value

clinical_v3_control = clinical_v3 %>% filter(T2D_status == "0")

# define the samples into their quartile for each RS
quartile_clinical_v3_PRS_control <- quantile(clinical_v3_control$PRS_Ge_et_al, q=4)
quartile_clinical_v3_MRS_control <- quantile(clinical_v3_control$MRS_lasso1se, q=4)


#Use case when to assign samples to different quartile

#classify RS
#RS reclassification 
clinical_v3 <- clinical_v3 %>%
  mutate(
    PRS_quartile = case_when(
      (clinical_v3$PRS_Ge_et_al > quartile_clinical_v3_PRS_control[4]  ~ 4),
      (clinical_v3$PRS_Ge_et_al > quartile_clinical_v3_PRS_control[3] & clinical_v3$PRS_Ge_et_al <= quartile_clinical_v3_PRS_control[4]  ~ 3),
      (clinical_v3$PRS_Ge_et_al > quartile_clinical_v3_PRS_control[2] & clinical_v3$PRS_Ge_et_al <= quartile_clinical_v3_PRS_control[3]  ~ 2),
      (clinical_v3$PRS_Ge_et_al <= quartile_clinical_v3_PRS_control[2]  ~ 1)
    ),
    MRS_quartile = case_when(
      (clinical_v3$MRS_lasso1se > quartile_clinical_v3_MRS_control[4]  ~ 4),
      (clinical_v3$MRS_lasso1se > quartile_clinical_v3_MRS_control[3] & clinical_v3$MRS_lasso1se <= quartile_clinical_v3_MRS_control[4]  ~ 3),
      (clinical_v3$MRS_lasso1se > quartile_clinical_v3_MRS_control[2] & clinical_v3$MRS_lasso1se <= quartile_clinical_v3_MRS_control[3]  ~ 2),
      (clinical_v3$MRS_lasso1se <= quartile_clinical_v3_MRS_control[2]  ~ 1)
    )
  )

#Check the number of case and control in each quartile
count(clinical_v3, T2D_status)
count(clinical_v3, T2D_status)

count(clinical_v3, PRS_quartile, T2D_status)
count(clinical_v3, MRS_quartile, T2D_status)
count_clinical_v3_quartile_summary <- data.frame(count(clinical_v3, PRS_quartile, T2D_status),count(clinical_v3, MRS_quartile, T2D_status))

#__________________

#as.factor the quartile and decile
#quartlie the RS
clinical_v3$PRS_quartile <- as.factor(clinical_v3$PRS_quartile)
clinical_v3$MRS_quartile <- as.factor(clinical_v3$MRS_quartile)

str(clinical_v3)

#Define prediabetic and nonprediabetic
#Define
#Pre-diabetic, glucose >=6.0mmol/L, HbA1c >= 6.0%
#High BMI (Asian BMI standard) >= 23 
#High waist >=90 for men, >=80 for woman
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3125562/

clinical_v3 <- clinical_v3 %>%
  mutate(
    prediabetic = case_when(
      Glu >= 6.0 | HbA1c >= 6.0 ~ 1,      TRUE ~ 0))

clinical_v3 %>% count(prediabetic)

clinical_v3 <- clinical_v3 %>%
  mutate(
    High_BMI_traditional = case_when(
      (BMI >= 30 ~ 2),
      (BMI >= 25.0 & BMI < 30 ~ 1),
      (TRUE ~ 0)))

#as.factor high bmi and prediabetic
clinical_v3$High_BMI_traditional <- as.factor(clinical_v3$High_BMI_traditional)
clinical_v3$prediabetic <- as.factor(clinical_v3$prediabetic)

#check number for prediabetic and non prediabetic
clinical_v3 %>% dplyr::count(prediabetic,T2D_status)
clinical_v3 %>% dplyr::count(prediabetic,High_BMI_traditional)

clinical_v3 %>% dplyr::count(T2D_status,PRS_quartile)
clinical_v3 %>% dplyr::count(T2D_status,MRS_quartile)

#________________

#basic
q_basic <- glm(formula = T2D_status ~ Age + Sex2
               family = binomial(link = "logit"), data = clinical_v3, na.action = na.omit)
summary(q_basic)$coefficients
q_basic_coef <- summary(q_basic)$coefficients
pR2(q_basic)

clinical_v3$phat0<-predict(q_basic, type='response')
q_basic.ROC<-roc(T2D_status ~ phat0, data=clinical_v3, print.auc=T)
q_basic.AUC<-as.numeric(q_basic.ROC$auc)
q_basic.AUC

ci.auc(clinical_v3$T2D_status, clinical_v3$phat0, conf.level = 0.95) 

#_________________

#PRS
q_basic_PRS <- glm(formula = T2D_status ~ Age + Sex2 + PRS_Ge_et_al, 
             family = binomial(link = "logit"), data = clinical_v3, na.action = na.omit)
summary(q_basic_PRS)$coefficients
q_basic_PRS_coef <- summary(q_basic_PRS)$coefficients
pR2(q_basic_PRS)

clinical_v3$phat0_PRS<-predict(q_basic_PRS, type='response')
q_basic_PRS.ROC<-roc(T2D_status ~ phat0_PRS, data=clinical_v3, print.auc=T)
q_basic_PRS.AUC<-as.numeric(q_basic_PRS.ROC$auc)
q_basic_PRS.AUC

ci.auc(clinical_v3$T2D_status, clinical_v3$phat0_PRS, conf.level = 0.95) 

#_____________________

#MRS_lasso1se 42 CpG
q_basic_MRS_lasso1se <- glm(formula = T2D_status ~ Age + Sex2 + scale(CD8T) + scale(CD4T) + scale(NK) + scale(Bcell) + scale(Mono) + scale(Gran)  + scale(PC1_cp) +scale(PC2_cp) +scale(PC3_cp) +scale(PC4_cp) +scale(PC5_cp) +scale(PC6_cp) +scale(PC7_cp) +scale(PC8_cp) +scale(PC9_cp) +scale(PC10_cp) +scale(PC11_cp) +scale(PC12_cp) +scale(PC13_cp) +scale(PC14_cp) +scale(PC15_cp) +scale(PC16_cp) +scale(PC17_cp) +scale(PC18_cp) +scale(PC19_cp) +scale(PC20_cp) +scale(PC21_cp) +scale(PC22_cp) +scale(PC23_cp) +scale(PC24_cp) +scale(PC25_cp) +scale(PC26_cp) +scale(PC27_cp) +scale(PC28_cp) +scale(PC29_cp) +scale(PC30_cp)+ MRS_lasso1se, 
             family = binomial(link = "logit"), data = clinical_v3, na.action = na.omit)
summary(q_basic_MRS_lasso1se)$coefficients
q_basic_MRS_lasso1se_coef <- summary(q_basic_MRS_lasso1se)$coefficients
pR2(q_basic_MRS_lasso1se)

clinical_v3$phat0_MRS_lasso1se<-predict(q_basic_MRS_lasso1se, type='response')
q_basic_MRS_lasso1se.ROC<-roc(T2D_status ~ phat0_MRS_lasso1se, data=clinical_v3, print.auc=T)
q_basic_MRS_lasso1se.AUC<-as.numeric(q_basic_MRS_lasso1se.ROC$auc)
q_basic_MRS_lasso1se.AUC

ci.auc(clinical_v3$T2D_status, clinical_v3$phat0_MRS_lasso1se, conf.level = 0.95) 

#_______________

#MRS + PRS
q_basic_PRSMRS_lasso1se <- glm(formula = T2D_status ~ Age + Sex2 + scale(CD8T) + scale(CD4T) + scale(NK) + scale(Bcell) + scale(Mono) + scale(Gran) + scale(PC1_cp) +scale(PC2_cp) +scale(PC3_cp) +scale(PC4_cp) +scale(PC5_cp) +scale(PC6_cp) +scale(PC7_cp) +scale(PC8_cp) +scale(PC9_cp) +scale(PC10_cp) +scale(PC11_cp) +scale(PC12_cp) +scale(PC13_cp) +scale(PC14_cp) +scale(PC15_cp) +scale(PC16_cp) +scale(PC17_cp) +scale(PC18_cp) +scale(PC19_cp) +scale(PC20_cp) +scale(PC21_cp) +scale(PC22_cp) +scale(PC23_cp) +scale(PC24_cp) +scale(PC25_cp) +scale(PC26_cp) +scale(PC27_cp) +scale(PC28_cp) +scale(PC29_cp) +scale(PC30_cp) + PRS_Ge_et_al + MRS_lasso1se, 
                            family = binomial(link = "logit"), data = clinical_v3, na.action = na.omit)
summary(q_basic_PRSMRS_lasso1se)$coefficients
q_basic_PRSMRS_lasso1se_coef <- summary(q_basic_PRSMRS_lasso1se)$coefficients
pR2(q_basic_PRSMRS_lasso1se)

clinical_v3$phat0<-predict(q_basic_PRSMRS_lasso1se, type='response')
q_basic_PRSMRS_lasso1se.ROC<-roc(T2D_status ~ phat0, data=clinical_v3, print.auc=T)
q_basic_PRSMRS_lasso1se.AUC<-as.numeric(q_basic_PRSMRS_lasso1se.ROC$auc)
q_basic_PRSMRS_lasso1se.AUC

ci.auc(clinical_v3$T2D_status, clinical_v3$phat0, conf.level = 0.95) 

#_______________

#Quartile analysis
#PRS
q_basic_PRS_quartile <- glm(formula = T2D_status ~ Age + Sex2 + scale(CD8T) + scale(CD4T) + scale(NK) + scale(Bcell) + scale(Mono) + scale(Gran) + PRS_quartile, 
                            family = binomial(link = "logit"), data = clinical_v3, na.action = na.omit)
summary(q_basic_PRS_quartile)$coefficients
q_basic_PRS_quartile_coef <- summary(q_basic_PRS_quartile)$coefficients
pR2(q_basic_PRS_quartile)

clinical_v3$phat_q_basic_PRS_quartile<-predict(q_basic_PRS_quartile, type='response')
q_basic_PRS_quartile.ROC<-roc(T2D_status ~ phat_q_basic_PRS_quartile, data=clinical_v3, print.auc=T)
q_basic_PRS_quartile.AUC<-as.numeric(q_basic_PRS_quartile.ROC$auc)
q_basic_PRS_quartile.AUC

ci.auc(clinical_v3$T2D_status, clinical_v3$phat_q_basic_PRS_quartile, conf.level = 0.95) 

#________________

#MRS_42cg_quartile
q_basic_MRS_quartile <- glm(formula = T2D_status ~ Age + Sex2 + scale(CD8T) + scale(CD4T) + scale(NK) + scale(Bcell) + scale(Mono) + scale(Gran) + scale(PC1_cp) +scale(PC2_cp) +scale(PC3_cp) +scale(PC4_cp) +scale(PC5_cp) +scale(PC6_cp) +scale(PC7_cp) +scale(PC8_cp) +scale(PC9_cp) +scale(PC10_cp) +scale(PC11_cp) +scale(PC12_cp) +scale(PC13_cp) +scale(PC14_cp) +scale(PC15_cp) +scale(PC16_cp) +scale(PC17_cp) +scale(PC18_cp) +scale(PC19_cp) +scale(PC20_cp) +scale(PC21_cp) +scale(PC22_cp) +scale(PC23_cp) +scale(PC24_cp) +scale(PC25_cp) +scale(PC26_cp) +scale(PC27_cp) +scale(PC28_cp) +scale(PC29_cp) +scale(PC30_cp) + MRS_quartile, 
                            family = binomial(link = "logit"), data = clinical_v3, na.action = na.omit)
summary(q_basic_MRS_quartile)$coefficients
q_basic_MRS_quartile_coef <- summary(q_basic_MRS_quartile)$coefficients
pR2(q_basic_MRS_quartile)

clinical_v3$phat_q_basic_MRS_quartile<-predict(q_basic_MRS_quartile, type='response')
q_basic_MRS_quartile.ROC<-roc(T2D_status ~ phat_q_basic_MRS_quartile, data=clinical_v3, print.auc=T)
q_basic_MRS_quartile.AUC<-as.numeric(q_basic_MRS_quartile.ROC$auc)
q_basic_MRS_quartile.AUC

ci.auc(clinical_v3$T2D_status, clinical_v3$phat_q_basic_MRS_quartile, conf.level = 0.95) 

#________________

#Filter for prediabetic
clinical_v3_prediabetic <- clinical_v3 %>%
  filter(prediabetic ==1)

#__________________

#
clinical_v3_prediabetic_control = clinical_v3_prediabetic %>% filter(T2D_status == "0")

# define the samples into their quartile for each RS
quartile_clinical_v3_prediabetic_PRS_control <- quantile(clinical_v3_prediabetic_control$PRS_Ge_et_al, q=4)
quartile_clinical_v3_prediabetic_MRS_control <- quantile(clinical_v3_prediabetic_control$MRS_lasso1se, q=4)

#Use case when to assign samples to different quartile

#classify RS
#RS reclassification 
clinical_v3_prediabetic <- clinical_v3_prediabetic %>%
  mutate(
    PRS_prediabetic_quartile = case_when(
      (clinical_v3_prediabetic$PRS_Ge_et_al > quartile_clinical_v3_prediabetic_PRS_control[4]  ~ 4),
      (clinical_v3_prediabetic$PRS_Ge_et_al > quartile_clinical_v3_prediabetic_PRS_control[3] & clinical_v3_prediabetic$PRS_Ge_et_al <= quartile_clinical_v3_prediabetic_PRS_control[4]  ~ 3),
      (clinical_v3_prediabetic$PRS_Ge_et_al > quartile_clinical_v3_prediabetic_PRS_control[2] & clinical_v3_prediabetic$PRS_Ge_et_al <= quartile_clinical_v3_prediabetic_PRS_control[3]  ~ 2),
      (clinical_v3_prediabetic$PRS_Ge_et_al <= quartile_clinical_v3_prediabetic_PRS_control[2]  ~ 1)
    ),
    MRS_prediabetic_quartile = case_when(
      (clinical_v3_prediabetic$MRS_lasso1se > quartile_clinical_v3_prediabetic_MRS_control[4]  ~ 4),
      (clinical_v3_prediabetic$MRS_lasso1se > quartile_clinical_v3_prediabetic_MRS_control[3] & clinical_v3_prediabetic$MRS_lasso1se <= quartile_clinical_v3_prediabetic_MRS_control[4]  ~ 3),
      (clinical_v3_prediabetic$MRS_lasso1se > quartile_clinical_v3_prediabetic_MRS_control[2] & clinical_v3_prediabetic$MRS_lasso1se <= quartile_clinical_v3_prediabetic_MRS_control[3]  ~ 2),
      (clinical_v3_prediabetic$MRS_lasso1se <= quartile_clinical_v3_prediabetic_MRS_control[2]  ~ 1)
    )
  )

#Check the number of case and control in each quartile
count(clinical_v3_prediabetic, T2D_status)
count(clinical_v3_prediabetic, T2D_status)

count(clinical_v3_prediabetic, PRS_prediabetic_quartile, T2D_status)
count(clinical_v3_prediabetic, MRS_prediabetic_quartile, T2D_status)

count_clinical_v3_prediabetic_quartile_summary <- data.frame(count(clinical_v3_prediabetic, PRS_prediabetic_quartile, T2D_status),count(clinical_v3_prediabetic, MRS_prediabetic_quartile, T2D_status))

#__________________

#as.factor the quartile and decile
#quartlie the RS
clinical_v3_prediabetic$PRS_prediabetic_quartile <- as.factor(clinical_v3_prediabetic$PRS_prediabetic_quartile)
clinical_v3_prediabetic$MRS_prediabetic_quartile <- as.factor(clinical_v3_prediabetic$MRS_prediabetic_quartile)

str(clinical_v3_prediabetic)

#________________________

#filter for nonprediabetic
clinical_v3_nonprediabetic <- clinical_v3 %>%
  filter(prediabetic ==0)

#select for nonprediabetic control
clinical_v3_nonprediabetic_control = clinical_v3_nonprediabetic %>% filter(T2D_status == "0")

# define the samples into their quartile for each RS
quartile_clinical_v3_nonprediabetic_PRS_control <- quantile(clinical_v3_nonprediabetic_control$PRS_Ge_et_al, q=4)
quartile_clinical_v3_nonprediabetic_MRS_control <- quantile(clinical_v3_nonprediabetic_control$MRS_lasso1se, q=4)

#Use case when to assign samples to different quartile

#classify RS
clinical_v3_nonprediabetic <- clinical_v3_nonprediabetic %>%
  mutate(
    PRS_nonprediabetic_quartile = case_when(
      (clinical_v3_nonprediabetic$PRS_Ge_et_al > quartile_clinical_v3_nonprediabetic_PRS_control[4]  ~ 4),
      (clinical_v3_nonprediabetic$PRS_Ge_et_al > quartile_clinical_v3_nonprediabetic_PRS_control[3] & clinical_v3_nonprediabetic$PRS_Ge_et_al <= quartile_clinical_v3_nonprediabetic_PRS_control[4]  ~ 3),
      (clinical_v3_nonprediabetic$PRS_Ge_et_al > quartile_clinical_v3_nonprediabetic_PRS_control[2] & clinical_v3_nonprediabetic$PRS_Ge_et_al <= quartile_clinical_v3_nonprediabetic_PRS_control[3]  ~ 2),
      (clinical_v3_nonprediabetic$PRS_Ge_et_al <= quartile_clinical_v3_nonprediabetic_PRS_control[2]  ~ 1)
    ),
    MRS_nonprediabetic_quartile = case_when(
      (clinical_v3_nonprediabetic$MRS_lasso1se > quartile_clinical_v3_nonprediabetic_MRS_control[4]  ~ 4),
      (clinical_v3_nonprediabetic$MRS_lasso1se > quartile_clinical_v3_nonprediabetic_MRS_control[3] & clinical_v3_nonprediabetic$MRS_lasso1se <= quartile_clinical_v3_nonprediabetic_MRS_control[4]  ~ 3),
      (clinical_v3_nonprediabetic$MRS_lasso1se > quartile_clinical_v3_nonprediabetic_MRS_control[2] & clinical_v3_nonprediabetic$MRS_lasso1se <= quartile_clinical_v3_nonprediabetic_MRS_control[3]  ~ 2),
      (clinical_v3_nonprediabetic$MRS_lasso1se <= quartile_clinical_v3_nonprediabetic_MRS_control[2]  ~ 1)
    )
  )

#Check the number of case and control in each quartile
count(clinical_v3_nonprediabetic, T2D_status)
count(clinical_v3_nonprediabetic, T2D_status)

count(clinical_v3_nonprediabetic, PRS_nonprediabetic_quartile, T2D_status)
count(clinical_v3_nonprediabetic, MRS_nonprediabetic_quartile, T2D_status)

count_clinical_v3_nonprediabetic_quartile_summary <- data.frame(count(clinical_v3_nonprediabetic, PRS_nonprediabetic_quartile, T2D_status),count(clinical_v3_nonprediabetic, MRS_nonprediabetic_quartile, T2D_status))

#__________________

#as.factor the quartile
#quartlie the RS
clinical_v3_nonprediabetic$PRS_nonprediabetic_quartile <- as.factor(clinical_v3_nonprediabetic$PRS_nonprediabetic_quartile)
clinical_v3_nonprediabetic$MRS_nonprediabetic_quartile <- as.factor(clinical_v3_nonprediabetic$MRS_nonprediabetic_quartile)

str(clinical_v3_nonprediabetic)

#______________

#logistic regression test for prediabeitc

#Prediabetic positive analysis

#count number of t2d case and controls
clinical_v3_prediabetic %>% count(T2D_status)
clinical_v3_nonprediabetic %>% count(T2D_status)

##Prediabetic analysis
#basic
q_prediabetic_basic_pos <- glm(formula = T2D_status ~ Age + Sex2, 
                               family = binomial(link = "logit"), data = clinical_v3_prediabetic, na.action = na.omit)
summary(q_prediabetic_basic_pos)$coefficients
q_prediabetic_basic_pos_coef <- summary(q_prediabetic_basic_pos)$coefficients
pR2(q_prediabetic_basic_pos)

clinical_v3_prediabetic$phat0<-predict(q_prediabetic_basic_pos, type='response')
q_prediabetic_basic_pos.ROC<-roc(T2D_status ~ phat0, data=clinical_v3_prediabetic, print.auc=T)
q_prediabetic_basic_pos.AUC<-as.numeric(q_prediabetic_basic_pos.ROC$auc)
q_prediabetic_basic_pos.AUC

ci.auc(clinical_v3_prediabetic$T2D_status, clinical_v3_prediabetic$phat0, conf.level = 0.95)

#_________________

#PRS
q_prediabetic_basic_pos_PRS <- glm(formula = T2D_status ~ Age + Sex2 + PRS_Ge_et_al, 
                                   family = binomial(link = "logit"), data = clinical_v3_prediabetic, na.action = na.omit)
summary(q_prediabetic_basic_pos_PRS)$coefficients
q_prediabetic_basic_pos_PRS_coef <- summary(q_prediabetic_basic_pos_PRS)$coefficients
pR2(q_prediabetic_basic_pos_PRS)

clinical_v3_prediabetic$phat0<-predict(q_prediabetic_basic_pos_PRS, type='response')
q_prediabetic_basic_pos_PRS.ROC<-roc(T2D_status ~ phat0, data=clinical_v3_prediabetic, print.auc=T)
q_prediabetic_basic_pos_PRS.AUC<-as.numeric(q_prediabetic_basic_pos_PRS.ROC$auc)
q_prediabetic_basic_pos_PRS.AUC

ci.auc(clinical_v3_prediabetic$T2D_status, clinical_v3_prediabetic$phat0, conf.level = 0.95) 

#_____________________

#MRS_lasso1se 42 CpG
q_prediabetic_basic_pos_MRS_lasso1se <- glm(formula = T2D_status ~ Age + Sex2 + scale(CD8T) + scale(CD4T) + scale(NK) + scale(Bcell) + scale(Mono) + scale(Gran) + scale(PC1_cp) +scale(PC2_cp) +scale(PC3_cp) +scale(PC4_cp) +scale(PC5_cp) +scale(PC6_cp) +scale(PC7_cp) +scale(PC8_cp) +scale(PC9_cp) +scale(PC10_cp) +scale(PC11_cp) +scale(PC12_cp) +scale(PC13_cp) +scale(PC14_cp) +scale(PC15_cp) +scale(PC16_cp) +scale(PC17_cp) +scale(PC18_cp) +scale(PC19_cp) +scale(PC20_cp) +scale(PC21_cp) +scale(PC22_cp) +scale(PC23_cp) +scale(PC24_cp) +scale(PC25_cp) +scale(PC26_cp) +scale(PC27_cp) +scale(PC28_cp) +scale(PC29_cp) +scale(PC30_cp) + MRS_lasso1se, 
                                                family = binomial(link = "logit"), data = clinical_v3_prediabetic, na.action = na.omit)
summary(q_prediabetic_basic_pos_MRS_lasso1se)$coefficients
q_prediabetic_basic_pos_MRS_lasso1se_coef <- summary(q_prediabetic_basic_pos_MRS_lasso1se)$coefficients
pR2(q_prediabetic_basic_pos_MRS_lasso1se)

clinical_v3_prediabetic$phat0<-predict(q_prediabetic_basic_pos_MRS_lasso1se, type='response')
q_prediabetic_basic_pos_MRS_lasso1se.ROC<-roc(T2D_status ~ phat0, data=clinical_v3_prediabetic, print.auc=T)
q_prediabetic_basic_pos_MRS_lasso1se.AUC<-as.numeric(q_prediabetic_basic_pos_MRS_lasso1se.ROC$auc)
q_prediabetic_basic_pos_MRS_lasso1se.AUC

ci.auc(clinical_v3_prediabetic$T2D_status, clinical_v3_prediabetic$phat0, conf.level = 0.95) 

#_______________

#PRS_MRS_lasso1se 42 CpG
q_prediabetic_basic_pos_PRS_MRS_lasso1se <- glm(formula = T2D_status ~ Age + Sex2 + scale(CD8T) + scale(CD4T) + scale(NK) + scale(Bcell) + scale(Mono) + scale(Gran) + scale(PC1_cp) +scale(PC2_cp) +scale(PC3_cp) +scale(PC4_cp) +scale(PC5_cp) +scale(PC6_cp) +scale(PC7_cp) +scale(PC8_cp) +scale(PC9_cp) +scale(PC10_cp) +scale(PC11_cp) +scale(PC12_cp) +scale(PC13_cp) +scale(PC14_cp) +scale(PC15_cp) +scale(PC16_cp) +scale(PC17_cp) +scale(PC18_cp) +scale(PC19_cp) +scale(PC20_cp) +scale(PC21_cp) +scale(PC22_cp) +scale(PC23_cp) +scale(PC24_cp) +scale(PC25_cp) +scale(PC26_cp) +scale(PC27_cp) +scale(PC28_cp) +scale(PC29_cp) +scale(PC30_cp) + PRS_Ge_et_al + MRS_lasso1se, 
                                            family = binomial(link = "logit"), data = clinical_v3_prediabetic, na.action = na.omit)
summary(q_prediabetic_basic_pos_PRS_MRS_lasso1se)$coefficients
q_prediabetic_basic_pos_PRS_MRS_lasso1se_coef <- summary(q_prediabetic_basic_pos_PRS_MRS_lasso1se)$coefficients
pR2(q_prediabetic_basic_pos_PRS_MRS_lasso1se)

clinical_v3_prediabetic$phat0<-predict(q_prediabetic_basic_pos_PRS_MRS_lasso1se, type='response')
q_prediabetic_basic_pos_PRS_MRS_lasso1se.ROC<-roc(T2D_status ~ phat0, data=clinical_v3_prediabetic, print.auc=T)
q_prediabetic_basic_pos_PRS_MRS_lasso1se.AUC<-as.numeric(q_prediabetic_basic_pos_PRS_MRS_lasso1se.ROC$auc)
q_prediabetic_basic_pos_PRS_MRS_lasso1se.AUC

ci.auc(clinical_v3_prediabetic$T2D_status, clinical_v3_prediabetic$phat0, conf.level = 0.95) 

#_______________
#Quartile analysis
#PRS
q_prediabetic_basic_PRS_prediabetic_quartile <- glm(formula = T2D_status ~ Age + Sex2 + PRS_prediabetic_quartile, 
                                        family = binomial(link = "logit"), data = clinical_v3_prediabetic, na.action = na.omit)
summary(q_prediabetic_basic_PRS_prediabetic_quartile)$coefficients
q_prediabetic_basic_PRS_prediabetic_quartile_coef <- summary(q_prediabetic_basic_PRS_prediabetic_quartile)$coefficients
pR2(q_prediabetic_basic_PRS_prediabetic_quartile)

clinical_v3_prediabetic$phat_q_prediabetic_basic_PRS_prediabetic_quartile<-predict(q_prediabetic_basic_PRS_prediabetic_quartile, type='response')
q_prediabetic_basic_PRS_prediabetic_quartile.ROC<-roc(T2D_status ~ phat_q_prediabetic_basic_PRS_prediabetic_quartile, data=clinical_v3_prediabetic, print.auc=T)
q_prediabetic_basic_PRS_prediabetic_quartile.AUC<-as.numeric(q_prediabetic_basic_PRS_prediabetic_quartile.ROC$auc)
q_prediabetic_basic_PRS_prediabetic_quartile.AUC

ci.auc(clinical_v3_prediabetic$T2D_status, clinical_v3_prediabetic$phat_q_prediabetic_basic_PRS_prediabetic_quartile, conf.level = 0.95)

#________________

#MRS_42cg
q_prediabetic_basic_MRS_prediabetic_quartile <- glm(formula = T2D_status ~ Age + Sex2 + scale(CD8T) + scale(CD4T) + scale(NK) + scale(Bcell) + scale(Mono) + scale(Gran) + scale(PC1_cp) +scale(PC2_cp) +scale(PC3_cp) +scale(PC4_cp) +scale(PC5_cp) +scale(PC6_cp) +scale(PC7_cp) +scale(PC8_cp) +scale(PC9_cp) +scale(PC10_cp) +scale(PC11_cp) +scale(PC12_cp) +scale(PC13_cp) +scale(PC14_cp) +scale(PC15_cp) +scale(PC16_cp) +scale(PC17_cp) +scale(PC18_cp) +scale(PC19_cp) +scale(PC20_cp) +scale(PC21_cp) +scale(PC22_cp) +scale(PC23_cp) +scale(PC24_cp) +scale(PC25_cp) +scale(PC26_cp) +scale(PC27_cp) +scale(PC28_cp) +scale(PC29_cp) +scale(PC30_cp) + MRS_prediabetic_quartile, 
                                        family = binomial(link = "logit"), data = clinical_v3_prediabetic, na.action = na.omit)
summary(q_prediabetic_basic_MRS_prediabetic_quartile)$coefficients
q_prediabetic_basic_MRS_prediabetic_quartile_coef <- summary(q_prediabetic_basic_MRS_prediabetic_quartile)$coefficients
pR2(q_prediabetic_basic_MRS_prediabetic_quartile)

clinical_v3_prediabetic$phat_q_prediabetic_basic_MRS_prediabetic_quartile<-predict(q_prediabetic_basic_MRS_prediabetic_quartile, type='response')
q_prediabetic_basic_MRS_prediabetic_quartile.ROC<-roc(T2D_status ~ phat_q_prediabetic_basic_MRS_prediabetic_quartile, data=clinical_v3_prediabetic, print.auc=T)
q_prediabetic_basic_MRS_prediabetic_quartile.AUC<-as.numeric(q_prediabetic_basic_MRS_prediabetic_quartile.ROC$auc)
q_prediabetic_basic_MRS_prediabetic_quartile.AUC

ci.auc(clinical_v3_prediabetic$T2D_status, clinical_v3_prediabetic$phat_q_prediabetic_basic_MRS_prediabetic_quartile, conf.level = 0.95)

#________________
