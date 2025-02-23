library(tidyverse) # to use mutate
library(glmnet) #lasso regression

## Lasso regression
#Read phenotype file that contain sample id and sentrix id
pheno <- read.table("TOAST_T2D_combine_old_newrerunsample_nogenderswap_SPHS.txt", sep = "\t", header = TRUE)
colnames(pheno)[2] <- "Sentrix.ID"
row.names(pheno) <- pheno$Sentrix.ID


#read csv file that contain scaled methylation value of MEC and EPI450K
beta_sig_SPHS_ALL_tv3 <- read.csv("beta_sig_SPHS_ALL_scaled.csv")

#recode gender to sex (1 for male, 2 for female so that it is the same as clinical data)
beta_sig_SPHS_ALL_tv3 <- beta_sig_SPHS_ALL_tv3 %>%
  mutate(Sex2 = case_when(
    Gender == "Male" ~ 0,
    Gender == "Female" ~ 1
  ))

#recode ethnic to sex (1 for chinese, 2 for indian, 3 for malay, alphabetical order)
beta_sig_SPHS_ALL_tv3 <- beta_sig_SPHS_ALL_tv3 %>%
  mutate(Ethnicity = case_when(
    Ethnic == "Chinese" ~ 1,
    Ethnic == "Indian" ~ 2,
    Ethnic == "Malay" ~ 3
  ))

#______________

#reformat data for lasso regression
xfactors <- model.matrix(T2D_status ~ Sex2 + Ethnicity, data=beta_sig_SPHS_ALL_tv3)[, -1]
x <- as.matrix(data.frame(xfactors,beta_sig_SPHS_ALL_tv3[,c(6:137)]))

beta_merge_meth = data.frame(x)
beta_merge_meth_tv3 = data.frame(beta_sig_SPHS_ALL_tv3$T2D_status,beta_sig_SPHS_ALL_tv3$Age,beta_merge_meth)
colnames(beta_merge_meth_tv3)[1] <- "T2D_status"
colnames(beta_merge_meth_tv3)[2] <- "Age"

beta_merge_meth_tv3$T2D_status <- as.factor(beta_merge_meth_tv3$T2D_status)
beta_merge_meth_tv3$Age <- as.numeric(beta_merge_meth_tv3$Age)
beta_merge_meth_tv3$Sex2 <- as.factor(beta_merge_meth_tv3$Sex2)
beta_merge_meth_tv3$Ethnicity <- as.factor(beta_merge_meth_tv3$Ethnicity)


#Lasso regression to calculate fit
fit <- glmnet(x, y=beta_sig_SPHS_ALL_tv3$T2D_status, family = "binomial")

#plot the L1 Norm against the different coefficient
#Each curve corresponds to a variable. It shows the path of its coefficient against the ℓ1-norm of the whole coefficient vector as λvaries. The axis above indicates the number of nonzero coefficients at the current 𝜆λ, which is the effective degrees of freedom (df) for the lasso.

plot(fit, label = TRUE)
print(fit)

cv.glmmod <- cv.glmnet(x, y=beta_sig_SPHS_ALL_tv3$T2D_status, family = "binomial")
plot(cv.glmmod)

coef_1selambda <- coef(cv.glmmod, s = "lambda.1se")
df_coef_1selambda <- (coef_1selambda)[,1]
write.csv(df_coef_1selambda, file="coef_1selambda.csv")

##End of lasso regression
