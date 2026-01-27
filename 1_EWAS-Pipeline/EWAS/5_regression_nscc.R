# setwd("/home/project/12000713/zhouli_data/20190225_EPIC_all/")

# regression model
lfla=as.formula('phe$t2d ~ beta[i, ] + phe$Age + phe$Sex + phe$CD8T + 
                phe$CD4T + phe$NK + phe$Bcell + phe$Mono + phe$Gran + phe$PC1_cp + phe$PC2_cp + 
                phe$PC3_cp + phe$PC4_cp + phe$PC5_cp + phe$PC6_cp + phe$PC7_cp + 
                phe$PC8_cp + phe$PC9_cp + phe$PC10_cp + phe$PC11_cp + phe$PC12_cp + 
                phe$PC13_cp + phe$PC14_cp + phe$PC15_cp + phe$PC16_cp + phe$PC17_cp + 
                phe$PC18_cp + phe$PC19_cp + phe$PC20_cp + phe$PC21_cp + phe$PC22_cp + 
                phe$PC23_cp + phe$PC24_cp + phe$PC25_cp + phe$PC26_cp + phe$PC27_cp + 
                phe$PC28_cp + phe$PC29_cp + phe$PC30_cp + phe$PC1 + phe$PC2 + phe$PC3 +
                phe$PC4 + phe$PC5')
logistic=TRUE   # linear or logistic regression

# regression
load("output/pheno.RData")
samples=rownames(phe)
load('output/beta_QN_rmGSwap.RData')
samples = intersect(samples, colnames(beta))
beta=beta[,as.character(samples)]
phe <- phe[rownames(phe) %in% samples, ]
nvar = nrow(beta)
nvar = nrow(beta)
res=matrix(ncol=4, nrow=nvar)
if(logistic){
    colnames(res) =c('Estimate','Std. Error', 'z value','Pr(>|z|)')  
    if(!is.factor(phe$t2d)){
        phe$t2d = as.factor(as.character(phe$t2d))
    }
    for(i in 1:nvar) {
        tryCatch({fit = summary(glm(lfla, family = binomial))}, error = function(error) {return(NA)})
        if(!exists("fit")){
            res[i,] = rep(NA,4)
        }else{
            res[i,] = tryCatch({fit$coefficients[2,]},error = function(error) {return(rep(NA,4))})
            rm(fit)
        }
    }
}else{
    colnames(res) =c('Estimate','Std. Error', 't value','Pr(>|t|)')  
    if(!is.numeric(phe$t2d)){
        phe$t2d = as.numeric(as.character(phe$t2d))
    }
    for(i in 1:nvar) {
        tryCatch({fit= summary(lm(lfla))}, error = function(error) {return(NA)})
        if(!exists("fit")){
            res[i,] = rep(NA,4)
        }else{
            res[i,] = tryCatch({fit$coefficients[2,]},error = function(error) {return(rep(NA,4))})
            rm(fit)
        }
    }
}
rownames(res) <- rownames(beta)
save(res, file = 'output/results_t2d.RData')


