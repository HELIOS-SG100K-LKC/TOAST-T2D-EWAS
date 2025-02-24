
outDir="output/"

load("output/estSex.RData")

estSex <- as.data.frame(estSex)

estSex$Sentrix_ID <- row.names(estSex)
## read in the phenotype information from the table
pheno <- readxl::read_xlsx("data/samplesheet_gender_samples_ID_DT.xlsx", skip=7)
pheno$Sentrix_ID <- paste(pheno$Sentrix_ID, pheno$Sentrix_Position, sep="_")
estSex_phe <- merge(pheno[, c("Sample_Name", "Sample_Well", "Sample_Plate", "Supplied_Gender", "Sentrix_ID")], estSex, by="Sentrix_ID")

estSex_phe$Gender_predictedSex_Match <- ifelse(estSex_phe$Supplied_Gender == estSex_phe$predictedSex, "pass", "fail")

write.table(estSex_phe, file=paste0(outDir, "check_gender_swap.txt"),
            row.names=F, col.names = T, sep="\t", quote=F)

########## gender swapping ##########
failed_genders <- estSex_phe[estSex_phe$Gender_predictedSex_Match == "fail",]
failed_genders <- failed_genders[grepl("_", failed_genders$Sentrix_ID), ] 
write.table(failed_genders, file=paste0(outDir, "gender_swap_list.txt"),
            row.names=F, col.names = T, sep="\t", quote=F)

dim(failed_genders)  ## 22   9
dim(estSex_phe)   ## 1680   99
# 22/1680=1.3%
    
estSex_phe <- estSex_phe[grepl("_", estSex_phe$Sentrix_ID), ]

require(ggplot2)

pdf(paste0(outDir, "predictedSex_HF.pdf"), width=8, height=6)
ggplot(estSex_phe, aes(x=xMed, y=yMed)) +
    geom_point(aes(color=Gender_predictedSex_Match, shape=Supplied_Gender), size=1) +
    scale_color_manual(values=c("red", "blue"),
                       name="Predicted_Gender", 
                       breaks=c("fail", "pass"),
                       labels=c("unexpected", "original")) +
    scale_shape_manual(values=c(6,8),
                       name="Supplied_Gender",
                       breaks=c("F", "M"),
                       labels=c("F", "M")) +
    geom_text(data=estSex_phe[estSex_phe$Gender_predictedSex_Match != "pass", ],
              aes(x=xMed, y=yMed, label=Sentrix_ID), size=1) +
    theme_classic() +
    theme(text=element_text(size=15),
          legend.title =element_text(size=10) )
dev.off()

