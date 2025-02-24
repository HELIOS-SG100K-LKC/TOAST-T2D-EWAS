
require(minfi)
require(limma)

anno <- read.csv("data/MethylationEPIC_15073387_v-1-0.csv",as.is=TRUE, skip = 7)
anno=anno[,c('Infinium_Design_Type','Color_Channel', 'CHR', 'MAPINFO', 'Name')]
cas=anno[substr(anno$Name, 1,3)=='ch.' & !(anno$CHR %in% c('X','Y')),]
cgs=anno[substr(anno$Name, 1,2)=='cg'& !(anno$CHR %in% c('X','Y')),]
auto = c(cgs$Name, cas$Name)
auto=as.matrix(auto)

#detection p-value
thres=0.01   ## detPrelax


load('output/intensities.RData')
load('output/detectionPvalue.RData')
d=dp.all[rownames(TypeII.Green.All),colnames(TypeII.Green.All)]
TypeII.Green.All.d = ifelse(d<thres,TypeII.Green.All,NA)
TypeII.Red.All.d = ifelse(d<thres,TypeII.Red.All,NA)

d=dp.all[rownames(TypeI.Green.M.All),colnames(TypeI.Green.M.All)]
TypeI.Green.M.All.d = ifelse(d<thres,TypeI.Green.M.All,NA)
TypeI.Green.U.All.d = ifelse(d<thres,TypeI.Green.U.All,NA)
d=dp.all[rownames(TypeI.Red.M.All),colnames(TypeI.Red.M.All)]
TypeI.Red.M.All.d = ifelse(d<thres,TypeI.Red.M.All,NA)
TypeI.Red.U.All.d = ifelse(d<thres,TypeI.Red.U.All,NA)
rm(dp.all,d)

# autosomes ------------------------------------------------------------------
samples=colnames(TypeI.Red.M.All)
category=auto
markers=as.matrix(intersect(rownames(TypeII.Green.All.d), category))
TypeII.Green = TypeII.Green.All.d[markers,samples]
TypeII.Red = TypeII.Red.All.d[markers,samples]
markers=intersect(rownames(TypeI.Green.M.All.d), category)
TypeI.Green.M = TypeI.Green.M.All.d[markers,samples]
TypeI.Green.U = TypeI.Green.U.All.d[markers,samples]
markers=intersect(rownames(TypeI.Red.M.All.d), category)
TypeI.Red.M = TypeI.Red.M.All.d[markers,samples]
TypeI.Red.U = TypeI.Red.U.All.d[markers,samples]

#raw betas
TypeII.betas = TypeII.Green/(TypeII.Red+TypeII.Green+100)
TypeI.Green.betas = TypeI.Green.M/(TypeI.Green.M+TypeI.Green.U+100)
TypeI.Red.betas = TypeI.Red.M/(TypeI.Red.M+TypeI.Red.U+100)
beta = as.matrix(rbind(TypeII.betas,TypeI.Green.betas,TypeI.Red.betas))
sample.call=colSums(!is.na(beta))/nrow(beta)
marker.call=rowSums(!is.na(beta))/ncol(beta)
save(sample.call, marker.call, file='output/callRates.RData')
save(beta, file='output/beta_raw.RData')

#call-rate filtering
callrate.thres=0.95
samples=names(sample.call[sample.call>=callrate.thres])

## gender_swap samples filtering
failed_genders <- read.csv("output/gender_swap_list.txt", sep="\t", header=T)
samples2=samples[!(samples %in% failed_genders$Sentrix_ID)]


markers=as.matrix(intersect(rownames(TypeII.Green.All.d), category))
markers=intersect(markers, names(marker.call[marker.call>=callrate.thres]))
TypeII.Green = TypeII.Green.All.d[markers,samples2]
TypeII.Red = TypeII.Red.All.d[markers,samples2]

markers=intersect(rownames(TypeI.Green.M.All.d), category)
markers=intersect(markers, names(marker.call[marker.call>=callrate.thres]))
TypeI.Green.M = TypeI.Green.M.All.d[markers,samples2]
TypeI.Green.U = TypeI.Green.U.All.d[markers,samples2]

markers=intersect(rownames(TypeI.Red.M.All.d), category)
markers=intersect(markers, names(marker.call[marker.call>=callrate.thres]))
TypeI.Red.M = TypeI.Red.M.All.d[markers,samples2]
TypeI.Red.U = TypeI.Red.U.All.d[markers,samples2]

#QN
TypeII.Green=normalizeQuantiles(TypeII.Green)
TypeII.Red = normalizeQuantiles(TypeII.Red)
TypeI.Green.M = normalizeQuantiles(TypeI.Green.M)
TypeI.Green.U = normalizeQuantiles(TypeI.Green.U)
TypeI.Red.M = normalizeQuantiles(TypeI.Red.M)
TypeI.Red.U = normalizeQuantiles(TypeI.Red.U)
TypeII.betas = TypeII.Green/(TypeII.Red+TypeII.Green+100)
TypeI.Green.betas = TypeI.Green.M/(TypeI.Green.M+TypeI.Green.U+100)
TypeI.Red.betas = TypeI.Red.M/(TypeI.Red.M+TypeI.Red.U+100)
beta = as.matrix(rbind(TypeII.betas,TypeI.Green.betas,TypeI.Red.betas))

save(beta, file="output/beta_QN_rmGSwap.RData")


rm(TypeII.Green.All.d,TypeII.Red.All.d,TypeI.Green.M.All.d,TypeI.Green.U.All.d,TypeI.Red.M.All.d,TypeI.Red.U.All.d)
rm(TypeII.Green,TypeII.Red,TypeI.Green.M,TypeI.Green.U,TypeI.Red.M,TypeI.Red.U,TypeII.betas,TypeI.Green.betas,TypeI.Red.betas)

sessionInfo()
