require(minfi)
require(IlluminaHumanMethylationEPICmanifest)
require(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
require(S4Vectors)

namelist <- read.csv("filename_list.txt", header=T)
filenames <- namelist$SentrixID

RGset <- read.metharray(file.path(paste0("IDAT_DT/", filenames)), force=TRUE, verbose=TRUE)
RGset <- bgcorrect.illumina(RGset)  # Illumina background subtraction 

GMsetEx <- mapToGenome(RGset)

estSex <- getSex(GMsetEx)

save(estSex, file="output/estSex.RData")


sessionInfo()



