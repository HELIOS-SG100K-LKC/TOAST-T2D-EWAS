  rm(list=ls())

  cat('\f')

  graphics.off()

  args=(commandArgs(trailingOnly = TRUE))

  job_number = as.integer( args[1] )
  total_jobs = as.integer( args[2] )

  cwd = paste0(as.character( args[4] ), "/")

  setwd(cwd)

loadRData <- function(fileName){
   #loads an RData file, and returns it
   load(fileName)
   get(ls()[ls() != "fileName"])
}


  MainEffect_Results = loadRData('sg10k_MainEffect_Results.RData') # Full_cis_DF

  chromosome_number = as.integer( args[3] )

  Chr_Indices = c()

  Chr_Indices = which( MainEffect_Results$CpG_chr == chromosome_number )

  Chr_CpGs_trans_meQTLs_DF = MainEffect_Results[Chr_Indices,]

  Unique_CpG_List = unique( as.character( Chr_CpGs_trans_meQTLs_DF$CpG ) )

  total_unique_CpGs = 0

  total_unique_CpGs = length( Unique_CpG_List )

  subset_size = ceiling(total_unique_CpGs/total_jobs)

  start_position = 0;         start_position = ( ( job_number - 1 ) * subset_size ) + 1

  stop_position = 0;          stop_position = ( ( job_number ) * subset_size )


  if ( start_position > total_unique_CpGs ){
      quit(save = "no", status = 0, runLast = F)
  }


  max_position = total_unique_CpGs

  if( stop_position > max_position ) {

      stop_position = max_position

  }


    
  # other MLMA specific directories and files
  # -----------------------------------------

  gcta = '/home/project/12000713/darwin/opt/gcta/gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static'

  grm_full = '/home/project/12002085/sg10k_imputed/v3/grm/maf0001_Rsq030/sg10k_maf0.001_hwe1e-03_Rsq0.30'

  cov_file = '/data/projects/12000713/darwin/meQTL/mainEffect_model_sg10k/data/covariates/sg10k_no_gPCs.cov.txt'

  qcov_file = '/data/projects/12000713/darwin/meQTL/mainEffect_model_sg10k/data/covariates/sg10k_no_gPCs.qcov.txt'

  phenotype_dir = paste0( "/data/projects/12000713/darwin/meQTL/mainEffect_model_sg10k/data/generate_phenotype/chr", as.character(chromosome_number),  "/Phenotype/")


  # performing MLMA LOCO analysis for the chosen subset of CpGs 
  # ------------------------------------------------------------


  for( cpg_number in c(start_position:stop_position) ){

       cpg_marker = '';           cpg_marker = Unique_CpG_List[cpg_number]

       # phenotype file ( methylation data ) 
       # -----------------------------------

       phenotype_file = '';       phenotype_file = paste( phenotype_dir, cpg_marker, "_phenotype.phen", sep='' )


       # RSIDs segregated by individual chromosomes
       # -------------------------------------------  

       CpG_Indices = c();         CpG_Indices = which( Chr_CpGs_trans_meQTLs_DF$CpG == cpg_marker )

       Current_CpG_DF = data.frame()

       Current_CpG_DF = Chr_CpGs_trans_meQTLs_DF[CpG_Indices,]

       # identify whether SNPs associated with the current CpG are spread over multiple chromosomes

       Unique_SNP_Chromosomes = c()

       Unique_SNP_Chromosomes = unique( as.integer( Current_CpG_DF$SNP_chr ) )

       
       # loop through individual SNP chromosomes and perform MLMA LOCO 

       for( SNP_chromosome in Unique_SNP_Chromosomes ){

            # genotype file for the current chromosome

            Chr_bfile = ''

            Chr_bfile = paste( "/home/project/12002085/sg10k_imputed/v3/qcd/maf0001_Rsq030/chr", SNP_chromosome, "/chr", SNP_chromosome, "_sg10k_maf0.001_hwe1e-3_Rsq0.30.final", sep='' )
          

            # grm for the current chromosome

            grm_chr = '' 

            grm_chr = paste( "/home/project/12002085/sg10k_imputed/v3/grm/maf0001_Rsq030/byChr/chr", SNP_chromosome, ".sg10k_maf0.001_hwe1e-03_Rsq0.30", sep='' )
            
 
            # results file and log file names

            outDir = paste0(cwd, "mlma_output/chr", as.character(chromosome_number), "/")
	    system(paste0("mkdir -p ", outDir))


            out_file = '';             out_file = paste0( outDir, cpg_marker, "_Chr", SNP_chromosome )

            log_file = '';             log_file = paste0( out_file, ".log" )
                         

            # MLMA LOCO

            MLMA_command = ''

            MLMA_command = paste( gcta, " --bfile ", Chr_bfile, " --grm ", grm_full, " --mlma-subtract-grm ", grm_chr, " --pheno ", phenotype_file, " --qcovar ", qcov_file, " --covar ", cov_file, " --mlma --reml-maxit 1000 --thread-num 4", " --out ", out_file, " > ", log_file, sep='' )          

            system( MLMA_command )
                                

       } 
          

  } 

 
