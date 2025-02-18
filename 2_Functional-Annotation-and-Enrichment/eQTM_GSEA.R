## Gene Set Enrichment Analysis ###

GSEA was performed using the gprofiler online tool. 

### cis-eQTMs

cis_results_all <- read.csv("cis-eQTM-Significant-GSEA.csv",header = T)
cis_results_all <- cis_results_all[,c(1,2,6,7,8,9,10)]
cis_results_all <- cis_results_all %>% mutate(
  FE = ((intersection_size/query_size)/(term_size/effective_domain_size)),
  P_Value = phyper(intersection_size-1,term_size,effective_domain_size-term_size,query_size,lower=F)
)

cis_results_all <- cis_results_all %>%
  group_by(source) %>%
  mutate(Adj_P = p.adjust(P_Value, method = "fdr")) %>% # You can choose the method you prefer
  ungroup()

cis_results_all <- cis_results_all[,c(1,2,3,5,8,9,10,7)]


#######################

Trans_eQTM

trans_results_all <- read.csv("trans-eQTM-GSEA-Results-new.csv",header = T)
trans_results_all <- trans_results_all[,c(1,2,6,7,8,9,10)]
trans_results_all <- trans_results_all %>% mutate(
  FE = ((intersection_size/query_size)/(term_size/effective_domain_size)),
  P_Value = phyper(intersection_size-1,term_size,effective_domain_size-term_size,query_size,lower=F)
)

trans_results_all <- trans_results_all %>%
  group_by(source) %>%
  mutate(Adj_P = p.adjust(P_Value, method = "fdr")) %>% # You can choose the method you prefer
  ungroup()

trans_results_all <- trans_results_all[,c(1,2,3,5,8,9,10,7)]





