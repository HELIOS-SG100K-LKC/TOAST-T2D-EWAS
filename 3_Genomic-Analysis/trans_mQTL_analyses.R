library(readxl)
library(dplyr)
library(igraph)
library(RColorBrewer)
library(circlize)
library(networkD3)
library(htmlwidgets)
library(webshot)
library(gprofiler2)
library(ieugwasr)
library(ggnewscale)
library(ggrepel)
library(CMplot)
library(dendextend)



sent_count_mlma <- read.table("SNP_trans_cpg.txt",header = T)
sent_count_mlma <- sent_count_mlma %>% group_by(RSID) %>%summarise(count=length(unique(CpG)))
names(sent_count_mlma) <- c("SNP","Sent_Count")


background_count <- read.table("Sentinel_background_count.txt",
                               col.names = c("SNP","Background_Count"))

All_count <- merge.data.frame(sent_count_mlma,background_count,by="SNP")

All_count$Background_Count <- ifelse(All_count$Background_Count==0,
                                     All_count$Background_Count + All_count$Sent_Count,
                                     All_count$Background_Count)



All_count$FE <- (All_count$Sent_Count/323)/(All_count$Background_Count/837722)
All_count$P_hyper <- phyper(All_count$Sent_Count-1,All_count$Background_Count,837722 - All_count$Background_Count,323,lower.tail = F)
All_count$P_adj <- p.adjust(All_count$P_hyper,method = "bonferroni")

Signif <- All_count %>% filter(P_adj<=0.05)

write.table(All_count,file = "~/Desktop/Toast_sentinel_trans_mqtl_enrichment.txt",row.names = F,quote = F,sep = "\t")



Coloc_Signif <- read.delim("~/Desktop/TOAST-transmQTL/TransmQTL_Coloc_signif_raw_table.txt",header = T)

snp_gene_assoc <- Coloc_Signif[,c(3,1)]
write.table(snp_gene_assoc,file = "~/Desktop/coloc_snp_gene_pair.txt",row.names = F,quote = F,sep = " ")

unique_genes <- unique(Coloc_Signif$name)
unique_cpgs <- unique(Coloc_Signif$CpG)
unique_snps <- unique(Coloc_Signif$SNP)

genes <- unique(Coloc_Signif$Gene)

#### All eqtl enrichment analysis ###########
enrichment_all_genes <- gost(genes,organism = "hsapiens",significant = T,
                            sources = c("GO:BP","KEGG","REAC","WP"))
enrichment_all_genes <- enrichment_all_genes$result


snp_gene <- unique(paste0("chr",Coloc_Signif$SNP_chr_hg38,":",Coloc_Signif$SNP_pos_hg38,"_",Coloc_Signif$Gene))



### HELIOS-Replicated Significant-ciseQTLs

helios_snp_gene_assoc <- read.table("~/Desktop/TOAST-transmQTL/snp_gene_helios_assoc.txt",header = T)
helios_snp_gene_assoc <- unique(helios_snp_gene_assoc)

tmp <- helios_snp_gene_assoc  %>% filter(p.value<=0.05)
length(unique(tmp$Gene))

tmp <- helios_snp_gene_assoc %>% filter(!(is.na(p.value)))
length(unique(tmp$Gene))


eqtlgen_signif <- a1_signif %>% filter(RSID %in% Coloc_Signif$RSID & Gene %in% Coloc_Signif$Gene)
eqtlgen_signif <- eqtlgen_signif[,c(3,10,11,14)]

eqtlgen_helios_signif <- merge.data.frame(eqtlgen_signif,tmp,by.x = c("RSID","Gene"),all.x = T)
eqtlgen_helios_signif$beta <- NULL
names(eqtlgen_helios_signif) <- c("SNP","Gene","GeneSymbol","eqtlgen_pvalue","helios_pvalue")

write.table(eqtlgen_helios_signif,file="eqtl_comparison_eqtlgen_helios.txt",row.names = F,quote = F,sep = " ")



###

library(dplyr)
library(tidyr)

df <- unique(Coloc_Signif[,c(2,13)])
snp_cpg_list <- with(df, split(CpG, SNP))


names_list <- names(snp_cpg_list)
intersection_matrix <- matrix(nrow = length(names_list), ncol = length(names_list), dimnames = list(names_list, names_list))

for (i in names_list) {
  for (j in names_list) {
    intersection_matrix[i, j] <- length(intersect(snp_cpg_list[[i]], snp_cpg_list[[j]]))
  }
}

count_matrix <- matrix(nrow = length(names_list), ncol = length(names_list), dimnames = list(names_list, names_list))

for (i in names_list) {
  for (j in names_list) {
    count_matrix[i, j] <- (length(unique(c(snp_cpg_list[[i]],snp_cpg_list[[j]]))))
  }
}

distance_matrix <- (1 - intersection_matrix/count_matrix )

# Perform hierarchical clustering
hc <- hclust(as.dist(distance_matrix),method = "average")
png("~/Desktop/clustering_plot_n33.png",height = 20,width = 20,units = "cm",res = 360)
plot(hc, main="Hierarchical Clustering Dendrogram", xlab="", sub="", cex=.9)
dev.off()

# Function to calculate average silhouette width for different k
library(cluster)

# Initialize a vector to store average silhouette widths
avg_silhouette_widths <- numeric()

# Loop over the range of possible k values
for (k in 2:15) {
  # Cut the tree to get k clusters
  cluster_labels <- cutree(hc, k)
  
  # Compute the silhouette widths
  sil <- silhouette(cluster_labels, dist(distance_matrix))
  
  # Compute the average silhouette width and store it
  avg_silhouette_widths <- c(avg_silhouette_widths, mean(sil[, 3]))
}

# Plot average silhouette widths
png("~/Desktop/Silhoutte.png",height = 12,width = 12,units = "cm",res = 480)
plot(2:15, avg_silhouette_widths, type = "b", xlab = "Number of Clusters", ylab = "Average Silhouette Width")
dev.off()

groups <- cutree(hc, k=7)

dend <- as.dendrogram(hc)
dend <- dend %>% color_branches(dend, k = 7) %>% color_labels(k=7)
dend <- hang.dendrogram(dend,hang = 0.35)


# Plot the dendrogram with colored branches
png("~/Desktop/clustering_plot_k7.png", height = 25, width = 40, units = "cm", res = 360)
par(mar = c(15,5,2,2))
plot(dend, main = "Hierarchical Clustering Dendrogram (K=7)", xlab = "", sub = "", cex = .75)
rect.dendrogram(dend,k=7,border = "black",lty=2,lwd=0.75)
dev.off()




# Add the group information back to the original data
df$Group <- groups[as.character(df$SNP)]
df %>% group_by(Group) %>% summarise(Count = length(unique(SNP)))
groups
names(df)[1] <- "SNP"

snp_groups <- unique(df[,c(1,3)])
write.table(snp_groups,file = "~/Desktop/TOAST-transmQTL/SNP_groups.txt",row.names = F,quote = F,sep = " ")


## sankey plots by clusters

group_1 <- df %>% dplyr::filter(Group==1)
group_2 <- df %>% dplyr::filter(Group==2)
group_3 <- df %>% dplyr::filter(Group==3)
group_4 <- df %>% dplyr::filter(Group==4)
group_5 <- df %>% dplyr::filter(Group==5)
group_6 <- df %>% dplyr::filter(Group==6)
group_7 <- df %>% dplyr::filter(Group==7)


############




## sankey plot code (sample code for 1 cluster)


links_tmp <- unique(Coloc_Signif[,c(5,17)])

links_snp_gene <- Coloc_Signif %>% dplyr::select(SNP,name) %>% group_by(SNP,name) %>% 
  mutate(value = n()) %>% slice_head()

names(links_snp_gene) <- c("source", "target", "value")

# Adjusting links for Gene to CpG
links_gene_cpg <- unique(Coloc_Signif[,c(11,13,16,25,26)])
links_gene_cpg$value <- ifelse(links_gene_cpg$Coloc_H4 >= 0.9, links_gene_cpg$Coloc_H4, links_gene_cpg$Coloc_H3)


links_gene_cpg$CpG_Annot <- paste0(links_gene_cpg$CpG," (",links_gene_cpg$GeneSymbol_nearest_to_CpG_bumphunter_Gencode,")")
links_gene_cpg <- links_gene_cpg[,c(1,7,6)]
names(links_gene_cpg) <- c("source", "target", "value")

# Combine links
links <- as.data.frame(rbind(links_snp_gene, links_gene_cpg))

# Creating a unique list of names for nodes
unique_snps <- unique(links_snp_gene$source)
unique_genes <- unique(c(links_snp_gene$target, links_gene_cpg$source))
unique_cpgs <- unique(links_gene_cpg$target)

# Create node list
node_names <- unique(c(unique_snps, unique_genes, unique_cpgs))
nodes <- data.frame(name = node_names)

######



# Filtering and preparing data as per your existing code
links_snp_gene_1 <- links_snp_gene %>% dplyr::filter(source %in% group_1$SNP)
links_gene_cpg_1 <- links_gene_cpg %>% dplyr::filter(source %in% links_snp_gene_1$target)
links_1 <- as.data.frame(rbind(links_snp_gene_1, links_gene_cpg_1))

# Generate a list of unique SNPs, Genes, and CpGs
unique_snps_1 <- unique(links_snp_gene_1$source)
unique_genes_1 <- unique(c(links_snp_gene_1$target, links_gene_cpg_1$source))
unique_cpgs_1 <- unique(links_gene_cpg_1$target)

# Create a color palette for the unique SNPs
# You can use RColorBrewer or manually specify colors if you have more SNPs
colors <- RColorBrewer::brewer.pal(min(8, length(unique_snps_1)), "Pastel1")
if (length(unique_snps_1) > 8) {colors <- colorRampPalette(colors)(length(unique_snps_1))}


# Map SNPs to colors
color_map <- setNames(colors, unique_snps_1)

# Assign colors to links based on the SNP
links_1$color <- sapply(links_1$source, function(snp) color_map[snp])

gene_colors <- na.omit(links_1[,c(2,4)])
color_map_gene <- gene_colors$color
color_map_gene <- setNames(color_map_gene,gene_colors$target)

colormap_full <- c(color_map,color_map_gene)

links_1$color <- sapply(links_1$source, function(x) colormap_full[x])
node_names_1 <- unique(c(unique_snps_1, unique_genes_1, unique(links_gene_cpg_1$target)))
nodes_1 <- data.frame(name = node_names_1)
nodes_1$color <- sapply(nodes_1$name, function(x) colormap_full[x])
nodes_1$color <- ifelse(is.na(nodes_1$color),"grey",nodes_1$color)
links_1$source <- match(links_1$source, node_names_1) - 1  # zero-based index
links_1$target <- match(links_1$target, node_names_1) - 1  # zero-based index

# Generate the Sankey plot
sankey_1 <- sankeyNetwork(Links = links_1, Nodes = nodes_1, Source = "source", Target = "target", 
                          Value = "value", NodeID = "name", fontSize = 12.5, nodePadding = 10, 
                          nodeWidth = 7.5, height = 1200, width = 1300, iterations = 50, 
                          sinksRight = FALSE, LinkGroup = "color",NodeGroup = "color")

sankey_1 <- htmlwidgets::onRender(
  sankey_1,
  'function(el) {
     d3.select(el).selectAll(".node text")
       .style("font-weight", "bold");
  }'
)

sankey_1 <- htmlwidgets::onRender(
  sankey_1,
  'function(el) {
     d3.select(el).selectAll(".link")
       .style("stroke-opacity", 0.4);  // Set to 0.5 for 50% transparency
  }'
)


# Print the sankey plot
print(sankey_1)

sankey_1$x$nodes <- sankey_1$x$nodes %>% mutate(is_snp_node = name %in% links_snp_gene$source)
sankey_1$x$nodes <- sankey_1$x$nodes %>% mutate(is_gene_node = name %in% links_gene_cpg$source)
sankey_1$x$nodes <- sankey_1$x$nodes %>% mutate(is_cpg_node = name %in% links_gene_cpg$target)
sankey_1 <- htmlwidgets::onRender(sankey_1,customJS)
htmlwidgets::saveWidget(sankey_1, "sankeyPlot_1_new.html", selfcontained = TRUE)
webshot("sankeyPlot_1_new.html", "sankeyPlot_1.png", delay = 3,zoom=3)  # Delay to ensure the JS renders


                        
####
## GSEA - gene clusters (sample code for cluster 1)

Pathways_gene_set1 <- gost(genes_set_1,organism = "hsapiens",ordered_query = F,significant = T,
                           sources = c("GO:BP","KEGG","REAC","WP"))
plot_gs1 <- gostplot(Pathways_gene_set1,capped = F,interactive = T)
Pathways_gene_set1 <- Pathways_gene_set1$result  
