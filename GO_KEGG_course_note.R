## use base R to draw heatmap
# Generate a matrix of 20 columns and 20 rows with normally distributed random values.
set.seed(1234)                                                     # Set seed for reproducibility
data <- matrix(rnorm(400, 0, 10), nrow = 20, ncol = 20)         # Create example data
colnames(data) <- paste0("col", 1:20)                             # Column names
rownames(data) <- paste0("row", 1:20)                             # Row names
# draw heatmap
heatmap(data)                                                     # Use ?heatmap() to check function arguments
heatmap(data, Rowv = NA, Colv = NA)                               # Remove dendogram
my_colors <- colorRampPalette(c("cyan", "deeppink3"))       # Manual color range
heatmap(data, col = my_colors(100))                               # Heatmap with manual colors

## use ggplot2 to draw heatmap
library(ggplot2)                                                   # load package. If need to install, run { install.packages("reshape")  }
library(reshape)
data_melt <- melt(data)                                           # Reorder data
head(data_melt)                                                   # First six rows of data
ggp <- ggplot(data_melt, aes(X1, X2)) +                           # Create heatmap with ggplot2
geom_tile(aes(fill = value))
print(ggp)                                                               # Print heatmap
ggp + scale_fill_gradient(low = "red", high = "black")          # Manual colors of heatmap

##=======================================================

# GO enrichment analysis
## Input RNA-Seq data set for one stats comparison between two sample groups (T2 and C1)
ENTREZID        REFSEQ  Symbol  log2FC_2_T2_C1  padj_2_T2_C1    FPKM_C1_R1      FPKM_C1_R2      FPKM_C1_R3      FPKM_T2_R1      FPKM_T2_R2      FPKM_T2_R3
79854   NR_024321       LINC00115       2.56117687410538        0.000980061044085233    0.293626        0.031558        0.313155        0.865013        1.24959 1.031439
148398  NM_152486       SAMD11  1.45018074678287        0.00684867349315434     3.520096        2.078321        4.672326        4.209266        5.661716        12.274103
9636    NM_005101       ISG15   1.23928784853599        0.0301917585476808      9.257691        8.80193 9.06228 6.830353        11.724176       35.494343

# Load packages
library(org.Hs.eg.db)
library(clusterProfiler)

# read in files to R
df1 <- read.delim("Decoy_human_RNASeq_DEG_grp1.tab", header=T)      # Read in RNA-Seq data set
df2 <- read.delim("Decoy_human_RNASeq_DEG_grp2.tab", header=T)

# GO enrichment analysis
gene1 <- as.character(df1$ENTREZID)          # extract gene entrez ID 

# groupGO classify genes by GO category distribution.
ggo1 <- groupGO(gene = gene1, OrgDb=org.Hs.eg.db, ont="CC", level=3, readable = TRUE)      # ont: One of "BP", "MF", and "CC"

# enrichGO performs Over-representation test(Boyle et al. 2004).
ego1 <- enrichGO(gene = gene1, OrgDb= org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)  # ont: one of three or "ALL".
write.table(ego1, file="grp1_go_enrich.tab", col.names=T, row.names=F, sep="\t", quote=F)

# gseGO performs Gene Set Enrichment Analysis (GSEA)(Subramanian et al. 2005). 
FC <- df1$log2FC_1_T1_C1                      # select gene fold change data as input
names(FC) <- gene1                            # add gene entrez ID as name for each fold change data point. 
FC <- rev(sort(FC))                           # input needs to be a degreasingly sorted vector
gsea1 <- gseGO(geneList = FC, OrgDb = org.Hs.eg.db, ont = "ALL", nPerm = 1000, minGSSize = 100, maxGSSize = 500, pvalueCutoff = 0.1, verbose = FALSE)
write.table(gsea1, file="grp1_go_GSEA.tab", col.names=T, row.names=F, sep="\t", quote=F)

#plotting GO enrichment results
library(ggplot2)
barplot(ego1, showCategory=10)
dotplot(ego1, showCategory=10) + ggtitle("dotplot for GO enrichment")
cnetplot(ego1, foldChange=FC )
heatplot(ego1, foldChange=FC)
library(enrichplot)
upsetplot(ego1)
##======================================================================

## KEGG enrichment
search_kegg_organism('Homo sapiens', by='scientific_name')
kegg1 <- enrichKEGG(gene = gene1, organism = 'hsa', pAdjustMethod = "BH", pvalueCutoff = 0.2)
library("pathview")
> head(kegg1$ID)
[1] "hsa04610" "hsa05150" "hsa00980" "hsa04640" "hsa00140" "hsa05204"
pathview(gene.data =FC, pathway.id = "hsa00140", species = "hsa", limit = list(gene=5, cpd=1))

# view all kegg pathways
for ( kegg_id in kegg$ID) { pathview(gene.data =FC, pathway.id = kegg_id, species = "hsa", limit = list(gene=5, cpd=1)) }



