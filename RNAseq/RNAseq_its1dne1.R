library(data.table)
library(openxlsx)
library(ggplot2)
library(biomaRt)
library(ggrepel)
library(DESeq2)
library(tidyr)
library(plyr)

TAIR10_GENES_MODELS <- fread("../ressources/TAIR10_annotations/TAIR10_representative_gene_models.txt", header = F, col.names = "AGI")

# RNAseq sample table
sampleTable <- fread("SampleTable.txt")

# Load featureCounts table
featureCounts <- fread("01_FeatureCounts/RNAseq_TAIR10.gene_models.featureCounts.txt", 
                       colClasses=c(rep("character", 6), rep("integer", 21)))

# Extract only the major isoform
featureCounts <- featureCounts[Geneid %in% TAIR10_GENES_MODELS$AGI]
colnames(featureCounts) <- gsub("\\.sorted.bam$", "", colnames(featureCounts))
setnames(featureCounts, c("Geneid", sampleTable$id), c("mRNA", paste(sampleTable$condition, sampleTable$replicate, sep = "_")))

# Count as matrix for DESeqDataSetFromMatrix
countdata <- as.matrix(featureCounts[, c("Chr", "Start", "End", "Strand", "Length") := NULL], rownames = TRUE)

# Coldata for DESeqDataSetFromMatrix
coldata <- data.frame(row.names=colnames(countdata), 
                      Condition = as.factor(sampleTable$condition), 
                      Replicate = as.factor(sampleTable$replicate))

# Load everything
dds <- DESeqDataSetFromMatrix(countData = countdata, 
                              colData = coldata, 
                              design= ~ Condition + Replicate)

# Minimal filtering step to remove very low covered genes
dds <- dds[ rowSums(counts(dds)) > 10, ]

## DESeq analysis
dds <- DESeq(dds, betaPrior=TRUE)

# Save as RDS
saveRDS(dds, "/02_DESeq2/RNAseq_dds.rds")
dds <- readRDS("02_DESeq2/RNAseq_dds.rds")

# Plots dispersion
plotDispEsts(dds, main="Dispersion plot")

# PCA to visualize variability between libraries
# VST (Variance Stabilizing Transformation)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE) #, fitType = "local")
pcaData <- plotPCA(vsd, intgroup=c("Condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# create pdf PCA.pdf"
my_colors = c("red", "green", "blue", "black", "gray80", "purple", "darkblue")
ggplot(pcaData, aes(PC1, PC2, color=Condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() +
  theme(plot.title = element_text(hjust = 0.5, size=22)) +
  ggtitle("PCA analysis_its1dne1_22") +
  scale_color_manual(values = my_colors) 

ggplot(pcaData, aes(PC1, PC2, color=Condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() +
  ggtitle("PCA analysis_its1dne1_22") +
  theme(plot.title = element_text(hjust = 0.5, size=22)) +
  geom_label_repel(aes(label = name, size = 2), show.legend = FALSE) + theme_bw()

## Statistical analysis
# The statistical results are then extracted for **pairwise comparison**. Desired comparison is specified by using the `contrast` argument. The ` results` function of DEseq2 generates a results table with log2 fold changes, p-values and adjusted p-values.

## Result table for pairwise comparison. Comparison between all genotypes and WT (Col0)
res_dne1_2_Col      <- results(dds, contrast=c("Condition", "dne1_2", "WT"), alpha = 0.05)
res_dne1_2_dcp2_Col <- results(dds, contrast=c("Condition", "dne1_2_dcp2", "WT"), alpha = 0.05)
res_xrn4_Col        <- results(dds, contrast=c("Condition", "xrn4", "WT"), alpha = 0.05)
res_dne1_3_Col      <- results(dds, contrast=c("Condition", "dne1_3", "WT"), alpha = 0.05)
res_dne1_3_dcp2_Col <- results(dds, contrast=c("Condition", "dne1_3_dcp2", "WT"), alpha = 0.05)
res_dcp2_Col        <- results(dds, contrast=c("Condition", "dcp2", "WT"), alpha = 0.05)

# Conversion to data.table
res_dne1_2_Col.dt      <- as.data.table(res_dne1_2_Col, keep.rownames = "mRNA")
res_dne1_2_dcp2_Col.dt <- as.data.table(res_dne1_2_dcp2_Col, keep.rownames = "mRNA")
res_xrn4_Col.dt        <- as.data.table(res_xrn4_Col, keep.rownames = "mRNA")
res_dne1_3_Col.dt      <- as.data.table(res_dne1_3_Col, keep.rownames = "mRNA")
res_dne1_3_dcp2_Col.dt <- as.data.table(res_dne1_3_dcp2_Col, keep.rownames = "mRNA")
res_dcp2_Col.dt        <- as.data.table(res_dcp2_Col, keep.rownames = "mRNA")

my_list <- list("xrn4vsCol" = res_xrn4_Col.dt,
               "dcp2vsCol" = res_dcp2_Col.dt,
               "dne1_2vsCol" = res_dne1_2_Col.dt,
               "dne1_3vsCol" = res_dne1_3_Col.dt, 
               "dne1_3_dcp2vsCol" = res_dne1_3_dcp2_Col.dt,
               "dne1_2_dcp2vsCol" = res_dne1_2_dcp2_Col.dt)         

# Remove ".1" from isoform name to get AGI annotation
llply(my_list, function(dt){ dt[, AGI := strsplit(mRNA, '[.]')[[1]][1], by = mRNA] })

# Function to annotate
addBiomaRtAnnotation <- function(input_dt, annotation_col = "ID", biomart = "plants_mart", dataset = "athaliana_eg_gene", host = "https://plants.ensembl.org", 
                                 features = c("ensembl_peptide_id", "ensembl_transcript_id", "ensembl_gene_id", "external_gene_name")) {
  
  stopifnot(!missing(input_dt))
  filter <- match.arg(features)
  mart <- useMart(biomart, dataset, host = host)
  atr <- c(features, "description")
  geneInfo <- data.table(getBM(values = get(annotation_col, input_dt), mart = mart, attributes = atr, filters = filter))
  if (length(geneInfo$description) == 0) {
    stop(sprintf("rowIDs from count table are different than IDs from \"%s\" from ensembl\n
                 Choose another feature ID, or another dataset from biomaRt, or check that your IDs are in the format as the one from ensembl", filter)
    )
  }
  message("Input # of elements: ", nrow(input_dt))
  message("Match found: ", nrow(geneInfo))
  ouput_dt <- merge(input_dt, geneInfo, by.x = annotation_col, by.y = features, all.x = T)
  return(ouput_dt)
}

# Add annotation based on AGI 
my_list.annotated <- llply(my_list, addBiomaRtAnnotation, annotation_col = "AGI", features = "ensembl_gene_id")

all_results.dt <- rbindlist(my_list.annotated, idcol = "Versus")
fwrite(all_results.dt, "02_DESeq2/All_versus.DESeq2.results.csv", sep = ",")

# Load back later
all_results.dt <- fread("02_DESeq2/All_versus.DESeq2.results.csv")

# Save in excel format with one tab per versus
raw_counts <- fread("01_FeatureCounts/RNAseq_TAIR10.gene_models.featureCounts.txt")
all_lists <- append(my_list.annotated, list(featureCounts=raw_counts))
write.xlsx(all_lists, file = "02_DESeq2/All_versus.DESeq2.results.xlsx")


up_regulated = all_results.dt[padj <= 0.05 & log2FoldChange > 0, .(up_regulated = .N), by = Versus]
down_regulated = all_results.dt[padj <= 0.05 & log2FoldChange < 0, .(down_regulated = .N), by = Versus]

fwrite(up_regulated, "02_DESeq2/up_regulated.log2FC>0&padj<=5%.stats.results.csv")
fwrite(down_regulated, "02_DESeq2/down_regulated.log2FC<0&padj<=5%.stats.results.csv")

# Save only log2FC and padj
summary_results <- rbindlist(llply(my_list.annotated, function(dt) { 
  setnames(dt, "log2FoldChange", "log2FC", skip_absent = T)
  dt[, c("mRNA", "AGI", "log2FC", "padj")]
}), idcol = "Versus")

summary_results.dt <- setDT(pivot_wider(summary_results, names_from = "Versus", values_from = c("log2FC", "padj"), names_sep = "."))
fwrite(summary_results.dt, "02_DESeq2/All_versus_log2FC_padj_only.results.csv", sep = ",")

# Add raw counts
res_dne1_2_Col.dt      <- merge(res_dne1_2_Col.dt, featureCounts[, c("mRNA", paste0("WT_rep", 1:3), paste0("dne1_2_rep", 1:3))], by = "mRNA")
res_dne1_2_dcp2_Col.dt <- merge(res_dne1_2_dcp2_Col.dt, featureCounts[, c("mRNA", paste0("WT_rep", 1:3), paste0("dne1_2_dcp2_rep", 1:3))], by = "mRNA")
res_xrn4_Col.dt        <- merge(res_xrn4_Col.dt, featureCounts[, c("mRNA", paste0("WT_rep", 1:3), paste0("xrn4_rep", 1:3))], by = "mRNA")
res_dne1_3_Col.dt      <- merge(res_dne1_3_Col.dt, featureCounts[, c("mRNA", paste0("WT_rep", 1:3), paste0("dne1_3_rep", 1:3))], by = "mRNA")
res_dne1_3_dcp2_Col.dt <- merge(res_dne1_3_dcp2_Col.dt, featureCounts[, c("mRNA", paste0("WT_rep", 1:3), paste0("dne1_3_dcp2_rep", 1:3))], by = "mRNA")
res_dcp2_Col.dt        <- merge(res_dcp2_Col.dt, featureCounts[, c("mRNA", paste0("WT_rep", 1:3), paste0("dcp2_rep", 1:3))], by = "mRNA")

my_list <-list("xrn4vsCol" = res_xrn4_Col.dt,
               "dcp2vsCol" = res_dcp2_Col.dt,
               "dne1_2vsCol" = res_dne1_2_Col.dt,
               "dne1_3vsCol" = res_dne1_3_Col.dt, 
               "dne1_3_dcp2vsCol" = res_dne1_3_dcp2_Col.dt,
               "dne1_2_dcp2vsCol" = res_dne1_2_dcp2_Col.dt)       

my_list.annotated <- llply(my_list, addBiomaRtAnnotation, annotation_col = "AGI", features = "ensembl_gene_id")

## SAVE RAW RESULTS IN CSV
for(i in 1:length(my_list.annotated)){
  fwrite(x = my_list.annotated[[i]], 
         file = file.path("02_DESeq2", paste0(names(my_list.annotated[i]), ".csv")), 
         sep = ","
  )
}

# MA-plot for paper
pdf("02_DESeq2/MA_plot_analyse_dne1_2_dcp2.pdf", 5, 10)
par(mfrow=c(2,1)) #here the three plots will be drawn next to each other
outdne1_2_Col <- plotMA(res_dne1_2_Col, main = "dne1_2vsCol", ylim=c(-3,3), alpha=0.05, xlab = "mean of normalized counts")
outdne1_2_dcp2_Col <- plotMA(res_dne1_2_dcp2_Col, main = "dne1_2_dcp2vsCol", ylim=c(-3,3), alpha=0.05, xlab = "mean of normalized counts")
outdne1_3_Col <- plotMA(res_dne1_3_Col, main = "dne1_3vsCol", ylim=c(-3,3), alpha=0.05, xlab = "mean of normalized counts")
outdne1_3_dcp2_Col <- plotMA(res_dne1_3_dcp2_Col, main = "dne1_3_dcp2vsCol", ylim=c(-3,3), alpha=0.05, xlab = "mean of normalized counts")
outxrn4_Col <- plotMA(res_xrn4_Col, main = "xrn4vsCol", ylim=c(-3,3), alpha=0.05, xlab = "mean of normalized counts")
outdcp2_Col <- plotMA(res_dcp2_Col, main = "dcp2vsCol", ylim=c(-3,3), alpha=0.05, xlab = "mean of normalized counts")
dev.off()

# Normalize count
table_norm <- counts(dds, normalized=TRUE)

up_regulated = all_results.dt[padj <= 0.05 & log2FoldChange > 0, .(up_regulated = .N), by = Versus]
down_regulated = all_results.dt[padj <= 0.05 & log2FoldChange < 0, .(down_regulated = .N), by = Versus]