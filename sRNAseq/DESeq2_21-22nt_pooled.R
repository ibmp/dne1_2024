library(rtracklayer)
library(data.table)
library(openxlsx)
library(pheatmap)
library(biomaRt)
library(ggrepel)
library(svglite)
library(ggplot2)
library(DESeq2)
library(scales)
library(tidyr)
library(dplyr)
library(knitr)
library(plyr)
library(cli)
library(DT)

TAIR10_GENES_MODELS <- fread("../ressources/TAIR10_annotations/TAIR10_representative_gene_models.txt", header = F, col.names = "AGI")

# Shortstack file after mapping with bowtie1 0 mismatch 18-26nt with unique mapping (-u)
shortstack_files <- list.files(path = "01_shortstack", "Results.txt", full.names = T, recursive = T)
names(shortstack_files) <- gsub("_shortstack", "", basename(dirname(shortstack_files)))

# Extract RowSum of specific columns from a shortstack Results.txt file
getRowSum <- function(shortstack_file, new_colname = "Count", columns_to_sum = NULL, new_extension = NULL) {
  outfile <- gsub(".txt", new_extension, shortstack_file)
  mydata <- fread(shortstack_file)
  cli_alert_info("Analysing {shortstack_file}")
  mydata[ ,(new_colname) := rowSums(.SD), .SDcols = columns_to_sum]
  mydata.subset <- mydata[, c("Name", new_colname), with = FALSE]
  fwrite(mydata.subset, file = outfile, col.names = FALSE, sep = "\t")
  cli_alert_success("New file written: {outfile}")
}

l_ply(shortstack_files, getRowSum, new_colname = "Sum21-22nt", columns_to_sum = c("21", "22"), new_extension = ".21-22nt.txt")

samples <- fread("SampleTable.txt") 
samples[, Condition := as.factor(Condition)]
samples[, Batch := as.factor(Batch)]
unique(samples$Condition)

# shortstack
dds <- DESeqDataSetFromHTSeqCount(directory = "01_shortstack", 
                                  sampleTable = samples, 
                                  design = ~ Condition + Batch)
counts(dds)

## Filter TAIR10_GENES_MODELS AGI only
dds <- dds[rownames(dds) %in% TAIR10_GENES_MODELS$AGI, ]

# Run DESeq2
dds <- DESeq(dds)

saveRDS(dds, "dds.21_22nt_counts.mRNA_only.shortstack.rds")

# read back once created
dds <- readRDS("dds.21_22nt_counts.mRNA_only.shortstack.rds")

# Plots dispersion
plotDispEsts(dds, main="Dispersion plot")

### PCA rapick check of the variability between the samples
## 1. Variance Stabilizing Transformation
vsd <- varianceStabilizingTransformation(dds, blind=FALSE) #, fitType = "local")
sampleDists <- dist(t(assay(vsd)))
plotPCA(vsd, intgroup=c("Condition")) + theme_bw()

# customize PCA
pcaData <- plotPCA(vsd, intgroup=c("Condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Condition)) +
  geom_point(size=3) + geom_label_repel(aes(label = Condition, size = 2), show.legend = FALSE) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme_bw(base_size = 16) + guides(fill = "none", color = "none")

ggsave("results/sRNAseq/PCA_vsd_method.mRNA_only.21_22nt_counts.pdf", width = 12, height = 9)
ggsave("results/sRNAseq/PCA_vsd_method.mRNA_only.21-22nt_counts.svg", width = 12, height = 9)

boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
raw_table <- as.data.table(counts(dds), keep.rownames = T)

dds@colData

# DEseq contrast
res0 <- results(dds, contrast = c("Condition", "dcp2", "WT"), alpha = 0.05)
res1 <- results(dds, contrast = c("Condition", "dne1_3", "WT" ), alpha = 0.05)
res2 <- results(dds, contrast = c("Condition", "dne1_2", "WT"), alpha = 0.05)
res3 <- results(dds, contrast = c("Condition", "dne1_3_dcp2", "WT"), alpha = 0.05)
res4 <- results(dds, contrast = c("Condition", "dne1_2_dcp2", "WT"), alpha = 0.05)
res5 <- results(dds, contrast = c("Condition", "xrn4", "WT"), alpha = 0.05)
res6 <- results(dds, contrast = c("Condition", "dne1_3_dcp2", "dcp2"), alpha = 0.05)
res7 <- results(dds, contrast = c("Condition", "dne1_2_dcp2", "dcp2"), alpha = 0.05)

# Conversion to data.table for speed and convenience
res0.dt <- as.data.table(as.data.frame(res0), keep.rownames = T)
res1.dt <- as.data.table(as.data.frame(res1), keep.rownames = T)
res2.dt <- as.data.table(as.data.frame(res2), keep.rownames = T)
res3.dt <- as.data.table(as.data.frame(res3), keep.rownames = T)
res4.dt <- as.data.table(as.data.frame(res4), keep.rownames = T)
res5.dt <- as.data.table(as.data.frame(res5), keep.rownames = T)
res6.dt <- as.data.table(as.data.frame(res6), keep.rownames = T)
res7.dt <- as.data.table(as.data.frame(res7), keep.rownames = T)

## ADD RAW COUNT
res0.dt <- merge(res0.dt, raw_table[, c("rn", "Col_0_R1", "Col_0_R2", "Col_0_R3", "its1_R1", "its1_R2", "its1_R3")], by="rn", sort=FALSE)
res1.dt <- merge(res1.dt, raw_table[, c("rn", "Col_0_R1", "Col_0_R2", "Col_0_R3", "dne1_3_11_R1", "dne1_3_11_R2", "dne1_3_11_R3")], by="rn", sort=FALSE)
res2.dt <- merge(res2.dt, raw_table[, c("rn", "Col_0_R1", "Col_0_R2", "Col_0_R3", "dne1_3_80_R1", "dne1_3_80_R2", "dne1_3_80_R3")], by="rn", sort=FALSE)
res3.dt <- merge(res3.dt, raw_table[, c("rn", "Col_0_R1", "Col_0_R2", "Col_0_R3", "its1xdne1_3_11_R1", "its1xdne1_3_11_R2", "its1xdne1_3_11_R3")], by="rn", sort=FALSE)
res4.dt <- merge(res4.dt, raw_table[, c("rn", "Col_0_R1", "Col_0_R2", "Col_0_R3", "its1xdne1_3_80_R1", "its1xdne1_3_80_R2", "its1xdne1_3_80_R3")], by="rn", sort=FALSE)
res5.dt <- merge(res5.dt, raw_table[, c("rn", "Col_0_R1", "Col_0_R2", "Col_0_R3", "xrn4_3_R1", "xrn4_3_R2", "xrn4_3_R3")], by="rn", sort=FALSE)
res6.dt <- merge(res6.dt, raw_table[, c("rn", "its1_R1", "its1_R2", "its1_R3", "its1xdne1_3_11_R1", "its1xdne1_3_11_R2", "its1xdne1_3_11_R3")], by="rn", sort=FALSE)
res7.dt <- merge(res7.dt, raw_table[, c("rn", "its1_R1", "its1_R2", "its1_R3", "its1xdne1_3_80_R1", "its1xdne1_3_80_R2", "its1xdne1_3_80_R3")], by="rn", sort=FALSE)

my_list <- list("dcp2_vs_col0" = res0.dt, 
                "dne1_3_vs_col0" = res1.dt, 
                "dne1_2_vs_col0" = res2.dt, 
                "dne1_3_dcp2_vs_col" = res3.dt, 
                "dne1_2_dcp2_vs_col" = res4.dt,
                "xrn4_vs_col0" = res5.dt,
                "dne1_3_dcp2_vs_dcp2" = res6.dt,
                "dne1_2_dcp2_vs_dcp2" = res7.dt)

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

# Remove ".1" from isoform name to get AGI annotation
llply(my_list, function(dt){ dt[, AGI := strsplit(rn, '[.]')[[1]][1], by = rn] })

# Query biomart with the gene AGI
my_list.annotated <- llply(my_list, addBiomaRtAnnotation, annotation_col = "AGI", features = "ensembl_gene_id")
l_ply(my_list.annotated, setnames, old = "rn", new = "mRNA")

all_lists <- append(my_list.annotated, list(featureCounts=raw_table))
write.xlsx(all_lists, file = "analyses/sRNAseq/02_DESeq2/All_versus.DESeq2.results.xlsx")

## SAVE RAW RESULTS IN CSV
for(i in 1:length(my_list.annotated)){
  fwrite(x = my_list.annotated[[i]], 
         file = file.path("analyses/sRNAseq/02_DESeq2", paste0(names(my_list.annotated[i]), ".csv")), 
         sep = ","
         )
}

summary_results <- rbindlist(llply(my_list, function(dt) { 
  setnames(dt, "log2FoldChange", "log2FC", skip_absent = T)
  dt[, c("rn", "AGI", "log2FC", "padj")]
  }), idcol = "Versus")

summary_results.dt <- setDT(pivot_wider(summary_results, names_from = "Versus", values_from = c("log2FC", "padj"), names_sep = ".", values_fill = NA))
setnames(summary_results.dt, "rn", "mRNA")
fwrite(summary_results.dt, "analyses/sRNAseq/02_DESeq2/All_versus.log2FC&padj_only.results.csv", sep = ",")

sRNA_files <- list.files("analyses/sRNAseq/02_DESeq2", "*.csv", full.names = T)[-1]
names(sRNA_files) <- gsub(".csv", "", basename(sRNA_files))
sRNA.dt <- rbindlist(lapply(sRNA_files, fread, select = 1:8), idcol = "Sample")

sRNA.dt[log2FoldChange >= 0 & padj < 0.05][, .N, by = Sample]


# DicerCall21 sites -------------------------------------------------------
shortstack_merged <- list.files("01_shortstack_merged", "Results.txt$", full.names = T, recursive = T) 
names(shortstack_merged) <- gsub("_shortstack", "", basename(dirname(dc21)))

shortstack_merged.dt <- rbindlist(lapply(shortstack_merged, fread), idcol = "Sample")

dc21.dt <- shortstack_merged.dt[DicerCall == 21 & grepl("AT[12345]G", Name)]
dc21.dt[, gene_id := tstrsplit(Name, ".", fixed = TRUE, keep = 1)]

dc21.dt[, .N, by = Sample]

fwrite(dc21.dt, "01_shortstack_merged/ShortStack_all_merged_DicerCall21.csv", sep = ",")

dc21_selected_name <- dc21.dt[, .N, by = Name][N >= 3][, c("Name")]
fwrite(dc21_selected_name, "01_shortstack_merged/DicerCall21.min3.mRNA.txt")

# Gene_id
dc21_selected <- dc21.dt[, .N, by = gene_id][N >= 3]$gene_id
