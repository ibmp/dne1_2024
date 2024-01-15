library(DEXSeq)
library(ggplot2)
library(openxlsx)
library(data.table)
library(rtracklayer)
library(BSgenome.Athaliana.TAIR.TAIR9)

theme_set(theme_bw()) # change ggplot default theme
BSGENOME = getSeq(BSgenome.Athaliana.TAIR.TAIR9) # Genome as BSgenome object

# Select the samples to compare
SAMPLES <- c("xrn4_R1", "xrn4_R2", "xrn4_R3", "xrn4dne1_R1", "xrn4dne1_R2", "xrn4dne1_R3")

# Load GMUCT counts table with the annotations
gmuct_mRNA_counts <- fread("GMUCT_TAIR10_mRNA.min1RPM.csv")
gmuct_mRNA_counts <- gmuct_mRNA_counts[rpmMeans_xrn4 >= 1 | rpmMeans_xrn4dne1 >= 1]

count.mat <- as.matrix(gmuct_mRNA_counts[, .SD, .SDcols = c("ID", paste0(SAMPLES, ".count"))], rownames = "ID")
colnames(count.mat) <- SAMPLES

# Sample table
sample.table <- data.frame(
  row.names = SAMPLES,
  condition = as.factor(c(rep("xrn4", 3), rep("xrn4dne1", 3)))
)

# Transform to matrix for DEXSeq
dxd <- DEXSeqDataSet(countData = count.mat, 
                     sampleData = sample.table,
                     design = ~ sample + exon + condition:exon,
                     featureID = rownames(count.mat), # Position as "exon"
                     groupID = gmuct_mRNA_counts$gene_id,  # GeneID
                     featureRanges = GRanges(seqnames = gmuct_mRNA_counts$Chr, 
                                             IRanges(start = gmuct_mRNA_counts$Position, 
                                                     end = gmuct_mRNA_counts$Position) , 
                                             strand = gmuct_mRNA_counts$Strand)
)

# Check colData
colData(dxd)
counts(dxd)

# ❗️ It needs a lot of RAM (~40G on my computer...) !!!
BPPARAM = MulticoreParam(workers = 8)

# DEXSeq step by step -----------------------------------------------------

## 1. Estimate the size factors
dxd <- estimateSizeFactors(dxd) 

# Check the sizeFactor
colData(dxd)$sizeFactor

## 2. Fit function to compute the variance in dependence on the mean
dxd <- estimateDispersions(dxd) 

# Having the dispersion estimates and the size factors, we can now test for differential exon usage. 
# For each gene, DEXSeq fits a generalized linear model with the formula ~sample + exon + condition:exon 
# and compare it to the smaller model (the null model) ~ sample + exon
#plotDispEsts(dxd)

## 3. Run the test for differential exon/position usage
dxd <- testForDEU(dxd, BPPARAM = BPPARAM) 

## 4. Estimate exon/position fold changes
dxd <- estimateExonFoldChanges(dxd , BPPARAM = BPPARAM)

# Save the dexseq object
saveRDS(dxd, "analyses/GMUCT/DEXSeq/dxd_xrn4_vs_xrn4dne1.rds")

# read back
#dxd <- readRDS("analyses/GMUCT/DEXSeq/dxd_xrn4_vs_xrn4dne1.rds")

# 5. Get DEXseq results 
dxd.res <- DEXSeqResults(dxd)


# Convert dxd to a datatable and remove the "genomicData" column
res.dt <- dplyr::select(as.data.table(dxd.res), -dplyr::starts_with("genomicData")) 
setkey(res.dt, groupID)

# Nr of unique geneID (groupID) with nr of position
res.dt[, .N, by = groupID] # 8453 unique geneID

# Add informations 
gmuct_mRNA_counts[, featureID := gsub(":", "", ID)]
res.dt.annot <- merge.data.table(res.dt, gmuct_mRNA_counts, by.x = c("groupID", "featureID"), by.y = c("gene_id", "featureID"))


#plotMA(res.dt, cex=0.5, alpha=0.2)
#mcols(dxd.res)

# Filtered positions by our thresholds
res.dt.annot.diff <- res.dt.annot[padj <= 0.05 & abs(log2fold_xrn4dne1_xrn4) >= 1]

res.dt.annot[padj <= 0.05 & log2fold_xrn4dne1_xrn4 >= 1][, .N, groupID]
res.dt.annot[padj <= 0.05 & log2fold_xrn4dne1_xrn4 <= 1][, .N, groupID]

res.dt.annot.diff[, .N, by = groupID][order(N)] #  1475 unique geneID

# Save results in excel sheet
fwrite(res.dt, "analyses/GMUCT/DEXSeq/DEXSeq_results.xrn4_vs_xrn4dne1.csv", sep = ";")
fwrite(res.dt.annot, "analyses/GMUCT/DEXSeq/DEXSeq_results.xrn4_vs_xrn4dne1.annnotated.csv", sep = ";")
fwrite(res.dt.annot.diff, "analyses/GMUCT/DEXSeq/DEXSeq_results.xrn4_vs_xrn4dne1.padj5%_log2FC>=1.annotated.csv", sep = ";")

# Plotting geneID with position
dxd.res[dxd.res$groupID == "AT1G01040",]
dxd.res[which(dxd.res$log2fold_xrn4dne1_xrn4 >= 2), ]
f.dxd.res <- dxd.res[which(dxd.res$log2fold_xrn4dne1_xrn4 >= 2), ] 


plotDEXSeq(dxd.res, geneID = "AT1G22190", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
plotDEXSeq(dxd.res, geneID = "AT1G78080", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)

# Graph -------------------------------------------------------------------
sub <- melt(res.dt.annot[groupID == "AT2G43020", 
                 c("transcript_id", "ID", "Chr", "Position", "Strand", "xrn4dne1_R1.rpm", 
                   "xrn4dne1_R2.rpm", "xrn4dne1_R3.rpm", "xrn4_R1.rpm", "xrn4_R2.rpm", "xrn4_R3.rpm", "log2fold_xrn4dne1_xrn4")], 
            id.vars = c("transcript_id", "ID", "Chr", "Position", "Strand", "log2fold_xrn4dne1_xrn4"), value.name = "RPM")


res.dt.annot[groupID == "AT1G22190"]

sub <- melt(res.dt.annot[groupID == "AT2G22430", 
                         c("transcript_id", "ID", "Chr", "Position", "Strand", 
                           "rpmMeans_xrn4", "rpmMeans_xrn4dne1", "log2fold_xrn4dne1_xrn4")], 
            id.vars = c("transcript_id", "ID", "Chr", "Position", "Strand", "log2fold_xrn4dne1_xrn4"), value.name = "RPM")

sub[, Diff := "no diff"]
sub[log2fold_xrn4dne1_xrn4 >= 2, Diff := "log2FC>=2"]
sub[log2fold_xrn4dne1_xrn4 <= -2, Diff := "log2FC<=-2"]

#ggplot(sub, aes(Position, log2(RPM + 1), fill = Diff)) + geom_col() + facet_wrap(~ variable) 

A = ggplot(sub, aes(Position, RPM, color = Diff)) + geom_point(size = 0.5) + facet_wrap(~ variable, ncol = 1) 
B = ggplot(sub, aes(Position, log2fold_xrn4dne1_xrn4, color = Diff)) + geom_point(alpha = 0.5) + facet_wrap(~ variable, ncol = 6) 


C = ggplot(sub, aes(Position, RPM, fill = Diff)) + geom_col() + facet_wrap(~ variable, ncol = 1) + scale_fill_manual(values = c("no diff" = "grey", "log2FC>=2" = "red", "log2FC<=-2"= "steelblue"))
D = ggplot(sub[variable=="rpmMeans_xrn4"], aes(Position, log2fold_xrn4dne1_xrn4, fill = Diff)) + geom_col() + scale_fill_manual(values = c("no diff" = "grey", "log2FC>=2" = "red", "log2FC<=-2"= "steelblue"))

C = ggplot(sub[variable=="rpmMeans_xrn4"], aes(Position, RPM, fill = Diff)) + geom_col() + ylim(c(0, 850)) + scale_fill_manual(values = c("no diff" = "grey", "log2FC>=2" = "red", "log2FC<=2"= "steelblue"))
D = ggplot(sub[variable=="rpmMeans_xrn4dne1"], aes(Position, RPM, fill = Diff)) + geom_col() + ylim(c(0, 850))
E = ggplot(sub[variable=="rpmMeans_xrn4"], aes(Position, log2fold_xrn4dne1_xrn4, fill = Diff)) + geom_col()

library(patchwork)
A / B

C / D + plot_layout(heights = c(0.66, 0.33))


# Missing genes -----------------------------------------------------------

# Fit for gene/exon AT3G29240 threw the next warning(s): the matrix is either rank-deficient or indefinite
# Fit for gene/exon AT4G01250 threw the next warning(s): the matrix is either rank-deficient or indefinite
# Fit for gene/exon AT4G03510 threw the next warning(s): the matrix is either rank-deficient or indefinite
# Fit for gene/exon AT1G14870 threw the next warning(s): the matrix is either rank-deficient or indefinite
# Fit for gene/exon AT1G36060 threw the next warning(s): the matrix is either rank-deficient or indefinite
# Fit for gene/exon AT5G05250 threw the next warning(s): the matrix is either rank-deficient or indefinite
# Fit for gene/exon AT2G21660 threw the next warning(s): the matrix is either rank-deficient or indefinite
# Fit for gene/exon AT2G23340 threw the next warning(s): the matrix is either rank-deficient or indefinite
# Fit for gene/exon AT2G26980 threw the next warning(s): the matrix is either rank-deficient or indefinite
# Fit for gene/exon AT2G28400 threw the next warning(s): the matrix is either rank-deficient or indefinite
