library(BSgenome.Athaliana.TAIR.TAIR9)
library(rtracklayer)
library(data.table)
library(openxlsx)
library(ggplot2)
library(ggvenn)
library(DEXSeq)
library(tidyr)

theme_set(theme_bw()) # change ggplot default theme

TAIR10_GFF_PATH = "../ressources/TAIR10_annotations/TAIR10_GFF3_genes.clean.gff3"
TAIR10_GENE_MODELS_PATH = "../ressources/TAIR10_annotations/TAIR10_representative_gene_models.txt"

# GFF and GENE_MODEL importation  -----------------------------------------
GFF <- as.data.table(import(TAIR10_GFF_PATH))
mRNA <- GFF[type %in% c("mRNA")]
GENE_MODELS <- fread(TAIR10_GENE_MODELS_PATH, col.names = "AGI")

# Select the gene_model (gives the most represented transcript/isoform per gene) 
mRNA <- mRNA[transcript_id %in% GENE_MODELS$AGI] 
setnames(mRNA, c("seqnames", "start", "end", "strand"), c("Chr", "Start", "End", "Strand"))
setkey(mRNA, Chr, Start, End)

BSGENOME = getSeq(BSgenome.Athaliana.TAIR.TAIR9) # Genome as BSgenome object

# ------------------------------------------------------------------------
SAMPLES = c("uniq-Col0-0-R1" = "WT_R1", "uniq-Col0-0-R2" = "WT_R2", "uniq-Col0-0-R3" = "WT_R3", 
            "uniq-xrn4-0-R1" = "xrn4_R1", "uniq-xrn4-0-R2" = "xrn4_R2", "uniq-xrn4-0-R3" = "xrn4_R3", 
            "uniq-xrn4dne1-0-R1" = "xrn4dne1_R1", "uniq-xrn4dne1-0-R2" = "xrn4dne1_R2", "uniq-xrn4dne1-0-R3" = "xrn4dne1_R3")

# Load the 5' coverage count table from all the GMUCT samples
gmuct_counts <- fread("GMUCT_all_samples.counts.txt.gz", 
                 colClasses = c(list(character = c("Sample","Chr", "Strand", "Nucleotide", "ID"), 
                                     numeric = c("Position", "Coverage", "RPM")))
                 )

gmuct_counts[, Name := SAMPLES[Sample], by = Sample]
gmuct_counts[Coverage >= 10, .N, by = Name]

# 21,769,823 rows
gmuct_counts[, ":=" (Start = Position, End = Position)] # Add Start & End for data.table::foverlaps function 
setkey(gmuct_counts, Chr, Start, End)

# Match 5' position to mRNA annotation in GTF
all_overlaps <- foverlaps(mRNA[, c("Chr", "Start","End", "Strand", "width", "type", "transcript_id", "gene_id")], gmuct_counts, mult = "all")
# 21,132,112 Positions 
setnames(all_overlaps, c("i.Start", "i.End", "i.Strand"), c("gene.Start", "gene.End", "gene.Strand")) 

# Get everything that matches in the mRNA annotations and also where the strands matches.
mRNA_overlaps <- all_overlaps[!is.na(Sample) & Strand == gene.Strand]

# 1. Identification of XRN4 sites -----------------------------------------
xrn4_vs_WT <- setDT(pivot_wider(mRNA_overlaps[Name %in% c("WT_R1", "WT_R2", "WT_R3", "xrn4_R1", "xrn4_R2", "xrn4_R3")],
            id_cols = c("Chr", "Position", "Strand", "Nucleotide", "ID", "transcript_id", "gene_id", "width", "gene.Start", "gene.End", "gene.Strand"), 
            names_from = "Name", 
            values_from = c("Coverage", "RPM"), 
            values_fill = 0,
            names_glue = "{Name}.{.value}", 
            names_sort = T
            ))

ggvenn(data = list("xrn4_R1" = xrn4_vs_WT[xrn4_R1.RPM >= 1]$ID, 
                   "xrn4_R2" = xrn4_vs_WT[xrn4_R2.RPM >= 1]$ID, 
                   "xrn4_R3" = xrn4_vs_WT[xrn4_R3.RPM >= 1]$ID),
       fill_color = c("#0073C2FF", "#EFC000FF", "#ABC000FF"),
       stroke_size = 0.5, set_name_size = 4
)

# 49960 sites ≥ 1 RPM in all xrn4 biological replicates
# In 4321 unique transcripts
xrn4_vs_WT_1CPM <- xrn4_vs_WT[xrn4_R1.RPM >= 1 & xrn4_R2.RPM >= 1 & xrn4_R3.RPM >= 1] # & `RPM_uniq-Col0-0-R1` >= 1 & `RPM_uniq-Col0-0-R2` >= 1 & `RPM_uniq-Col0-0-R3` >= 1]
length(unique(xrn4_vs_WT_1CPM$ID))
xrn4_vs_WT_1CPM[, .N, gene_id]

fwrite(x = xrn4_vs_WT_1CPM, file = "01_GMUCT_counts/01_xrn4_vs_WT_all1CPM.csv", sep = ",")


# 2. Sites enriched in xrn4 compared to Col0 ------------------------------

# Log2 FC xrn4/Col0 ≥ 1 in all biological replicates
ZERO_RPM = 0.01 # PseudoCount, we use 0.01 RPM when RPM is 0 to avoid zero division 

xrn4_vs_WT_1CPM[, log2FC_xrn4_vs_WT_R1 := round(log2(xrn4_R1.RPM / WT_R1.RPM), 3)]
xrn4_vs_WT_1CPM[log2FC_xrn4_vs_WT_R1 == Inf, log2FC_xrn4_vs_WT_R1 := round(log2(xrn4_R1.RPM / ZERO_RPM), 3)]
xrn4_vs_WT_1CPM[, log2FC_xrn4_vs_WT_R2 := round(log2(xrn4_R2.RPM / WT_R2.RPM), 3)]
xrn4_vs_WT_1CPM[log2FC_xrn4_vs_WT_R2 == Inf, log2FC_xrn4_vs_WT_R2 := round(log2(xrn4_R2.RPM / ZERO_RPM), 3)]
xrn4_vs_WT_1CPM[, log2FC_xrn4_vs_WT_R3 := round(log2(xrn4_R3.RPM / WT_R3.RPM), 3)]
xrn4_vs_WT_1CPM[log2FC_xrn4_vs_WT_R3 == Inf, log2FC_xrn4_vs_WT_R3 := round(log2(xrn4_R3.RPM / ZERO_RPM), 3)]

ggvenn(data = list("log2FC_xrn4_vs_WT_R1" = xrn4_vs_WT_1CPM[log2FC_xrn4_vs_WT_R1 >= 1]$ID, 
                   "log2FC_xrn4_vs_WT_R2" = xrn4_vs_WT_1CPM[log2FC_xrn4_vs_WT_R2 >= 1]$ID, 
                   "log2FC_xrn4_vs_WT_R3" = xrn4_vs_WT_1CPM[log2FC_xrn4_vs_WT_R3 >= 1]$ID),
       fill_color = c("#0073C2FF", "#EFC000FF", "#ABC000FF"),
       stroke_size = 0.5, set_name_size = 4
)

# 12856 common enriched sites with log2foldchange ≥ 1 
xrn4_WT_log2FC1 <- xrn4_vs_WT_1CPM[log2FC_xrn4_vs_WT_R1 >= 1 & log2FC_xrn4_vs_WT_R2 >= 1 & log2FC_xrn4_vs_WT_R3 >= 1]
fwrite(x = xrn4_WT_log2FC1, file = "02_xrn4_WT.log2FC>=1.csv", sep = ",")

xrn4_WT_log2FC2 <- xrn4_vs_WT_1CPM[log2FC_xrn4_vs_WT_R1 >= 2 & log2FC_xrn4_vs_WT_R2 >= 2 & log2FC_xrn4_vs_WT_R3 >= 2][, .N, gene_id]
fwrite(x = xrn4_WT_log2FC2, file = "02_xrn4_WT.log2FC>=2.csv", sep = ",")

length(unique(xrn4_WT_log2FC1$ID))
xrn4_WT_log2FC1[, .N, gene_id]

# Check position in multiple transcripts (overlapping transcripts)
xrn4_WT_log2FC1$ID[duplicated(xrn4_WT_log2FC1$ID)]
xrn4_WT_log2FC1[ID == "Chr1_6709786_+:A"]

# 3. Sites enriched in xrn4 compared to xrn4dne1: DNE1 targets ------------
xrn4_vs_xrn4dne1 <- setDT(pivot_wider(mRNA_overlaps[Name %in% c("xrn4_R1", "xrn4_R2", "xrn4_R3", "xrn4dne1_R1", "xrn4dne1_R2", "xrn4dne1_R3")],
                                id_cols = c("Chr", "Position", "Strand", "Nucleotide", "ID", "transcript_id", "gene_id", "width", "gene.Start", "gene.End", "gene.Strand"), 
                                names_from = "Name", 
                                values_from = c("Coverage", "RPM"), 
                                values_fill = 0,
                                names_glue = "{Name}.{.value}", 
                                names_sort = T
))

subset <- xrn4_vs_xrn4dne1[ID %in% xrn4_WT_log2FC1$ID]

# Log2 FC xrn4/xrn4dne1 ≥ 1
subset[, log2FC_xrn4_vs_xrn4dne1_R1 := round(log2(xrn4_R1.RPM / xrn4dne1_R1.RPM), 3)]
subset[log2FC_xrn4_vs_xrn4dne1_R1 == Inf, log2FC_xrn4_vs_xrn4dne1_R1 := round(log2(xrn4_R1.RPM / ZERO_RPM), 3)]
subset[, log2FC_xrn4_vs_xrn4dne1_R2 := round(log2(xrn4_R2.RPM/ xrn4dne1_R2.RPM), 3)]
subset[log2FC_xrn4_vs_xrn4dne1_R2 == Inf, log2FC_xrn4_vs_xrn4dne1_R2 := round(log2(xrn4_R2.RPM / ZERO_RPM), 3)]
subset[, log2FC_xrn4_vs_xrn4dne1_R3 := round(log2(xrn4_R3.RPM / xrn4dne1_R3.RPM), 3)]
subset[log2FC_xrn4_vs_xrn4dne1_R3 == Inf, log2FC_xrn4_vs_xrn4dne1_R3 := round(log2(xrn4_R3.RPM / ZERO_RPM), 3)]

ggvenn(data = list("log2FC_xrn4_vs_xrn4dne1_R1" = subset[log2FC_xrn4_vs_xrn4dne1_R1 >= 1]$ID, 
                   "log2FC_xrn4_vs_xrn4dne1_R2" = subset[log2FC_xrn4_vs_xrn4dne1_R2 >= 1]$ID, 
                   "log2FC_xrn4_vs_xrn4dne1_R3" = subset[log2FC_xrn4_vs_xrn4dne1_R3 >= 1]$ID),
       fill_color = c("#0073C2FF", "#EFC000FF", "#ABC000FF"),
       stroke_size = 0.5, set_name_size = 4
)

subset_xrn4_xrn4dne1_log2FC1 <- subset[log2FC_xrn4_vs_xrn4dne1_R1 >= 1 & log2FC_xrn4_vs_xrn4dne1_R2 >= 1 & log2FC_xrn4_vs_xrn4dne1_R3 >= 1]
subset_xrn4_xrn4dne1_log2FC2 <- subset[log2FC_xrn4_vs_xrn4dne1_R1 >= 2 & log2FC_xrn4_vs_xrn4dne1_R2 >= 2 & log2FC_xrn4_vs_xrn4dne1_R3 >= 2]

subset[log2FC_xrn4_vs_xrn4dne1_R1 <= -1 & log2FC_xrn4_vs_xrn4dne1_R2 <= -1 & log2FC_xrn4_vs_xrn4dne1_R3 <= -1]

length(unique(subset_xrn4_xrn4dne1_log2FC1$ID))
subset_xrn4_xrn4dne1_log2FC1[, .N, by = gene_id][order(N)]
 
# TOP 3 genes

# 1. AT1G22190 93
# 2. AT1G78080 72
# 3. AT3G02470 25

AT1G22190_sub <- melt(
  subset[gene_id == "AT1G22190", c("gene_id", "Position", "ID",  
                                   "xrn4_R1.RPM", "xrn4_R2.RPM", "xrn4_R3.RPM", 
                                   "xrn4dne1_R1.RPM", "xrn4dne1_R2.RPM", "xrn4dne1_R3.RPM", 
                                   "log2FC_xrn4_vs_xrn4dne1_R1", "log2FC_xrn4_vs_xrn4dne1_R2", "log2FC_xrn4_vs_xrn4dne1_R3")], 
  id.vars = c("gene_id", "Position", "ID"), 
  measure.vars = list(c("xrn4_R1.RPM", "xrn4_R2.RPM", "xrn4_R3.RPM", 
                        "xrn4dne1_R1.RPM", "xrn4dne1_R2.RPM", "xrn4dne1_R3.RPM"), 
                      c("log2FC_xrn4_vs_xrn4dne1_R1", "log2FC_xrn4_vs_xrn4dne1_R2", "log2FC_xrn4_vs_xrn4dne1_R3")),
  value.name = c("RPM", "log2FC"))

ggplot(AT1G22190_sub, aes(x = Position, y = log2FC, color = Sample)) + geom_point()
ggplot(mRNA_overlaps[gene_id == "AT1G22190"], aes(Position, RPM, Name)) + geom_col() + facet_wrap(~Name)


# Filter data ----------------------------------------------------------
wt_samples <- c("WT_R1", "WT_R2", "WT_R3")
xrna_samples <- c("xrn4_R1", "xrn4_R2", "xrn4_R3")
xrn4dne1_samples <- c("xrn4dne1_R1", "xrn4dne1_R2", "xrn4dne1_R3")

# Add rowSums/Means in selected SAMPLES
mRNA_count.dc[, ":=" (countSums = rowSums(.SD), countMeans = round(rowMeans(.SD), 3)), .SDcols = paste0(SAMPLES, ".count")]
mRNA_rpm.dc[, ":=" (rpmSums  = rowSums(.SD), rpmMeans = round(rowMeans(.SD), 3)), .SDcols = paste0(names(SAMPLES), ".rpm")]
mRNA_rpm.dc[, ':='(rpmMeans_WT = round(rowMeans(.SD), 3)), .SDcols = paste0(wt_samples, ".rpm")]
mRNA_rpm.dc[, ':='(rpmMeans_xrn4 = round(rowMeans(.SD), 3)), .SDcols = paste0(xrna_samples, ".rpm")]
mRNA_rpm.dc[, ':='(rpmMeans_xrn4dne1 = round(rowMeans(.SD), 3)), .SDcols = paste0(xrn4dne1_samples, ".rpm")]

mRNA_count_rpm <- merge(mRNA_count.dc, mRNA_rpm.dc)

# Filter position where rmpMeans is lower than 1 RPM in the 3 conditions (WT, xrn4, xrn4dne1)
mRNA_count_rpm_min1RPM <- mRNA_count_rpm[rpmMeans_xrn4 >= 1 | rpmMeans_WT >= 1 | rpmMeans_xrn4dne1 >= 1] # 226'464
mRNA_count_rpm_min1RPM[, .N, by = transcript_id][order(N)]

fwrite(mRNA_count_rpm_min1RPM, "GMUCT_TAIR10_mRNA.min1RPM.csv", sep = ";")
gmuct.dt <- fread("analyses/GMUCT/GMUCT_COUNTS/GMUCT_TAIR10_mRNA.min1RPM.csv")

# Descriptive statistics
ggvenn(data = list("Col0" = gmuct.dt[rpmMeans_WT >= 1]$ID, 
                   "xrn4" = gmuct.dt[rpmMeans_xrn4 >= 1]$ID, 
                   "xrn4dne1" = gmuct.dt[rpmMeans_xrn4dne1 >= 1]$ID), 
       stroke_size = 0.5, set_name_size = 4, text_size = 3)

gmuct.dt[rpmMeans_xrn4 >= 1]
gmuct.dt[rpmMeans_xrn4dne1 >= 1]

# Periodicity analysis ----------------------------------------------------

GFF <- as.data.table(import(TAIR10_GFF_PATH))
setkey(GFF, seqnames, start, end)
GENE_MODELS <- fread(TAIR10_GENE_MODELS_PATH, col.names = "AGI")

# Subset with only major mRNA isoform
gff_subset <- GFF[type %in% c("CDS") & transcript_id %in% GENE_MODELS$AGI, 
                  c("seqnames", "start", "end", "width", "strand", "type", "gene_id", "transcript_id")]
setkey(gff_subset, seqnames, start, end)

# ATG start
start_fwd <- gff_subset[strand == "+", .SD[1], transcript_id] # First CDS
start_rev <- gff_subset[strand == "-", .SD[.N], transcript_id] # Last CDS 
ATG <- rbind(start_fwd, start_rev)
ATG[strand == "+", ATG := start]
ATG[strand == "-", ATG := end]
setkey(ATG, seqnames, start, end)

# STOP codon
stop_fwd <- gff_subset[strand == "+", .SD[.N], transcript_id] # Last CDS
stop_rev <- gff_subset[strand == "-", .SD[1], transcript_id] # First CDS
stop_codon[strand == "+", stop_codon := end]
stop_codon[strand == "-", stop_codon := start]
stop_codon <- rbind(stop_fwd, stop_rev)
setkey(stop_codon, seqnames, start, end)

test <- merge(mRNA_overlaps, ATG[, c("transcript_id", "ATG")], by = "transcript_id", all.x = T)
test <- merge(test, stop_codon[, c("transcript_id", "stop_codon")], by = "transcript_id", all.x = T)

test[Strand == "-", ATG_relative := ATG - Position]
test[Strand == "-", stop_codon_relative := stop_codon - Position]

test[Strand == "+", ATG_relative := Position - ATG]
test[Strand == "+", stop_codon_relative := Position - stop_codon]

test[, Condition := tstrsplit(Name, "_", keep = 1)]

ggplot(test[ATG_relative <= 50 & ATG_relative >= -100][, list(RPM = sum(RPM)), by = c("Sample", "ATG_relative")], aes(ATG_relative, RPM, color = Sample)) + 
  geom_line()

ggplot(test[stop_codon_relative <= 50 & stop_codon_relative >= -100][, list(RPM = sum(RPM)), by = c("Sample", "stop_codon_relative")], 
       aes(stop_codon_relative, RPM, color = Sample)) + 
  geom_line() + geom_vline(xintercept = -18)

# By condition with sum of RPM in the position
ggplot(test[ATG_relative <= 50 & ATG_relative >= -100][, list("RPM" = sum(RPM)), by = c("Condition", "ATG_relative")], 
       aes(ATG_relative, RPM, color = Condition)) + 
  geom_line()

ggplot(test[stop_codon_relative <= 50 & stop_codon_relative >= -100][, list("RPM" = sum(RPM)), by = c("Condition", "stop_codon_relative")], 
       aes(stop_codon_relative, RPM, color = Condition)) + 
  geom_line() + geom_vline(xintercept = -18, linetype="dashed")

# By condition (Count of position) 
ggplot(test[ATG_relative <= 50 & ATG_relative >= -100][, list("N" = sum(.N)), by = c("Condition", "ATG_relative")], 
       aes(ATG_relative, N, color = Condition)) + 
  geom_line()

ggplot(test[stop_codon_relative <= 50 & stop_codon_relative >= -100][, list("N" = sum(.N)), by = c("Condition", "stop_codon_relative")], 
       aes(stop_codon_relative, N, color = Condition)) + 
  geom_line() + geom_vline(xintercept = -18, linetype="dashed")

GOI <- test[, .N, by = c("gene_id", "Sample")][N >= 200][order(N)]$gene_id

fwrite(x = test[gene_id %in% GOI][, c("transcript_id", "Chr", "Position", 
                                      "Coverage", "Strand", "RPM", "Name", "gene_id", 
                                      "ATG", "stop_codon", "ATG_relative", "stop_codon_relative", "Condition")], file = "~/Desktop/GMUCT_cleavage_position.csv.gz", sep = ",")

subtest <- fread("~/Desktop/GMUCT_cleavage_position.csv.gz")
# Global view all samples 
ggplot(subtest[ATG_relative <= 50 & ATG_relative >= -100][, list(RPM = sum(RPM)), by = c("Name", "ATG_relative")], 
       aes(ATG_relative, RPM, color = Name)) + 
  geom_line()

# Global view per conditon (3 replicates per condition)
ggplot(subtest[stop_codon_relative <= 50 & stop_codon_relative >= -100][, list("N" = sum(.N)), by = c("Condition", "stop_codon_relative")], 
       aes(stop_codon_relative, N, color = Condition)) + 
  geom_line() + geom_vline(xintercept = -18, linetype="dashed")

# Profile for one transcript
ggplot(subtest[gene_id == "AT1G01960"],
       aes(Position, RPM, color = Condition)) + 
  geom_line() 

# Focus on region around the start & stop codons
ggplot(subtest[stop_codon_relative <= 50 & stop_codon_relative >= -100][gene_id == "AT5G20250"],
       aes(stop_codon_relative, RPM, color = Condition)) + 
  geom_line() + geom_vline(xintercept = -18, linetype="dashed")

ggplot(subtest[ATG_relative <= 50 & ATG_relative >= -100][gene_id == "AT5G20250"],
       aes(ATG_relative, RPM, color = Condition)) + 
  geom_line() 

# By gene
ggplot(test[stop_codon_relative <= 50 & stop_codon_relative >= -100][gene_id == "AT5G20250"],  aes(stop_codon_relative, RPM, color = Condition)) + 
  geom_line() + geom_vline(xintercept = -18, linetype="dashed")

ggplot(test[ATG_relative <= 50 & ATG_relative >= -100][gene_id == "AT1G01960"],
       aes(ATG_relative, RPM, color = Condition)) + 
  geom_line() 

test[ATG_relative <= 50 & ATG_relative >= -100][gene_id == "AT1G01960"]

ggplot(test[gene_id == "AT1G01960"],
       aes(Position, RPM, color = Condition)) + 
  geom_line() 



