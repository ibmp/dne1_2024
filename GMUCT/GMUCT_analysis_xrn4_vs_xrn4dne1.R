library(rtracklayer)
library(data.table)
library(ggplot2)
library(ggrepel)
library(ggvenn)
library(viridis)
library(gridExtra)
library(patchwork)
library(cli)

# Venn Diagram UP and DOWN  -----------------------------------------------
dexseq_xrn4_vs_xrn4dne1 <- fread("analyses/GMUCT/02_DEXSeq_analysis/DEXSeq_results.xrn4_vs_xrn4dne1.csv") # 155100

#dexseq_xrn4_vs_xrn4dne1[groupID %in% pgreen.dt$Gene][log2fold_xrn4dne1_xrn4 >= 2 & padj <= 0.05][, .N, by = groupID]
#dexseq_xrn4_vs_xrn4dne1[groupID %in% pgreen.dt$Gene][, .N, by = groupID]

#dexseq_xrn4_vs_xrn4dne1[featureID %like% "181383"]


# Number of GMUCT positions
dexseq_xrn4_vs_xrn4dne1[, .N, by = groupID][order(N)]
dexseq_xrn4_vs_xrn4dne1[, Diff := "No diff"]
dexseq_xrn4_vs_xrn4dne1[log2fold_xrn4dne1_xrn4 >= 1 & padj <= 0.05, Diff := "Up"]
dexseq_xrn4_vs_xrn4dne1[log2fold_xrn4dne1_xrn4 <= -1 & padj <= 0.05, Diff := "Down"]
dexseq_xrn4_vs_xrn4dne1[, .N, by = Diff]

View(dexseq_xrn4_vs_xrn4dne1[Diff %in% c("Up", "Down"), .N, by = groupID][order(N)])

# Number of genes with differential positions
dexseq_xrn4_vs_xrn4dne1[, .N, by = c("groupID", "Diff")]

# Venn Diagram Genes with UP and DOWN position  ---------------------------
UP = unique(dexseq_xrn4_vs_xrn4dne1[Diff == "Up"]$groupID)  # 575  (UNIQUE AGI) UP
DOWN = unique(dexseq_xrn4_vs_xrn4dne1[Diff == "Down"]$groupID)  # 1296 (UNIQUE AGI) DOWN
UP_DOWN = intersect(UP, DOWN)
ONLY_DOWN = setdiff(DOWN, UP)
ONLY_UP = setdiff(UP, DOWN)



# Excel fives subtables
library(openxlsx)

wb <- createWorkbook()
addWorksheet(wb,sheetName = "UP"); writeData(wb, sheet = "UP", x = full_dexseq_xrn4_vs_xrn4dne1_annotated[gene_id %in% UP])
addWorksheet(wb, sheetName = "DOWN"); writeData(wb, sheet = "DOWN", x = full_dexseq_xrn4_vs_xrn4dne1_annotated[gene_id %in% DOWN])
addWorksheet(wb, sheetName = "UP_DOWN"); writeData(wb, sheet = "UP_DOWN", x = full_dexseq_xrn4_vs_xrn4dne1_annotated[gene_id %in% UP_DOWN])
addWorksheet(wb, sheetName = "ONLY_DOWN"); writeData(wb, sheet = "ONLY_DOWN", x = full_dexseq_xrn4_vs_xrn4dne1_annotated[gene_id %in% ONLY_DOWN])
addWorksheet(wb, sheetName = "ONLY_UP"); writeData(wb, sheet = "ONLY_UP", x = full_dexseq_xrn4_vs_xrn4dne1_annotated[gene_id %in% ONLY_UP])
saveWorkbook(wb, "results/GMUCT/GMUCT_Venn_genes_classification.xlsx", overwrite = TRUE)

# Stats ONLY UP
only_up_stats <- full_dexseq_xrn4_vs_xrn4dne1_annotated[gene_id %in% ONLY_UP][, .N, by = c("type", "Diff")]
# Stats ONLY DOWN
only_down_stats <- full_dexseq_xrn4_vs_xrn4dne1_annotated[gene_id %in% ONLY_DOWN][, .N, by = c("type", "Diff")]

fwrite(only_up_stats, "results/GMUCT/only_up_stats.csv")
fwrite(only_down_stats, "results/GMUCT/only_down_stats.csv")

length(unique(c(UP, DOWN)))

ggvenn(data = list("UP" = unique(UP), "DOWN" = unique(DOWN)), 
       fill_color = c("#0073C2FF", "#EFC000FF"),
       stroke_size = 0.5, set_name_size = 4) + ggtitle("xrn4dne1_vs_xrn4")

# 3'UTR 5'UTR & CDS -------------------------------------------------------
TAIR10_GFF_PATH = "ressources/TAIR10_annotations/TAIR10_GFF3_genes.clean.gff3"
TAIR10_GENE_MODELS_PATH = "ressources/TAIR10_annotations/TAIR10_representative_gene_models.txt"

# GFF, GENE_MODEL and uORF importation  -----------------------------------------
GFF <- as.data.table(import(TAIR10_GFF_PATH))

setkey(GFF, seqnames, start, end)
GENE_MODELS <- fread(TAIR10_GENE_MODELS_PATH, col.names = "AGI")
GFF <- GFF[transcript_id %in% GENE_MODELS$AGI]
# Subset with only major mRNA isoform
gff_subset <- GFF[type %in% c("five_prime_UTR", "CDS", "three_prime_UTR") & transcript_id %in% GENE_MODELS$AGI, 
                  c("seqnames", "start", "end", "width", "strand", "type", "gene_id", "transcript_id")]
setkey(gff_subset, seqnames, start, end)

# Load the gmuct result table to add RPM to DEXSeq result table
gmuct_counts <- fread("analyses/GMUCT/01_GMUCT_counts/GMUCT_TAIR10_mRNA.min1RPM.csv")
gmuct_counts[, featureID := gsub(":", "", ID)]
# 226464 GMUCT sites in mRNA annotations 

# Merge GMUCT table with DEXseq table                         
full_dexseq_xrn4_vs_xrn4dne1 <- merge(dexseq_xrn4_vs_xrn4dne1[, c("featureID", "exonBaseMean", "pvalue", "padj", "log2fold_xrn4dne1_xrn4", "Diff")],
                                      gmuct_counts, 
                                      by = "featureID", all.x = T)
setkey(full_dexseq_xrn4_vs_xrn4dne1, Chr, Position)
full_dexseq_xrn4_vs_xrn4dne1[, .N, by = gene_id][order(N)]
full_dexseq_xrn4_vs_xrn4dne1[, .N, by = Diff]

# Overlap now the annotation (CDS, UTR)
full_dexseq_xrn4_vs_xrn4dne1[, `:=` (start=Position, end=Position)]
setkey(full_dexseq_xrn4_vs_xrn4dne1, Chr, start, end)
# 157642
full_dexseq_xrn4_vs_xrn4dne1_annotated <- foverlaps(full_dexseq_xrn4_vs_xrn4dne1, gff_subset)
full_dexseq_xrn4_vs_xrn4dne1_annotated[gene_id %in% UP_DOWN]

pdf("results/GMUCT/Repartition_UTR_CDS.pdf", width = 8, height = 8)


# Venn diagram
ggvenn(data = list("UP" = unique(UP), "DOWN" = unique(DOWN)), 
       fill_color = c("#0073C2FF", "#EFC000FF"),
       stroke_size = 0.5, set_name_size = 4) + ggtitle("xrn4dne1_vs_xrn4")

# Draw barplot for gene with UP and DOWN
#full_dexseq_xrn4_vs_xrn4dne1_annotated[!is.na(type)][, .N, by = c("type", "Diff")][order(type)]
A_tbl <- tableGrob(full_dexseq_xrn4_vs_xrn4dne1_annotated[!is.na(type)][, .N, by = c("type", "Diff")][order(type)], rows=NULL)
A = ggplot(full_dexseq_xrn4_vs_xrn4dne1_annotated[!is.na(type)], 
       aes(x = type, fill = Diff)) + geom_bar() + 
  ggtitle(label = "Repartition of all gmuct positions detected") +
  xlab("Annotations")
A + A_tbl

#grid.arrange(plot1, plot2, tbl, nrow = 3, heights = c(2, 2, 0.5))
B_tbl <- tableGrob(full_dexseq_xrn4_vs_xrn4dne1_annotated[gene_id %in% UP_DOWN][, .N, by = c("type", "Diff")][order(type)], rows = NULL)
B <- ggplot(full_dexseq_xrn4_vs_xrn4dne1_annotated[gene_id %in% UP_DOWN], 
       aes(x = type, fill = Diff)) + geom_bar() + 
  ggtitle(label = "Repartition of all gmuct positions detected", subtitle = "in the 396 transcripts") +
  xlab("Annotations")
B + B_tbl

C_tbl <- tableGrob(full_dexseq_xrn4_vs_xrn4dne1_annotated[gene_id %in% UP_DOWN & Diff %in% c("Up", "Down")][, .N, by = c("type", "Diff")][order(type)], rows = NULL)
C <- ggplot(full_dexseq_xrn4_vs_xrn4dne1_annotated[gene_id %in% UP_DOWN & Diff %in% c("Up", "Down")], 
       aes(x = type, fill = Diff)) + geom_bar() + 
  ggtitle(label = "Repartition of differential gmuct positions", subtitle = "in the 396 transcripts (showing up & down patterns)") +
  xlab("Annotations")
C + C_tbl

D_tbl <- tableGrob(full_dexseq_xrn4_vs_xrn4dne1_annotated[gene_id %in% UP_DOWN & Diff %in% c("Up", "Down")][, .N, by = c("type", "Diff")][order(type)], rows = NULL)
D <- ggplot(full_dexseq_xrn4_vs_xrn4dne1_annotated[gene_id %in% UP_DOWN & Diff %in% c("Up", "Down")], 
       aes(x = type, fill = Diff)) + geom_bar(position = "dodge") + 
  ggtitle(label = "Repartition of differential gmuct positions", subtitle = "in the 396 transcripts (showing up & down patterns)") +
  xlab("Annotations")
D + D_tbl

# For all transcripts
E_tbl <- tableGrob(full_dexseq_xrn4_vs_xrn4dne1_annotated[!is.na(type) & Diff %in% c("Up", "Down")][, .N, by = c("type", "Diff")][order(type)], rows = NULL)
E <- ggplot(full_dexseq_xrn4_vs_xrn4dne1_annotated[!is.na(type) & Diff %in% c("Up", "Down")], 
       aes(x = type, fill = Diff)) + geom_bar() + 
  ggtitle("Repartition of differential gmuct positions", subtitle =  "within all transcripts (1475)") +
  xlab("Annotations")
E + E_tbl

length(unique(full_dexseq_xrn4_vs_xrn4dne1_annotated[!is.na(type) & Diff %in% c("Up", "Down")]$gene_id)) # 1487

G_tbl <- tableGrob(full_dexseq_xrn4_vs_xrn4dne1_annotated[!is.na(type) & Diff %in% c("Up", "Down")][, .N, by = c("type", "Diff")][order(type)], rows = NULL)
G <- ggplot(full_dexseq_xrn4_vs_xrn4dne1_annotated[!is.na(type) & Diff %in% c("Up", "Down")], 
       aes(x = type, fill = Diff)) + geom_bar(position = "dodge") + 
  ggtitle("Repartition of differential gmuct positions", subtitle =  "within all transcripts (1475)") +
  xlab("Annotations")
G + G_tbl
dev.off()

# Which Position is first (UP or DOWN) ------------------------------------
# Brin + extract first position 
A_positive <- full_dexseq_xrn4_vs_xrn4dne1[Diff == "Up" & Strand == "+", .SD[1], by = gene_id][, c("transcript_id", "featureID", "Position", "Strand", "log2fold_xrn4dne1_xrn4")]
B_positive <- full_dexseq_xrn4_vs_xrn4dne1[Diff == "Down" & Strand == "+", .SD[1], by = gene_id][, c("transcript_id", "featureID", "Position", "Strand", "log2fold_xrn4dne1_xrn4")]
C_positive <- merge(A_positive, B_positive, by = "transcript_id", suffixes = c(".UP", ".DOWN"))
C_positive[, First := ifelse(Position.DOWN < Position.UP, "DownFirst", "UpFirst")]
C_positive[, .N, by = First ]
# First N
# UpFirst 128
# DownFirst  61

# Brin - extracting first position is the last since gene in the other direction 
A_negative <- full_dexseq_xrn4_vs_xrn4dne1[Diff == "Up" & Strand == "-", .SD[.N], by = gene_id][, c("transcript_id", "featureID", "Position", "Strand", "log2fold_xrn4dne1_xrn4")]
B_negative <- full_dexseq_xrn4_vs_xrn4dne1[Diff == "Down" & Strand == "-", .SD[.N], by = gene_id][, c("transcript_id", "featureID", "Position", "Strand", "log2fold_xrn4dne1_xrn4")]
C_negative <- merge(A_negative, B_negative, by = "transcript_id", suffixes = c(".UP", ".DOWN"))
C_negative[, First := ifelse(Position.DOWN < Position.UP, "UpFirst", "DownFirst")]
C_negative[, .N, by = First ]
# First N
# UpFirst 130
# DownFirst  77


# Graphs publications -----------------------------------------------------

# Graph by genes
# AT1G22190 --> OK 
# AT1G78080 --> OK
# AT1G52342 --> OK
sub_xrn4_vs_xrn4dne1 <- melt(full_dexseq_xrn4_vs_xrn4dne1[gene_id == "AT3G22120", c("transcript_id", "ID", "Chr", "Position", "Strand", "xrn4dne1_R1.rpm",
                           "xrn4dne1_R2.rpm", "xrn4dne1_R3.rpm", "xrn4_R1.rpm", "xrn4_R2.rpm", "xrn4_R3.rpm", "log2fold_xrn4dne1_xrn4", "padj", "Diff")],
            id.vars = c("transcript_id", "ID", "Chr", "Position", "Strand", "log2fold_xrn4dne1_xrn4", "padj", "Diff"), value.name = "RPM")

sub_xrn4_vs_xrn4dne1[, Diff := "no diff"]

sub_xrn4_vs_xrn4dne1[log2fold_xrn4dne1_xrn4 >= 1 & padj <= 0.05, Diff := "log2FC>=1"]
sub_xrn4_vs_xrn4dne1[log2fold_xrn4dne1_xrn4 <= -1 & padj <= 0.05, Diff := "log2FC<=-1"]
sub_xrn4_vs_xrn4dne1[, variable := factor(sub_xrn4_vs_xrn4dne1$variable, levels = c("rpmMeans_xrn4", "rpmMeans_xrn4dne1"))]
sub_xrn4_vs_xrn4dne1[, .N, Diff]
                     
library(patchwork)
#ggplot(sub, aes(Position, log2(RPM + 1), fill = Diff)) + geom_col() + facet_wrap(~ variable) 
A = ggplot(sub_xrn4_vs_xrn4dne1, aes(Position, RPM, color = Diff)) + geom_point(size = 0.5) + facet_wrap(~ variable, ncol = 1) 
B = ggplot(sub_xrn4_vs_xrn4dne1, aes(Position, log2fold_xrn4dne1_xrn4, color = Diff)) + geom_point(alpha = 0.5) + facet_wrap(~ variable, ncol = 6) 
A / B

C = ggplot(sub_xrn4_vs_xrn4dne1, aes(Position, log2(RPM + 1), fill = Diff)) + geom_col() + facet_wrap(~ variable, ncol = 1) + 
  scale_fill_manual(values = c("no diff" = "grey", "log2FC>=1" = "red", "log2FC<=-1"= "steelblue")) + theme_bw(base_size = 14) +
  xlab("") + ylab("log2(RPM + 1)")

sub[variable=="rpmMeans_xrn4" & Diff == "log2FC<=-1"]

D = ggplot(sub[variable=="rpmMeans_xrn4"], aes(Position, log2fold_xrn4dne1_xrn4, fill = Diff)) + geom_col() + 
  scale_fill_manual(values = c("no diff" = "grey", "log2FC>=1" = "red", "log2FC<=-1"= "steelblue")) + theme_bw(base_size = 14) +
  xlab("") + ylab("log2FC (xrn4dne1/xrn4)")

E = ggplot(gff_subset[gene_id == "AT1G78080"]) +
  # geom_segment(aes(x = ifelse(strand == "-", end, start), xend = ifelse(strand == "-", start - 80, end + 80),
  #                  y = transcript_id, yend = transcript_id), color = "black", arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +  
  geom_segment(data = GFF[gene_id == "AT1G78080" & type == "mRNA"], aes(x = ifelse(strand == "-", end, start), xend = ifelse(strand == "-", start - 80, end + 80),
                   y = transcript_id, yend = transcript_id), color = "black", arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  geom_tile(aes(x = start + (width / 2), y = transcript_id, width = width, height = 0.2, fill = type)) +
  labs(y = "Isoform", x = "Genomic coordinates", fill= "Region") + scale_fill_brewer(palette = "Paired") +
  theme_minimal() + 
  theme(panel.grid.major = element_line(linewidth = rel(0.5)), panel.grid.minor = element_blank(), axis.text = element_text(size = 7)) 
E

# Combined them 
p_combined = C / D / E
p_ranges_x <- c(ggplot_build(E)$layout$panel_scales_x[[1]]$range$range) # extract limits 
myplot <- p_combined + plot_layout(guides = "collect", heights = c(0.60, 0.3, 0.1)) + theme_bw(base_size = 14) & xlim(min(p_ranges_x), max(p_ranges_x)) 
myplot

#ggsave(plot = myplot, path = "results/GMUCT/", filename = "AT1G22190.GMUCT_DEXSeq_log2FC>=1.pdf", width = 12, height = 8, dpi = 300)
ggsave(plot = myplot, path = "results/GMUCT/", filename = "AT1G78080.GMUCT_DEXSeq_log2FC>=1.pdf", width = 12, height = 8, dpi = 300)
#ggsave(plot = myplot, path = "results/GMUCT/", filename = "AT2G22430.GMUCT_DEXSeq_log2FC>=1.pdf", width = 12, height = 8, dpi = 300)

# Stats UTR CDS -----------------------------------------------------------
ggplot(full_res, aes(type, fill = type)) + geom_bar() + scale_fill_brewer(palette = "Paired") + theme_bw() +
  #geom_text(stat = "count", aes(label = after_stat(count)), vjust = -1)
  geom_text(
    aes(label=paste(after_stat(round(count / sum(count) * 100, 2)), "%")),
    stat='count', 
    vjust = -0.5
  )
ggsave("results/GMUCT/Stats_UTR_CDS.xrn4dne1_vs_xrn4.pdf")


# log2FC>=2 ---------------------------------------------------------------
# log2fold_xrn4dne1_xrn4
makeGraph_log2FC2 <- function(dt, gene, path, filename){
  cli_alert_info("Analysing {gene}")
  subset <- melt(dt[gene_id == gene, c("transcript_id", "ID", "Chr", "Position", "Strand", "rpmMeans_xrn4dne1", "rpmMeans_xrn4", "log2fold_xrn4dne1_xrn4", "padj")], 
              id.vars = c("transcript_id", "ID", "Chr", "Position", "Strand", "log2fold_xrn4dne1_xrn4", "padj"), value.name = "RPM")
  subset[, Diff := "no diff"]
  subset[log2fold_xrn4dne1_xrn4 >= 2 & padj <= 0.05, Diff := "log2FC>=2"]
  subset[log2fold_xrn4dne1_xrn4 <= -2 & padj <= 0.05, Diff := "log2FC<=-2"]
  print(subset[, .N, Diff])
  subset[, variable := factor(subset$variable, levels = c("rpmMeans_xrn4", "rpmMeans_xrn4dne1"))]
  subset[, Diff := factor(subset$Diff, levels = c("no diff", "log2FC>=2", "log2FC<=-2"))]
  print(subset)
  # RPM means 
  A <- ggplot(subset[order(Diff)], aes(Position, log2(RPM + 1), color = Diff, shape = Diff, fill = Diff)) + geom_point() + #position=position_jitter(h=0.1, w=0.1)
    scale_fill_manual(values = c("log2FC>=2" = "darkgreen", "log2FC<=-2"= "darkred", "no diff" = "lightgrey")) + 
    scale_color_manual(values = c("log2FC>=2" = "darkgreen", "log2FC<=-2"= "darkred", "no diff" = "lightgrey")) +
    scale_shape_manual(values = c("log2FC>=2" = 24, "log2FC<=-2"= 25, "no diff" = 21)) +
    theme_bw(base_size = 14) + xlab("") + ylab("log2(RPM + 1)") + facet_wrap(~ variable, ncol = 1) + theme(panel.grid.minor = element_blank())
  
  # DEXSeq results
  B <- ggplot(subset[variable=="rpmMeans_xrn4"][order(Diff)], aes(Position, log2fold_xrn4dne1_xrn4, color = Diff, fill = Diff, shape = Diff)) + 
    geom_point(size = 1.5) + 
    scale_fill_manual(values = c("log2FC>=2" = "darkgreen", "log2FC<=-2"= "darkred", "no diff" = "lightgrey")) + 
    scale_color_manual(values = c("log2FC>=2" = "darkgreen", "log2FC<=-2"= "darkred", "no diff" = "lightgrey")) +
    scale_shape_manual(values = c("log2FC>=2" = 24, "log2FC<=-2"= 25, "no diff" = 21)) +
    theme_bw(base_size = 14) +xlab("") + ylab("log2FC (xrn4dne1/xrn4)") + theme(panel.grid.minor = element_blank())
  
  # Draw isoform from gff 
  C <- ggplot(gff_subset[gene_id == gene]) +
    geom_segment(data = GFF[gene_id == gene & type == "mRNA"], 
                 aes(x = ifelse(strand == "-", end + 20, start - 20), 
                     xend = ifelse(strand == "-", start - 80, end + 80), 
                     y = transcript_id, 
                     yend = transcript_id), 
                 color = "black", arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
    geom_tile(aes(x = start + (width / 2), y = transcript_id, width = width, height = 0.2, fill = type)) +
    labs(y = "Isoform", x = "Genomic coordinates", fill= "Region") + scale_fill_brewer(palette = "Paired") +
    theme_minimal() + 
    theme(panel.grid.major = element_line(linewidth = rel(0.5)), panel.grid.minor = element_blank(), axis.text = element_text(size = 7)) 
  
  # Combined them 
  p_combined = A / B / C
  p_ranges_x <- c(ggplot_build(C)$layout$panel_scales_x[[1]]$range$range) # extract limits 
  # Match x axis to all graphs
  myplot <- p_combined + plot_layout(guides = "collect", heights = c(0.6, 0.3, 0.1)) + theme_bw(base_size = 14) & xlim(min(p_ranges_x), max(p_ranges_x)) 
  print(myplot)
  cli_alert_info("Saving into {path}{filename}")
  ggsave(plot = myplot, path = path, filename = filename, width = 12, height = 8, dpi = 300)
  cli_alert_success("Done!")
}

makeGraph_log2FC2(dt = full_dexseq_xrn4_vs_xrn4dne1, gene = "AT1G78080", path = "results/GMUCT/", filename = "AT1G78080.xrn4dne1_vs_xrn4.GMUCT_DEXSeq_log2FC>=2.pdf")
makeGraph_log2FC2(dt = full_dexseq_xrn4_vs_xrn4dne1, gene = "AT2G22430", path = "results/GMUCT/", filename = "AT2G22430.xrn4dne1_vs_xrn4.GMUCT_DEXSeq_log2FC>=2.pdf")
makeGraph_log2FC2(dt = full_dexseq_xrn4_vs_xrn4dne1, gene = "AT1G22190", path = "results/GMUCT/", filename = "AT1G22190.xrn4dne1_vs_xrn4.GMUCT_DEXSeq_log2FC>=2.pdf")


 # Log2FC>=1 ---------------------------------------------------------------
# log2fold_xrn4dne1_xrn4
makeGraph_log2FC1 <- function(dt, gene, path, filename){
  cli_alert_info("Analysing {gene}")
  subset <- melt(dt[gene_id == gene, c("transcript_id", "ID", "Chr", "Position", "Strand", "rpmMeans_xrn4dne1", "rpmMeans_xrn4", "log2fold_xrn4dne1_xrn4", "padj")], 
                 id.vars = c("transcript_id", "ID", "Chr", "Position", "Strand", "log2fold_xrn4dne1_xrn4", "padj"), value.name = "RPM")
  subset[, Diff := "no diff"]
  subset[log2fold_xrn4dne1_xrn4 >= 1 & padj <= 0.05, Diff := "log2FC>=1"]
  subset[log2fold_xrn4dne1_xrn4 <= -1 & padj <= 0.05, Diff := "log2FC<=-1"]
  print(subset[, .N, Diff])
  subset[, variable := factor(subset$variable, levels = c("rpmMeans_xrn4", "rpmMeans_xrn4dne1"))]
  subset[, Diff := factor(subset$Diff, levels = c("no diff", "log2FC>=1", "log2FC<=-1"))]
  print(subset)
  # RPM means 
  # A <- ggplot(subset[order(Diff)], aes(Position, log2(RPM + 1), color = Diff, shape = Diff, fill = Diff)) + geom_point() + #position=position_jitter(h=0.1, w=0.1)
  #   scale_fill_manual(values = c("log2FC>=1" = "darkgreen", "log2FC<=-1"= "darkred", "no diff" = "lightgrey")) + 
  #   scale_color_manual(values = c("log2FC>=1" = "darkgreen", "log2FC<=-1"= "darkred", "no diff" = "lightgrey")) +
  #   scale_shape_manual(values = c("log2FC>=1" = 24, "log2FC<=-1"= 25, "no diff" = 21)) +
  #   theme_bw(base_size = 14) + xlab("") + ylab("log2(RPM + 1)") + facet_wrap(~ variable, ncol = 1) + theme(panel.grid.minor = element_blank())
  # 
  A <- ggplot(subset[order(Diff)], aes(Position, RPM, fill = Diff)) + geom_col() + #position=position_jitter(h=0.1, w=0.1)
    theme_bw(base_size = 14) + xlab("") + ylab("RPM") + facet_wrap(~ variable, ncol = 1) + 
    theme(panel.grid.minor = element_blank()) + 
    scale_fill_manual(values = c("log2FC>=1" = "darkred", "log2FC<=-1"= "blue", "no diff" = "lightgrey")) 
  
  # DEXSeq results
  B <- ggplot(subset[variable=="rpmMeans_xrn4"][order(Diff)], aes(Position, log2fold_xrn4dne1_xrn4, color = Diff, fill = Diff, shape = Diff)) + 
    geom_point(size = 1.5) + 
    scale_fill_manual(values = c("log2FC>=1" = "darkred", "log2FC<=-1"= "blue", "no diff" = "lightgrey")) + 
    scale_color_manual(values = c("log2FC>=1" = "darkred", "log2FC<=-1"= "blue", "no diff" = "lightgrey")) +
    scale_shape_manual(values = c("log2FC>=1" = 24, "log2FC<=-1"= 25, "no diff" = 21)) +
    theme_bw(base_size = 14) +xlab("") + ylab("log2FC (xrn4dne1/xrn4)") + theme(panel.grid.minor = element_blank())
  
  # Draw isoform from gff 
  C <- ggplot(gff_subset[gene_id == gene]) +
    geom_segment(data = GFF[gene_id == gene & type == "mRNA"], 
                 aes(x = ifelse(strand == "-", end + 20, start - 20), 
                     xend = ifelse(strand == "-", start - 80, end + 80), 
                     y = transcript_id, 
                     yend = transcript_id), 
                 color = "black", arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
    geom_tile(aes(x = start + (width / 2), y = transcript_id, width = width, height = 0.2, fill = type)) +
    labs(y = "Isoform", x = "Genomic coordinates", fill= "Region") + scale_fill_brewer(palette = "Paired") +
    theme_minimal() + 
    theme(panel.grid.major = element_line(linewidth = rel(0.5)), panel.grid.minor = element_blank(), axis.text = element_text(size = 7)) 
  
  # Combined them 
  p_combined = A / B / C
  p_ranges_x <- c(ggplot_build(C)$layout$panel_scales_x[[1]]$range$range) # extract limits 
  # Match x axis to all graphs
  myplot <- p_combined + plot_layout(guides = "collect", heights = c(0.6, 0.3, 0.1)) + theme_bw(base_size = 14) & xlim(min(p_ranges_x), max(p_ranges_x)) 
  print(myplot)
  cli_alert_info("Saving into {path}{filename}")
  ggsave(plot = myplot, path = path, filename = filename, width = 12, height = 8, dpi = 300)
  cli_alert_success("Done!")
}

makeGraph_log2FC1(dt = full_dexseq_xrn4_vs_xrn4dne1, gene = "AT1G78080", path = "results/GMUCT/tmp/", filename = "AT1G78080.xrn4dne1_vs_xrn4.GMUCT_DEXSeq_log2FC1.pdf")
makeGraph_log2FC1(dt = full_dexseq_xrn4_vs_xrn4dne1, gene = "AT2G22430", path = "results/GMUCT/tmp/", filename = "AT2G22430.xrn4dne1_vs_xrn4.GMUCT_DEXSeq_log2FC1.pdf")
makeGraph_log2FC1(dt = full_dexseq_xrn4_vs_xrn4dne1, gene = "AT1G22190", path = "results/GMUCT/tmp/", filename = "AT1G22190.xrn4dne1_vs_xrn4.GMUCT_DEXSeq_log2FC1.pdf")

makeGraph_log2FC1(dt = full_dexseq_xrn4_vs_xrn4dne1, gene = "AT3G02470", path = "results/GMUCT/tmp/", filename = "AT3G02470.xrn4dne1_vs_xrn4.GMUCT_DEXSeq_log2FC1.pdf")
makeGraph_log2FC1(dt = full_dexseq_xrn4_vs_xrn4dne1, gene = "AT5G25280", path = "results/GMUCT/tmp/", filename = "AT5G25280.xrn4dne1_vs_xrn4.GMUCT_DEXSeq_log2FC1.pdf")
makeGraph_log2FC1(dt = full_dexseq_xrn4_vs_xrn4dne1, gene = "AT5G63160", path = "results/GMUCT/tmp/", filename = "AT5G63160.xrn4dne1_vs_xrn4.GMUCT_DEXSeq_log2FC1.pdf")
makeGraph_log2FC1(dt = full_dexseq_xrn4_vs_xrn4dne1, gene = "AT1G57680", path = "results/GMUCT/tmp/", filename = "AT1G57680.xrn4dne1_vs_xrn4.GMUCT_DEXSeq_log2FC1.pdf")
makeGraph_log2FC1(dt = full_dexseq_xrn4_vs_xrn4dne1, gene = "AT1G72150", path = "results/GMUCT/tmp/", filename = "AT1G72150.xrn4dne1_vs_xrn4.GMUCT_DEXSeq_log2FC1.pdf")
makeGraph_log2FC1(dt = full_dexseq_xrn4_vs_xrn4dne1, gene = "AT5G11090", path = "results/GMUCT/tmp/", filename = "AT5G11090.xrn4dne1_vs_xrn4.GMUCT_DEXSeq_log2FC1.pdf")
makeGraph_log2FC1(dt = full_dexseq_xrn4_vs_xrn4dne1, gene = "AT2G22430", path = "results/GMUCT/tmp/", filename = "AT2G22430.xrn4dne1_vs_xrn4.GMUCT_DEXSeq_log2FC1.pdf")

# Genes of interest
GOI <- c("AT1G78080", "AT2G43020", "AT3G18780", "AT5G11580", "AT5G64260", "AT5G19120", 
         "AT1G08630", "AT1G22360", "AT4G34138", "AT2G22010", "AT5G21940", "AT3G59940", 
         "AT1G32170", "AT3G16150", "AT1G22190", "AT1G57680")

# 
GOI_DOWN <- unique(c("AT5G67300", "AT5G64260", "AT2G44670", "AT5G67300", "AT3G04070", 
        "AT1G21910", "AT2G44670", "AT2G44670", "AT5G64260", "AT2G46220", 
        "AT2G47400", "AT5G67300", "AT3G59940", "AT3G48990", "AT2G44670", 
        "AT5G25350", "AT5G64260"))

GOI_UP <- unique(c("AT2G02060", "AT1G76410", "AT2G23430", "AT4G13940", "AT2G36950",
            "AT1G03080", "AT2G31090", "AT1G80920", "AT5G43580", "AT5G02380"))

plyr::l_ply(GOI, function(AGI) {
  #try(makeGraph_log2FC1_withWT(dt = full_dexseq_xrn4_vs_xrn4dne1, gene = AGI, path = "results/GMUCT/tmp/", filename = paste0(AGI, ".all_samples.RPM_means.GMUCT_DEXSeq_log2FC1.pdf")))
  try(makeGraph_log2FC1(dt = full_dexseq_xrn4_vs_xrn4dne1, gene = AGI, path = "results/GMUCT/tmp/", filename = paste0(AGI, ".xrn4dne1_vs_xrn4.GMUCT_DEXSeq_log2FC1.pdf")))
})

plyr::l_ply(GOI_DOWN, function(AGI) {
  try(makeGraph_log2FC1_withWT(dt = full_dexseq_xrn4_vs_xrn4dne1, gene = AGI, path = "results/GMUCT/tmp/", filename = paste0(AGI, ".xrn4dne1_vs_xrn4.ONLY_DOWN.GMUCT_DEXSeq_log2FC1.pdf")))
})

plyr::l_ply(GOI_UP, function(AGI) {
  try(makeGraph_log2FC1_withWT(dt = full_dexseq_xrn4_vs_xrn4dne1, gene = AGI, path = "results/GMUCT/tmp/", filename = paste0(AGI, ".xrn4dne1_vs_xrn4.ONLY_UP.GMUCT_DEXSeq_log2FC1.pdf")))
})


full_dexseq_xrn4_vs_xrn4dne1[gene_id == "AT2G22430" & rpmMeans_xrn4dne1 > 100]

makeGraph_log2FC1_withWT <- function(dt, gene, path, filename){
  cli_alert_info("Analysing {gene}")
  subset <- melt(dt[gene_id == gene, c("transcript_id", "ID", "Chr", "Position", "Strand", "rpmMeans_WT", "rpmMeans_xrn4dne1", "rpmMeans_xrn4", "log2fold_xrn4dne1_xrn4", "padj")], 
                 id.vars = c("transcript_id", "ID", "Chr", "Position", "Strand", "log2fold_xrn4dne1_xrn4", "padj"), value.name = "RPM")
  subset[, Diff := "no diff"]
  subset[log2fold_xrn4dne1_xrn4 >= 1 & padj <= 0.05, Diff := "log2FC>=1"]
  subset[log2fold_xrn4dne1_xrn4 <= -1 & padj <= 0.05, Diff := "log2FC<=-1"]
  print(subset[, .N, Diff])
  subset[, variable := factor(subset$variable, levels = c("rpmMeans_WT", "rpmMeans_xrn4", "rpmMeans_xrn4dne1"))]
  subset[, Diff := factor(subset$Diff, levels = c("no diff", "log2FC>=1", "log2FC<=-1"))]
  print(subset)
 
  A <- ggplot(subset[order(Diff)], aes(Position, RPM, fill = Diff)) + geom_col() + #position=position_jitter(h=0.1, w=0.1)
    theme_bw(base_size = 14) + xlab("") + ylab("RPM") + facet_wrap(~ variable, ncol = 1) + 
    theme(panel.grid.minor = element_blank()) + 
    scale_fill_manual(values = c("log2FC>=1" = "darkred", "log2FC<=-1"= "blue", "no diff" = "lightgrey")) 
  
  # DEXSeq results
  B <- ggplot(subset[variable=="rpmMeans_xrn4"][order(Diff)], aes(Position, log2fold_xrn4dne1_xrn4, color = Diff, fill = Diff, shape = Diff)) + 
    geom_point(size = 1.5) + 
    scale_fill_manual(values = c("log2FC>=1" = "darkred", "log2FC<=-1"= "blue", "no diff" = "lightgrey")) + 
    scale_color_manual(values = c("log2FC>=1" = "darkred", "log2FC<=-1"= "blue", "no diff" = "lightgrey")) +
    scale_shape_manual(values = c("log2FC>=1" = 24, "log2FC<=-1"= 25, "no diff" = 21)) +
    theme_bw(base_size = 14) +xlab("") + ylab("log2FC (xrn4dne1/xrn4)") + theme(panel.grid.minor = element_blank())
  
  # Draw isoform from gff 
  C <- ggplot(gff_subset[gene_id == gene]) +
    geom_segment(data = GFF[gene_id == gene & type == "mRNA"], 
                 aes(x = ifelse(strand == "-", end + 20, start - 20), 
                     xend = ifelse(strand == "-", start - 80, end + 80), 
                     y = transcript_id, 
                     yend = transcript_id), 
                 color = "black", arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
    geom_tile(aes(x = start + (width / 2), y = transcript_id, width = width, height = 0.2, fill = type)) +
    labs(y = "Isoform", x = "Genomic coordinates", fill= "Region") + scale_fill_brewer(palette = "Paired") +
    theme_minimal() + 
    theme(panel.grid.major = element_line(linewidth = rel(0.5)), panel.grid.minor = element_blank(), axis.text = element_text(size = 7)) 
  
  # Combined them 
  p_combined = A / B / C
  p_ranges_x <- c(ggplot_build(C)$layout$panel_scales_x[[1]]$range$range) # extract limits 
  # Match x axis to all graphs
  myplot <- p_combined + plot_layout(guides = "collect", heights = c(0.6, 0.3, 0.1)) + theme_bw(base_size = 14) & xlim(min(p_ranges_x), max(p_ranges_x)) 
  print(myplot)
  cli_alert_info("Saving into {path}{filename}")
  ggsave(plot = myplot, path = path, filename = filename, width = 12, height = 8, dpi = 300)
  cli_alert_success("Done!")
}

makeGraph_log2FC1_withWT(dt = full_dexseq_xrn4_vs_xrn4dne1, gene = "AT1G78080", path = "results/GMUCT/tmp/", filename = "AT1G78080.all_samples_RPM_means.GMUCT_DEXSeq_log2FC1.pdf")


makeGraph_log2FC1_test <- function(dt, gene, path, filename){
  cli_alert_info("Analysing {gene}")
  subset <- melt(dt[gene_id == gene, c("transcript_id", "ID", "Chr", "Position", "Strand", "rpmMeans_xrn4dne1", "rpmMeans_xrn4", "log2fold_xrn4dne1_xrn4", "padj")], 
                 id.vars = c("transcript_id", "ID", "Chr", "Position", "Strand", "log2fold_xrn4dne1_xrn4", "padj"), value.name = "RPM")
  subset[, Diff := "no diff"]
  subset[log2fold_xrn4dne1_xrn4 >= 1 & padj <= 0.05, Diff := "log2FC>=1"]
  subset[log2fold_xrn4dne1_xrn4 <= -1 & padj <= 0.05, Diff := "log2FC<=-1"]
  print(subset[, .N, Diff])
  subset[, variable := factor(subset$variable, levels = c("rpmMeans_xrn4", "rpmMeans_xrn4dne1"))]
  subset[, Diff := factor(subset$Diff, levels = c("no diff", "log2FC>=1", "log2FC<=-1"))]
  print(subset)
  # RPM means 
  # A <- ggplot(subset[order(Diff)], aes(Position, log2(RPM + 1), color = Diff, shape = Diff, fill = Diff)) + geom_point() + #position=position_jitter(h=0.1, w=0.1)
  #   scale_fill_manual(values = c("log2FC>=1" = "darkgreen", "log2FC<=-1"= "darkred", "no diff" = "lightgrey")) + 
  #   scale_color_manual(values = c("log2FC>=1" = "darkgreen", "log2FC<=-1"= "darkred", "no diff" = "lightgrey")) +
  #   scale_shape_manual(values = c("log2FC>=1" = 24, "log2FC<=-1"= 25, "no diff" = 21)) +
  #   theme_bw(base_size = 14) + xlab("") + ylab("log2(RPM + 1)") + facet_wrap(~ variable, ncol = 1) + theme(panel.grid.minor = element_blank())
  # 
  A <- ggplot(subset[order(Diff)], aes(Position, log2(RPM + 1), fill = Diff)) + geom_col() + #position=position_jitter(h=0.1, w=0.1)
    theme_bw(base_size = 14) + xlab("") + ylab("log2(RPM + 1)") + facet_wrap(~ variable, ncol = 1) + 
    theme(panel.grid.minor = element_blank()) + 
    scale_fill_manual(values = c("log2FC>=1" = "darkgreen", "log2FC<=-1"= "darkred", "no diff" = "lightgrey")) 
  
  # DEXSeq results
  B <- ggplot(subset[variable=="rpmMeans_xrn4"][order(Diff)], aes(Position, log2fold_xrn4dne1_xrn4, color = Diff, fill = Diff, shape = Diff)) + 
    geom_point(size = 1.5) + 
    scale_fill_manual(values = c("log2FC>=1" = "darkgreen", "log2FC<=-1"= "darkred", "no diff" = "lightgrey")) + 
    scale_color_manual(values = c("log2FC>=1" = "darkgreen", "log2FC<=-1"= "darkred", "no diff" = "lightgrey")) +
    scale_shape_manual(values = c("log2FC>=1" = 24, "log2FC<=-1"= 25, "no diff" = 21)) +
    theme_bw(base_size = 14) +xlab("") + ylab("log2FC (xrn4dne1/xrn4)") + theme(panel.grid.minor = element_blank())
  
  # Draw isoform from gff 
  C <- ggplot(gff_subset[gene_id == gene]) +
    geom_segment(data = GFF[gene_id == gene & type == "mRNA"], 
                 aes(x = ifelse(strand == "-", end + 20, start - 20), 
                     xend = ifelse(strand == "-", start - 80, end + 80), 
                     y = transcript_id, 
                     yend = transcript_id), 
                 color = "black", arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
    geom_tile(aes(x = start + (width / 2), y = transcript_id, width = width, height = 0.2, fill = type)) +
    labs(y = "Isoform", x = "Genomic coordinates", fill= "Region") + scale_fill_brewer(palette = "Paired") +
    theme_minimal() + 
    theme(panel.grid.major = element_line(linewidth = rel(0.5)), panel.grid.minor = element_blank(), axis.text = element_text(size = 7)) 
  
  # Combined them 
  p_combined = A / B / C
  p_ranges_x <- c(ggplot_build(C)$layout$panel_scales_x[[1]]$range$range) # extract limits 
  # Match x axis to all graphs
  myplot <- p_combined + plot_layout(guides = "collect", heights = c(0.6, 0.3, 0.1)) + theme_bw(base_size = 14) & xlim(min(p_ranges_x), max(p_ranges_x)) 
  print(myplot)
  cli_alert_info("Saving into {path}{filename}")
  ggsave(plot = myplot, path = path, filename = filename, width = 12, height = 8, dpi = 300)
  cli_alert_success("Done!")
}

makeGraph_log2FC1(dt = full_dexseq_xrn4_vs_xrn4dne1, gene = "AT1G78080", path = "results/GMUCT/tmp", filename = "AT1G78080.xrn4dne1_vs_xrn4.GMUCT_DEXSeq_log2FC1.pdf")
makeGraph_log2FC1(dt = full_dexseq_xrn4_vs_xrn4dne1, gene = "AT2G22430", path = "results/GMUCT/tmp", filename = "AT2G22430.xrn4dne1_vs_xrn4.GMUCT_DEXSeq_log2FC1.pdf")
makeGraph_log2FC1(dt = full_dexseq_xrn4_vs_xrn4dne1, gene = "AT1G22190", path = "results/GMUCT/tmp", filename = "AT1G22190.xrn4dne1_vs_xrn4.GMUCT_DEXSeq_log2FC1.pdf")




# -------------------------------------------------------------------------
# full_dexseq_col0_vs_xrn4
makeGraph_log2FC2_col0_vs_xrn4 <- function(dt, gene, path, filename){
  cli_alert_info("Analysing {gene}")
  subset <- melt(dt[gene_id == gene, c("transcript_id", "ID", "Chr", "Position", "Strand", "rpmMeans_xrn4", "rpmMeans_WT", "log2fold_xrn4dne1_xrn4", "padj")], 
                 id.vars = c("transcript_id", "ID", "Chr", "Position", "Strand", "log2fold_xrn4dne1_xrn4", "padj"), value.name = "RPM")
  subset[, Diff := "no diff"]
  subset[log2fold_xrn4dne1_xrn4 >= 2 & padj <= 0.05, Diff := "log2FC>=2"]
  subset[log2fold_xrn4dne1_xrn4 <= -2 & padj <= 0.05, Diff := "log2FC<=-2"]
  print(subset[, .N, Diff])
  subset[, variable := factor(subset$variable, levels = c("rpmMeans_xrn4", "rpmMeans_WT"))]
  subset[, Diff := factor(subset$Diff, levels = c("no diff", "log2FC>=2", "log2FC<=-2"))]
  print(subset)
  # RPM means 
  A <- ggplot(subset[order(Diff)], aes(Position, log2(RPM + 1), color = Diff, shape = Diff, fill = Diff)) + geom_point() + #position=position_jitter(h=0.1, w=0.1)
    scale_fill_manual(values = c("log2FC>=2" = "darkgreen", "log2FC<=-2"= "darkred", "no diff" = "lightgrey")) + 
    scale_color_manual(values = c("log2FC>=2" = "darkgreen", "log2FC<=-2"= "darkred", "no diff" = "lightgrey")) +
    scale_shape_manual(values = c("log2FC>=2" = 24, "log2FC<=-2"= 25, "no diff" = 21)) +
    theme_bw(base_size = 14) + xlab("") + ylab("log2(RPM + 1)") + facet_wrap(~ variable, ncol = 1) + theme(panel.grid.minor = element_blank())
  
  # DEXSeq results
  B <- ggplot(subset[variable=="rpmMeans_xrn4"][order(Diff)], aes(Position, log2fold_xrn4dne1_xrn4, color = Diff, fill = Diff, shape = Diff)) + 
    geom_point(size = 1.5) + 
    scale_fill_manual(values = c("log2FC>=2" = "darkgreen", "log2FC<=-2"= "darkred", "no diff" = "lightgrey")) + 
    scale_color_manual(values = c("log2FC>=2" = "darkgreen", "log2FC<=-2"= "darkred", "no diff" = "lightgrey")) +
    scale_shape_manual(values = c("log2FC>=2" = 24, "log2FC<=-2"= 25, "no diff" = 21)) +
    theme_bw(base_size = 14) +xlab("") + ylab("log2FC (xrn4/WT)") + theme(panel.grid.minor = element_blank())
  
  # Draw isoform from gff 
  C <- ggplot(gff_subset[gene_id == gene]) +
    geom_segment(data = GFF[gene_id == gene & type == "mRNA"], 
                 aes(x = ifelse(strand == "-", end + 20, start - 20), 
                     xend = ifelse(strand == "-", start - 80, end + 80), 
                     y = transcript_id, 
                     yend = transcript_id), 
                 color = "black", arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
    geom_tile(aes(x = start + (width / 2), y = transcript_id, width = width, height = 0.2, fill = type)) +
    labs(y = "Isoform", x = "Genomic coordinates", fill= "Region") + scale_fill_brewer(palette = "Paired") +
    theme_minimal() + 
    theme(panel.grid.major = element_line(linewidth = rel(0.5)), panel.grid.minor = element_blank(), axis.text = element_text(size = 7)) 
  
  # Combined them 
  p_combined = A / B / C
  p_ranges_x <- c(ggplot_build(C)$layout$panel_scales_x[[1]]$range$range) # extract limits 
  # Match x axis to all graphs
  myplot <- p_combined + plot_layout(guides = "collect", heights = c(0.6, 0.3, 0.1)) + theme_bw(base_size = 14) & xlim(min(p_ranges_x), max(p_ranges_x)) 
  print(myplot)
  cli_alert_info("Saving into {path}{filename}")
  ggsave(plot = myplot, path = path, filename = filename, width = 12, height = 8, dpi = 300)
  cli_alert_success("Done!")
}

makeGraph_log2FC2_col0_vs_xrn4(dt = full_dexseq_col0_vs_xrn4, gene = "AT1G78080", path = "results/GMUCT/", filename = "AT1G78080.col0_vs_xrn4.GMUCT_DEXSeq_log2FC>=2.pdf")
makeGraph_log2FC2_col0_vs_xrn4(dt = full_dexseq_col0_vs_xrn4, gene = "AT2G22430", path = "results/GMUCT/", filename = "AT2G22430.col0_vs_xrn4.GMUCT_DEXSeq_log2FC>=2.pdf")
makeGraph_log2FC2_col0_vs_xrn4(dt = full_dexseq_col0_vs_xrn4, gene = "AT1G22190", path = "results/GMUCT/", filename = "AT1G22190.col0_vs_xrn4.GMUCT_DEXSeq_log2FC>=2.pdf")


# log2FC>=1 col0_vs_xrn4 --------------------------------------------------
makeGraph_log2FC1_col0_vs_xrn4 <- function(dt, gene, path, filename){
  cli_alert_info("Analysing {gene}")
  subset <- melt(dt[gene_id == gene, c("transcript_id", "ID", "Chr", "Position", "Strand", "rpmMeans_xrn4", "rpmMeans_WT", "log2fold_xrn4dne1_xrn4", "padj")], 
                 id.vars = c("transcript_id", "ID", "Chr", "Position", "Strand", "log2fold_xrn4dne1_xrn4", "padj"), value.name = "RPM")
  subset[, Diff := "no diff"]
  subset[log2fold_xrn4dne1_xrn4 >= 1 & padj <= 0.05, Diff := "log2FC>=1"]
  subset[log2fold_xrn4dne1_xrn4 <= -1 & padj <= 0.05, Diff := "log2FC<=-1"]
  print(subset[, .N, Diff])
  subset[, variable := factor(subset$variable, levels = c("rpmMeans_WT", "rpmMeans_xrn4"))]
  subset[, Diff := factor(subset$Diff, levels = c("no diff", "log2FC>=1", "log2FC<=-1"))]
  print(subset)
  # RPM means 
  A <- ggplot(subset[order(Diff)], aes(Position, log2(RPM + 1), color = Diff, shape = Diff, fill = Diff)) + geom_point() + #position=position_jitter(h=0.1, w=0.1)
    scale_fill_manual(values = c("log2FC>=1" = "darkgreen", "log2FC<=-1"= "darkred", "no diff" = "lightgrey")) + 
    scale_color_manual(values = c("log2FC>=1" = "darkgreen", "log2FC<=-1"= "darkred", "no diff" = "lightgrey")) +
    scale_shape_manual(values = c("log2FC>=1" = 24, "log2FC<=-1"= 25, "no diff" = 21)) +
    theme_bw(base_size = 14) + xlab("") + ylab("log2(RPM + 1)") + facet_wrap(~ variable, ncol = 1) + theme(panel.grid.minor = element_blank())
  
  # DEXSeq results
  B <- ggplot(subset[variable=="rpmMeans_WT"][order(Diff)], aes(Position, log2fold_xrn4dne1_xrn4, color = Diff, fill = Diff, shape = Diff)) + 
    geom_point(size = 1.5) + 
    scale_fill_manual(values = c("log2FC>=1" = "darkgreen", "log2FC<=-1"= "darkred", "no diff" = "lightgrey")) + 
    scale_color_manual(values = c("log2FC>=1" = "darkgreen", "log2FC<=-1"= "darkred", "no diff" = "lightgrey")) +
    scale_shape_manual(values = c("log2FC>=1" = 24, "log2FC<=-1"= 25, "no diff" = 21)) +
    theme_bw(base_size = 14) +xlab("") + ylab("log2FC (xrn4/WT)") + theme(panel.grid.minor = element_blank())
  
  # Draw isoform from gff 
  C <- ggplot(gff_subset[gene_id == gene]) +
    geom_segment(data = GFF[gene_id == gene & type == "mRNA"], 
                 aes(x = ifelse(strand == "-", end + 20, start - 20), 
                     xend = ifelse(strand == "-", start - 80, end + 80), 
                     y = transcript_id, 
                     yend = transcript_id), 
                 color = "black", arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
    geom_tile(aes(x = start + (width / 2), y = transcript_id, width = width, height = 0.2, fill = type)) +
    labs(y = "Isoform", x = "Genomic coordinates", fill= "Region") + scale_fill_brewer(palette = "Paired") +
    theme_minimal() + 
    theme(panel.grid.major = element_line(linewidth = rel(0.5)), panel.grid.minor = element_blank(), axis.text = element_text(size = 7)) 
  
  # Combined them 
  p_combined = A / B / C
  p_ranges_x <- c(ggplot_build(C)$layout$panel_scales_x[[1]]$range$range) # extract limits 
  # Match x axis to all graphs
  myplot <- p_combined + plot_layout(guides = "collect", heights = c(0.6, 0.3, 0.1)) + theme_bw(base_size = 14) & xlim(min(p_ranges_x), max(p_ranges_x)) 
  print(myplot)
  cli_alert_info("Saving into {path}{filename}")
  ggsave(plot = myplot, path = path, filename = filename, width = 12, height = 8, dpi = 300)
  cli_alert_success("Done!")
}

makeGraph_log2FC1_col0_vs_xrn4(dt = full_dexseq_col0_vs_xrn4, gene = "AT1G78080", path = "results/GMUCT/", filename = "AT1G78080.col0_vs_xrn4.GMUCT_DEXSeq_log2FC>=1.pdf")
makeGraph_log2FC1_col0_vs_xrn4(dt = full_dexseq_col0_vs_xrn4, gene = "AT2G22430", path = "results/GMUCT/", filename = "AT2G22430.col0_vs_xrn4.GMUCT_DEXSeq_log2FC>=1.pdf")
makeGraph_log2FC1_col0_vs_xrn4(dt = full_dexseq_col0_vs_xrn4, gene = "AT1G22190", path = "results/GMUCT/", filename = "AT1G22190.col0_vs_xrn4.GMUCT_DEXSeq_log2FC>=1.pdf")


# Periodicity analysis ----------------------------------------------------
test <- gmuct_counts[rpmSums >= 10][, diff := Position - data.table::shift(Position), by = Chr]

# Add the rleid to search back after rle 
test[, rleid := rleid(diff)]

# Compute the Run Length Encoding
dt_rle <- rle(test$rleid)

test[diff > 1, .N, by = c("diff", "rleid")][N >= 2][order(N)]

test[rleid == "38733"]
test[transcript_id == "AT5G08450.1"]

library(gghighlight)
ggplot(test[transcript_id == "AT5G08450.1"], aes(x = Position, y = rpmMeans)) + geom_point() +
  gghighlight(rleid == "38733")

