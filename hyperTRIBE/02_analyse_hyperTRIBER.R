library(rtracklayer)
library(data.table)
library(patchwork)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(scales)
library(cli)

# ANNOTATION FILEPATHS
TAIR10_GFF_PATH         = "ressources/TAIR10_annotations/TAIR10_GFF3_genes.clean.gff3"
TAIR10_GENE_MODELS_PATH = "ressources/TAIR10_annotations/TAIR10_representative_gene_models.txt"

# Load gff and gene_models
GFF <- as.data.table(import(TAIR10_GFF_PATH))
setkey(GFF, seqnames, start, end)
GENE_MODELS <- fread(TAIR10_GENE_MODELS_PATH, col.names = "AGI")

# Table with 3'UTR, 5'UTR and CDS information
gff_subset <- GFF[type %in% c("five_prime_UTR", "CDS", "three_prime_UTR") & transcript_id %in% GENE_MODELS$AGI, 
                  c("seqnames", "start", "end", "width", "strand", "type", "gene_id", "transcript_id")]

# Load all your hyperTRIBER results (.csv format) in one big data.table
hypertribe_results <- list.files("analyses/HyperTRIBE/01_HyperTRIBER_TAIR10", "*.with_AGI_description.csv$", recursive = T, full.names = T)
names(hypertribe_results) <- basename(dirname(hypertribe_results))

# Save into one excel file with one tab per result
res_list <- lapply(hypertribe_results, fread)
write.xlsx(res_list, file = "analyses/HyperTRIBE/01_HyperTRIBER_TAIR10/All_versus.HyperTRIBER.results.xlsx")

# Put everyting in a table with a new column versus
res <- rbindlist(res, idcol = "Versus")

# Create a new column for the position of the editing site from the name column --> Chr1_32658,-
res[, c("gene_id", "start", "end", "width", "strand") := NULL]
res[, c("Position") := tstrsplit(name, "_", keep = 2)]
res[, c("Position", "strand") := tstrsplit(Position, ",")]
res[, Position := as.numeric(Position)]
res[, start := Position]
res[, end := Position]
setkey(res, seqnames, start, end)

# Filtering options to reduce false positive
res_subset <- res[tags_treat >= 10 & padj <= 0.01 & ref == "A" & targ == "G" & abs(fold_change >= 1)]

# Overlap GFF annotations with editing positions
res_subset_annot <- foverlaps(res_subset, gff_subset)
res_subset_annot[, `:=`(i.start=NULL, i.end = NULL)]
setkey(res_subset_annot, seqnames, start, end)

# INFO: an edit position can be found in multiple annotations, since annotations can overlap each other.
# like this position Chr1_1399716 found in AT1G04940.1 and AT1G04945.3
res_subset_annot[transcript_id != i.transcript_id]

# If not in 3/5'UTR or CDS, it's located in the intron
res_subset_annot[is.na(type), type := "Intron"]

# Get stats of the number of hypertribe position found in the different type
ggplot(res_subset_annot[type != "Intron"][, .N, by = c("Versus", "type")][Versus %in% c("ADAR_vs_DNE1-D153N", "ADAR_vs_DNE1")], aes(x = Versus, y = N, fill = type)) + 
  geom_col(position = position_fill(reverse = TRUE)) + 
  coord_flip() + theme_bw() + 
  scale_fill_brewer(labels = c("five_prime_UTR" = "5'UTR", "CDS" = "CDS", "three_prime_UTR" = "3'UTR"), palette = "Paired") +
  scale_y_continuous(labels = scales::percent) + ylab("Proportion of editing sites")

ggplot(res_subset_annot[Versus %in% c("ADAR_vs_DNE1-D153N", "ADAR_vs_DNE1")], aes(x = Versus, fill = type)) + geom_bar() + coord_flip() + theme_bw() + scale_fill_brewer(palette = "Dark2")

makeGraph <- function(dt, gene, path, filename){
  cli_alert_info("Analysing {gene}")
  subset <- dt[transcript_id == gene]
  subset[, Diff := "no diff"]
  subset[fold_change >= 1 & padj <= 0.05, Diff := "log2FC>=1"]
  subset[fold_change <= -1 & padj <= 0.05, Diff := "log2FC<=-1"]
  print(subset[, .N, Diff])
  subset[, Diff := factor(subset$Diff, levels = c("no diff", "log2FC>=1", "log2FC<=-1"))]
  print(subset)
  # RPM means 
  A <- ggplot(subset[order(Diff)], aes(Position, fold_change, color = Diff, shape = Diff, fill = Diff)) + geom_point() + #position=position_jitter(h=0.1, w=0.1)
    scale_fill_manual(values = c("log2FC>=1" = "darkgreen", "log2FC<=-1"= "darkred", "no diff" = "lightgrey")) + 
    scale_color_manual(values = c("log2FC>=1" = "darkgreen", "log2FC<=-1"= "darkred", "no diff" = "lightgrey")) +
    scale_shape_manual(values = c("log2FC>=1" = 24, "log2FC<=-1"= 25, "no diff" = 21)) +
    theme_bw(base_size = 14) + xlab("") + ylab("log2FC") + facet_wrap(~ Versus, ncol = 1) + theme(panel.grid.minor = element_blank())
  
  # Draw isoform from gff 
  C <- ggplot(gff_subset[transcript_id == gene]) +
    geom_segment(data = GFF[transcript_id == gene & type == "mRNA"], 
                 aes(x = ifelse(strand == "-", end + 20, start - 20), 
                     xend = ifelse(strand == "-", start - 80, end + 80), 
                     y = transcript_id, 
                     yend = transcript_id), 
                 color = "black", arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
    geom_tile(aes(x = start + (width / 2), y = transcript_id, width = width, height = 0.2, fill = type)) +
    labs(y = "Isoform", x = "Genomic coordinates", fill= "Region") + scale_fill_brewer(palette = "Paired") +
    theme_minimal() + 
    theme(panel.grid.major = element_line(size = rel(0.5)), panel.grid.minor = element_blank(), axis.text = element_text(size = 7)) 
  
  # Combined them 
  p_combined = A / C
  p_ranges_x <- c(ggplot_build(C)$layout$panel_scales_x[[1]]$range$range) # extract limits 
  # Match x axis to all graphs
  myplot <- p_combined + plot_layout(guides = "collect", heights = c(0.6, 0.3, 0.1)) + theme_bw(base_size = 14) & xlim(min(p_ranges_x), max(p_ranges_x)) 
  print(myplot)
  cli_alert_info("Saving into {path}{filename}")
  ggsave(plot = myplot, path = path, filename = filename, width = 12, height = 8, dpi = 300)
  cli_alert_success("Done!")
}

res_subset_annot[gene_id == "AT1G78080"]
makeGraph(dt = res_subset_annot, gene = "AT1G78080.1", path = "results/HyperTRIBE", filename = "AT1G78080.1.pdf")
res_subset_annot[gene_id == "AT1G22190"]
makeGraph(dt = res_subset_annot, gene = "AT1G22190.1", path = "results/HyperTRIBE", filename = "AT1G22190.1.pdf")

res[transcript_id == "AT1G78080.1"]
res[transcript_id == "AT1G22190.1"]

# Graphs ------------------------------------------------------------------

# Plot showing editing positions with labels 
ggplot(res_subset_annot[gene_id == "AT1G78080"], aes(x = start, y = gene_id, label = name)) + 
  geom_tile(data = gff_subset[gene_id == "AT1G78080"], 
            aes(x = start + (width / 2), y = gene_id, width = width, height = 0.2, fill = type, label = gene_id)) +
  geom_text_repel(point.padding = 0.2, 
                   nudge_x = .15,
                   nudge_y = .5,
                   segment.curvature = -1e-20) +
  labs(y = "Genes", x = "Genomic coordinates", fill= "Region") +
  scale_fill_brewer(palette = "Paired") +
  theme_minimal() + facet_wrap(~ Versus, ncol = 1, drop = FALSE, strip.position = "right") +
  theme(panel.grid.major = element_line(size = rel(0.5)), panel.grid.minor = element_blank(), axis.text = element_text(size = 7))


# Plot showing editing positions without label
ggplot(res_subset_annot[gene_id == "AT1G78080"], aes(x = start, y = gene_id)) + 
  geom_tile(data = GFF[gene_id == "AT1G78080" & type %in% c("five_prime_UTR", "CDS", "three_prime_UTR")], aes(x = start + (width / 2), y = gene_id, width = width, height = 0.2, fill = type)) +
  geom_tile(data = res_subset_annot[gene_id == "AT1G78080" & Type %in% c("five_prime_UTR & uORF")], aes(x = start + (width / 2), y = uORF_ID, width = width, height = 0.2), fill = "grey") +
  labs(y = "Genes", x = "Genomic coordinates", fill= "Region") +
  scale_fill_brewer(palette = "Paired") +
  theme(panel.grid.major = element_line(size = rel(0.5)), panel.grid.minor = element_blank(), axis.text = element_text(size = 7)) +
  theme_minimal() + facet_wrap(~ Versus, ncol = 1, drop = FALSE, strip.position = "right") +
  geom_tile(aes(x = Position, y = gene_id, width = 1, height = 0.4), fill = "red") #+ 
  #geom_point(aes(x = Position), color = "red", size = 0.2)

ggplot(res_subset_annot[, .N, by = c("Versus", "type")], aes(x = Versus, y = N, fill = type)) + geom_col() + theme_bw() +
  scale_fill_brewer(palette = "Set2") + 
  labs(y = "Number of editing sites", x = "Analyzes", fill= "Regions") + coord_flip()
