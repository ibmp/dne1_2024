# BiocManager::install("rtracklayer")
library(rtracklayer)
library(data.table)

# TAIR10 GENES MODELS file. Available on arabidopsis.org:
# https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FTAIR10_genome_release%2FTAIR10_gene_lists
TAIR10_GENES_MODELS <- fread("ressources/TAIR10_annotations/TAIR10_representative_gene_models.txt", header = F, col.names = "AGI")

WANTED <- c("miRNA", "tRNA", "ncRNA", "pseudogene", "snoRNA", "snRNA", "rRNA", "transposable_element")

# Load the annotations with the "import" function of rtracklayer package and convert to data.table 
gff <- as.data.table(import(con = "ressources/TAIR10_annotations/TAIR10_GFF3_genes_transposons.gff"))

# Select mRNA based on TAIR gene models file
mRNA_gene_models_gff <- gff[type == "mRNA" & Name %in% TAIR10_GENES_MODELS$AGI][, c("seqnames", "start", "end", "Name", "strand", "type", "Parent", "Alias")]
# Select other annotation from WANTED list
other_annotations_gff <- gff[type %in% WANTED][, c("seqnames", "start", "end", "Name", "strand", "type", "Parent", "Alias")]

# Regroup in one data.table
gff.dt <- rbindlist(list(mRNA_gene_models_gff, other_annotations_gff))
gff.dt[, score := "."]
setkey(gff.dt, seqnames, start, end)

# Create new column for ShortStack. ShortStack needs exactly this format --> Gene:start-end
gff.dt[, shortstack := paste0(seqnames, ":", start, "-", end)]

# Save the new annotation txt file for ShortStack
fwrite(gff.dt[, .(shortstack, Name)], "ressources/TAIR10_annotations/TAIR10_GFF3_genes_transposons.ShortStack.txt", sep = "\t", col.names = F)
# Save as BED FORMAT
fwrite(gff.dt[, .(seqnames, start, end, Name, score, strand)], "ressources/TAIR10_annotations/TAIR10_GFF3_genes_transposons.ShortStack.bed", sep = "\t", col.names = F)
