library(GenomicRanges)
library(rtracklayer)
library(devtools)
#install_github("sarah-ku/hyperTRIBER") # See https://github.com/sarah-ku/hyperTRIBER 
library(hyperTRIBER) 
library(data.table)
library(doParallel)
library(reshape2)
library(DEXSeq)

# Load data count ---------------------------------------------------------
data <- fread("all_samples.stranded.mpileup.basecounts.txt.gz", nThread = 4)
dim(data)

# Sample names definition. Must be exactly in same order as the "mpileup *.bam" command that we run 
# to generate the basecounts.txt file since there is no header in it.
# Order is ADAR, DNE1-D153N, DNE1, NYN-D153N.
# We have 5 replicates per condition
samp.names <- c(paste0("ADAR_REP", 1:5), paste0("DNE1-D153N_REP", 1:5), paste0("DNE1_REP", 1:5), paste0("NYN-D153N_REP", 1:5))
data_list <- extractCountData(data, samp.names, strand = TRUE)

# Create a GRanges object of the reference genome base for each of the considered positions.
locsGR <- GRanges(Rle(data$V1), IRanges(data$V2,width=1), ref=data$V3, names=paste(data$V1,data$V2,sep="_"))
rm(data)

# SETUP -------------------------------------------------------------------

# 1. ADAR vs DNE1 
# 2. ADAR vs DNE1-D153N
# 3. ADAR vs NYN-D153N
# 4. DNE1-D153N vs DNE1
# 5. DNE1-D153N vs NYN-D153N
# 6. NYN-D153N vs DNE1
# 7. NYN-D153N vs DNE1-D153N # here it needs all the edits otherwise it fails

# Change here the conditions and re-run the script
condition1 = "DNE1-D153N" # control 
condition2 = "DNE1"       # treat
nr_rep = 5                # Number of replicates
main_dir <- "analyses/HyperTRIBE/"         # Main output folder
my_edits <- rbind(c("A","G"), c("T", "C")) # Choose your edits

# Create the output directories: {results_dir}/{model_dir}
model_dir <- paste0(condition1, "_vs_", condition2)
results_dir <- file.path(main_dir, model_dir)
dir.create(results_dir, recursive = TRUE)

# WARNING: Sometimes the following error arises:  * > { : task 511 failed - "subscript contains NAs"
# It is because there were two possible edits here, and the sub-function to handle this threw the error 
# because none of those possibilities were in the allowed edits list. 
# It will improve in the next version of the package, but we can always just give 
# the getHits command the full list (my_edits below)and then subset for 
# the ref-target pairs you need after.

## For 7. NYN-D153N vs DNE1-D153N, it needs all edits

# my_edits <- rbind(c("A","G"), c("G","A"), c("T","C"), c("C","T"),
#                   c("A","T"), c("T","A"), c("G","T"), c("T","G"),
#                   c("G","C"), c("C","G"), c("A","C"), c("C","A"))
# -------------------------------------------------------------------------

# Now produce one design vector per experiment
design_vector <- setNames(c(rep("control", nr_rep), rep("treat", nr_rep)),
                          c(paste0(condition1, "_REP", 1:nr_rep), paste0(condition2, "_REP", 1:nr_rep)))

print(design_vector)

data_list <- data_list[names(design_vector)]

data_list_restricted <- restrict_data(data_list = data_list,
                                      design_vector = design_vector,
                                      min_samp_control = 4,
                                      min_samp_treat = 4,
                                      min_count = 2,
                                      edits_of_interest = my_edits)

generateCountFiles(data_list = data_list_restricted,
                   stranded = TRUE, 
                   names_vec = row.names(data_list_restricted[[1]]),
                   out_dir = results_dir, 
                   design_vector = design_vector)

dxd.res <- make_test(out_dir = results_dir, 
                     design_vector = design_vector, 
                     ncores = 12)

posGR <- getHits(res = dxd.res,
                 stranded = TRUE, 
                 fdr = 0.1,
                 fold = 1,
                 addMeta = TRUE,
                 ncore = 1,
                 include_ref = TRUE,
                 refGR = locsGR,
                 edits_of_interest = my_edits,
                 design_vector = design_vector,
                 data_list = data_list_restricted)

saveRDS(dxd.res, file.path(results_dir, paste0(model_dir, ".dxd_res.rds")))
saveRDS(posGR, file.path(results_dir, paste0(model_dir, ".posGR.rds")))

results <- as.data.table(posGR)
setkey(results, seqnames, start, end)

# Keep the A->G & T->C edits and save
fwrite(results[(ref=="A" & targ=="G") | (ref == "T" & targ == "C")], 
       file.path(results_dir, paste0(model_dir, ".DEXSeq.results.csv")))

# Add TAIR10 genes annotation to the result file
TAIR10_GENES_MODELS <- fread("ressources/TAIR10_annotations/TAIR10_representative_gene_models.txt", header = F, col.names = "AGI")
GFF <- as.data.table(import("ressources/TAIR10_annotations/TAIR10_GFF3_genes.clean.gff3"))
setkey(GFF, seqnames, start, end)

# Keep the most represented isoform
mRNA <- GFF[type == "mRNA" & transcript_id %in% TAIR10_GENES_MODELS$AGI]

# overlap annotation with position of the edit
results.annotated <- foverlaps(results, mRNA[, c("seqnames", "start", "end", "width", "gene_id", "transcript_id")])

# Save results with annotation
fwrite(results.annotated[, -c("i.start", "i.end", "i.width")], 
       file.path(results_dir, paste0(model_dir, ".DEXSeq.results.annotated.csv")))

