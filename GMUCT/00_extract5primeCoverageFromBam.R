library(BSgenome.Athaliana.TAIR.TAIR9)
library(future.apply)
library(data.table)
library(openxlsx)
library(ggplot2)
library(cli)

# Settings ----------------------------------------------------------------
if (Sys.getenv("RSTUDIO")=="1") options(cli.palette = "dichro") # rstudio cli color palette
theme_set(theme_bw()) # set ggplot default theme
# -------------------------------------------------------------------------

# Global variables --------------------------------------------------------
BSGENOME = getSeq(BSgenome.Athaliana.TAIR.TAIR9) # Genome as BSgenome object
BEDTOOLS = Sys.which("bedtools") # Search bedtools in the PATH. 

# If not present in PATH, you should put the path of the binary
# BEDTOOLS = "/usr/bin/bedtools"
# -------------------------------------------------------------------------

try(if(BEDTOOLS == "") stop("bedtools was not found in the PATH. You need to specify the binary path in the BEDTOOLS variable or add bedtools to your PATH"))
try(system(paste("bedtools", "--version"), intern = TRUE, ignore.stderr = TRUE))

# Functions ---------------------------------------------------------------

#' Parse bedtools genomecov function to extract the coverage of each mapped reads 5' positions with strand information
#'
#' @param bamfile The path of a bam file. (Must be indexed .bam.bai)
#' @param bedtools The path of bedtools binary. Default is Sys.which("bedtools"). If not in your PATH you should specify the fullpath here
#' @param genomecov_header The header of the bedtools genomecov. Default is c("Chr", "Position", "Coverage") from the "-dz -5" option
#' @param save.stats Boolean to save stats in a {bamfile}.stats.txt file. Default is FALSE. 
#' 
#' @return A data.table with Sample, Chr, Position, Coverage, Strand, RPM & Nucleotide information per 5' positions
#' @export
#'
#' @examples
#' get_5prime_coverage_stranded(bamfile = "xrn4.bam", bedtools = "/usr/bin/bedtools", save.stats = T)

get_5prime_coverage_stranded <- function(bamfile, bedtools = BEDTOOLS, genomecov_header = c("Chr", "Position", "Coverage"), save.stats = F) {
  
  cli_alert_success("Analysing {.path {bamfile}}") 
  # We need per-base report so we use -dz option in bedtools instead of the
  # -bg option (used in Fivepseq) for bedgraph output that reports position within regions >=1bp. 
  # The '-5' option gives the coverage of the 5' position of each mapped read. 
  # The '-strand' gives the strand information
  
  # Extract minus strand 5' coverage
  rev_coverage <- fread(cmd = paste(bedtools, "genomecov -dz -5 -strand - -ibam", bamfile), col.names = genomecov_header)
  rev_coverage[, Strand := "-"]
  # Extract plus strand 5' coverage
  fwd_coverage <- fread(cmd = paste(bedtools, "genomecov -dz -5 -strand + -ibam", bamfile), col.names = genomecov_header)
  fwd_coverage[, Strand := "+"]
  # Concat rev & fwd in same table
  coverage_data <- rbindlist(list(rev_coverage, fwd_coverage))  
  setkey(coverage_data, Chr, Position) # Index by Chr and Position
  
  # Number of reads in the dataset
  rev_reads <- sum(rev_coverage$Coverage)
  fwd_reads <- sum(fwd_coverage$Coverage)
  total_reads <- sum(coverage_data$Coverage)
  cli_alert_info("Number of 5' positions: {total_reads} (+: {fwd_reads} / -: {rev_reads})")
  
  # Add a column for RPM normalization
  coverage_data[, RPM := round((Coverage / total_reads) * 1e6, 3)]
  
  # Since the "bedtools genomecov" does not output the 5' nucleotide
  # we add it ourself by extracting the nucleotide from the genome (using here BSgenome lib) 
  # We need to add +1 to each position because bedtools genomecov positions are 0-based. 
  # We use the BSgenome::getSeq function to extract the nucleotide with the positions, as a GRanges object. 
  
  coverage_data[, Position := Position + 1] # TODO: Changes position name to "0_based_pos" or "1_based_pos" to avoid confusion
  coverage_data[, Nucleotide := as.character(
    getSeq(BSGENOME, 
           makeGRangesFromDataFrame(coverage_data, 
                                    start.field = "Position", 
                                    end.field = "Position", 
                                    strand.field = "Strand"))
  )]
  
  # Add the position ID to be like this: Chr1_1008_+:T (for DEXSeq only)
  coverage_data[, ID := paste0(Chr, "_", Position, "_", Strand, ":", Nucleotide)]
  
  # Save the position ID with the coverage in a count txt file. It adds ".count.txt" to the given bamfile path)
  fwrite(coverage_data, paste0(bamfile, ".coverage.txt"), col.names = T, sep = "\t")
  
  # Just print some stats per chromosome
  stats = coverage_data[, list(Nr_of_positions = .N, 
                               Min = min(Coverage), 
                               Max = max(Coverage), 
                               Mean = round(mean(Coverage), 2), 
                               Median = round(median(Coverage)),
                               Std = round(sd(Coverage), 2)), 
                        by = Chr]
  print(stats)
  if(save.stats) fwrite(stats, paste0(bamfile, ".coverage.stats.txt"), col.names = T, sep = "\t")
  return(coverage_data)
}


# Analysis ----------------------------------------------------------------

# GMUCT bamfiles to analyse
mybamfiles <- list.files(path = "analyses/GMUCT/GMUCT_BAMFILES", 
                         pattern = ".sort.bam$", 
                         full.names = T) 

# Speed up things, one core per bamfile
plan(multisession, workers = length(mybamfiles)) 

# Run everything in parallel and save it in a list 
mycount <- future_lapply(mybamfiles, get_5prime_coverage_stranded, save.stats = T)

# Keep the names clean
sample_names <- gsub(pattern = "-vs_TAIR10.sort.bam", "", basename(mybamfiles))
# Add names for rbindlist 
names(mycount) <- sample_names

# Rbind all the results 
mycount.dt <- rbindlist(mycount, idcol = "Sample")

# Save it
fwrite(mycount.dt, "GMUCT_all_samples.counts.txt.gz", sep = "\t")

# Load it back if already computed
mycount.dt <- fread("GMUCT_all_samples.counts.txt.gz")

# Brief description of data
ggplot(mycount.dt[RPM > 1], aes(Sample, log2(Coverage), fill = Chr)) + geom_boxplot()
ggplot(mycount.dt[RPM > 1], aes(Sample, log2(RPM), fill = Chr)) + geom_boxplot()
ggplot(mycount.dt[RPM > 1], aes(log2(RPM), fill = Sample)) + geom_histogram(binwidth = 0.01) + facet_wrap(~ Sample, ncol = 3)
summary(mycount.dt[RPM > 1]$RPM)
