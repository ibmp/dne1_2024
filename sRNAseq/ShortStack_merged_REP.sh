#!/bin/bash
genome="../ressources/TAIR10_genome/TAIR10_Chr.all.fasta"
locifile="../ressources/TAIR10_annotations/TAIR10_GFF3_genes_transposons.ShortStack.txt"

bamfolder="../analyses/sRNAseq/00_sRNAseq_mapping_TAIR10"
outdir="../analyses/sRNAseq/01_shortstack"

mkdir -p $outdir

## SINGLE BAM ANALYSIS
for bam in "$bamfolder"/*.R1_sorted.bam
do
    dirname=$(basename ${bam/.R1_sorted.bam/_shortstack})
    echo $dirname
    ShortStack --threads 6 --nohp --dicermin 15 --dicermax 30 --locifile $locifile --bamfile $bam --genomefile $genome --outdir ${outdir}/${dirname}
done

## MERGE BAM ANALYSIS
## The replicates were merged together in a single bam file
bamdir_merged="../../analyses/sRNAseq/00_sRNAseq_mapping_TAIR10/Merged_rep/Merged_bam"
outdir_merged="../../analyses/sRNAseq/01_shortstack_merged"

mkdir -p $outdir_merged

for bam in "$bamdir_merged"/merged*.bam
do
    dirname=$(basename ${bam/.bam/_shortstack})
    echo $dirname
    ShortStack --threads 8 --nohp --dicermin 15 --dicermax 30 --locifile $locifile --bamfile $bam --genomefile $genome --outdir ${outdir_merged}/${dirname}
done







