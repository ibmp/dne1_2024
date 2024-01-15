# Snakemake HyperTRIBEseq

Snakemake pipeline for HyperTRIBEseq data analysis

## Description

The workflow automates all the steps described here: 

- https://github.com/sarah-ku/hyperTRIBER

## Installation

Just clone the repository. All the tools are installed in a conda environment when running the snakemake with the ```--use-conda``` option 
See here */envs/hypertribeseq.yaml* for all the dependencies.

## Getting started

You need to create a new config file (ex: config/my-config-tribeseq-June2022.yaml) based on the default config file present in the config folder ("config/config.yaml")

## Settings in the config file

```yaml 
sample_table: Sample_Names_TRIBE.xlsx

outdir: "TRIBEseq_analysis/"
rawdata_dir: "/data/TRIBEseq/0_RawData/NGS287-TRIBEseq/rawdata/220701"

runs_dir: "1_Runs"
mapping_dir: "2_Mapping"
results_dir: "3_Results"
genome_dir: "ReferenceGenome"

reference:
    genome: "/Genomes/Arabidopsis_thaliana/Araport11/TAIR10_Chr.all.fasta" # Reference genome in fasta

# R1 extension format with file extension. 
# generally "_R1.fastq.gz" or "_R1.fastq". It will be used to find R2 files.
# We will replace 1 by 2 in the extension. 
extension: ".R1.fastq.gz"

# hisat2 options
hisat2_options: "-t -k 50 --max-intronlen 2000 --rna-strandness RF --no-unal --rg LB:2x100bp_RF_stranded --rg PL:Illumina_HiSeq_4000 --rg PU:Arabidopsis_thaliana"
```

## Usage

Just run the snakemake command using the new config created. (Snakemake version >=7 is required)
It needs the **--configfile**, the **--use-conda** as well as the **--profile** option if running on slurm.

```bash
conda activate snakemake7
snakemake --configfile config/config_projectXX_01022022.yaml --use-conda --profile slurm
```

## Problems

The creation of the mpileup file can be probematic when running with slurm on nodes since the output file is very large (>600Go)
It is best to redirect the ```tmp_dir``` to a folder the main node or to a ssd storage space.

Here is a command to speed up the computation. It is best to write the mpileup on a ssd.

```sh
parallel --tmpdir /ssd_workspace/tmp -k --colsep '\t' samtools mpileup --max-depth 20000 -Q 40 --skip-indels -f ReferenceGenome/TAIR10_Chr.all.fasta -r {1} 2_Mapping/STRANDED/*.bam :::: ReferenceGenome/TAIR10_Chr.all.fasta.fai > /ssd_workspace/ANALYSIS/all_samples.mpileup
```

Best way to launch the perl script in parallel on the large mpileup file. We make use of the ```parallel``` command.

```sh
srun -p gpu -c 40 parallel -j 40 --block-size 500M --keep --pipepart TRIBEseq/.snakemake/conda/cdfe291fb5757b4cfc46abd106fcef5c/bin/perl ~/Snakemake_hyperTRIBE/workflow/scripts/RNAeditR_mpileup2bases.pl :::: all_samples.mpileup > all_samples.mpileup.stranded.basecounts.txt
```

## Output

The workflow maps the reads with hisat2 on the reference genome. 
Then it creates a single mpileup file from all the bams.   
**Note:** This step generates a very big mpileup file, ~500Go for 20 bams.  

Finally it runs a perl script (see in the *scripts/* folder) using GNU parallel on the big mpileup file. 
It creates a final count file called **'all.pileup.count.sorted.txt'**, which is the input required for the hyperTRIBE R package.

The rest of the analysis is made in Rstudio.
