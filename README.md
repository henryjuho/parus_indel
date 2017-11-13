# Great tit (*Parus major*) INDEL analysis
Henry Juho Barton  
Department of Animal and Plant Sciences, The University of Sheffield  

# Introduction

This repository outlines the pipeline used to generate and analyse an INDEL dataset from 10 high coverage (mean coverage = 44X) great tit (*Parus major*) genomes (described here: [Corcoran et al. 2017](https://academic.oup.com/gbe/article-lookup/doi/10.1093/gbe/evx213)). The repository is subdivided by processing steps.

## Programs required

  * Python 2.7.2
  * GATK version 3.4-46-gbc02625 available from: <https://software.broadinstitute.org/gatk/download/archive>
  * VCFtools version 0.1.12b available from: <https://sourceforge.net/projects/vcftools/files/>
  * SAMtools version 1.2 available from: <https://sourceforge.net/projects/samtools/files/samtools/>
  * BCFtools version 1.3 
  * bedtools version 2.23.0
  * anavar version 1.2.2
  * q_sub.py and qsub_gen.py available from <https://github.com/henryjuho/python_qsub_wrapper>
  * pysam version 0.11.2.1 available from <https://github.com/pysam-developers/pysam>

\* Note \* that most scripts make use of the script 'qsub_gen.py' which is designed to submit jobs in the form of shell scripts to the 'Sun Grid Engine', if shell scripts only are required the '-OM' option in the 'qsub_gen.py' command line within the scripts can be changed from 'q' to 'w'. Alternatively some scripts make use of the python qsub wrapper module ```qsub.py``` described here: <https://github.com/henryjuho/python_qsub_wrapper>.

## Pre-prepared files required for analysis

  * Reference genome: **/fastdata/bop15hjb/GT_ref/Parus_major_1.04.rename.fa**
  * Reference genome index file: **/fastdata/bop15hjb/GT_ref/Parus_major_1.04.rename.fa.fai**
  * GFF annotation file: **/fastdata/bop15hjb/GT_ref/GCF_001522545.1_Parus_major1.0.3_genomic.gff.gz**
  * All sites VCF: **/fastdata/bop15hjb/GT_data/BGI/bgi_10birds.raw.snps.indels.all_sites.vcf**
  * Repeat masker bed file: **/fastdata/bop15hjb/GT_data/BGI_10_repeats/ParusMajorBuild1_v24032014_reps.bed**
  * BAM files for SAMtools calling: **/fastdata/bop15hjb/GT_data/BGI_10_BAM/\*.bam**

# Pipeline

## Generating the dataset

The variant calling and filtering pipeline for both SNPs and INDELs is described here: [variant_calling/](variant_calling).

## Multispecies alignment and INDEL polarisation

The generation of a multiple species alignment between zebra finch, great tit and fly catcher and its use in polarisating variants and identifying ancestral repeats is described here: [alignment_and_polarisation/](alignment_and_polarisation).

## Annotating the data 

Variant annotation using the NCBI ```GFF``` file is described here: [annotation/](annotation).

## Summary statistics and analyses

The calculation of summary statistics and other data summary analyses are documented here: [summary_analyses/](summary_analyses).

## Anavar analyses

Analysis of the INDEL data with the ```anavar``` package is described here: [anavar_analyses/](anavar_analyses).

## Length analyses

Analysis of the amount of sequence inserted and deleted in various genomic contexts is documented here: [length_analyses/](length_analyses).
