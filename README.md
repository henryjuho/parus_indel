# Great tit (*Parus major*) INDEL analysis pipeline
Henry Juho Barton  
Department of Animal and Plant Sciences, The University of Sheffield  

# Introduction

This document outlines the pipeline used to generate and analyse an INDEL dataset from 10 high coverage (mean coverage = 44X) great tit (*Parus major*) genomes. The document is subdivided by processing steps.

## Programs required

  * Python 2.7.2
  * GATK version 3.4-46-gbc02625
  * VCFtools version 0.1.12b
  * SAMtools version 1.2
  * BCFtools version 1.3
  
## Python scripts used in this pipeline

\* Note \* that most scripts make use of the script 'qsub_gen.py' which is designed to submit jobs in the form of shell scripts to the 'Sun Grid Engine', if shell scripts only are required the '-OM' option in the 'qsub_gen.py' command line within the scripts can be changed from 'q' to 'w'.

  * qsub_gen.py
  * SAMtools_calling_v2.py
  * CatVariants.py
  * get_consensus_vcf.py
  * hardfilter_indels.py
  * depth_filter.py
  * repeat_filter.py
  * filter_length_biallelic.py
  * run_VQSR.py
  * trancheSTATS.py 
  

## Pre-prepared files required for analysis

  * Reference genome: **/fastdata/bop15hjb/GT_ref/Parus_major_1.04.rename.fa**
  * Reference genome index file: **/fastdata/bop15hjb/GT_ref/Parus_major_1.04.rename.fa.fai**
  * All sites VCF: **/fastdata/bop15hjb/GT_data/BGI/bgi_10birds.raw.snps.indels.all_sites.vcf**
  * Repeat masker bed file: **/fastdata/bop15hjb/GT_data/BGI_10_repeats/ParusMajorBuild1_v24032014_reps.bed**
  * BAM files for SAMtools calling: **/fastdata/bop15hjb/GT_data/BGI_10_BAM/\*.bam**
  
# Generating the dataset

## Creating a 'set of known variants' or 'truth set'

### 1) GATK - SAMtools consensus set

This step takes the intersection of INDELs discovered with GATK and SAMtools variant discovery pipelines. 
Firstly raw INDELs were extracted from an all sites VCF file produced by GATK using the following command:

```
java -Xmx6g -jar /usr/local/packages6/apps/binapps/GATK/3.4-46/GenomeAnalysisTK.jar -T SelectVariants -R /fastdata/bop15hjb/GT_ref/Parus_major_1.04.rename.fa -V /fastdata/bop15hjb/GT_data/BGI/bgi_10birds.raw.snps.indels.all_sites.vcf -selectType INDEL -trimAlternates -env -o /fastdata/bop15hjb/GT_data/BGI/bgi_10birds.raw.snps.indels.all_sites.rawindels.vcf
```

Secondly SAMtools was used to call INDELs and SNPs from the BAM files using the python script 'SAMtools_calling_v2.py' which implements SAMtools and BCFtools calling pipeline per chromosome, before using GATK CatVariants (from within the script 'CatVariants.py') to produce a genome wide file. This step uses the following command:

```
python SAMtools_calling_v2.py -bam_list /fastdata/bop15hjb/GT_data/BGI_10_BAM/bgi10_bam.list -ref /fastdata/bop15hjb/GT_ref/Parus_major_1.04.rename.fa -out /fastdata/bop15hjb/GT_data/BGI/Consensus/SAMtools_run2/BGI_10.dedup.real.recal
```

The consensus set of the SAMtools and GATK INDELs was then obtained using GATK's SelectVariants with '--concordance', implemented in the script 'get_consensus_vcf.py', used as follows:

```
python get_consensus_vcf.py -vcf_I /fastdata/bop15hjb/GT_data/BGI/bgi_10birds.raw.snps.indels.all_sites.rawindels.vcf -vcf_II /fastdata/bop15hjb/GT_data/BGI/Consensus/SAMtools_run2/BGI_10.dedup.real.recal.allsites.vcf -ref /fastdata/bop15hjb/GT_ref/Parus_major_1.04.rename.fa -out /fastdata/bop15hjb/GT_data/BGI/Consensus/
```

### 2) Hardfiltering 

The consensus set was filtered using the GATK best practice hard filters for INDELs of "QD<2.0", "FS>200.0" and "ReadPosRankSum<-20.0" (see <https://www.broadinstitute.org/gatk/guide/tagged?tag=filtering>). This was implemented in the python script 'hardfilter_indels.py' as follows:

```
python hardfilter_indels.py -vcf /fastdata/bop15hjb/GT_data/BGI/Consensus/bgi_10birds.raw.snps.indels.all_sites.rawindels.consensus.rawindels.vcf
```

### 3) Coverage filtering

Next, INDELs with coverage more than twice, or less than half, the mean coverage were excluded with the in-house script 'depth_filter.py'. Usage follows:

```
python depth_filter.py -vcf /fastdata/bop15hjb/GT_data/BGI/Consensus/bgi_10birds.raw.snps.indels.all_sites.rawindels.consensus.rawindels.hardfiltered.pass.vcf
```

### 4) Repeat filtering

Repeat regions identified by a previous repeat masker analysis were excluded using GATK's VariantFiltration implemented in the script 'repeat_filter.py':

```
python repeat_filter.py -vcf /fastdata/bop15hjb/GT_data/BGI/Consensus/bgi_10birds.raw.snps.indels.all_sites.rawindels.consensus.rawindels.hardfiltered.pass.coveragefiltered.pass.vcf -bed /fastdata/bop15hjb/GT_data/BGI_10_repeats/ParusMajorBuild1_v24032014_reps.bed -ref /fastdata/bop15hjb/GT_ref/Parus_major_1.04.rename.fa
```

### 5) Filter by length and number of alleles

INDELs that were longer than 50 bp or were not biallelic were excluded using GATK SelectVariants implemented in the script 'filter_length_biallelic.py':

```
python filter_length_biallelic.py -vcf /fastdata/bop15hjb/GT_data/BGI/Consensus/bgi_10birds.raw.snps.indels.all_sites.rawindels.consensus.rawindels.hardfiltered.pass.coveragefiltered.pass.repeatfilter.pass.vcf -ref /fastdata/bop15hjb/GT_ref/Parus_major_1.04.rename.fa
```

### Truth set step summary

The number of INDELs remove at each step of filtering during creation of the 'truth set' are listed in the table below: 

|Step                              |INDELs passed                     |INDELs failed      |
|----------------------------------|:--------------------------------:|:-----------------:|
|Raw INDELs                        |1869484                           |NA                 |
|Consensus calling                 |1408029                           |461455             |
|Hardfiltering                     |1406784                           |1245               |
|Coverage filtering                |1320289                           |86495              |
|Repeat filtering                  |1156233                           |164056             |
|Length and allele number filtering|1047463                           |108770             |

## Variant quality score recalibration (VQSR)

### 1) VQSR

VQSR was performed on raw INDELs called with the GATK HaplotypeCaller, using the INDEL set descirbed above as the 'truth set' to train the VQSR model. VQSR was implement with the python script 'run_VQSR.py' as follows:

```
python run_VQSR.py -vcf /fastdata/bop15hjb/GT_data/BGI/bgi_10birds.raw.snps.indels.all_sites.rawindels.vcf -ref /fastdata/bop15hjb/GT_ref/Parus_major_1.04.rename.fa -truth_set /fastdata/bop15hjb/GT_data/BGI/Consensus/bgi_10birds.raw.snps.indels.all_sites.rawindels.consensus.rawindels.hardfiltered.pass.coveragefiltered.pass.repeatfilter.pass.maxlength50.biallelic.vcf -out /fastdata/bop15hjb/GT_data/BGI/VQSR_consensus/
```

This produced output VCFs for six tranche levels, as show in the table below, along with the number of variants retained at each level.

|Tranche level | INDELs retained | INDELs excluded |
|:------------:|:---------------:|:---------------:|
|90.0          |1420423          |449061           |
|98.0          |1584287          |285197           |
|99.0          |1612948          |256536           |
|99.5          |1632060          |237424           |
|99.9          |1817437          |52047            |
|100.0         |1869484          |0                |

The custom script 'trancheSTATS.py' was used to evaluate which tranche level was most suitable by obtaining data on the number of novel INDELs discovered at each tranche level. The script uses a list of vcfs generated using ls. The command used is as follows:

```
python trancheSTATS.py -vcf /fastdata/bop15hjb/GT_data/BGI/VQSR_consensus/post_vqsr_vcfs.list -out /fastdata/bop15hjb/GT_data/BGI/VQSR_consensus/bgi_10_tranchelevel_stats.txt
```

The numbers of novel variants in each tranche are shown in the table below.

|Tranche  |All sites	|Positive training sites|	Negative training sites|	Novel sites|
|:-------:|:---------:|:---------------------:|:----------------------:|:-----------:|
|100.0	  |1869484	  |1047463	              |321122	                 |500899       |
|99.9     |1817437	  |1046415                |302038 	               |468984       |
|99.5     |1632060	  |1042225                |133408                  |456427       |
|99.0     |1612948	  |1036992                |116234                  |459722       |
|98.0	    |1584287	  |1026513                |91570	                 |466204       |
|90.0     |1420423	  |942720	                |14487	                 |463216       |

### 2) Post VQSR filters

As conducted for the truth set generating step, the output vcf from the 99.0 tranche level was filtered to include only biallelic sites, INDELs 50 bp long or less, non-repetative regions, and sites with a depth between half and twice the mean depth of 44X. The command lines are as follows:

```
python filter_length_biallelic.py -vcf /fastdata/bop15hjb/GT_data/BGI/VQSR_consensus/bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.vcf -ref /fastdata/bop15hjb/GT_ref/Parus_major_1.04.rename.fa
python depth_filter.py -vcf /fastdata/bop15hjb/GT_data/BGI/VQSR_consensus/bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.vcf
python repeat_filter.py -vcf /fastdata/bop15hjb/GT_data/BGI/VQSR_consensus/bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.vcf -ref /fastdata/bop15hjb/GT_ref/Parus_major_1.04.rename.fa -bed /fastdata/bop15hjb/GT_data/BGI_10_repeats/ParusMajorBuild1_v24032014_reps.bed
```
The number of INDELs removed at each step are listed below:

|Filter             |INDELs retained |INDELs removed|
|:------------------|:--------------:|:------------:|
|None               |1612948         |NA            |
|Allele no. & length|1431198         |181750        |
|Depth              |1373399         |57799         |
|Repeat             |1240366         |133033        |

# Multispecies alignment and INDEL polarisation

TODO

# Performing the INDEL analysis

## Repeat region analysis

The occurance of INDELs in repeat regions was investigated using GATK's 'VariantAnnotator' to annotate homoploymer runs and tandem repeats in the INDEL dataset. This used a python wrapper as follows:

```
python annotate_hr_tr.py -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Analysis_ready_data/bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.vcf -ref /fastdata/bop15hjb/GT_ref/Parus_major_1.04.rename.fa -out /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Repeat_analysis/
```
Summary data was the generated using the script ```indel_repeat_stats.py``` as follows:

```
python indel_repeat_stats.py -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Repeat_analysis/bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.hr.tr.vcf
```

This yielded the following results:

|Class	          |Number of INDELs|
|:----------------|:--------------:|
|Hompolymer runs	|241528          |
|Tandem repeats	  |668234          |
|Overlap	        |241528          |
|Total repetative	|668234          |
|Total no. INDELs	|1240366         | 
|% INDEL in repeat|53              |

## Binning data by recombination rate

### 1) Predicting recombination rate for each INDEL

Linkage map data was used to estimate recombination rate at each INDEL. First, 3rd order polynomials were fitted to plots of physical position versus map length. Second, the derivative of each chromosome's polynomial was used to estimate recombination rate for each INDEL start position. This predicition was implemented in the following python script:

```
python predict_recomb.py -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Analysis_ready_data/bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.vcf -out /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Recomb_data/indel_recomb_data.txt -poly /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Recomb_data/3rd_polynom.txt
```

The file specified by ```-poly``` is a list of variables for each chromosome's polynomial.

### 2) Binning the data

The output from the previous step was used to bin the data by recombination yielding x new vcf files, where x = the required number of recombination bins. Bins had equal numbers of variants. The script used for binning is as follows:

```
python recomBINdels.py -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Analysis_ready_data/bgi_10birds.raw.snps.indels.all_sites.rawindels.recalibrated.filtered_t99.0.pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.vcf -recomb /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Recomb_data/indel_recomb_data.txt -out /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/SFS/5_recomb_bins/
```
