# Great tit (*Parus major*) SNP analysis pipeline
Henry Juho Barton  
Department of Animal and Plant Sciences, The University of Sheffield  

# Introduction

This document outlines the pipeline used to generate and analyse a SNP dataset from 10 high coverage (mean coverage = 44X) great tit (*Parus major*) genomes. This dataset is used in a comparative approach in the INDEL analysis described in this repository. The document is subdivided by processing steps.

# Analysis
## SNP calling and VQSR

SNP calling and VQSR were performed with GATK - see Toni

## Post VQSR filters

The output vcf from the 99.0 tranche level was filtered to include only biallelic sites, non-repetative regions, and sites with a depth between half and twice the mean depth of 44X. The command lines are as follows:

```
python ~/qsub_gen.py -cmd "cd /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/SNP_data/" -cmd "java -Xmx6g -jar /usr/local/packages6/apps/binapps/GATK/3.4-46/GenomeAnalysisTK.jar -T SelectVariants -R ../../../GT_ref/Parus_major_1.04.rename.fa -V gt_10birds_recalibrated_snps_99.vcf.gz -selectType SNP -trimAlternates -env -o gt_10birds_recalibrated_snps_only_99.vcf.gz" -mem 10 -rmem 10 -o /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/SNP_data/snp_extract -OM q -evolgen
python ~/qsub_gen.py -cmd "cd /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/SNP_data/" -cmd "java -Xmx6g -jar \$GATKHOME/GenomeAnalysisTK.jar -T SelectVariants -R ../../../GT_ref/Parus_major_1.04.rename.fa -V  gt_10birds_recalibrated_snps_only_99.vcf.gz -o gt_10birds_recalibrated_snps_only_99pass.vcf.gz --excludeFiltered" -mem 10 -rmem 10 -o /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/SNP_data/vqsr_pass -evolgen -OM q
gunzip gt_10birds_recalibrated_snps_only_99pass.vcf.gz
python filter_length_biallelic.py -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/SNP_data/gt_10birds_recalibrated_snps_only_99pass.vcf -ref /fastdata/bop15hjb/GT_ref/Parus_major_1.04.rename.fa
python ~/depth_filter.py -vcf gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.vcf
python repeat_filter.py -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/SNP_data/gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.vcf -ref /fastdata/bop15hjb/GT_ref/Parus_major_1.04.rename.fa -bed /fastdata/bop15hjb/GT_ref/ParusMajorBuild1_v24032014_reps.bed 
```

The number of SNPs removed at each step are listed below:

|Filter             |SNPs retained |SNPs removed  |
|:------------------|:------------:|:------------:|
|None               |11784675      |NA            |
|Allele no          |11648889      |135786        |
|Depth              |11644621      |4268          |
|Repeat             |10772087      |872534        |

# Multispecies alignment and polarisation

## Alignment

Multispecies alignment was performed between zebra finch, flycatcher and great tit as described here: <https://github.com/henryjuho/bird_alignments/tree/master/Zebrafinch_Flycatcher_Greattit>.


## SNP polarisation

Firstly the sequence alignments across species for each SNP in the dataset were pulled out of the multiple alignment file with the following script:

```
./VARfromMAF.py  -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/SNP_data/gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.vcf -maf /fastdata/bop15hjb/bird_alignments/UCSC_pipeline/multiple_zhang_param/Zebrafinch.Flycatcher.Greattit.maf -target_spp Greattit -out /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/SNP_data/polarisation/ -no_jobs 500
```

Secondly the output of the above script was used to annotate the ancestral state for each variant in the INDEL vcf as follows:

```
./polarise_vcf_snps.py -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/SNP_data/gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.vcf -align_data /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/SNP_data/polarisation/all_variants.alignment_states.txt -target_spp Greattit -sub
```

|Category         | Number of SNPs |
|:----------------|:--------------:|
|Total no SNPs    | 10772087       |
|SNPs polarised   | 4430706        |
|Low spp coverage | 446296         |
|Ambiguous        | 2399524        |
|Not in alignment | 2515988        |
|Total unpolarised| 6341381        |