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

# Annotating the data

Firstly SNPs were annotated in the vcf file as belonging to either 'CDS_non_frameshift', 'intron' or 'intergenic' as follows:

```
./annotate_all_vcf_chr.py -gff /data/bop15hjb/annotating_gt_snps/GCF_001522545.1_Parus_major1.0.3_genomic.rename.gff.gz -vcf /data/bop15hjb/annotating_gt_snps/gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.polarised.vcf -evolgen
```

A breakdown of this annotation is shown below:

|Category             | Number SNPs  |
|:--------------------|:------------:|
|All                  | 10643865     |
|CDS                  | 155463       |
|Intron               | 5489556      |
|Intergenic           | 4718020      |
|Unannotated          | 280826       |

Secondly the degeneracy of SNPs in coding seqences was annotated as follows:

```
./annotate_degeneracy.py -gff /fastdata/bop15hjb/GT_ref/GCF_001522545.1_Parus_major1.0.3_genomic.rename.gff.gz -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/SNP_data/gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.polarised.annotated.vcf -ref /fastdata/bop15hjb/GT_ref/Parus_major_1.04.rename.fa -db_dir /data/bop15hjb/databases/greattit/ -out /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/SNP_data/ degeneracy_annotation/ -evolgen
```

Annotation summary below:

|Category          | Number SNPs     |
|:-----------------|:---------------:|
|All               | 10643865        |
|CDS               | 155463          |
|0fold             | 46792           |
|2fold             | 46786           |
|3fold             | 4240            |
|4fold             | 57169           |

# Site frequency spectrum analysis

Firstly the unfolded SFS for W<->W and S<->S snps was obtained from the snp data as follows:

```  
./snpSFS.py -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/SNP_data/gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.polarised.annotated.vcf -folded N -bin intergenic_ww_ss -sfs_out /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/SFS/snp_ww_ss_spectra -sub -evolgen
```

Secondly the unfolded SFS for zerofold coding snps was generated as follows:

```
./snpSFS.py -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/SNP_data/gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.polarised.annotated.degen.vcf -folded N -bin zerofold -sfs_out /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/SFS/zerofold_snps -evolgen -sub
```

Finally both these SFS were also generated folded:

```
./snpSFS.py -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/SNP_data/gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.polarised.annotated.degen.vcf -folded Y -bin intergenic_ww_ss -bin zerofold -sfs_out /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/SFS/folded_snps -sub -evolgen
```

# pi theta_w and Tajima's D

In order to calculate  pi, theta_w and Tajima's D, first a fasta of callable sites with the bases coded as 0s (Ns), 1s (not callable) or 2s (callable) was generated and indexed as follows:

```
./callable_sites_from_vcf.py -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/bgi_10birds.raw.snps.indels.all_sites.vcf.bgz -chr ALL -bed /fastdata/bop15hjb/GT_ref/ParusMajorBuild1_v24032014_reps.bed -out /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Callable_sites/bgi_10birds_callable -evolgen -sub
samtools faidx bgi_10birds_callable.ALL.fa
```

The statistics were then calculated with the following script:

```
./summarise_vcf.py -vcf ../debug_data/gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.polarised.annotated.degen.vcf.gz -call_fa ../debug_data/bgi_10birds_callable.ALL.fa -mode SNP > /Users/henryjuho/genomics/indel_pi_theta_tajd/gt_snp_sum_stats.txt
```
