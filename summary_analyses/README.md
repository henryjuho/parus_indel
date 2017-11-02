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


## pi theta_w and Tajima's D

Fasta files of callable sites were created and summarised using the following codes:

| Case            | code  |
|:----------------|:-----:|
| N               | 0     |
| Filtered        | 1     |
| Pass polarised  | K     |
| Pass unpolarised| k     |
| AR polarised    | R     |
| AR unpolarised  | r     |

```
$ ~/parus_indel/summary_analyses/callable_sites_parallel.py -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/bgi_10birds.raw.snps.indels.all_sites.vcf.bgz -bed /fastdata/bop15hjb/GT_ref/ParusMajorBuild1_v24032014_reps.bed -ar_bed /fastdata/bop15hjb/GT_ref/Zebrafinch.Flycatcher.Greattit.ancLINEs.sorted.bed.gz -chr_bed /fastdata/bop15hjb/GT_ref/chromosome_list.bed -pol /fastdata/bop15hjb/bird_alignments/UCSC_pipeline/multiple_zhang_param/Zebrafinch.Flycatcher.Greattit.wga.bed.gz -out_pre /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Callable_sites/bgi_10birds.callable
$ ~/parus_indel/summary_analyses/callable_sites_summary_nogff.py -call_fa /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Callable_sites/bgi_10birds.callable.fa -chr_list /fastdata/bop15hjb/GT_ref/gt_autosomes.txt -opt_bed /fastdata/bop15hjb/GT_ref/gt_cds.bed.gz,CDS -opt_bed /fastdata/bop15hjb/GT_ref/gt_introns.bed.gz,intron -opt_bed /fastdata/bop15hjb/GT_ref/gt_intergenic.bed.gz,intergenic > /fastdata/bop15hjb/GT_ref/gt_callable_summary.csv
```

The statistics were then calculated with the following script:

```
$ ~/parus_indel/summary_analyses/summary_stats_gt.py -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Analysis_ready_data/final/bgi_10birds.filtered_indels.pol.anno.recomb.line.vcf.gz -call_csv /fastdata/bop15hjb/GT_ref/gt_callable_summary.csv -mode INDEL -no_sex -sub -out /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Summary_stats/bgi_10birds.indels.summary_stats.txt
$ ~/parus_indel/summary_analyses/summary_stats_gt.py -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Analysis_ready_data/final/bgi_10birds.filtered_indels.pol.anno.recomb.line.vcf.gz -call_csv /fastdata/bop15hjb/GT_ref/gt_callable_summary.csv -mode INDEL -no_sex -bootstrap 10000 -sub -out /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Summary_stats/bgi_10birds.indels.summary_stats.bs10000.txt
$ ~/parus_indel/summary_analyses/summary_stats_gt.py -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Analysis_ready_data/final/bgi_10birds.filtered_snps.pol.anno.degen.vcf.gz -call_csv /fastdata/bop15hjb/GT_ref/gt_callable_summary.csv -mode SNP -no_sex -sub -out /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Summary_stats/bgi_10birds.snps.summary_stats.txt
$ ~/parus_indel/summary_analyses/summary_stats_gt.py -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/SNP_data/gt_10birds_recalibrated_snps_only_99pass.maxlength50.biallelic.coveragefiltered.pass.repeatfilter.pass.polarised.annotated.degen.vcf.gz -call_csv /fastdata/bop15hjb/GT_ref/gt_callable_summary.csv -mode SNP -no_sex -bootstrap 10000 -sub -out /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Summary_stats/bgi_10birds.snps.summary_stats.bs10000.txt
```

And plotted:

```
$ cd ~/parus_indel/summary_analyses/
$ Rscript summary_stats.R
```

![stats_plot](gt_summary_stats.png)

## INDEL divergence and alpha

Simple INDEL divergence estimates were obtained from the whole genome alignment for coding and non-coding regions and plotted as follows.

```
$ ~/parus_indel/summary_analyses/indel_divergence.py -wga /fastdata/bop15hjb/bird_alignments/UCSC_pipeline/multiple_zhang_param/Zebrafinch.Flycatcher.Greattit.wga.bed.gz -bed /fastdata/bop15hjb/GT_ref/gt_cds.bed.gz -out /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Summary_stats/gt_cds_indel_divergence.txt
$ ~/parus_indel/summary_analyses/indel_divergence.py -wga /fastdata/bop15hjb/bird_alignments/UCSC_pipeline/multiple_zhang_param/Zebrafinch.Flycatcher.Greattit.wga.bed.gz -bed /fastdata/bop15hjb/GT_ref/gt_noncoding.bed.gz -out /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Summary_stats/gt_noncoding_indel_divergence.txt

$ Rscript collate_indel_divergence.R 
```

Alpha was calculated (see Equation 1 Eyre-walker 2006) for INDELs:

```
$ Rscript alpha.R 
```

This yields an alpha estimate of **0.1524897**
