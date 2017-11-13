# Performing the analyses of sequence length deleted and inserted

## Regional length analysis

Firstly the total amount of sequence deleted and inserted per individual was calculated from the vcf file, both genome 
wide and within coding, intronic and intergenic regions.

```
$ cd parus_indel/length_analyses/
$ ~/parus_indel/length_analyses/region_length_summary.py -vcf /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Analysis_ready_data/final/bgi_10birds.filtered_indels.pol.anno.recomb.line.vcf.gz -call_fa /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Callable_sites/bgi_10birds.callable.fa -chr_bed /fastdata/bop15hjb/GT_ref/gt_chromosomes.bed.gz -opt_bed /fastdata/bop15hjb/GT_ref/gt_cds.bed.gz,CDS -opt_bed /fastdata/bop15hjb/GT_ref/gt_introns.bed.gz,intron -opt_bed /fastdata/bop15hjb/GT_ref/gt_intergenic.bed.gz,intergenic > bgi_10birds_bp_loss_gain.csv
```

Length data can be found [here](bgi_10birds_bp_loss_gain.csv).

From this data per haplotype, per base insertion and deletion rates were calculated, we define these rates as the 
total amount of sequence either deleted or inserted in our vcf file divided by the number of haplotypes (2n) and divided
by the number of callable sites. Additionally the length adjusted rDI (lrDI), the ratio of base pairs deleted to 
inserted was calculated. This was undertaken with the following script:

```
$ Rscript regional_length_analysis.R
```

![regional_len](regional_length_plots.png)
