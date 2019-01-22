# Investigating correlations between INDEL length and Tajima's D and pi

We looked into the relationship between INDEL length and pi and Tajima's D.

```bash
bedtools intersect -header -a /fastdata/bop15hjb/h_j_b/GT_data/BGI_BWA_GATK/Analysis_ready_data/final/bgi_10birds.filtered_indels.pol.anno.recomb.line.vcf.gz -b /fastdata/bop15hjb/h_j_b/GT_ref/gt_cds.bed.gz | ./length_summary_stats.py -region CDS > gt_cds_indels_length_sum.csv
bedtools intersect -header -a /fastdata/bop15hjb/h_j_b/GT_data/BGI_BWA_GATK/Analysis_ready_data/final/bgi_10birds.filtered_indels.pol.anno.recomb.line.vcf.gz -b /fastdata/bop15hjb/h_j_b/GT_ref/gt_noncoding.bed.gz | ./length_summary_stats.py -region noncoding > gt_nc_indels_length_sum.csv

Rscript length_sum_stats.R 
```

<img src='indel_length_summary_stats.png' width=700 height=700>