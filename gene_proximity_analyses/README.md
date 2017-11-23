# Analysis of neutral INDEL data different distances from genes 

In order to investigate the possibility of linked selection influencing the ratio of insertions to deletions (rDI) in
regions of low recombination, INDELs where binned into a number of 50Kb windows according to gene distance.

```
$ cd ~/parus_indel/gene_proximity_analyses
$ mkdir distance_bin_beds
$ zcat /fastdata/bop15hjb/GT_ref/gt_intergenic.bed.gz | ./create_gene_proximity_bins.py -bin_size 50000 -out_prefix distance_bin_beds/gt_intergenic_gene_proximity_50kbwindows
```

