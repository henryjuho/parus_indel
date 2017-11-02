#!/usr/bin/env bash

source ~/.bash_profile

#$-l h_rt=8:00:00
#$-l mem=6G
#$-l rmem=2G


#$-N prepare_windows.sh
#$-o /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/anavar_analysis/window_analysis/prepare_windows.out
#$-e /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/anavar_analysis/window_analysis/prepare_windows.error

#$-V

cd /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/anavar_analysis/window_analysis/

# Create genome file for each species
# Format is as follows
# <chromName><TAB><chromSize>

faCount < /fastdata/bop15hjb/GT_ref/Parus_major_1.04.rename.fa | grep -v ^# | cut -f1,2 | grep ^chr > gt_genome.txt

# Create the 1Mb windows for the zebra finch genome and create a numerical id for each window.
# Format is as follows;
# <chromName><TAB><StartPos><TAB><EndPos><TAB><WindowID>

bedtools makewindows -g gt_genome.txt -w 2000000 | ~/biased_gene_conversion/create_windows/add_window_id.py | bgzip -c > gt_windows.2Mb.bed.gz
tabix -p bed gt_windows.2Mb.bed.gz

zgrep ^# /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Analysis_ready_data/final/bgi_10birds.filtered_indels.pol.anno.recomb.line.vcf.gz > gt_indels_non_coding.vcf
zgrep -v ^# /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Analysis_ready_data/final/bgi_10birds.filtered_indels.pol.anno.recomb.line.vcf.gz | grep ANNO=int >> gt_indels_non_coding.vcf

bgzip gt_indels_non_coding.vcf
tabix -p vcf gt_indels_non_coding.vcf.gz

python ~/biased_gene_conversion/create_windows/calc_window_rec_rates.py --infile ~/biased_gene_conversion/create_windows/rec_rate_files/gt_genetic_physical_pos.txt \
	--genome gt_genome.txt --bed gt_windows.2Mb.bed.gz \
	--vcf gt_indels_non_coding.vcf.gz \
> gt_2Mb_window_rec_rates.txt