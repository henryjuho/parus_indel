#!/bin/bash

source ~/.bash_profile

#$-l h_rt=72:00:00
#$-l mem=6G
#$-l rmem=2G

#$-pe openmp 3
#$-P evolgen
#$-q evolgen.q

#$-N wga_subregions.sh
#$-o /fastdata/bop15hjb/hjb/GT_data/BGI_BWA_GATK/snp_divergence/wga_subregions.out
#$-e /fastdata/bop15hjb/hjb/GT_data/BGI_BWA_GATK/snp_divergence/wga_subregions.error

#$-V

cd /fastdata/bop15hjb/hjb/GT_data/BGI_BWA_GATK/snp_divergence/

bedtools intersect -a /fastdata/bop15hjb/hjb/GT_ref/Greattit.Zebrafinch.Flycatcher.wga.bed.gz -b /fastdata/bop15hjb/hjb/GT_ref/gt_introns.bed.gz | bgzip -c > gt_introns.wga.bed.gz &
bedtools intersect -a /fastdata/bop15hjb/hjb/GT_ref/Greattit.Zebrafinch.Flycatcher.wga.bed.gz -b /fastdata/bop15hjb/hjb/GT_ref/gt_intergenic.bed.gz | bgzip -c > gt_intergenic.wga.bed.gz &
bedtools intersect -a /fastdata/bop15hjb/hjb/GT_ref/Greattit.Zebrafinch.Flycatcher.wga.bed.gz -b /fastdata/bop15hjb/hjb/GT_ref/gt_cds.bed.gz | bgzip -c > gt_cds.wga.bed.gz &

wait


