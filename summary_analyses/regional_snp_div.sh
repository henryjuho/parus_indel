#!/bin/bash

source ~/.bash_profile

#$-l h_rt=72:00:00
#$-l rmem=6G

#$-pe openmp 5
#$-P evolgen
#$-q evolgen.q

#$-N snp_div
#$-o /fastdata/bop15hjb/hjb/GT_data/BGI_BWA_GATK/snp_divergence/snp_div.out
#$-e /fastdata/bop15hjb/hjb/GT_data/BGI_BWA_GATK/snp_divergence/snp_div.error

#$-V

cd /fastdata/bop15hjb/hjb/GT_data/BGI_BWA_GATK/snp_divergence/

Rscript ~/parus_indel/summary_analyses/k80_div_est.R gt_introns.wga.fa introns &
Rscript ~/parus_indel/summary_analyses/k80_div_est.R gt_intergenic.wga.fa intergenic &
Rscript ~/parus_indel/summary_analyses/k80_div_est.R gt_cds.wga.fa CDS &
Rscript ~/parus_indel/summary_analyses/k80_div_est.R gt_ar.wga.fa AR &
Rscript ~/parus_indel/summary_analyses/k80_div_est.R gt_all.wga.fa ALL &

wait


