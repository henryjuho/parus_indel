#!/bin/bash

source ~/.bash_profile

#$-l h_rt=72:00:00
#$-l mem=6G
#$-l rmem=2G

#$-pe openmp 5
#$-P evolgen
#$-q evolgen.q

#$-N wga_subregions2fa.sh
#$-o /fastdata/bop15hjb/hjb/GT_data/BGI_BWA_GATK/snp_divergence/wga_subregions2fa.out
#$-e /fastdata/bop15hjb/hjb/GT_data/BGI_BWA_GATK/snp_divergence/wga_subregions2fa.error

#$-V

cd /fastdata/bop15hjb/hjb/GT_data/BGI_BWA_GATK/snp_divergence/

zcat gt_introns.wga.bed.gz | ~/parus_indel/summary_analyses/wga2fa.py -out_stem gt_introns.wga &
zcat gt_intergenic.wga.bed.gz | ~/parus_indel/summary_analyses/wga2fa.py -out_stem gt_intergenic.wga &
zcat gt_cds.wga.bed.gz | ~/parus_indel/summary_analyses/wga2fa.py -out_stem gt_cds.wga &
zcat /fastdata/bop15hjb/hjb/GT_data/BGI_BWA_GATK/divergence_data/repeat_coordinates/LINE_intersect/Greattit.Zebrafinch.Flycatcher.LINEs.sorted.wga.bed.gz | ~/parus_indel/summary_analyses/wga2fa.py -out_stem gt_ar.wga &
zcat /fastdata/bop15hjb/hjb/GT_ref/Greattit.Zebrafinch.Flycatcher.wga.bed.gz | ~/parus_indel/summary_analyses/wga2fa.py -out_stem gt_all.wga &

wait


