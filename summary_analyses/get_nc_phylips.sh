#!/bin/bash

source ~/.bash_profile

#$-l h_rt=72:00:00
#$-l mem=6G
#$-l rmem=2G

#$-pe openmp 4
#$-P evolgen
#$-q evolgen.q

#$-N gt_non-coding_phylips.sh
#$-o /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/snp_divergence/non-coding/phylips/nc_phylip_gen.out
#$-e /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/snp_divergence/non-coding/phylips/nc_phylip_gen.error

#$-V

cd /fastdata/bop15hjb/GT_data/BGI_BWA_GATK/snp_divergence/non-coding/phylips/

bedtools intersect -a /fastdata/bop15hjb/GT_ref/Greattit.Zebrafinch.Flycatcher.wga.bed.gz -b /fastdata/bop15hjb/GT_ref/Greattit.Zebrafinch.Flycatcher.ancLINEs.sorted.bed.gz | ~/parus_indel/summary_analyses/wga2phy_nc.py > gt_ar_lines.phy &
bedtools intersect -a /fastdata/bop15hjb/GT_ref/Greattit.Zebrafinch.Flycatcher.wga.bed.gz -b /fastdata/bop15hjb/GT_ref/gt_introns.bed.gz | ~/parus_indel/summary_analyses/wga2phy_nc.py > gt_introns.phy &
bedtools intersect -a /fastdata/bop15hjb/GT_ref/Greattit.Zebrafinch.Flycatcher.wga.bed.gz -b /fastdata/bop15hjb/GT_ref/gt_intergenic.bed.gz | ~/parus_indel/summary_analyses/wga2phy_nc.py > gt_intergenic.phy &
bedtools intersect -a /fastdata/bop15hjb/GT_ref/Greattit.Zebrafinch.Flycatcher.wga.bed.gz -b /fastdata/bop15hjb/GT_ref/gt_ucne.filtered2.sorted.bed.gz | ~/parus_indel/summary_analyses/wga2phy_nc.py > gt_ucne.phy &

wait


