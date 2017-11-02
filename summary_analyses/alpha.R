# Title     : alpha calulation for dmel INDELs
# Objective : calculates alpha
# Created by: henryjuho
# Created on: 09/10/2017

calc_alpha <- function(dn, ds, pn, ps) {
  alpha_val = 1 - ((ds * pn) / (dn * ps))
  return(alpha_val)
}

cds_indel_subs = read.delim('/Users/henryjuho/sharc_fastdata/GT_data/BGI_BWA_GATK/Summary_stats/gt_cds_indel_divergence.txt')
non_coding_indel_subs = read.delim('/Users/henryjuho/sharc_fastdata/GT_data/BGI_BWA_GATK/Summary_stats/gt_noncoding_indel_divergence.txt')

indel_poly = read.delim('/Users/henryjuho/sharc_fastdata/GT_data/BGI_BWA_GATK/Summary_stats/bgi_10birds.indels.summary_stats.txt')

dn_data = subset(cds_indel_subs, chromo!='XHet' & chromo!='X' & chromo!='YHet' &
    chromo!='Y' & chromo!='U' & chromo!='dmel_mitochondrion_genome' & chromo != 'chrZ')

ds_data = subset(non_coding_indel_subs, chromo!='XHet' & chromo!='X' & chromo!='YHet' &
    chromo!='Y' & chromo!='U' & chromo!='dmel_mitochondrion_genome' & chromo != 'chrZ')

cds_poly = subset(indel_poly, bin=='CDS' & type=='indel')
non_coding_poly = subset(indel_poly, bin=='non-coding' & type=='indel')

d_n = sum(dn_data$indels)
d_s = sum(ds_data$indels)

p_n =  as.numeric(cds_poly$seg_sites)
p_s = as.numeric(non_coding_poly$seg_sites)

#print(c(d_n, d_s, p_n, p_s))

indel_alpha = calc_alpha(d_n, d_s, p_n, p_s)

print(indel_alpha)