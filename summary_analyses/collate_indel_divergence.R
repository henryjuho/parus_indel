library(ggplot2)
library(dplyr)

cds_div = read.delim('/Users/henryjuho/sharc_fastdata/GT_data/BGI_BWA_GATK/Summary_stats/gt_cds_indel_divergence.txt')
non_coding_div = read.delim('/Users/henryjuho/sharc_fastdata/GT_data/BGI_BWA_GATK/Summary_stats/gt_noncoding_indel_divergence.txt')
line_div = read.delim('/Users/henryjuho/sharc_fastdata/GT_data/BGI_BWA_GATK/Summary_stats/gt_ancLINEs_indel_divergence.txt')

cds_div$type = 'cds'
non_coding_div$type = 'non-coding'
line_div$type = 'ancestral LINEs'

div = rbind(cds_div, non_coding_div, line_div)
div$type = factor(div$type, levels = c('ancestral LINEs', 'non-coding', 'cds'))

all_chr = summarise(group_by(subset(div, chromo != 'X' & chromo != 'XHet' & chromo != 'YHet' & chromo != 'chrZ'), type),
    chromo='autosomes', indels=sum(indels), callable=sum(callable))

all_chr$divergence = all_chr$indels / all_chr$callable

div_plot = ggplot(all_chr, aes(x=type, y=divergence)) +
    geom_bar(stat='identity', position='dodge') +
    theme_bw(base_size = 10) +
    xlab('') + ylab('divergence')

pdf(file='indel_divergence.pdf', width=3, height=3)
div_plot
dev.off()

png(file='indel_divergence.png', width=3, height=3, units='in', res=320)
div_plot
dev.off()

out_csv = rbind(div, all_chr)
write.csv(out_csv, file='indel_divergence.csv', row.names=FALSE)