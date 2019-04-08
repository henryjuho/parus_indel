library(ggplot2)
library(dplyr)
library(gridExtra)

length_data = read.delim('../summary_analyses/gt_indel_lengths.txt')
gwide = subset(length_data, region=='gwide' & variant=='indel')

# prop length distro
non_coding_pol = subset(length_data, variant!='indel' & region=='non-coding')
ins_data = subset(non_coding_pol, variant=='ins')
no_ins = sum(ins_data$count)
ins_data$prop = ins_data$count/no_ins
del_data = subset(non_coding_pol, variant=='del')
no_del = sum(del_data$count)
del_data$prop = del_data$count/no_del
prop_data = rbind(ins_data, del_data)

# length prop plot
# prop_plot = ggplot(prop_data, aes(x=length, y=prop, fill=variant))+
#     geom_bar(stat='identity', position='dodge') +
#     theme_bw(base_size = 11) +
#     labs(x='INDEL length (bp)', y='Proportion of variants')

length_data$region = factor(length_data$region, levels=c('gwide', 'non-coding', 'CDS'),
    labels=c('genome-wide', 'non-coding', 'coding sequence'))

len_plot = ggplot(subset(length_data, variant!='indel' & region!='non-coding'), aes(x=length, y=count, fill=variant))+
    geom_bar(stat='identity', position='dodge') +
    theme_bw(base_size = 11) +
    labs(x='INDEL length (bp)', y='Number of variants') +
    facet_wrap(~region, ncol=1, scales='free') +
    theme(legend.title=element_blank(), legend.position=c(0.8, 0.8))

# pdf('gt_indel_lengths.pdf', height=6, width=6)
# len_plot
# dev.off()

# percents by length
total_indel = sum(gwide$count)
gwide$prop = gwide$count / total_indel
gwide$cum = cumsum(gwide$prop)
print(gwide)