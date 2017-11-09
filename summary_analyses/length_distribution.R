# Title     : gt indel lengths
# Objective : generates summary plot of indel lengths
# Created by: henryjuho
# Created on: 28/09/2017

library(ggplot2)
library(dplyr)
library(gridExtra)

cbPalette <- c("#E69F00", 'tomato 3', 'steel blue', "#999999", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

length_data = read.delim('gt_indel_lengths.txt')

# prop length distro
non_coding_pol = subset(length_data, variant!='indel' & region=='non-coding')
ins_data = subset(non_coding_pol, variant='ins')
no_ins = sum(ins_data$count)
ins_data$prop = ins_data$count/no_ins
del_data = subset(non_coding_pol, variant='del')
no_del = sum(del_data$count)
del_data$prop = del_data$count/no_del
prop_data = rbind(ins_data, del_data)

# length prop plot
prop_plot = ggplot(prop_data, aes(x=length, y=prop, fill=variant))+
    geom_bar(stat='identity', position='dodge') +
    theme_bw(base_size = 11) +
    labs(x='INDEL length (bp)', y='Proportion of variants') +
    scale_fill_manual(values=cbPalette)

len_plot = ggplot(subset(length_data, variant!='indel' & region!='gwide'), aes(x=length, y=count, fill=variant))+
    geom_bar(stat='identity', position='dodge') +
    theme_bw(base_size = 11) +
    labs(x='INDEL length (bp)', y='Number of variants') +
    scale_fill_manual(values=cbPalette) +
    facet_wrap(~region, ncol=1, scales='free')

png('gt_lengths.png', height=5, width=10, units='in', res=320)
grid.arrange(len_plot, prop_plot, nrow=1)
dev.off()
