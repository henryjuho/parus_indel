# Title     : regional length analysis
# Objective : plot INDEL length stats in different regions
# Created by: henryjuho
# Created on: 13/11/2017

library(dplyr)
library(ggplot2)
library(gridExtra)

cbPalette <- c("#E69F00", 'tomato 3', 'steel blue', "#999999", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

bgi_lengths = read.csv('bgi_10birds_bp_loss_gain.csv')

gwide_lens = summarise(group_by(subset(bgi_lengths, chr!='chrZ'), indel_type, region),
    call = sum(callable), total_bp = sum(total))

gwide_lens$indel_per_base = (gwide_lens$total_bp / 20.0) / gwide_lens$call
gwide_lens$region = factor(gwide_lens$region,
    levels=c('gwide', 'intergenic', 'intron', 'CDS'))

rdi_data = subset(gwide_lens, indel_type=='del', select=c(region, total_bp))
rdi_data$ins_bp = subset(gwide_lens, indel_type=='ins', select=c(region, total_bp))$total_bp
rdi_data$rdi = rdi_data$total_bp / rdi_data$ins_bp

length_plot = ggplot(gwide_lens, aes(x=region, y=total_bp/20.0, fill=indel_type)) +
    geom_bar(stat='identity', position='dodge') +
    theme_bw() +
    scale_fill_manual(values=cbPalette) +
    ylab('INDEL length per haplotype') + xlab('') +
    theme(legend.title=element_blank(), legend.position=c(0.8, 0.8))

rate_plot = ggplot(gwide_lens, aes(x=region, y=indel_per_base, colour=indel_type)) +
    geom_point(stat='identity', position=position_dodge(width=.5)) +
    theme_bw() +
    scale_colour_manual(values=cbPalette) +
    ylab('Per base INDEL rate') + xlab('') +
    theme(legend.position='none')

rdi_plot = ggplot(rdi_data, aes(x=region, y=rdi,)) +
    geom_point(stat='identity') +
    theme_bw() +
    scale_fill_manual(values=cbPalette) +
    ylab('Ratio of deleted bp to inserted bp') + xlab('') +
    theme(legend.position='none')


png('regional_length_plots.png', height=3, width=9, units='in', res=320)
grid.arrange(length_plot, rate_plot, rdi_plot, nrow=1)
dev.off()
