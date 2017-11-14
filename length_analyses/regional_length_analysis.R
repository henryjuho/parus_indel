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
    call=sum(callable), total_bp=sum(total),
    n_indel=sum(n_indel), total_event_length=sum(total_event_length),
    freq_sum=sum(freq_sum))

gwide_lens$indel_per_base = (gwide_lens$total_bp / 20.0) / gwide_lens$call
gwide_lens$region = factor(gwide_lens$region,
    levels=c('gwide', 'intergenic', 'intron', 'CDS'))

gwide_lens$mean_length = gwide_lens$total_event_length / gwide_lens$n_indel
gwide_lens$mean_freq = gwide_lens$freq_sum / gwide_lens$n_indel

len_rdi_data = subset(gwide_lens, indel_type=='del', select=c(region, total_event_length))
len_rdi_data$ins_bp = subset(gwide_lens, indel_type=='ins', select=c(region, total_event_length))$total_event_length
len_rdi_data$rdi = len_rdi_data$total_event_length / len_rdi_data$ins_bp
len_rdi_data$type = 'Length ratio'

colnames(len_rdi_data) = c('region', 'del', 'ins', 'rdi', 'type')

mute_rdi_data = subset(gwide_lens, indel_type=='del', select=c(region, n_indel))
mute_rdi_data$ins_bp = subset(gwide_lens, indel_type=='ins', select=c(region, n_indel))$n_indel
mute_rdi_data$rdi = mute_rdi_data$n_indel / mute_rdi_data$ins_bp
mute_rdi_data$type = 'Mutation ratio'

colnames(mute_rdi_data) = c('region', 'del', 'ins', 'rdi', 'type')

rdi_data = rbind(len_rdi_data, mute_rdi_data)

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

freq_plot = ggplot(gwide_lens, aes(x=region, y=mean_freq, colour=indel_type)) +
    geom_point(stat='identity', position=position_dodge(width=.5)) +
    theme_bw() +
    scale_colour_manual(values=cbPalette) +
    ylab('Mean derived allele frequency') + xlab('') +
    theme(legend.position='none')

rdi_plot = ggplot(rdi_data, aes(x=region, y=rdi, colour=type)) +
    geom_point(stat='identity') +
    theme_bw() +
    scale_colour_manual(values=c('steel blue', "#999999")) +
    ylab('Ratio of deletions to insertions') + xlab('') +
    theme(legend.position=c(0.3, 0.85), legend.title=element_blank())

mean_len = ggplot(gwide_lens, aes(x=region, y=mean_length, fill=indel_type)) +
    geom_bar(stat='identity', position='dodge') +
    theme_bw() +
    scale_fill_manual(values=cbPalette) +
    ylab('Mean INDEL length') + xlab('') +
    theme(legend.position='none')

png('regional_length_plots.png', height=6, width=9, units='in', res=320)
grid.arrange(length_plot, mean_len, freq_plot, rate_plot, rdi_plot, nrow=2)
dev.off()
