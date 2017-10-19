library(ggplot2)
library(gridExtra)
library(dplyr)

# The palette with grey:
cbPalette <- c("#E69F00", 'tomato 3', 'steel blue', "#999999", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

snp_sum_stats <- read.delim('~/sharc_fastdata/GT_data/BGI_BWA_GATK/Summary_stats/bgi_10birds.snps.summary_stats.txt')
indel_sum_stats <- read.delim('~/sharc_fastdata/GT_data/BGI_BWA_GATK/Summary_stats/bgi_10birds.indels.summary_stats.txt')

snp_no_ar = subset(snp_sum_stats, bin!='AR')

plot_stats = subset(rbind(indel_sum_stats, snp_no_ar), bin!='CDS_frameshift' & bin!='CDS_non_frameshift' & bin!='intron' & bin!='intergenic')

#snp_sum_stats$bin <- as.character(snp_sum_stats$bin)
#for(x in 1:length(snp_sum_stats$bin)){
#  line = snp_sum_stats[x,]
#  if(line$bin == 'intergenic'){snp_sum_stats[x,]$bin <- 'neutral'}
#}

#all_sum_stats = rbind(snp_sum_stats, indel_sum_stats, LINE_sum_stats)

#all_sum_stats$bin = as.character(all_sum_stats$bin)
#for(x in 1:length(all_sum_stats$bin)){
#  line = all_sum_stats[x,]
#  if(line$bin == 'CDS'){all_sum_stats[x,]$bin <- 'coding'}
#}

#all_sum_stats$bin <- factor(all_sum_stats$bin, levels=c('neutral', 'intergenic', 'intron', 'coding'))

# theta_w
tw = ggplot(plot_stats, aes(x=bin, y=theta_w, colour=type))+
  geom_point(stat='identity', position = position_dodge(width=0.9), size = 2) +
  #geom_errorbar(aes(ymin = t_lwr, ymax = t_upr), stat = 'identity', position = position_dodge(width=0.9), width=0.2) +
  theme_bw() +
  scale_color_manual(values=cbPalette) +
  xlab('')  + ylab(expression(theta[w])) +
  theme(legend.title=element_blank(), legend.position='none')

# pi
pi = ggplot(plot_stats, aes(x=bin, y=pi, colour=type))+
  geom_point(stat='identity', position = position_dodge(width=0.9), size = 2) +
  #geom_errorbar(aes(ymin = pi_lwr, ymax = pi_upr), stat = 'identity', position = position_dodge(width=0.9), width=0.2) +
  theme_bw() +
  scale_color_manual(values=cbPalette) +
  xlab('')  + ylab(expression(pi)) +
  theme(legend.title=element_blank(), legend.position='none')

# tajD
tajd = ggplot(plot_stats, aes(x=bin, y=tajD, colour=type))+
  geom_point(stat='identity', position = position_dodge(width=0.9), size = 2) +
  #geom_errorbar(aes(ymin = tajD_lwr, ymax = tajD_upr), stat = 'identity', position = position_dodge(width=0.9), width=0.2) +
  theme_bw() +
  scale_color_manual(values=cbPalette) +
  xlab('')  + ylab("Tajima's D") +
  theme(legend.title=element_blank(), legend.position=c(0.16, 0.4))

pdf(file='gt_summary_stats.pdf', width=9, height=3)

grid.arrange(tajd, tw, pi, nrow=1)

dev.off()
