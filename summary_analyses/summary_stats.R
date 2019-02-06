library(ggplot2)
library(gridExtra)
library(dplyr)
library(viridis)

summary_data = read.delim('bgi10_stats.txt', na=0)
summary_data_sfs_correct = read.delim('bgi10_stats_sfs_corrected.txt', na=0)

summary_data$sfs = 'raw'
summary_data_sfs_correct$sfs = 'corrected'

all_sum = rbind(summary_data, summary_data_sfs_correct)

call_data = read.delim('bgi10_call.txt')

all_data = dplyr::left_join(all_sum, call_data, by=c('category', 'variation'))

all_data$tw_per_site = all_data$theta_w / all_data$call
all_data$pi_per_site = all_data$pi / all_data$call

str(all_data)

write.csv(all_data, 'bgi10_summary_stats.csv', row.names=FALSE)

all_data = subset(all_data, category != 'noncoding_noUCNEs'
  & category != 'UCNE' & category != 'ALL'
  & category != 'noncoding')
#
# all_data$category = factor(all_data$category,
#   levels=rev(c('cds_frameshift', 'cds_non_frameshift', 'CDS', '0fold', 'nonsense', '4fold', 'introns', 'intergenic', 'AR')))
#

indel_data = subset(all_data, variation!='SNP' & category!='0fold' & category!='4fold' & category!='nonsense')
snp_data = subset(all_data, variation=='SNP' & category!='cds_frameshift' & category!='cds_non_frameshift' &
                  sfs=='raw')
str(snp_data)
filtered_data = rbind(indel_data, snp_data)

indel_data$category = factor(indel_data$category, levels=c('cds_frameshift', 'cds_non_frameshift', 'CDS',
                                                           'introns', 'intergenic', 'AR'))

snp_data$category = factor(snp_data$category, levels=c('0fold', 'nonsense', 'CDS', '4fold',
                                                         'introns', 'intergenic', 'AR'))

# theta_w
tw = ggplot(indel_data, aes(x=category, y=tw_per_site, fill=variation))+
  geom_bar(stat='identity', position = position_dodge(width=0.9), size = 2) +
  theme_bw() +
  facet_wrap(~sfs, ncol=1) +
  xlab('')  + ylab(expression(theta[w])) +
  theme(legend.title=element_blank(), legend.position='none',
  axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_viridis(discrete=T)

# pi
pi = ggplot(indel_data, aes(x=category, y=pi_per_site, fill=variation))+
  geom_bar(stat='identity', position = position_dodge(width=0.9), size = 2) +
  theme_bw() +
  facet_wrap(~sfs, ncol=1) +
  xlab('')  + ylab(expression(pi)) +
  theme(legend.title=element_blank(), legend.position='none',
  axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_viridis(discrete=T)

# tajD
tajd = ggplot(indel_data, aes(x=category, y=tajD, fill=variation))+
  geom_bar(stat='identity', position = position_dodge(width=0.9), size = 2) +
  theme_bw() +
  facet_wrap(~sfs, ncol=1) +
  xlab('')  + ylab("Tajima's D") +
  #ylim(-1.3, -0.2) +
  theme(legend.title=element_blank(), legend.position=c(0.83, 0.65),
  axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_viridis(discrete=T)
#
# # pdf(file='gt_summary_stats.pdf', width=9, height=3)
# #
# # grid.arrange(tajd, tw, pi, nrow=1)
# #
# # dev.off()
#
# # pdf(file='gt_tajd.pdf', width=3, height=3)
# #
# # tajd
# #
# # dev.off()
#
png(file='gt_summary_stats.png', width=9, height=6, units='in', res=360)

grid.arrange(tajd, tw, pi, nrow=1)

dev.off()

# SNP sum stats #

# theta_w
tw_snp = ggplot(snp_data, aes(x=category, y=tw_per_site))+
  geom_bar(stat='identity', position = position_dodge(width=0.9), size = 2) +
  theme_bw() +
  xlab('')  + ylab(expression(theta[w])) +
  theme(legend.title=element_blank(), legend.position='none',
  axis.text.x=element_text(angle=45, hjust=1))

# pi
pi_snp = ggplot(snp_data, aes(x=category, y=pi_per_site))+
  geom_bar(stat='identity', position = position_dodge(width=0.9), size = 2) +
  theme_bw() +
  xlab('')  + ylab(expression(pi)) +
  theme(legend.title=element_blank(), legend.position='none',
  axis.text.x=element_text(angle=45, hjust=1))

# tajD
tajd_snp = ggplot(snp_data, aes(x=category, y=tajD))+
  geom_bar(stat='identity', position = position_dodge(width=0.9), size = 2) +
  theme_bw() +
  xlab('')  + ylab("Tajima's D") +
  #ylim(-1.3, -0.2) +
  theme(legend.title=element_blank(), legend.position=c(0.83, 0.65),
  axis.text.x=element_text(angle=45, hjust=1))

png(file='gt_summary_stats_snps.png', width=9, height=3, units='in', res=360)

grid.arrange(tajd_snp, tw_snp, pi_snp, nrow=1)

dev.off()