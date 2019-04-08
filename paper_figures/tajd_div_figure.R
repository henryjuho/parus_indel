library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)

pal = viridis(n=5)
## divergence
div = read.delim('../summary_analyses/gt_indel_div_frameshifts.txt')

div = subset(div, category!='noncoding' & category!='ALL')
div$category = as.character(div$category)

snp_div = read.delim('../summary_analyses/snp_div.txt', header=FALSE)
snp_div$variation = 'SNP'
snp_div$V1 = as.character(snp_div$V1)
# snp_div[2,]$V1 = 'cds_non_frameshift'
snp_missing = data.frame(category=c('cds_frameshift', 'cds_non_frameshift'),
                         variation=c('SNP', 'SNP'),
                         divergence=c(0, 0))

snp_div = data.frame(category=as.character(snp_div$V1), variation=snp_div$variation,
                     divergence=as.numeric(snp_div$V2))
snp_div = rbind(snp_div, snp_missing)
snp_div$divergence = snp_div$divergence/10

div = rbind(subset(div, select=c('category', 'variation', 'divergence')), snp_div)

div$category = factor(div$category,
    levels=rev(c('cds_frameshift', 'cds_non_frameshift', 'CDS', 'introns', 'intergenic', 'AR')))

div_plot = ggplot(div, aes(x=category, y=divergence, fill=variation)) +
    geom_bar(stat='identity', size=2, position = position_dodge(width=0.9)) +
    theme_bw() +
    xlab('') + ylab('divergence') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title=element_blank(), legend.position='None',
    legend.background=element_blank(), plot.title=element_text(hjust=-0.325, vjust=-3)) +
    guides(colour = guide_legend(nrow = 1)) + ggtitle('(b)') +
    scale_x_discrete(labels=rev(c('Frame-shift', 'In-frame', 'CDS', 'introns', 'intergenic', 'AR'))) +
    #annotate("text", x = 3.5, y = 0.0065,
    #          label = '(SNPs: divergence x 10)', size=2.9) +
    scale_fill_manual(values=pal)

## tajd
summary_data = read.delim('../summary_analyses/bgi10_stats_nocorrection.txt', na=0)

summary_data = subset(summary_data, category!='noncoding' & category!='ALL'
                      & category != 'AR_frameshift' & category != 'AR_non_frameshift' & category != '4fold',
                      select=c('category', 'variation', 'tajD'))


# reset SNP cat indel res to NA
# 0fold
summary_data[2,]$tajD = NA
summary_data[3,]$tajD = NA
summary_data[4,]$tajD = NA

# 4fold
# summary_data[6,]$tajD = NA
# summary_data[7,]$tajD = NA
# summary_data[8,]$tajD = NA

# reset SNP in in frame cat to value from CDS snp then remove CDS cat
# summary_data[17,]$tajD = summary_data[21,]$tajD
# summary_data = subset(summary_data, category!='CDS')

# relevel
summary_data$category = factor(summary_data$category,
  levels=rev(c('cds_frameshift', 'cds_non_frameshift', 'CDS', '0fold', 'nonsense',
               'introns', 'intergenic', 'AR')))

lables= rev(c('Frameshift', 'In-frame', 'CDS', '0fold', 'nonsense', 'introns', 'intergenic', 'AR'))
# tajD
tajd = ggplot(summary_data, aes(x=category, y=tajD, fill=variation))+
  geom_bar(stat='identity', position = position_dodge(width=0.9), size = 2) +
  theme_bw() +
  xlab('')  + ylab("Tajima's D") +
  #ylim(-1.3, -0.2) +
  theme(legend.title=element_blank(), legend.position=c(0.2, 0.3), legend.key.size = unit(.4, "cm"),
  axis.text.x=element_text(angle=45, hjust=1),
  plot.title=element_text(hjust=-0.28, vjust=-3)) + ggtitle('(a)') +
  guides(fill=guide_legend(ncol=1))+
  scale_x_discrete(labels=rev(c('Frame-shift', 'In-frame', 'CDS', '0fold',
                                'nonsense', 'introns', 'intergenic', 'AR'))) +
  scale_fill_manual(values=pal)

pdf('tajd_div_plot.pdf', width=6, height=3)

grid.arrange(tajd, div_plot, nrow=1)

dev.off()

