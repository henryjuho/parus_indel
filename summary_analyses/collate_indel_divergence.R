library(ggplot2)
library(dplyr)

div = read.delim('gt_indel_div.txt')

div = subset(div, category!='noncoding')

div$category = factor(div$category,
  levels=rev(c('CDS', 'UCNE', 'ALL', 'intergenic', 'introns', 'AR')))

div_plot = ggplot(div, aes(x=category, y=divergence, colour=variation)) +
    geom_point(stat='identity') +
    theme_bw(base_size = 10) +
    xlab('') + ylab('divergence') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title=element_blank(), legend.position=c(0.4, 0.1),
    legend.background=element_blank()) +
    guides(colour = guide_legend(nrow = 1))

# pdf(file='indel_divergence.pdf', width=3, height=3)
# div_plot
# dev.off()

png(file='indel_divergence.png', width=3, height=3, units='in', res=320)
div_plot
dev.off()
