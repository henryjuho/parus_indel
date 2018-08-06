library(ggplot2)
library(dplyr)
library(viridis)

div = read.delim('gt_indel_div.txt')

div = subset(div, category!='noncoding')
snp_div = read.delim('snp_div.txt', header=FALSE)
snp_div$variation = 'SNP'
snp_div = data.frame(cbind(as.character(snp_div$V1),snp_div$variation, snp_div$V2))
colnames(snp_div) = c('category', 'variation', 'divergence')
snp_div$divergence = as.numeric(as.character(snp_div$divergence))/10
div = rbind(subset(div, select=c('category', 'variation', 'divergence')), snp_div)

div$category = factor(div$category,
  levels=rev(c('CDS', 'UCNE', 'ALL', 'intergenic', 'introns', 'AR')))

div_plot = ggplot(div, aes(x=category, y=divergence, fill=variation)) +
    geom_bar(stat='identity', position='dodge') +
    theme_bw(base_size = 10) +
    xlab('') + ylab('divergence') +
    theme(#axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title=element_blank(), legend.position=c(0.85, 0.8),
    legend.background=element_blank()) +
    scale_fill_viridis(discrete=TRUE)+
    guides(fill = guide_legend(nrow = 2))

# pdf(file='indel_divergence.pdf', width=3, height=3)
# div_plot
# dev.off()

png(file='indel_divergence.png', width=6, height=3, units='in', res=320)
div_plot
dev.off()
