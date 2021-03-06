library(ggplot2)
library(dplyr)
library(viridis)

div = read.delim('gt_indel_div_frameshifts.txt')

div = subset(div, category!='noncoding' & category!='CDS' & category!='ALL')
div$category = as.character(div$category)

snp_div = read.delim('snp_div.txt', header=FALSE)
snp_div$variation = 'SNP'
snp_div$V1 = as.character(snp_div$V1)
snp_div[2,]$V1 = 'cds_non_frameshift'

snp_missing = data.frame(category=c('cds_frameshift'), variation=c('SNP'), divergence=c(0))

snp_div = data.frame(category=as.character(snp_div$V1), variation=snp_div$variation,
                     divergence=as.numeric(snp_div$V2))

snp_div = rbind(snp_div, snp_missing)

snp_div$divergence = snp_div$divergence/10

div = rbind(subset(div, select=c('category', 'variation', 'divergence')), snp_div)

div$category = factor(div$category, levels=rev(c('cds_frameshift', 'cds_non_frameshift', 'introns', 'intergenic', 'AR')))

div_plot = ggplot(div, aes(x=category, y=divergence, fill=variation)) +
    geom_bar(stat='identity', position='dodge') +
    theme_bw(base_size = 10) +
    xlab('') + ylab('divergence') +
    theme(#axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title=element_blank(), legend.position=c(0.85, 0.8),
    legend.background=element_blank()) +
    scale_fill_viridis(discrete=TRUE)+
    guides(fill = guide_legend(nrow = 2)) +
    annotate("text", x = 3.5, y = 0.0065, label = '(SNPs: divergence x 10)', size=2.9)

# # pdf(file='indel_divergence.pdf', width=3, height=3)
# # div_plot
# # dev.off()
#
png(file='indel_divergence.png', width=6, height=3, units='in', res=320)
div_plot
dev.off()
