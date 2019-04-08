library(ggplot2)
library(viridis)

non_r = read.csv('polarisation_summary_indels_nonrep.csv')
ar = subset(read.csv('polarisation_summary_indels_ar.csv'), cat!='not_aligned' & cat!='low_coverage')
#ar$percent = (ar$count / (19851 * (1.24 + 0.3))) * 100
all_sum = rbind(non_r, ar)

print(all_sum)

# make plot of pol success
png('pol_success_indels.png', width=4, height=3, res=320, units='in')

ggplot(subset(all_sum, cat!='total' & cat!='unpolarised'), aes(x=region, y=percent, fill=cat)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=viridis(6),
    labels = c('Ambiguous', 'Hotspot', 'Low coverage', 'Not aligned', 'Polarised')) +
    theme_bw() + theme(legend.title=element_blank()) +
    labs(y='Percentage of INDELs', x='') +
    theme(axis.text.x=element_text(angle=45, hjust=1))

dev.off()

pdf('pol_success_indels.pdf', width=4, height=3)

ggplot(subset(all_sum, cat!='total' & cat!='unpolarised'), aes(x=region, y=percent, fill=cat)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=viridis(6),
    labels = c('Ambiguous', 'Hotspot', 'Low coverage', 'Not aligned', 'Polarised')) +
    theme_bw() + theme(legend.title=element_blank()) +
    labs(y='Percentage of INDELs', x='') +
    theme(axis.text.x=element_text(angle=45, hjust=1))

dev.off()