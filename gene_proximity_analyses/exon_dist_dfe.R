library(ggplot2)
library(viridis)

dfe = read.delim('nc_dfe_exon_dist.txt')

dfe$cat_theta = dfe$theta * dfe$prop

pdf('nc_exons_dist_dfe_prop.pdf', height=3, width=6)

ggplot(dfe, aes(x=dist*2, y=prop, fill=cat)) +
    geom_bar(stat='identity') +
    # geom_smooth(method='lm', se=F) +
    theme_bw() +
    scale_fill_manual(values=viridis(3)) +
    facet_wrap(~var, nrow=1, ncol=2) +
    xlab('Distance from exons (kb)') + ylab('Density') +
    theme(legend.title=element_blank(), legend.position=c(0.08, 0.8),
    legend.key.size = unit(.4, "cm"), legend.text=element_text(size=8)) +
    guides(colour=guide_legend(ncol=2))

dev.off()


pdf('nc_exons_dist_dfe_theta.pdf', height=3, width=6)

ggplot(dfe, aes(x=dist*2, y=cat_theta*10000, colour=toupper(var))) +
    geom_point(stat='identity') +
    # geom_smooth(method='lm', se=F) +
    theme_bw() +
    scale_colour_manual(values=viridis(3)) +
    facet_wrap(~cat, nrow=1, ncol=2) + #, scales='free') +
    xlab('Distance from exons (kb)') + ylab(expression(theta*x*10^-4)) +
    theme(legend.title=element_blank(), legend.position=c(0.07, 0.15),
    legend.key.size = unit(.4, "cm"), legend.text=element_text(size=8)) +
    guides(fill=guide_legend(ncol=2))

dev.off()

neu_del = subset(dfe, var=='del' & cat=='0 - 1')
neu_del
neu_ins = subset(dfe, var=='ins' & cat=='0 - 1')

sel_del = subset(dfe, var=='del' & cat=='>1')
sel_ins = subset(dfe, var=='ins' & cat=='>1')

cor.test(neu_del$cat_theta, neu_del$dist*2000, method='spearman')
cor.test(neu_ins$cat_theta, neu_ins$dist*2000, method='spearman')
cor.test(sel_del$cat_theta, sel_del$dist*2000, method='spearman')
cor.test(sel_ins$cat_theta, sel_ins$dist*2000, method='spearman')

write.csv(dfe, 'split_dfe_dist.csv', row.names=F)