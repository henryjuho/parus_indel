library(ggplot2)
library(viridis)
library(gridExtra)

cds_dat_raw = read.csv('gt_cds_indels_length_sum_raw.csv')
cds_dat_raw$region = 'CDS'
cds_dat_raw$sfs = 'raw'
nc_dat_raw = read.csv('gt_nc_indels_length_sum_raw.csv')
nc_dat_raw$sfs = 'raw'
nc_dat_raw$region = 'Non-coding'
cds_dat = read.csv('gt_cds_indels_length_sum.csv')
cds_dat$region = 'CDS'
cds_dat$sfs = 'corrected'
nc_dat = read.csv('gt_nc_indels_length_sum.csv')
nc_dat$sfs = 'corrected'
nc_dat$region = 'Non-coding'

all_dat = rbind(cds_dat, nc_dat, cds_dat_raw, nc_dat_raw)

str(all_dat)

# pi plots
pi_plot = ggplot(all_dat, aes(x=length, y=pi, colour=var_type)) +
    geom_point(stat='identity') +
    theme_bw(base_size=10) +
    scale_colour_manual(values=viridis(3)) +
    facet_grid(sfs~region, scales='free_y') +
    labs(x='INDEL length (bp)', y=expression(pi)) +
    theme(legend.title=element_blank(), legend.position=c(0.9, 0.8))

# taj cor
nc_ins = subset(all_dat, region=='Non-coding' & var_type=='ins' & sfs=='raw')
ncins_cor = cor.test(nc_ins$tajd, nc_ins$length, method='spearman')

nc_del = subset(all_dat, region=='Non-coding' & var_type=='del' & sfs=='raw')
ncdel_cor = cor.test(nc_del$tajd, nc_del$length, method='spearman')

print('insertions')
print(ncins_cor)
print('deletions')
print(ncdel_cor)

avg_dat = data.frame(tajd=c(-0.972881511536, -1.26992351684, -0.325084323868, -0.463188762236),
                     var_type=c('ins', 'del', 'ins', 'del'),
                     region=c('CDS', 'CDS', 'Non-coding', 'Non-coding'))


taj_plot = ggplot(subset(all_dat, sfs=='raw'), aes(x=length, y=tajd, colour=var_type, fill=var_type)) +
    geom_hline(data = avg_dat, aes(yintercept=tajd, colour=var_type)) +
    geom_point(stat='identity') +
    theme_bw(base_size=10) +
    scale_colour_manual(values=viridis(3)) +
    facet_grid(sfs~region) +
    labs(x='INDEL length (bp)', y="Tajima's D")+
    theme(legend.title=element_blank(), legend.position=c(0.9, 0.8),
    plot.title = element_text(size=10, colour='steel blue', face='bold')) +
    ggtitle(paste("non coding insertions: Spearman's", expression(rho), "=", signif(ncins_cor$estimate, 3),
           'p =', signif(ncins_cor$p.value, 6),
           "\nnon coding deletions: Spearman's", expression(rho), "=", signif(ncdel_cor$estimate, 3),
           'p =', signif(ncdel_cor$p.value, 3),sep=' '))

png('indel_length_summary_stats.png', width=12, height=6, res=320, units='in')

grid.arrange(pi_plot, taj_plot, ncol=2)

dev.off()

# supp mat plot

pdf('indel_length_tajd.pdf', width=6, height=3)

ggplot(all_dat, aes(x=length, y=tajd, colour=var_type, fill=var_type)) +
    geom_hline(data = avg_dat, aes(yintercept=tajd, colour=var_type)) +
    geom_point(stat='identity') +
    theme_bw(base_size=10) +
    scale_colour_manual(values=viridis(3)) +
    facet_wrap(~region, nrow=1) +
    labs(x='INDEL length (bp)', y="Tajima's D")+
    theme(legend.title=element_blank(), legend.position=c(0.4, 0.8), legend.background=element_blank())

dev.off()