library(ggplot2)
library(gridExtra)
library(gtools)

sum_file = read.delim('gt_prox_summary_stats_raw.txt')
call_file = read.delim('gt_prox_call.txt')

all_data = cbind(sum_file, call_file)
all_data[7] = NULL
all_data[7] = NULL

all_data$dist = as.numeric(sapply(strsplit(as.character(all_data$category), 'bin'), "[[", 2)) * 2

all_data$pi_per_site = all_data$pi / all_data$callable

ins_theta = subset(all_data, variation == 'INS', select=c(category, pi_per_site, dist))
colnames(ins_theta) = c('bin', 'ins_theta', 'distance')

cor.test(as.numeric(ins_theta$ins_theta), ins_theta$distance, method='spearman', exact=NULL)

del_theta = subset(all_data, variation == 'DEL', select=c(category, pi_per_site, dist))
colnames(del_theta) = c('bin', 'del_theta', 'distance')

cor.test(as.numeric(del_theta$del_theta), del_theta$distance, method='spearman', exact=NULL)

rdi_data = as.data.frame(cbind(ins_theta, del_theta$del_theta))
rdi_data$rdi = rdi_data$del_theta / rdi_data$ins_theta

cor.test(as.numeric(rdi_data$rdi), as.numeric(rdi_data$bin), method='spearman', exact=NULL)

theta_plot = ggplot(subset(all_data, variation!='SNP' & variation!='INDEL'), aes(x=dist, y=pi_per_site, colour=variation)) +
      geom_point(stat='identity', size = 2) +
      theme_bw() +
      geom_smooth(method='lm')+
      #xlim(0,10000)+
      ylab(expression(pi))  + xlab('Distance from CDS (kb)') +
      theme(legend.title = element_blank(), legend.position = c(0.15, 0.8), legend.background = element_blank())

gamma_plot = ggplot(subset(all_data, variation!='SNP' & variation!='INDEL'), aes(x=dist, y=tajD, colour=variation)) +
      geom_point(stat='identity', size = 2) +
      theme_bw() +
      #xlim(0,10000)+
      ylab("Tajima's D")  + xlab('Distance from CDS (kb)') +
      theme(legend.title = element_blank(), legend.position = 'none', legend.background = element_blank())

rdi_plot = ggplot(rdi_data, aes(x=distance, y=rdi)) +
      geom_point(stat='identity', size = 2) +
      theme_bw() +
      #xlim(0,10000)+
      geom_smooth(method='lm') +
      ylab('rDI')  + xlab('Distance from CDS (kb)')

png(file='2kb_nc_cds_dist_sum.png', width=9, height=3, units='in', res=360)

grid.arrange(theta_plot, gamma_plot, rdi_plot, nrow=1)

dev.off()
