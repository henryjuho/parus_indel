library(ggplot2)
library(gridExtra)
library(gtools)
library(viridis)


in_file = '../gene_proximity_analyses/gt_nc_v_ar_2kb_wind_cds.results.csv'
window_size = 2

prox_data = read.csv(in_file)

prox_data$bin = factor(prox_data$bin, levels=mixedsort(as.character(subset(prox_data, sel_type=='sel' & var_type=='ins')$bin)))
prox_data$mean_gamma = as.numeric(as.character(prox_data$shape)) * as.numeric(as.character(prox_data$scale)) * -1
prox_data$distance = as.numeric(prox_data$bin) * window_size

ins_theta = subset(prox_data, sel_type == 'sel' & var_type == 'ins', select=c(bin, theta, distance))
colnames(ins_theta) = c('bin', 'ins_theta', 'distance')

cor.test(as.numeric(ins_theta$ins_theta), ins_theta$distance, method='spearman', exact=NULL)
length(ins_theta$distance)

del_theta = subset(prox_data, sel_type == 'sel' & var_type == 'del', select=c(bin, theta, distance))
colnames(del_theta) = c('bin', 'del_theta', 'distance')

cor.test(as.numeric(del_theta$del_theta), del_theta$distance, method='spearman', exact=NULL)

rdi_data = as.data.frame(cbind(ins_theta, del_theta$del_theta))
rdi_data$rdi = rdi_data$del_theta / rdi_data$ins_theta

cor.test(as.numeric(rdi_data$rdi), as.numeric(rdi_data$bin), method='spearman', exact=NULL)

pal = viridis(n=3)[1:3]

theta_plot = ggplot(subset(prox_data, sel_type=='sel'), aes(x=distance, y=theta*10000, colour=toupper(var_type))) +
      geom_hline(yintercept=2.93, colour=pal[1], size=1, linetype=2) +
      geom_hline(yintercept=1.7, colour=pal[2], size=1, linetype=2) +
      geom_point(stat='identity', size = 2) +
      theme_bw() +
      # geom_smooth(method='lm', se=FALSE)+
      #xlim(0,10000)+
      ylab(expression(theta * x *10^-4))  + xlab('Distance from exons (kb)') +
      theme(legend.title = element_blank(), legend.position = c(0.15, 0.87), legend.background = element_blank()) +
      scale_colour_manual(values=pal)

gamma_plot = ggplot(subset(prox_data, sel_type=='sel'), aes(x=distance, y=mean_gamma, colour=toupper(var_type))) +
      geom_point(stat='identity', size = 2) +
      # geom_smooth(method='lm', se=FALSE) +
      theme_bw() +
      ylim(-40, 0)+
      ylab(expression(gamma))  + xlab('Distance from exons (bp)') +
      theme(legend.position = c(0.85, 0.15), legend.background = element_blank(), legend.title=element_blank())+
      scale_colour_manual(values=pal)
#
# rdi_plot = ggplot(rdi_data, aes(x=distance, y=rdi)) +
#       geom_point(stat='identity', size = 2) +
#       theme_bw() +
#       #xlim(0,10000)+
#       geom_smooth(method='lm') +
#       ylab('rDI')  + xlab('Distance from conserved region (bp)')

pdf(file='linked_sel_plot.pdf', width=3, height=3)

theta_plot

dev.off()

#sum stats plot
sum_file = read.delim('../gene_proximity_analyses/gt_prox_summary_stats.txt')
call_file = read.delim('../gene_proximity_analyses/gt_prox_call.txt')

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
      # geom_smooth(method='lm', se=FALSE)+
      #xlim(0,10000)+
      ylab(expression(pi))  + xlab('Distance from exons (kb)') +
      theme(legend.title = element_blank(), legend.position = c(0.15, 0.8), legend.background = element_blank()) +
      scale_colour_manual(values=pal)

# gamma_plot = ggplot(subset(all_data, variation!='SNP' & variation!='INDEL'), aes(x=dist, y=tajD, colour=variation)) +
#       geom_point(stat='identity', size = 2) +
#       theme_bw() +
#       #xlim(0,10000)+
#       ylab("Tajima's D")  + xlab('Distance from CDS (kb)') +
#       theme(legend.title = element_blank(), legend.position = 'none', legend.background = element_blank())
#
# rdi_plot = ggplot(rdi_data, aes(x=distance, y=rdi)) +
#       geom_point(stat='identity', size = 2) +
#       theme_bw() +
#       #xlim(0,10000)+
#       geom_smooth(method='lm') +
#       ylab('rDI')  + xlab('Distance from CDS (kb)')

pdf('linked_sel_pi.pdf', width=3, height=3)

theta_plot

dev.off()

pdf('dist_gamma_plot.pdf', width=3, height=3)

gamma_plot

dev.off()
