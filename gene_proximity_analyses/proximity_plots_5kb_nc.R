# Title     : distance anavar 5kb nc
# Objective : plots anavar results for 5kb distance windows using non coding data
# Created by: henryjuho
# Created on: 28/11/2017

library(ggplot2)
library(gridExtra)
library(gtools)

cbPalette <- c("#E69F00", 'tomato 3', 'steel blue', "#999999", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

prox_data = subset(read.csv('gt_sel_neu_ref_cdsdist_5kb.results.csv'), bin!='47-157')

prox_data$bin = factor(prox_data$bin, levels=mixedsort(as.character(subset(prox_data, sel_type=='sel' & var_type=='ins')$bin)))
prox_data$mean_gamma = as.numeric(as.character(prox_data$shape)) * as.numeric(as.character(prox_data$scale)) * -1
prox_data$distance = as.numeric(prox_data$bin) * 5000

ins_theta = subset(prox_data, sel_type == 'sel' & var_type == 'ins', select=c(bin, theta, distance))
colnames(ins_theta) = c('bin', 'ins_theta', 'distance')

cor.test(as.numeric(ins_theta$ins_theta), ins_theta$distance, method='spearman', exact=NULL)

del_theta = subset(prox_data, sel_type == 'sel' & var_type == 'del', select=c(bin, theta, distance))
colnames(del_theta) = c('bin', 'del_theta', 'distance')

cor.test(as.numeric(del_theta$del_theta), del_theta$distance, method='spearman', exact=NULL)

rdi_data = as.data.frame(cbind(ins_theta, del_theta$del_theta))
rdi_data$rdi = rdi_data$del_theta / rdi_data$ins_theta

cor.test(as.numeric(rdi_data$rdi), as.numeric(rdi_data$bin), method='spearman', exact=NULL)

theta_plot = ggplot(subset(prox_data, sel_type=='sel'), aes(x=distance, y=theta, colour=var_type)) +
      geom_point(stat='identity', size = 2) +
      theme_bw() +
      geom_smooth(method='lm')+
      scale_colour_manual(values=cbPalette) +
      ylab(expression(theta))  + xlab('Distance from CDS (bp)') +
      theme(legend.title = element_blank(), legend.position = c(0.15, 0.8), legend.background = element_blank())

gamma_plot = ggplot(subset(prox_data, sel_type=='sel'), aes(x=distance, y=mean_gamma, colour=var_type)) +
      geom_point(stat='identity', size = 2) +
      theme_bw() +
      scale_colour_manual(values=cbPalette) +
      ylab(expression(gamma))  + xlab('Distance from CDS (bp)') +
      theme(legend.title = element_blank(), legend.position = 'none', legend.background = element_blank())

rdi_plot = ggplot(rdi_data, aes(x=distance, y=rdi)) +
      geom_point(stat='identity', size = 2) +
      theme_bw() +
      geom_smooth(method='lm') +
      ylab('rDI')  + xlab('Distance from CDS (bp)')

png(file='distance_estimates_5kb_nc.png', width=9, height=3, units='in', res=360)

grid.arrange(theta_plot, gamma_plot, rdi_plot, nrow=1)

dev.off()
