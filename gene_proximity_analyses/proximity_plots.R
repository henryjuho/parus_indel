# Title     : TODO
# Objective : TODO
# Created by: henryjuho
# Created on: 28/11/2017

library(ggplot2)
library(gridExtra)
library(gtools)

cbPalette <- c("#E69F00", 'tomato 3', 'steel blue', "#999999", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

prox_data = read.csv('gt_sel_neu_ref_genedist.results.csv')

prox_data$bin = factor(prox_data$bin, levels=mixedsort(as.character(subset(prox_data, sel_type=='sel' & var_type=='ins')$bin)))
prox_data$mean_gamma = as.numeric(as.character(prox_data$shape)) * as.numeric(as.character(prox_data$scale)) * -1

ins_theta = subset(prox_data, sel_type == 'sel' & var_type == 'ins', select=c(bin, theta))
colnames(ins_theta) = c('bin', 'ins_theta')

del_theta = subset(prox_data, sel_type == 'sel' & var_type == 'del', select=c(bin, theta))
colnames(del_theta) = c('bin', 'del_theta')

rdi_data = as.data.frame(cbind(ins_theta, del_theta$del_theta))
rdi_data$rdi = rdi_data$del_theta / rdi_data$ins_theta

theta_plot = ggplot(subset(prox_data, sel_type=='sel'), aes(x=bin, y=theta, colour=var_type)) +
      geom_point(stat='identity', size = 2) +
      theme_bw() +
      scale_colour_manual(values=cbPalette) +
      ylab(expression(theta))  + xlab('Distance from gene / 50kb') +
      theme(legend.title = element_blank(), legend.position = c(0.8, 0.5), legend.background = element_blank())

gamma_plot = ggplot(subset(prox_data, sel_type=='sel'), aes(x=bin, y=mean_gamma, colour=var_type)) +
      geom_point(stat='identity', size = 2) +
      theme_bw() +
      scale_colour_manual(values=cbPalette) +
      ylab(expression(gamma))  + xlab('Distance from gene / 50kb') +
      theme(legend.title = element_blank(), legend.position = c(0.8, 0.5), legend.background = element_blank())

rdi_plot = ggplot(rdi_data, aes(x=bin, y=rdi)) +
      geom_point(stat='identity', size = 2) +
      theme_bw() +
      ylab('rDI')  + xlab('Distance from gene / 50kb')

png(file='distance_estimates.png', width=9, height=3, units='in', res=360)

grid.arrange(theta_plot, gamma_plot, rdi_plot, nrow=1)

dev.off()
