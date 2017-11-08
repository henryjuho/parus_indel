# Title     : anavar results
# Objective : analyses anavar results and plot
# Created by: henryjuho
# Created on: 08/11/2017

library(ggplot2)
library(gridExtra)
library(reshape2)

cbPalette <- c("#E69F00", 'tomato 3', 'steel blue', "#999999", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

window_data = read.csv('gt_windows_noncoding_continuous_full_results.windowdata.csv')
window_data$mean_gamma = as.numeric(window_data$scale) * as.numeric(window_data$shape) * -1

window_data_sel = subset(window_data, sel_type=='sel')

ins_data = subset(window_data_sel, var_type=='ins')
del_data = subset(window_data_sel, var_type=='del')

# theta correlations
theta_ins_test = cor.test(ins_data$theta, ins_data$rec_rate, method='spearman', exact=NULL)
theta_del_test = cor.test(del_data$theta, del_data$rec_rate, method='spearman', exact=NULL)

# gamma correlations
gamma_ins_test = cor.test(ins_data$mean_gamma, ins_data$rec_rate, method='spearman', exact=NULL)
gamma_del_test = cor.test(del_data$mean_gamma, del_data$rec_rate, method='spearman', exact=NULL)

# error correlations
error_ins_test = cor.test(ins_data$e, ins_data$rec_rate, method='spearman', exact=NULL)
error_del_test = cor.test(del_data$e, del_data$rec_rate, method='spearman', exact=NULL)

# theta plot

theta_plot = ggplot(window_data_sel, aes(x=rec_rate, y=theta, colour=var_type)) +
      geom_point(stat='identity', size = 2) +
      geom_smooth(method='lm') +
      theme_bw() +
      scale_colour_manual(values=cbPalette) +
      xlab('Recombination rate')  + ylab(expression(theta[w])) +
      ggtitle(paste('Ins: ', expression(rho), '=', round(theta_ins_test$estimate, digits=2),
      'p =', round(theta_ins_test$p.value, digits=7),
      'Del: ', expression(rho), '=', round(theta_del_test$estimate, digits=2),
      'p =', round(theta_del_test$p.value, digits=7)))

# rdi plot

rdi_data = as.data.frame(cbind(ins_data$rec_rate, del_data$rec_rate, ins_data$theta, del_data$theta))
colnames(rdi_data) <- c('ins_rec_rate', 'del_rec_rate', 'ins_theta', 'del_theta')
rdi_data$rdi = rdi_data$del_theta / rdi_data$ins_theta

rdi_test = cor.test(rdi_data$ins_rec_rate, rdi_data$rdi, method='spearman', exact=NULL)

rdi_plot = ggplot(rdi_data, aes(x=ins_rec_rate, y=rdi)) +
      geom_point(stat='identity', size = 2) +
      geom_smooth(method='lm') +
      theme_bw() +
      xlab('Recombination rate')  + ylab('rDI') +
      ggtitle(paste(expression(rho), '=', round(rdi_test$estimate, digits=2),
      'p =', round(rdi_test$p.value, digits=7)))

# gamma plot

gamma_plot = ggplot(window_data_sel, aes(x=rec_rate, y=mean_gamma, colour=var_type)) +
      geom_point(stat='identity', size = 2) +
      geom_smooth(method='lm') +
      theme_bw() +
      scale_colour_manual(values=cbPalette) +
      xlab('Recombination rate')  + ylab("Gamma") +
      ggtitle(paste('Ins: ', expression(rho), '=', round(gamma_ins_test$estimate, digits=2),
      'p =', round(gamma_ins_test$p.value, digits=7),
      'Del: ', expression(rho), '=', round(gamma_del_test$estimate, digits=2),
      'p =', round(gamma_del_test$p.value, digits=7)))

# error plot

error_plot = ggplot(window_data_sel, aes(x=rec_rate, y=e, colour=var_type)) +
      geom_point(stat='identity', size = 2) +
      geom_smooth(method='lm') +
      theme_bw() +
      scale_colour_manual(values=cbPalette) +
      xlab('Recombination rate')  + ylab("Polarisation error") +
      ggtitle(paste('Ins: ', expression(rho), '=', round(error_ins_test$estimate, digits=2),
      'p =', round(error_ins_test$p.value, digits=7),
      'Del: ', expression(rho), '=', round(error_del_test$estimate, digits=2),
      'p =', round(error_del_test$p.value, digits=7)))

# saving
png(file='window_anavar.png', width=9, height=6, units='in', res=360)

grid.arrange(theta_plot, gamma_plot, nrow=2)

dev.off()

png(file='rdi_anavar.png', width=9, height=3, units='in', res=360)

rdi_plot

dev.off()

png(file='error_anavar.png', width=9, height=3, units='in', res=360)

error_plot

dev.off()