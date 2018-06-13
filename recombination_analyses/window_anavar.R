# Title     : anavar results
# Objective : analyses anavar results and plot
# Created by: henryjuho
# Created on: 08/11/2017

library(ggplot2)
library(gridExtra)
library(reshape2)

in_file = 'gt_windows_noncoding_continuous_full_results.windowdata.csv'
out_stem = 'continuous'

cbPalette <- c("#E69F00", 'tomato 3', 'steel blue', "#999999", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

window_data = read.csv(in_file)
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

theta_plot = ggplot(window_data_sel, aes(x=log(rec_rate + 1), y=theta, colour=var_type)) +
      geom_point(stat='identity', size = 2) +
      geom_smooth(method='lm') +
      theme_bw() +
      scale_colour_manual(values=cbPalette) +
      xlab('Recombination rate (log)')  + ylab(expression(theta)) +
      ggtitle(paste('Ins: ', expression(rho), '=', round(theta_ins_test$estimate, digits=2),
      'p =', round(theta_ins_test$p.value, digits=7),
      '\nDel: ', expression(rho), '=', round(theta_del_test$estimate, digits=2),
      'p =', round(theta_del_test$p.value, digits=7))) +
      theme(legend.title = element_blank(), legend.position=c(0.9, 0.8))

# rdi plot

rdi_data = as.data.frame(cbind(ins_data$rec_rate, del_data$rec_rate, ins_data$theta, del_data$theta))
colnames(rdi_data) <- c('ins_rec_rate', 'del_rec_rate', 'ins_theta', 'del_theta')
rdi_data$rdi = rdi_data$del_theta / rdi_data$ins_theta

rdi_test = cor.test(rdi_data$ins_rec_rate, rdi_data$rdi, method='spearman', exact=NULL)

rdi_plot = ggplot(rdi_data, aes(x=log(ins_rec_rate + 1), y=rdi)) +
      geom_point(stat='identity', size = 2) +
      geom_smooth(method='lm') +
      theme_bw() +
      xlab('Recombination rate (log)')  + ylab('\nrDI') +
      ggtitle(paste('\n', expression(rho), '=', round(rdi_test$estimate, digits=2),
      'p =', round(rdi_test$p.value, digits=7)))

# gamma plot

gamma_plot = ggplot(window_data_sel, aes(x=log(rec_rate + 1), y=mean_gamma, colour=var_type)) +
      geom_point(stat='identity', size = 2) +
      geom_smooth(method='lm') +
      theme_bw() +
      scale_colour_manual(values=cbPalette) +
      xlab('Recombination rate (log)')  + ylab("Gamma") +
      ggtitle(paste('Ins: ', expression(rho), '=', round(gamma_ins_test$estimate, digits=2),
      'p =', round(gamma_ins_test$p.value, digits=7),
      '\nDel: ', expression(rho), '=', round(gamma_del_test$estimate, digits=2),
      'p =', round(gamma_del_test$p.value, digits=7))) + theme(legend.position='none')

# error plot

error_plot = ggplot(window_data_sel, aes(x=log(rec_rate + 1), y=e, colour=var_type)) +
      geom_point(stat='identity', size = 2) +
      geom_smooth(method='lm') +
      theme_bw() +
      scale_colour_manual(values=cbPalette) +
      xlab('Recombination rate (log)')  + ylab("Polarisation error") +
      ggtitle(paste('Ins: ', expression(rho), '=', round(error_ins_test$estimate, digits=2),
      'p =', round(error_ins_test$p.value, digits=7),
      '\nDel: ', expression(rho), '=', round(error_del_test$estimate, digits=2),
      'p =', round(error_del_test$p.value, digits=7))) + theme(legend.position='none')

# saving
png(file=paste(out_stem, 'window_anavar.png', sep='_'), width=9, height=6, units='in', res=360)

grid.arrange(theta_plot, gamma_plot, rdi_plot, error_plot, nrow=2, ncol=2)

dev.off()
