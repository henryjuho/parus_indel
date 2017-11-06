# Title     : rdi plotting
# Objective : plot rdi against recomb per window
# Created by: henryjuho
# Created on: 06/11/2017

library(ggplot2)
library(gridExtra)

window_data = read.delim('filtered_2Mb_windows_with_subs.txt')
window_data$rdi = as.numeric(as.character(window_data$n_del)) / as.numeric(as.character(window_data$n_ins))
window_data$sub_rdi = as.numeric(as.character(window_data$n_del_sub)) / as.numeric(as.character(window_data$n_ins_sub))

cor_res = cor.test(window_data$rec_rate, window_data$rdi, method='spearman', exact=NULL)
sub_cor_res = cor.test(window_data$rec_rate, window_data$sub_rdi, method='spearman', exact=NULL)
rdi_test = cor.test(window_data$rdi, window_data$sub_rdi, method='spearman', exact=NULL)

rdi_plot = ggplot(window_data, aes(x=rec_rate, y=rdi)) +
      geom_point(stat='identity', size = 2) +
      theme_bw() +
      xlab('Recombination rate')  + ylab('rDI') +
      ggtitle(paste('Polymorphisms: rho =', round(cor_res$estimate, digits=3), 'p =', round(cor_res$p.value, digits=5)))

sub_rdi_plot = ggplot(window_data, aes(x=rec_rate, y=sub_rdi)) +
      geom_point(stat='identity', size = 2) +
      theme_bw() +
      xlab('Recombination rate')  + ylab('rDI') +
      ggtitle(paste('Substitutions: rho =', round(sub_cor_res$estimate, digits=3), 'p =', round(sub_cor_res$p.value, digits=5)))

rdi_rdi_plot = ggplot(window_data, aes(x=rdi, y=sub_rdi)) +
      geom_point(stat='identity', size = 2) +
      theme_bw() +
      xlab('rDI polymorphisms')  + ylab('rDI substitutions') +
      ggtitle(paste('rho =', round(rdi_test$estimate, digits=3), 'p =', round(rdi_test$p.value, digits=7)))


png(file='window_rdi.png', width=9, height=6, units='in', res=360)

rec_rates = grid.arrange(rdi_plot, sub_rdi_plot, nrow=2)
grid.arrange(rec_rates, rdi_rdi_plot, nrow=1)

dev.off()