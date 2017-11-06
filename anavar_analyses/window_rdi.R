# Title     : rdi plotting
# Objective : plot rdi against recomb per window
# Created by: henryjuho
# Created on: 06/11/2017

library(ggplot2)

window_data = read.delim('filtered_2Mb_windows.txt')
window_data$rdi = as.numeric(as.character(window_data$n_del)) / as.numeric(as.character(window_data$n_ins))

cor_res = cor.test(window_data$rec_rate, window_data$rdi, method='spearman', exact=NULL)

rdi_plot = ggplot(window_data, aes(x=rec_rate, y=rdi)) +
      geom_point(stat='identity', size = 2) +
      theme_bw() +
      xlab('Recombination rate')  + ylab('rDI') +
      ggtitle(paste('rho =', cor_res$estimate, 'p =', cor_res$p.value))

png(file='window_rdi.png', width=9, height=3, units='in', res=360)

rdi_plot

dev.off()