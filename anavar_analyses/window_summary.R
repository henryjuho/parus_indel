# Title     : TODO
# Objective : TODO
# Created by: henryjuho
# Created on: 07/11/2017

library(ggplot2)
library(gridExtra)
library(reshape2)

cbPalette <- c("#E69F00", 'tomato 3', 'steel blue', "#999999", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

window_data = read.delim('filtered_2Mb_windows.txt')

theta_ins_test = cor.test(window_data$theta_ins, window_data$rec_rate, method='spearman', exact=NULL)
tajd_ins_test = cor.test(window_data$tajd_ins, window_data$rec_rate, method='spearman', exact=NULL)

theta_del_test = cor.test(window_data$theta_del, window_data$rec_rate, method='spearman', exact=NULL)
tajd_del_test = cor.test(window_data$tajd_del, window_data$rec_rate, method='spearman', exact=NULL)

# theta plot
theta_data = subset(window_data, select=c(window, rec_rate, theta_ins, theta_del))
theta_data = melt(theta_data, id=c('window', 'rec_rate'))

theta_plot = ggplot(theta_data, aes(x=rec_rate, y=value, colour=variable)) +
      geom_point(stat='identity', size = 2) +
      geom_smooth(method='lm') +
      theme_bw() +
      scale_colour_manual(values=cbPalette) +
      xlab('Recombination rate')  + ylab(expression(theta[w])) +
      ggtitle(paste('Ins: ', expression(rho), '=', round(theta_ins_test$estimate, digits=2),
      'p =', round(theta_ins_test$p.value, digits=5),
      'Del: ', expression(rho), '=', round(theta_del_test$estimate, digits=2),
      'p =', round(theta_del_test$p.value, digits=5)))


# tajd plot
tajd_data = subset(window_data, select=c(window, rec_rate, tajd_ins, tajd_del))
tajd_data = melt(tajd_data, id=c('window', 'rec_rate'))

tajd_plot = ggplot(tajd_data, aes(x=rec_rate, y=value, colour=variable)) +
      geom_point(stat='identity', size = 2) +
      geom_smooth(method='lm') +
      theme_bw() +
      scale_colour_manual(values=cbPalette) +
      xlab('Recombination rate')  + ylab("Tajima's D") +
      ggtitle(paste('Ins: ', expression(rho), '=', round(tajd_ins_test$estimate, digits=2),
      'p =', round(tajd_ins_test$p.value, digits=5),
      'Del: ', expression(rho), '=', round(tajd_del_test$estimate, digits=2),
      'p =', round(tajd_del_test$p.value, digits=5)))

png(file='window_summary.png', width=9, height=6, units='in', res=360)

grid.arrange(theta_plot, tajd_plot, nrow=2)

dev.off()