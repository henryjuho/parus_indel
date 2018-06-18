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
colnames(theta_data) = c('window', 'rec_rate', 'ins', 'del')
theta_data = melt(theta_data, id=c('window', 'rec_rate'))
theta_data$variable = factor(theta_data$variable, levels=c('del', 'ins'))

theta_plot = ggplot(theta_data, aes(x=log(rec_rate + 1), y=value, colour=variable)) +
      geom_point(stat='identity', size = 2) +
      geom_smooth(method='lm') +
      theme_bw() +
      scale_colour_manual(values=cbPalette) +
      xlab('Recombination rate (log)')  + ylab(expression(theta[w])) +
      ggtitle(paste('Ins: ', expression(rho), '=', round(theta_ins_test$estimate, digits=2),
      'p =', round(theta_ins_test$p.value, digits=5),
      '\nDel: ', expression(rho), '=', round(theta_del_test$estimate, digits=2),
      'p =', round(theta_del_test$p.value, digits=5))) +
      theme(legend.position = 'none')


# tajd plot
tajd_data = subset(window_data, select=c(window, rec_rate, tajd_ins, tajd_del))
colnames(tajd_data) = c('window', 'rec_rate', 'ins', 'del')
tajd_data = melt(tajd_data, id=c('window', 'rec_rate'))
tajd_data$variable = factor(tajd_data$variable, levels=c('del', 'ins'))

tajd_plot = ggplot(tajd_data, aes(x=log(rec_rate + 1), y=value, colour=variable)) +
      geom_point(stat='identity', size = 2) +
      geom_smooth(method='lm') +
      theme_bw() +
      scale_colour_manual(values=cbPalette) +
      xlab('Recombination rate (log)')  + ylab("Tajima's D") +
      ggtitle(paste('Ins: ', expression(rho), '=', round(tajd_ins_test$estimate, digits=2),
      'p =', round(tajd_ins_test$p.value, digits=5),
      '\nDel: ', expression(rho), '=', round(tajd_del_test$estimate, digits=2),
      'p =', round(tajd_del_test$p.value, digits=5))) +
      theme(legend.position=c(0.9, 0.2), legend.title=element_blank(), legend.background=element_blank())

png(file='window_summary.png', width=6, height=3, units='in', res=360)

grid.arrange(theta_plot, tajd_plot, nrow=1)

dev.off()