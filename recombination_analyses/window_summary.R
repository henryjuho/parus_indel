library(ggplot2)
library(gridExtra)
library(reshape2)
library(viridis)
library(dplyr)
library(ppcor)

opts <- commandArgs(trailingOnly = TRUE)
in_dat = opts[1]
out_plot = opts[2]

window_data = read.delim(in_dat)

theta_ins_test = cor.test(window_data$theta_ins, window_data$rec_rate, method='spearman', exact=NULL)
tajd_ins_test = cor.test(window_data$tajd_ins, window_data$rec_rate, method='spearman', exact=NULL)
pi_ins_test = cor.test(window_data$pi_ins, window_data$rec_rate, method='spearman', exact=NULL)

theta_del_test = cor.test(window_data$theta_del, window_data$rec_rate, method='spearman', exact=NULL)
tajd_del_test = cor.test(window_data$tajd_del, window_data$rec_rate, method='spearman', exact=NULL)
pi_del_test = cor.test(window_data$pi_del, window_data$rec_rate, method='spearman', exact=NULL)

# theta plot
theta_data = subset(window_data, select=c(window, rec_rate, theta_ins, theta_del))
colnames(theta_data) = c('window', 'rec_rate', 'ins', 'del')
theta_data = melt(theta_data, id=c('window', 'rec_rate'))
theta_data$variable = factor(theta_data$variable, levels=c('del', 'ins'))

theta_plot = ggplot(theta_data, aes(x=log(rec_rate + 1), y=value, colour=variable)) +
      geom_point(stat='identity', size = 2) +
      theme_bw() +
      xlab('Recombination rate (log)')  + ylab(expression(theta[w])) +
      ggtitle(paste('Ins: ', expression(rho), '=', round(theta_ins_test$estimate, digits=2),
      'p =', round(theta_ins_test$p.value, digits=5),
      '\nDel: ', expression(rho), '=', round(theta_del_test$estimate, digits=2),
      'p =', round(theta_del_test$p.value, digits=5))) +
      theme(legend.position = 'none')+
      scale_colour_manual(values=viridis(3))

# pi plot
pi_data = subset(window_data, select=c(window, rec_rate, pi_ins, pi_del))
colnames(pi_data) = c('window', 'rec_rate', 'ins', 'del')
pi_data = melt(pi_data, id=c('window', 'rec_rate'))
pi_data$variable = factor(pi_data$variable, levels=c('del', 'ins'))

pi_plot = ggplot(pi_data, aes(x=log(rec_rate + 1), y=value, colour=variable)) +
      geom_point(stat='identity', size = 2) +
      theme_bw() +
      xlab('Recombination rate (log)')  + ylab(expression(pi)) +
      ggtitle(paste('Ins: ', expression(rho), '=', round(pi_ins_test$estimate, digits=2),
      'p =', round(pi_ins_test$p.value, digits=5),
      '\nDel: ', expression(rho), '=', round(pi_del_test$estimate, digits=2),
      'p =', round(pi_del_test$p.value, digits=5))) +
      theme(legend.position = 'none')+
      scale_colour_manual(values=viridis(3))


# rdi plot
rdi_data = subset(window_data, select=c(window, rec_rate, theta_ins, theta_del))
# rdi_data = subset(window_data, select=c(window, rec_rate, pi_ins, pi_del))
colnames(rdi_data) = c('window', 'rec_rate', 'ins', 'del')
rdi_data$rdi = rdi_data$del / rdi_data$ins
rdi_test = cor.test(rdi_data$rdi, rdi_data$rec_rate, method='spearman', exact=NULL)

rdi_plot = ggplot(rdi_data, aes(x=log(rec_rate + 1), y=rdi)) +
      geom_point(stat='identity', size = 2) +
      theme_bw() +
      xlab('Recombination rate (log)')  + ylab('rDI') +
      ggtitle(paste(expression(rho), '=', round(rdi_test$estimate, digits=2),
      'p =', round(rdi_test$p.value, digits=5)))+
      scale_colour_manual(values=viridis(3))


# tajd plot
tajd_data = subset(window_data, select=c(window, rec_rate, tajd_ins, tajd_del))
colnames(tajd_data) = c('window', 'rec_rate', 'ins', 'del')
tajd_data = melt(tajd_data, id=c('window', 'rec_rate'))
tajd_data$variable = factor(tajd_data$variable, levels=c('del', 'ins'))

tajd_plot = ggplot(tajd_data, aes(x=log(rec_rate + 1), y=value, colour=variable)) +
      geom_point(stat='identity', size = 2) +
      theme_bw() +
      xlab('Recombination rate (log)')  + ylab("Tajima's D") +
      ggtitle(paste('Ins: ', expression(rho), '=', round(tajd_ins_test$estimate, digits=2),
      'p =', round(tajd_ins_test$p.value, digits=5),
      '\nDel: ', expression(rho), '=', round(tajd_del_test$estimate, digits=2),
      'p =', round(tajd_del_test$p.value, digits=5))) +
      theme(legend.position=c(0.9, 0.2), legend.title=element_blank(), legend.background=element_blank()) +
      scale_colour_manual(values=viridis(3))

png(file=out_plot, width=9, height=3, units='in', res=360)

grid.arrange(theta_plot, pi_plot, tajd_plot, nrow=1)

dev.off()

# divergence plot for reviewer
div_data = read.delim('filtered_2Mb_windows_with_subs.txt')
div_data = subset(div_data, select=c('window', 'rec_rate', 'n_ins_sub', 'n_del_sub'))

call_sites = subset(window_data, select=c('window', 'callable'))
div_data = dplyr::left_join(div_data, call_sites, by='window')
div_data$ins = div_data$n_ins_sub / div_data$callable
div_data$del = div_data$n_del_sub / div_data$callable

div_ins_test = cor.test(div_data$rec_rate, div_data$ins, method='spearman')
div_del_test = cor.test(div_data$rec_rate, div_data$del, method='spearman')

div_data = melt(subset(div_data, select=c('window', 'rec_rate', 'ins', 'del')),
                id=c('window', 'rec_rate'), variable='var_type')

div_data$var_type = factor(div_data$var_type, levels=c("del", "ins"))
str(div_data)

div_plot = ggplot(div_data, aes(x=log(rec_rate+1), y=value, colour=var_type)) +
               geom_point(stat='identity', size = 2) +
               theme_bw() +
               labs(x='Recombination rate (log)', y='INDEL divergence') +
               theme(legend.title=element_blank(), legend.position='None') +
               scale_colour_manual(values=viridis(3)) +
               ggtitle(paste('Ins: ', expression(rho), '=', round(div_ins_test$estimate, digits=2),
              'p =', round(div_ins_test$p.value, digits=5),
              '\nDel: ', expression(rho), '=', round(div_del_test$estimate, digits=2),
              'p =', round(div_del_test$p.value, digits=5)))

# png('div_recomb_plot.png', width=6, height=3, res=320, units='in')
#
# grid.arrange(div_plot, pi_plot, nrow=1)
#
# dev.off()

# partial correlation

colnames(div_data) = c('window', 'rec_rate', 'var_type', 'div')
colnames(pi_data) = c('window', 'rec_rate', 'var_type', 'pi')
div_pi_dat = dplyr::left_join(div_data, pi_data, by=c('window', 'rec_rate', 'var_type'))

# ins
print('insertions')
ins_pi_div = subset(div_pi_dat, var_type=='ins')
partial_cor_ins = pcor.test(ins_pi_div$pi, ins_pi_div$rec_rate, ins_pi_div$div, method='spearman')
print(partial_cor_ins)

# del
print('deletions')
del_pi_div = subset(div_pi_dat, var_type=='del')
partial_cor_del = pcor.test(del_pi_div$pi, del_pi_div$rec_rate, del_pi_div$div, method='spearman')
print(partial_cor_del)