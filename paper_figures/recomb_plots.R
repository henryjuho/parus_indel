library(ggplot2)
library(gridExtra)
library(reshape2)
library(viridis)

# summary stat results
window_data = read.delim('../recombination_analyses/filtered_2Mb_windows.txt')

cor.test(window_data$tajd_ins, window_data$rec_rate, method='spearman', exact=NULL)
cor.test(window_data$tajd_del, window_data$rec_rate, method='spearman', exact=NULL)

# rdi data
rdi_data = subset(window_data, select=c(window, rec_rate, theta_ins, theta_del))
colnames(rdi_data) = c('window', 'rec_rate', 'ins', 'del')
rdi_data$rdi = rdi_data$del / rdi_data$ins
cor.test(rdi_data$rdi, rdi_data$rec_rate, method='spearman', exact=NULL)
rdi_data$filter = 'No filter'

pal = viridis(n=3)[1:3]

# pi plot
pi_data = subset(window_data, select=c(window, rec_rate, pi_ins, pi_del))
colnames(pi_data) = c('window', 'rec_rate', 'ins', 'del')
pi_data = melt(pi_data, id=c('window', 'rec_rate'))
pi_data$variable = factor(pi_data$variable, levels=c('del', 'ins'))

pi_plot = ggplot(pi_data, aes(x=log(rec_rate + 1), y=value*10000, colour=toupper(variable))) +
      geom_point(stat='identity', size = 2) +
      # geom_smooth(method='lm', se=FALSE) +
      theme_bw() +
      xlab('Recombination rate (log)')  + ylab(expression(pi * x *10^-4)) +
      theme(legend.position = 'none', plot.title=element_text(hjust=-0.2, vjust=-3)) +
      ggtitle('(a)')+ scale_colour_manual(values=pal)


# tajd plot
tajd_data = subset(window_data, select=c(window, rec_rate, tajd_ins, tajd_del))
colnames(tajd_data) = c('window', 'rec_rate', 'ins', 'del')
tajd_data = melt(tajd_data, id=c('window', 'rec_rate'))
tajd_data$variable = factor(tajd_data$variable, levels=c('del', 'ins'))

tajd_plot = ggplot(tajd_data, aes(x=log(rec_rate + 1), y=value, colour=toupper(variable))) +
      geom_point(stat='identity', size = 2) +
      # geom_smooth(method='lm', se=FALSE) +
      theme_bw() +
      xlab('Recombination rate (log)')  + ylab("Tajima's D") +
      theme(legend.position=c(0.85, 0.15), legend.title=element_blank(), legend.background=element_blank(),
      plot.title=element_text(hjust=-0.3, vjust=-3))+
      ggtitle('(b)') + scale_colour_manual(values=pal)

pdf('recomb_tajd_pi.pdf', width=6, height=3)

grid.arrange(pi_plot, tajd_plot, nrow=1)

dev.off()

# ## with filters
#
# # summary stat results
# window_data2 = read.delim('../recombination_analyses/filtered_2Mb_windows_percentpol.txt')
# window_data2 = subset(window_data2, pol_success > 0.6)
#
# cor.test(window_data2$tajd_ins, window_data2$rec_rate, method='spearman', exact=NULL)
# cor.test(window_data2$tajd_del, window_data2$rec_rate, method='spearman', exact=NULL)
#
# # rdi plot
# rdi_data2 = subset(window_data2, select=c(window, rec_rate, theta_ins, theta_del))
# colnames(rdi_data2) = c('window', 'rec_rate', 'ins', 'del')
# rdi_data2$rdi = rdi_data2$del / rdi_data2$ins
# cor.test(rdi_data2$rdi, rdi_data2$rec_rate, method='spearman', exact=NULL)
# rdi_data2$filter = '60% min polarised'
#
# rdi_all = rbind(rdi_data, rdi_data2)
#
#
# rdi_plot = ggplot(rdi_all, aes(x=log(rec_rate + 1), y=rdi)) +
#       geom_point(stat='identity', size = 2) +
#       # geom_smooth(method='lm', se=FALSE) +
#       theme_bw() +
#       xlab('Recombination rate (log)')  + ylab('rDI') +
#       theme(plot.title=element_text(hjust=-0.2, vjust=-3)) +
#       facet_wrap(~filter, scale = 'free_x')
#
# pdf('recomb_rdi_filter.pdf', width=6, height=3)
#
# rdi_plot
#
# dev.off()
#
#
#
# # pi plot
# pi_data2 = subset(window_data2, select=c(window, rec_rate, pi_ins, pi_del))
# colnames(pi_data2) = c('window', 'rec_rate', 'ins', 'del')
# pi_data2 = melt(pi_data2, id=c('window', 'rec_rate'))
# pi_data2$variable = factor(pi_data2$variable, levels=c('del', 'ins'))
#
# pi_plot2 = ggplot(pi_data2, aes(x=log(rec_rate + 1), y=value*10000, colour=toupper(variable))) +
#       geom_point(stat='identity', size = 2) +
#       # geom_smooth(method='lm', se=FALSE) +
#       theme_bw() +
#       xlab('Recombination rate (log)')  + ylab(expression(pi * x *10^-4)) +
#       theme(legend.position = 'none', plot.title=element_text(hjust=-0.2, vjust=-3)) +
#       ggtitle('(a)')+ scale_colour_manual(values=pal)
#
#
# # tajd plot
# tajd_data2 = subset(window_data2, select=c(window, rec_rate, tajd_ins, tajd_del))
# colnames(tajd_data2) = c('window', 'rec_rate', 'ins', 'del')
# tajd_data2 = melt(tajd_data2, id=c('window', 'rec_rate'))
# tajd_data2$variable = factor(tajd_data2$variable, levels=c('del', 'ins'))
#
# tajd_plot2 = ggplot(tajd_data2, aes(x=log(rec_rate + 1), y=value, colour=toupper(variable))) +
#       geom_point(stat='identity', size = 2) +
#       # geom_smooth(method='lm', se=FALSE) +
#       theme_bw() +
#       xlab('Recombination rate (log)')  + ylab("Tajima's D") +
#       theme(legend.position=c(0.85, 0.2), legend.title=element_blank(), legend.background=element_blank(),
#       plot.title=element_text(hjust=-0.3, vjust=-3))+
#       ggtitle('(b)') + scale_colour_manual(values=pal)
#
# pdf('recomb_tajd_pi_filter.pdf', width=6, height=3)
#
# grid.arrange(pi_plot2, tajd_plot2, nrow=1)
#
# dev.off()
