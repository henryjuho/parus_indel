library(ggplot2)
library(gridExtra)
library(gtools)

args = commandArgs(TRUE)

in_file = args[[1]]
out_file = args[[2]]
window_size = as.numeric(args[[3]])

prox_data = read.csv(in_file)

prox_data$bin = factor(prox_data$bin, levels=mixedsort(as.character(subset(prox_data, sel_type=='sel' & var_type=='ins')$bin)))
prox_data$mean_gamma = as.numeric(as.character(prox_data$shape)) * as.numeric(as.character(prox_data$scale)) * -1
prox_data$distance = as.numeric(prox_data$bin) * window_size

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
      #xlim(0,10000)+
      ylab(expression(theta))  + xlab('Distance from conserved region (bp)') +
      theme(legend.title = element_blank(), legend.position = c(0.15, 0.8), legend.background = element_blank())

gamma_plot = ggplot(subset(prox_data, sel_type=='sel'), aes(x=distance, y=mean_gamma, colour=var_type)) +
      geom_point(stat='identity', size = 2) +
      theme_bw() +
      #xlim(0,10000)+
      ylab(expression(gamma))  + xlab('Distance from conserved region (bp)') +
      theme(legend.title = element_blank(), legend.position = 'none', legend.background = element_blank())

rdi_plot = ggplot(rdi_data, aes(x=distance, y=rdi)) +
      geom_point(stat='identity', size = 2) +
      theme_bw() +
      #xlim(0,10000)+
      geom_smooth(method='lm') +
      ylab('rDI')  + xlab('Distance from conserved region (bp)')

png(file=out_file, width=9, height=3, units='in', res=360)

grid.arrange(theta_plot, gamma_plot, rdi_plot, nrow=1)

dev.off()
