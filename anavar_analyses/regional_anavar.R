# Title     : TODO
# Objective : TODO
# Created by: henryjuho
# Created on: 03/11/2017

library(gridExtra)
library(ggplot2)

cbPalette <- c("#E69F00", 'tomato 3', 'steel blue', "#999999", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

results = read.csv('gt_regional_anavar_results.csv')
results$mean_gamma = as.numeric(as.character(results$shape)) * as.numeric(as.character(results$scale)) * -1

results_no_neu = subset(results, sel_type!='neu')

# gamma plot
gamma = ggplot(results_no_neu, aes(x=region, y=mean_gamma, colour=var_type)) +
    geom_point(stat='identity', position = position_dodge(width=0.9), size = 2) +
    theme_bw() +
    scale_color_manual(values=cbPalette) +
    xlab('Genomic region')  + ylab(expression(gamma)) +
    theme(legend.title=element_blank(), legend.position=c(0.8, 0.2))

# theta plot
theta = ggplot(results_no_neu, aes(x=region, y=theta, colour=var_type)) +
    geom_point(stat='identity', position = position_dodge(width=0.9), size = 2) +
    theme_bw() +
    scale_color_manual(values=cbPalette) +
    xlab('Genomic region')  + ylab(expression(theta)) +
    theme(legend.title=element_blank(), legend.position=c(0.8, 0.2))

pdf('regional_anavar.pdf', width=6, height=3)
grid.arrange(gamma, theta, nrow=1)
dev.off()

png('regional_anavar.png', width=6, height=3, res=320, units='in')
grid.arrange(gamma, theta, nrow=1)
dev.off()