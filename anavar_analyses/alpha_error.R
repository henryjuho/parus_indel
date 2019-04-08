library(ggplot2)
library(viridis)

alpha_data = read.delim('alpha_change_with_error.txt')

alpha_plot = ggplot(alpha_data, aes(x=percent_bias_adjust, y=alpha, colour=toupper(var_type))) +
    geom_point(stat='identity') +
    theme_bw() +
    scale_colour_manual(values=viridis(n=3)) +
    theme(legend.title=element_blank(), legend.position=c(0.85, 0.15), legend.background=element_blank()) +
    xlab('Percentage shift to deletions') +
    ylab(expression(alpha))

pdf('alpha_error.pdf', width=3, height=3)

alpha_plot

dev.off()
