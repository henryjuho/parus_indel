library(ggplot2)

res = read.delim('all_models_lnl.txt')

png('all_mods_lnl.png', width=12, height=12, units='in', res=320)
ggplot(res, aes(x=rank, y=lnL)) +
    geom_line(stat='identity') +
    facet_wrap(~model, ncol=4, scales = "free") +
    xlab('searches')

dev.off()