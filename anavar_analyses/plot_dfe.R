library(ggplot2)
library(viridis)

dfe = read.delim('dfe_plotdata.txt')

dfe
dfe_coding = dfe
dfe_coding$prop = NULL
dfe_coding$prop = c(0, 0.03882046, 0, 0.9611795, 0, 0.03541607, 0, 0.9645839)

dfe$new_var = paste(dfe$var, 'NC', sep=' ')
dfe_coding$new_var = paste(dfe_coding$var, 'CDS', sep=' ')

all_dfe = rbind(dfe, dfe_coding)

all_dfe$new_var = factor(all_dfe$new_var, levels=c('DEL NC', 'INS NC', 'DEL CDS', 'INS CDS'))

pdf('nc_dfe.pdf', height=3, width=3)

ggplot(all_dfe, aes(x=cat, y=prop, fill=new_var)) +
    geom_bar(stat='identity', position='dodge') +
    theme_bw() +
    scale_fill_manual(values=viridis(5)) +
    xlab(expression(-gamma)) + ylab('Density') +
    theme(legend.title=element_blank(), legend.position=c(0.3, 0.9),
    legend.key.size = unit(.4, "cm"), legend.text=element_text(size=8)) +
    guides(fill=guide_legend(ncol=2))

dev.off()