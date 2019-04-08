
# |Genes     |CDS frameshift   |CDS non-frameshift |Ratio   |
# |:---------|:---------------:|:-----------------:|:------:|
# |All       | 973             | 1054              | 0.9    |
# |Conserved | 357             | 512               | 0.7    |

library(ggplot2)

  conserved <- vector()
  all_genes <- vector()
  stat_data <- data.frame()
  
  for(i in 1:357){
    conserved <- c(conserved, 1)
    line = c('con', 1)
    stat_data = rbind(stat_data, line)
  }
  
  for(i in 1:512){
    conserved <- c(conserved, 0)
    line = c('con', 0)
    stat_data = rbind(stat_data, line)
  }
  
  for(i in 1:973){
    all_genes <- c(all_genes, 1)
    line = c('all', 1)
    stat_data = rbind(stat_data, line)
  }
  
  for(i in 1:1054){
    all_genes <- c(all_genes, 0)
    line = c('all', 0)
    stat_data = rbind(stat_data, line)
  }
  
  # bootstrap 
  resam_all_data <- data.frame()
  resam_con_data <- data.frame()
  for(x in 1:10000){
    resampled_all = sample(all_genes, replace=TRUE)
    resampled_con = sample(conserved, replace=TRUE)
    all_counts <- table(resampled_all)
    cons_counts <- table(resampled_con)
    resam_all_data = rbind(resam_all_data, all_counts)
    resam_con_data = rbind(resam_con_data, cons_counts)
  }
  
colnames(resam_all_data) = c('inframe', 'shift')
colnames(resam_con_data) = c('inframe', 'shift')
resam_all_data$prop = resam_all_data$shift / (resam_all_data$shift + resam_all_data$inframe)
resam_con_data$prop = resam_con_data$shift / (resam_con_data$shift + resam_con_data$inframe)

plot_data = data.frame()
all_line = c(c('All genes', mean(resam_all_data$prop)), quantile(resam_all_data$prop, probs=c(0.025, 0.975)))
con_line = c(c('Conserved genes', mean(resam_con_data$prop)), quantile(resam_con_data$prop, probs=c(0.025, 0.975)))
plot_data = rbind(plot_data, all_line, con_line)

colnames(plot_data) = c('category', 'proportion_shift', 'lwr', 'upr')

pdf('conserved_gene_shift.pdf', height=3, width=3)

ggplot(plot_data, aes(x=category, y=1-round(as.numeric(as.character(proportion_shift)), 3))) +
  geom_errorbar(aes(ymin = 1-round(as.numeric(as.character(lwr)), 3), ymax = 1-round(as.numeric(as.character(upr)), 3)),
                stat = 'identity', position = 'dodge', width=0.25) +
  geom_point(stat='identity') +
  theme_bw() +
  labs(x='', y='Proportion of in-frame INDELs')

dev.off()